"""
Partial Mantel test: does geographic distance predict disease co-reporting strength
after controlling for direct air-travel connectivity between countries?

The plain Mantel test in mantel_test.py shows distance is significantly associated with
co-reporting, but does not rule out that this is really just a travel-connectivity effect
(close countries tend to also be well-connected by direct flights, and flight connectivity
is itself plausibly tied to shared trade, tourism, and surveillance-reporting practices).
The partial Mantel test (Smouse, Long, and Sokal, 1986) asks whether distance and
co-reporting remain correlated once their shared association with a third matrix, here
flight-route connectivity (built in build_connectivity_matrix.py), is removed.

Partial correlation: r_XY.Z = (r_XY - r_XZ * r_YZ) / sqrt((1 - r_XZ^2) * (1 - r_YZ^2))
where X = distance, Y = co-reporting strength, Z = connectivity, each vectorized over the
upper-triangular entries of their respective matrices.

Significance is assessed by permuting the row/column labels of the distance matrix (X)
1,999 times, holding Y and Z fixed, and recomputing the partial correlation each time —
the permutation scheme recommended for partial Mantel tests (Legendre, 2000), since
permuting X alone tests whether the X-Y relationship persists beyond what X shares with Z,
without disturbing the genuine Y-Z relationship being controlled for.

This script also repeats the same procedure within each disease transmission-mode stratum
(vector-borne, respiratory, waterborne/foodborne; see disease_classification.py), to test
whether the distance effect, net of connectivity, differs by how a disease is transmitted.
"""

import itertools
import json

import numpy as np
import pandas as pd

from disease_classification import classify

RNG_SEED = 42
N_PERMUTATIONS = 1999  # partial Mantel is more compute-heavy per permutation; fewer reps

EARTH_RADIUS_KM = 6371.0


def haversine_matrix(lat, lon):
    lat_r = np.radians(lat.to_numpy())[:, None]
    lon_r = np.radians(lon.to_numpy())[:, None]
    dlat = lat_r - lat_r.T
    dlon = lon_r - lon_r.T
    a = np.sin(dlat / 2) ** 2 + np.cos(lat_r) * np.cos(lat_r.T) * np.sin(dlon / 2) ** 2
    return 2 * EARTH_RADIUS_KM * np.arcsin(np.sqrt(np.clip(a, 0, 1)))


def co_report_matrix(outbreaks, idx, n):
    mat = np.zeros((n, n), dtype=int)
    for (_, _), group in outbreaks.groupby(["Disease", "Year"]):
        countries = sorted(set(group["iso3"]) & idx.keys())
        for a, b in itertools.combinations(countries, 2):
            i, j = idx[a], idx[b]
            mat[i, j] += 1
            mat[j, i] += 1
    return mat


def partial_correlation(x, y, z):
    r_xy = np.corrcoef(x, y)[0, 1]
    r_xz = np.corrcoef(x, z)[0, 1]
    r_yz = np.corrcoef(y, z)[0, 1]
    denom = np.sqrt((1 - r_xz**2) * (1 - r_yz**2))
    return (r_xy - r_xz * r_yz) / denom, r_xy, r_xz, r_yz


def partial_mantel(dist, co_report, connectivity, n, rng):
    triu = np.triu_indices(n, k=1)
    x, y, z = dist[triu], co_report[triu], connectivity[triu]
    observed_partial_r, r_xy, r_xz, r_yz = partial_correlation(x, y, z)

    perm_r = np.empty(N_PERMUTATIONS)
    for k in range(N_PERMUTATIONS):
        perm = rng.permutation(n)
        x_perm = dist[np.ix_(perm, perm)][triu]
        perm_r[k], _, _, _ = partial_correlation(x_perm, y, z)

    p_value = (np.sum(np.abs(perm_r) >= np.abs(observed_partial_r)) + 1) / (N_PERMUTATIONS + 1)
    return {
        "partial_r": float(observed_partial_r),
        "partial_mantel_p": float(p_value),
        "raw_r_distance_vs_coreport": float(r_xy),
        "raw_r_distance_vs_connectivity": float(r_xz),
        "raw_r_connectivity_vs_coreport": float(r_yz),
        "null_mean": float(perm_r.mean()),
        "null_sd": float(perm_r.std()),
    }


# ---------------------------------------------------------------------
# Load data, restricted to the same country set used in mantel_test.py and
# build_connectivity_matrix.py (loaded from file to guarantee identical ordering).
# ---------------------------------------------------------------------
with open("connectivity_country_order.txt") as f:
    common_iso3 = f.read().splitlines()
idx = {iso: i for i, iso in enumerate(common_iso3)}
n = len(common_iso3)

coords = pd.read_csv("country-coord.csv").dropna(subset=["lat", "lon"]).drop_duplicates(subset="iso3")
coords = coords.set_index("iso3").loc[common_iso3]
dist = haversine_matrix(coords["lat"], coords["lon"])

connectivity = np.load("connectivity_matrix.npy").astype(float)

outbreaks_all = pd.read_csv("Outbreaks.csv").dropna(subset=["iso3", "Disease", "Year"]).drop_duplicates(
    subset=["iso3", "Disease", "Year"]
)
outbreaks_all = outbreaks_all[outbreaks_all["iso3"].isin(idx)]
outbreaks_all["transmission_mode"] = outbreaks_all["Disease"].apply(classify)

rng = np.random.default_rng(RNG_SEED)
results = {}

# ---------------------------------------------------------------------
# Full network
# ---------------------------------------------------------------------
co_full = co_report_matrix(outbreaks_all, idx, n)
triu = np.triu_indices(n, k=1)
print(f"Full network: {n} countries, {(co_full[triu] > 0).sum()} pairs with shared disease-year")
print(f"Connectivity: {(connectivity[triu] > 0).sum()} pairs with at least one direct flight route")

res_full = partial_mantel(dist, co_full, connectivity, n, rng)
results["full_network"] = res_full
print("\n=== Full network, partial Mantel (controlling for flight connectivity) ===")
for k, v in res_full.items():
    print(f"  {k}: {v:.4f}" if isinstance(v, float) else f"  {k}: {v}")

# ---------------------------------------------------------------------
# Stratified by transmission mode
# ---------------------------------------------------------------------
for stratum in ["vector-borne", "respiratory", "waterborne/foodborne"]:
    sub = outbreaks_all[outbreaks_all["transmission_mode"] == stratum]
    co_sub = co_report_matrix(sub, idx, n)
    n_pairs = int((co_sub[triu] > 0).sum())
    n_countries_in_stratum = len(set(sub["iso3"]))
    print(f"\n=== Stratum: {stratum} ({len(sub)} records, "
          f"{n_countries_in_stratum} countries, {n_pairs} co-reporting pairs) ===")
    if n_pairs < 30:
        print("  Too few co-reporting pairs for a meaningful permutation test; skipping.")
        results[stratum] = {"skipped": True, "n_pairs": n_pairs}
        continue
    res = partial_mantel(dist, co_sub, connectivity, n, rng)
    res["n_records"] = len(sub)
    res["n_countries"] = n_countries_in_stratum
    res["n_pairs_with_shared_disease_year"] = n_pairs
    results[stratum] = res
    for k, v in res.items():
        print(f"  {k}: {v:.4f}" if isinstance(v, float) else f"  {k}: {v}")

with open("results_partial_mantel_test.json", "w") as f:
    json.dump(results, f, indent=2)
print("\nResults written to results_partial_mantel_test.json")
