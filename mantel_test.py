"""
Mantel test: does geographic distance between countries predict the strength
of their disease co-reporting relationship?

The original analysis observed that closer countries appear more strongly
connected in the weighted distance network and in the clustering-coefficient
map, but never tested this statistically (it was flagged explicitly as
future work). A naive row-wise correlation or regression between distance
and co-reporting count is not valid here because both quantities are defined
over the same set of countries: every country's distances and co-reporting
counts to all others are not independent observations, they are entries in
two matrices indexed by the same nodes. The Mantel test is the standard tool
for exactly this situation: it permutes the labeling of one matrix relative
to the other to build a null distribution for the correlation between them,
which respects the dependency structure that a naive test would violate.

Method: build (1) the symmetric matrix of pairwise great-circle distances
between all countries with coordinate data, and (2) the symmetric matrix of
co-reporting strength, the number of distinct disease-year combinations two
countries both reported. Compute the Pearson correlation between the
upper-triangular entries of the two matrices, then assess significance by
randomly permuting the row/column labels of one matrix 9,999 times and
recomputing the correlation each time, giving an empirical p-value.
"""

import itertools
import numpy as np
import pandas as pd

RNG_SEED = 42
N_PERMUTATIONS = 9999

# ---------------------------------------------------------------------
# 1. Geographic distance matrix (Haversine, matching the original paper)
# ---------------------------------------------------------------------
coords = pd.read_csv("country-coord.csv")
coords = coords.dropna(subset=["lat", "lon"]).drop_duplicates(subset="iso3")
coords = coords.set_index("iso3")[["lat", "lon"]]

EARTH_RADIUS_KM = 6371.0


def haversine_matrix(lat, lon):
    lat_r = np.radians(lat.to_numpy())[:, None]
    lon_r = np.radians(lon.to_numpy())[:, None]
    dlat = lat_r - lat_r.T
    dlon = lon_r - lon_r.T
    a = np.sin(dlat / 2) ** 2 + np.cos(lat_r) * np.cos(lat_r.T) * np.sin(dlon / 2) ** 2
    return 2 * EARTH_RADIUS_KM * np.arcsin(np.sqrt(np.clip(a, 0, 1)))


# ---------------------------------------------------------------------
# 2. Co-reporting strength matrix (shared disease-year combinations)
# ---------------------------------------------------------------------
outbreaks = pd.read_csv("Outbreaks.csv")
outbreaks = outbreaks.dropna(subset=["iso3", "Disease", "Year"]).drop_duplicates(
    subset=["iso3", "Disease", "Year"]
)

# Restrict to countries present in both datasets.
common_iso3 = sorted(set(outbreaks["iso3"]) & set(coords.index))
coords = coords.loc[common_iso3]
n = len(common_iso3)
idx = {iso: i for i, iso in enumerate(common_iso3)}

co_report = np.zeros((n, n), dtype=int)
for (_, _), group in outbreaks[outbreaks["iso3"].isin(common_iso3)].groupby(["Disease", "Year"]):
    countries = sorted(set(group["iso3"]))
    for a, b in itertools.combinations(countries, 2):
        i, j = idx[a], idx[b]
        co_report[i, j] += 1
        co_report[j, i] += 1

dist = haversine_matrix(coords["lat"], coords["lon"])

print(f"Countries included: {n}")
print(f"Total country pairs: {n * (n - 1) // 2}")
print(f"Country pairs with at least one shared disease-year: {(co_report[np.triu_indices(n, k=1)] > 0).sum()}")

# ---------------------------------------------------------------------
# 3. Mantel test
# ---------------------------------------------------------------------
triu = np.triu_indices(n, k=1)
dist_vec = dist[triu]
co_vec = co_report[triu]

observed_r = np.corrcoef(dist_vec, co_vec)[0, 1]
print(f"\nObserved Pearson correlation (distance vs. co-reporting strength): {observed_r:.4f}")

rng = np.random.default_rng(RNG_SEED)
perm_r = np.empty(N_PERMUTATIONS)
for k in range(N_PERMUTATIONS):
    perm = rng.permutation(n)
    dist_perm = dist[np.ix_(perm, perm)][triu]
    perm_r[k] = np.corrcoef(dist_perm, co_vec)[0, 1]

p_value = (np.sum(np.abs(perm_r) >= np.abs(observed_r)) + 1) / (N_PERMUTATIONS + 1)
print(f"Mantel permutation test ({N_PERMUTATIONS} permutations): p = {p_value:.5f}")
print(f"Null distribution of r: mean = {perm_r.mean():.4f}, sd = {perm_r.std():.4f}")

# ---------------------------------------------------------------------
# 4. Spearman version as a robustness check (monotonic, not just linear)
# ---------------------------------------------------------------------
from scipy.stats import rankdata

dist_rank = rankdata(dist_vec)
co_rank = rankdata(co_vec)
observed_rho = np.corrcoef(dist_rank, co_rank)[0, 1]

perm_rho = np.empty(N_PERMUTATIONS)
for k in range(N_PERMUTATIONS):
    perm = rng.permutation(n)
    dist_perm = dist[np.ix_(perm, perm)][triu]
    perm_rho[k] = np.corrcoef(rankdata(dist_perm), co_rank)[0, 1]

p_value_rho = (np.sum(np.abs(perm_rho) >= np.abs(observed_rho)) + 1) / (N_PERMUTATIONS + 1)
print(f"\nObserved Spearman correlation: {observed_rho:.4f}")
print(f"Mantel permutation test (Spearman, {N_PERMUTATIONS} permutations): p = {p_value_rho:.5f}")

# ---------------------------------------------------------------------
# Save results
# ---------------------------------------------------------------------
import json

results = {
    "n_countries": n,
    "n_pairs": int(n * (n - 1) // 2),
    "n_pairs_with_shared_disease_year": int((co_vec > 0).sum()),
    "pearson_r": float(observed_r),
    "pearson_mantel_p": float(p_value),
    "spearman_rho": float(observed_rho),
    "spearman_mantel_p": float(p_value_rho),
    "n_permutations": N_PERMUTATIONS,
}
with open("results_mantel_test.json", "w") as f:
    json.dump(results, f, indent=2)
print("\nResults written to results_mantel_test.json")
