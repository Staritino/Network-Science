"""
Builds a country-level air-travel connectivity matrix from OpenFlights data.

This matrix is the confound for the partial Mantel test: countries that are close
geographically are often also well connected by direct flights, and flight connectivity
is itself a plausible driver of shared disease reporting (independent of raw distance) via
shared trade, tourism, and surveillance-reporting ties between well-connected countries.
Controlling for this lets the partial Mantel test ask whether distance still predicts
co-reporting once travel connectivity is accounted for.

Source: OpenFlights (https://openflights.org/data.php), airports.dat and routes.dat,
fetched directly from the project's public GitHub mirror. Connectivity between two
countries is defined as the number of distinct scheduled routes (airline + source + dest
airport combinations) connecting any airport in one country to any airport in the other,
summed over both directions. This is a coarse proxy (it counts route legs, not passenger
volume, and OpenFlights' route coverage is known to be incomplete for some regions), which
is disclosed in the paper rather than treated as ground truth.
"""

import csv
import numpy as np
import pandas as pd
import pycountry

# ---------------------------------------------------------------------
# 1. Map OpenFlights country names to ISO3 (manual overrides for known
#    OpenFlights-specific naming that pycountry's fuzzy search cannot resolve)
# ---------------------------------------------------------------------
MANUAL_OVERRIDES = {
    "Burma": "MMR", "Congo (Brazzaville)": "COG", "Congo (Kinshasa)": "COD",
    "Cote d'Ivoire": "CIV", "East Timor": "TLS", "South Korea": "KOR",
    "North Korea": "PRK", "Vietnam": "VNM", "Laos": "LAO", "Macau": "MAC",
    "Reunion": "REU", "Virgin Islands": "VIR", "British Virgin Islands": "VGB",
    "Czech Republic": "CZE", "Swaziland": "SWZ", "Cape Verde": "CPV",
    "Saint Helena": "SHN", "Wallis and Futuna": "WLF",
    "Western Sahara": "ESH", "Kosovo": "XKX", "Brunei": "BRN",
    "Falkland Islands": "FLK", "Faroe Islands": "FRO", "Christmas Island": "CXR",
    "Cocos (Keeling) Islands": "CCK", "Heard Island and McDonald Islands": "HMD",
    "Saint Pierre and Miquelon": "SPM", "South Georgia and the Islands": "SGS",
    "Sint Maarten": "SXM", "Saint Martin": "MAF", "Saint Barthelemy": "BLM",
    "Turkey": "TUR", "West Bank": "PSE",
    # Not in any modern ISO3 list / no permanent population center; excluded.
    "Netherlands Antilles": None, "Midway Islands": None, "Wake Island": None,
    "Johnston Atoll": None,
}


def country_name_to_iso3(name):
    if name in MANUAL_OVERRIDES:
        return MANUAL_OVERRIDES[name]
    try:
        return pycountry.countries.search_fuzzy(name)[0].alpha_3
    except LookupError:
        return None


# ---------------------------------------------------------------------
# 2. Airport ID/IATA -> country ISO3
# ---------------------------------------------------------------------
airport_to_country = {}
with open("airports_raw.dat", encoding="utf-8", errors="replace") as f:
    for row in csv.reader(f):
        if len(row) < 5:
            continue
        airport_id, country_name, iata = row[0], row[3], row[4]
        iso3 = country_name_to_iso3(country_name)
        if iso3 is None:
            continue
        airport_to_country[airport_id] = iso3
        if iata and iata != "\\N":
            airport_to_country[iata] = iso3

# ---------------------------------------------------------------------
# 3. Routes -> country-pair counts
# ---------------------------------------------------------------------
outbreaks = pd.read_csv("Outbreaks.csv").dropna(subset=["iso3", "Disease", "Year"])
coords = pd.read_csv("country-coord.csv").dropna(subset=["lat", "lon"]).drop_duplicates(subset="iso3")
common_iso3 = sorted(set(outbreaks["iso3"]) & set(coords["iso3"]))
idx = {iso: i for i, iso in enumerate(common_iso3)}
n = len(common_iso3)

route_count = np.zeros((n, n), dtype=int)
unresolved_routes = 0
with open("routes_raw.dat", encoding="utf-8", errors="replace") as f:
    for row in csv.reader(f):
        if len(row) < 6:
            continue
        src_airport_id, dst_airport_id = row[3], row[5]
        src_country = airport_to_country.get(src_airport_id)
        dst_country = airport_to_country.get(dst_airport_id)
        if src_country is None or dst_country is None:
            unresolved_routes += 1
            continue
        if src_country not in idx or dst_country not in idx:
            continue
        if src_country == dst_country:
            continue
        i, j = idx[src_country], idx[dst_country]
        route_count[i, j] += 1
        route_count[j, i] += 1

print(f"Countries in connectivity matrix: {n}")
print(f"Routes with unresolved airport-to-country mapping: {unresolved_routes} of "
      f"{sum(1 for _ in open('routes_raw.dat'))}")
triu = np.triu_indices(n, k=1)
n_connected_pairs = (route_count[triu] > 0).sum()
print(f"Country pairs with at least one direct route: {n_connected_pairs} of {len(triu[0])}")

np.save("connectivity_matrix.npy", route_count)
with open("connectivity_country_order.txt", "w") as f:
    f.write("\n".join(common_iso3))
print("Saved connectivity_matrix.npy and connectivity_country_order.txt")
