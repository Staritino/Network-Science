# Analyzing the Network of Disease Distribution by Country and Distance

## Overview
This repository builds and analyzes a global network of infectious disease co-reporting
across 233 countries and seventy diseases (1996-2022), where two countries are linked if
they reported the same disease in the same year. The analysis covers temporal centrality
evolution, k-core decomposition, and geospatial structure, and statistically tests whether
geographic distance between countries predicts the strength of their disease co-reporting
relationship.

## Motivation
Closer countries appear, descriptively, to be more strongly connected in the disease
co-reporting network. This repository tests that directly with a Mantel permutation test,
the appropriate tool when both quantities being compared (distance and co-reporting
strength) are defined over the same set of nodes, which violates the independence
assumption of a naive correlation test, and then tests whether the result survives
controlling for air-travel connectivity with a partial Mantel test, since geographically
close countries are also often well connected by direct flights.

## Methods
- Network construction: undirected weighted graph over 233 countries, edges weighted by
  shared disease-year reports
- Degree, betweenness, closeness, and eigenvector centrality, computed per year (1996-2022)
  and in aggregate
- K-core decomposition of the aggregate network
- Haversine great-circle distance between all country pairs from latitude/longitude data
- Mantel permutation test (Pearson and Spearman, 9,999 permutations) of geographic distance
  against co-reporting strength, with a pre-2020 robustness check excluding the COVID-19
  pandemic years
- Partial Mantel test (1,999 permutations) of distance against co-reporting strength,
  controlling for an air-travel connectivity matrix built from OpenFlights route data,
  run on the full network and within three disease transmission-mode strata (vector-borne,
  respiratory, waterborne/foodborne)

## Requirements
Python 3.9+ with `pandas`, `numpy`, `scipy`, `matplotlib`, `networkx`, `pycountry`. Install
with:
```bash
python3 -m pip install pandas numpy scipy matplotlib networkx pycountry
```

## How to Run
The original temporal and geospatial network analysis is in `ANDS - Python Code.ipynb`.
The statistical validation added in this refinement is three standalone scripts, run in
order:
```bash
python3 mantel_test.py                  # plain Mantel test (distance vs. co-reporting)
python3 build_connectivity_matrix.py    # builds the air-travel connectivity confound matrix
python3 partial_mantel_test.py          # partial Mantel test, full network + by disease stratum
```
`disease_classification.py` is imported by `partial_mantel_test.py` and classifies each
disease by transmission mode; it can also be run directly to print the classification.
`build_connectivity_matrix.py` downloads no data itself — it expects `airports_raw.dat` and
`routes_raw.dat` (OpenFlights data) in the working directory. Each script writes its
results to a `results_*.json` file. Combined runtime: under two minutes.

## Data
- `Outbreaks.csv`: outbreak surveillance records, 1996-2022 (Torres Munguía et al.,
  *Scientific Data* 9, 683, 2022, https://doi.org/10.1038/s41597-022-01797-2)
- `country-coord.csv`: country latitude/longitude coordinates (metal3d, 2024,
  https://gist.github.com/metal3d/5b925077e66194551df949de64e910f6)
- `airports_raw.dat`, `routes_raw.dat`: airport and route records (OpenFlights, 2024,
  https://openflights.org/data.php)

## Results
- Aggregate network: 233 nodes, 26,203 edges; maximum k-core = 221
- Network connectivity surges sharply during the COVID-19 pandemic years (2020-2022) while
  betweenness centrality falls to near zero, indicating broadly distributed rather than
  bridge-dependent connectivity during that period
- Mantel test: geographic distance is significantly associated with co-reporting strength
  (Pearson r = -0.192, permutation p < 0.0001; Spearman rho = -0.145, p < 0.0001), and this
  association strengthens rather than weakens when the pandemic years are included
  (pre-2020 only: r = -0.107, p < 0.0001)
- Partial Mantel test: the distance effect survives controlling for air-travel connectivity
  (partial r = -0.176, p < 0.001 for the full network), and is weakest for vector-borne
  disease (partial r = -0.070) and strongest for respiratory disease (partial r = -0.144)
  when stratified by transmission mode; see `paper.tex` for the proposed climate-zone
  explanation for this pattern, presented from established ecological theory and not
  independently verified within this dataset
- See `paper.tex` / `paper.pdf` for the full write-up, all figures, and discussion

## Availability
Code and figures for this paper are publicly available on GitHub: [TBD - URL to be
added].

## Paper
`paper.tex` is the canonical source, formatted for submission to *Applied Network Science*
(Springer), a single-column journal format chosen as the best fit for a figure-heavy
network science paper. Numbered bracket citations throughout. Compile with:
```bash
pdflatex paper.tex
pdflatex paper.tex   # second pass resolves cross-references
```
`paper.pdf` is the compiled output (20 pages).

## License
MIT License (code) / CC BY 4.0 (paper and figures)
