# Reproducibility Notes

This note records the analysis settings behind the manuscript revision. Large raw datasets and full intermediate grids are not committed to this repository; the committed files are lightweight derived summaries, figures, scripts, and manuscripts.

## Version Control

- Repository: `https://github.com/SugiyamaMitsuki/jma-intensity-network-analysis`
- Branch used for the manuscript revision: `main`
- Parent commit before the review-response revision: `0093442`
- Exact revision after cloning: run `git rev-parse HEAD`

## Environment

The analysis was run in a conda environment named `seismo`.

```bash
conda env create -f environment.yml
conda activate seismo
```

Key packages observed in the local environment were:

| Package | Version |
| --- | --- |
| Python | 3.11 |
| numpy | 1.26.4 |
| pandas | 2.2.3 |
| scipy | 1.17.1 |
| matplotlib | 3.10.8 |
| xarray | 2026.2.0 |
| pygmt | 0.17.0 |
| shapely | 2.1.2 |
| cartopy | 0.25.0 |

## Raw Data Sources

Expected local files are described in [data/README.md](../data/README.md). Raw data should be obtained from the original providers and stored under a local `data/` directory.

| Dataset | Source | Notes |
| --- | --- | --- |
| JMA seismic-intensity annual files | JMA intensity archive | Annual files for 1980-2022 were parsed into event-level and station-level records. |
| JMA station-history metadata | Local `code_p.dat` archive | Used to reconstruct active station counts and historical station geometry. |
| JMA unified hypocenter catalog | https://www.data.jma.go.jp/eqev/data/bulletin/hypo.html | Used to replace coarse hypocenter information embedded in intensity records where possible. |
| J-SHIS surface-ground data | https://www.j-shis.bosai.go.jp/ | AVS30 and engineering-bedrock-to-surface PGV amplification factors. |
| Prefectural boundary GMT file | User-provided local boundary file | Used only for map rendering. |

Record raw-data checksums after download:

```bash
find data -type f \( -name '*.zip' -o -name '*.dat' -o -name '*.csv' -o -name '*.txt' \) \
  -print | sort | xargs shasum -a 256 > data/raw_checksums_sha256.txt
```

## Main Analysis Order

Run commands from the repository root. Paths can be changed if the raw data and output directories are elsewhere.

```bash
conda activate seismo

python src/analyze_jma_intensity.py

python src/prepare_hypocenter_intensity_catalog.py \
  --data-dir data \
  --event-summary outputs/csv/jma_intensity_event_summary.csv \
  --out-dir outputs/csv

python src/plot_jshis_surface_ground_pygmt.py \
  --data-dir data \
  --csv-dir outputs/csv \
  --png-dir outputs/png

python src/estimate_jma_intensity_distribution.py \
  --methods gmpe_raw,idw,kriging,gmpe_kriging \
  --min-station-intensity 1.0 \
  --region-padding-deg 0.55 \
  --kriging-nugget 0.02 \
  --gmpe-fault-type crustal

python src/analyze_station_thinning_interpolation.py \
  --event-id all-targets \
  --methods gmpe_raw,gmpe_calibrated,idw,kriging,gmpe_kriging \
  --keep-fractions 0.1,0.2,0.3,0.5,0.7,0.9 \
  --n-random 5 \
  --seed 20260509 \
  --min-station-intensity 1.0 \
  --region-padding-deg 0.55 \
  --kriging-nugget 0.02 \
  --gmpe-fault-type crustal

python src/analyze_intensity_dependent_predictability.py \
  --predictions outputs/csv/station_thinning_predictions_effective0.csv

python src/analyze_probabilistic_uncertainty.py \
  --predictions outputs/csv/station_thinning_predictions_effective0.csv \
  --amp-intensity-coef 1.72 \
  --sigma-amp-intensity-coef 0.20 \
  --sigma-log10-amplification 0.08 \
  --min-sigma 0.05

python src/analyze_osaka_2018_network_counterfactual.py \
  --event-id i2018_000842 \
  --counterfactual-year 1994 \
  --radius-km 100 \
  --map-radius-km 115 \
  --match-radius-km 10 \
  --methods gmpe_raw,idw,kriging,gmpe_kriging \
  --map-method gmpe_kriging \
  --min-station-intensity 1.0
```

## Key Analysis Settings

| Setting | Value |
| --- | --- |
| Target events for interpolation | Events with maximum JMA intensity 6 upper or 7 in the parsed catalog |
| Station-thinning seed | `20260509` |
| Retained fractions | `0.1,0.2,0.3,0.5,0.7,0.9` |
| Random draws per fraction | `5` |
| Minimum reported station intensity | `1.0` |
| Event-region padding | `0.55` degrees |
| Site amplification intensity coefficient | `1.72` |
| IDW neighbors / power | `16` / `2.0` |
| Kriging nugget | `0.02` |
| Default attenuation fault type | `crustal` |
| PGV-to-intensity intercept | `2.68` |
| Intensity clipping range | `0.0` to `7.2` |
| Osaka counterfactual year | `1994` |
| Osaka station assignment radius | `10 km` |

## Interpretation Constraints

The 57.7 stations per 10,000 km2 benchmark is defined as an effective retained-station density in event-specific evaluation regions. It is not a nationwide operating density. It is also an all-withheld-reporting-station point-validation result; it is not a sufficient condition for reproducing strong-motion footprints.

The `gmpe_*` method names are retained for code compatibility. In the manuscript, `gmpe_raw` is interpreted as a simplified Si and Midorikawa-type attenuation baseline because finite-fault distance, event-specific rupture geometry, and the published aleatory variability of a full GMPE are not yet represented.
