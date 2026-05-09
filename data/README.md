# Data Notes

This repository includes only lightweight derived CSV tables needed to inspect the manuscript results.

## Included

- `data/derived/yearly_magbin_summary_land_approx_depth20.csv`
- `data/derived/yearly_shortest_epicentral_distance_m4plus_land_approx_depth20.csv`
- `data/derived/station_network_density_selected_years.csv`
- `data/derived/station_network_density_period_summary.csv`
- `data/derived/station_network_year_end_comparison.csv`
- `data/derived/station_network_transition_analysis.csv`
- `data/derived/period_magbin_summary_land_approx_depth20.csv`
- `data/derived/jshis_ground_variability_metrics_0p02deg.csv.gz`
- `data/derived/observed_predicted_intensity_scatter_summary.csv`
- `data/derived/gmpe_vs_spatial_interpolation_effective0/`
- `data/derived/station_thinning_interpolation_6lower_plus_class/`
- `data/derived/osaka_2018_network_counterfactual/`
- `data/raw_data_manifest_sha256.csv`

These files are small derived summaries. They are not a substitute for the raw JMA intensity catalog, JMA hypocenter catalog, or J-SHIS surface-ground grid.

## Excluded

The following are not committed:

- Raw JMA annual seismic-intensity files.
- JMA `code_p.dat` station-history metadata.
- JMA unified hypocenter annual files.
- Raw or gridded J-SHIS surface-ground files.
- Full estimated intensity grids for all target events.
- Full station-thinning trial prediction tables.
- Monthly station-map animation frames and videos.

Reasons:

- Some source datasets have external terms of use and should be obtained from the original providers.
- Several generated products are hundreds of MB to multiple GB.
- Lightweight summaries are sufficient for manuscript inspection.

## Expected Local Reproduction Layout

For local reproduction, use a workspace layout like:

```text
data/
  i1980.zip ... i2022.zip
  code_p.zip or code_p.dat
  hypocenter/
  j-shis/
outputs/
```

Scripts write derived outputs to `outputs/csv/`, `outputs/png/`, and `outputs/manuscript/`.

## Reproducibility Metadata

The manuscript revision uses random station thinning with seed `20260509`, a minimum station intensity of `1.0`, event-region padding of `0.55` degrees, and smoothing kriging with nugget `0.02` unless otherwise noted. See [`../docs/reproducibility.md`](../docs/reproducibility.md) for the command order, key options, raw-data source URLs, and checksum workflow.

`data/raw_data_manifest_sha256.csv` records local raw-file names, byte counts, modification timestamps, source URLs, and SHA-256 checksums for the local archive used in this manuscript revision. The original download timestamps for older local files were not logged before this revision, so the manifest should be read as a fixed local-file inventory rather than a provider-side acquisition log.
