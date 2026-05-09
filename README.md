# JMA Intensity Network Analysis

This repository contains a manuscript, figures, analysis scripts, and lightweight derived tables for a study of how the evolution of the Japan Meteorological Agency (JMA) seismic-intensity observation network affects long-term intensity statistics and estimated intensity distributions.

## Read the Paper

- [English manuscript](paper/jma_intensity_network_interpolation_manuscript.md)
- [日本語版原稿](paper/jma_intensity_network_interpolation_manuscript_ja.md)

## Manuscript

- English manuscript: [`paper/jma_intensity_network_interpolation_manuscript.md`](paper/jma_intensity_network_interpolation_manuscript.md)
- Japanese manuscript: [`paper/jma_intensity_network_interpolation_manuscript_ja.md`](paper/jma_intensity_network_interpolation_manuscript_ja.md)
- Manuscript figures: [`paper/assets_en/`](paper/assets_en/) and [`paper/assets_ja/`](paper/assets_ja/)

## Main Findings

- The post-1995 expansion of the seismic-intensity network shortened nearest reporting-station distances; this is consistent with higher reported maximum intensity for comparable event groups.
- The within-event mean over reporting stations can decrease after network expansion because many more distant, low-intensity observations enter the event average.
- Treating zero retained stations as a simplified attenuation baseline gives an exact JMA intensity-class hit rate of 0.467 for all withheld reporting stations in the 63 target events with maximum intensity 6 lower or larger.
- Observation-constrained interpolation reaches about 70% exact class agreement only for all withheld reporting stations at roughly 57.4 retained stations per 10,000 km2.
- That 57.4 stations per 10,000 km2 benchmark is not a sufficient condition for strong-motion reconstruction: in the highest-density bin, exact class accuracy is only 0.45-0.46 for I>=5 lower and 0.37-0.39 for I>=6 lower.
- The 80% and 90% exact-class criteria are not reached within the tested density range, implying that density increase alone is insufficient under the present workflow.
- A counterfactual analysis of the 2018 northern Osaka earthquake shows that resampling the modern network to the 1994 active-site geometry reduces the 100-km station count from 465 to 34 and lowers exact class accuracy from 0.877 to 0.544.

## Repository Layout

```text
paper/                     Manuscripts and manuscript figures
src/                       Python analysis and plotting scripts
data/derived/              Lightweight derived CSV tables used in the manuscript
docs/reports/              Internal analysis reports retained as supplemental notes
```

Large raw and intermediate data are intentionally not committed. See [`data/README.md`](data/README.md).
Reproducibility commands, analysis options, random seeds, and environment notes are summarized in [`docs/reproducibility.md`](docs/reproducibility.md).

## Data Sources

The analysis uses the following external/public sources:

- JMA seismic-intensity annual files and station metadata.
- JMA unified hypocenter data from the JMA hypocenter data portal.
- J-SHIS surface-ground information, including AVS30 and site amplification factors.
- A prefectural boundary GMT file used only for map rendering.

The repository includes only lightweight derived summaries. Raw JMA/J-SHIS data and large generated grids should be downloaded or regenerated separately.

## Related Study

This work builds on the observation-network interpretation problem discussed by:

- Sugiyama, M., Yoshioka, Y., Hirai, T., and Fukuwa, N. (2020). 震度観測体制の年代差・地域差の定量評価と震度情報の解釈. Journal of Japan Association for Earthquake Engineering, 20(7), 7_101-7_119. https://doi.org/10.5610/jaee.20.7_101

## Reproduction Notes

The scripts were developed in a local conda environment named `seismo`.

```bash
conda env create -f environment.yml
conda activate seismo
```

The main workflow scripts are:

- `src/analyze_jma_intensity.py`: station-network and event-level intensity statistics.
- `src/prepare_hypocenter_intensity_catalog.py`: merges intensity events with JMA hypocenter data.
- `src/plot_jshis_surface_ground_pygmt.py`: maps J-SHIS AVS30 and site amplification.
- `src/estimate_jma_intensity_distribution.py`: generates JMA-style estimated intensity fields.
- `src/analyze_station_thinning_interpolation.py`: station-thinning cross validation.
- `src/resolve_final_analysis_concerns.py`: final sensitivity and uncertainty checks.
- `src/analyze_osaka_2018_network_counterfactual.py`: 2018 northern Osaka network-density counterfactual.

The exact local paths used during development are not required by the repository design; users should place raw data under a local `data/` directory and write regenerated outputs under `outputs/`.
See [`docs/reproducibility.md`](docs/reproducibility.md) for the command order and key options used for the manuscript revision.

## Authorship and AI Assistance

The current manuscript lists Mitsuki Sugiyama and Codex as authors. Codex/OpenAI was used for coding, analysis support, figure/report preparation, and manuscript drafting.

Recommended repository citation is provided in [`CITATION.cff`](CITATION.cff).

## License

No reuse license has been selected yet. Until a license is added, all rights are reserved by default. Before public release, choose an explicit license for code and manuscript materials, for example MIT/BSD for code and CC BY 4.0 for text and figures if open reuse is intended.
