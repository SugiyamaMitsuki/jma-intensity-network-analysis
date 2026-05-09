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

- The post-1995 expansion of the seismic-intensity network increases the probability of observing strong shaking near the epicenter, raising reported maximum intensity for comparable events.
- Event mean station intensity can decrease after network expansion because many more distant, low-intensity observations enter the event average.
- Treating zero retained stations as attenuation-only prediction gives an exact JMA intensity-class hit rate of 0.455 for the 22 target events with maximum intensity 6 upper or 7.
- Observation-constrained interpolation reaches about 70% exact class agreement only at roughly 57.7 retained stations per 10,000 km2.
- The 80% and 90% exact-class criteria are not reached within the tested density range, implying that density increase alone is insufficient under the present workflow.
- A counterfactual analysis of the 2018 northern Osaka earthquake shows that resampling the modern network to the 1995 active-site geometry reduces the 100-km station count from 465 to 59 and lowers exact class accuracy from 0.877 to 0.572.

## Repository Layout

```text
paper/                     Manuscripts and manuscript figures
src/                       Python analysis and plotting scripts
data/derived/              Lightweight derived CSV tables used in the manuscript
docs/reports/              Internal analysis reports retained as supplemental notes
```

Large raw and intermediate data are intentionally not committed. See [`data/README.md`](data/README.md).

## Data Sources

The analysis uses the following external/public sources:

- JMA seismic-intensity annual files and station metadata.
- JMA unified hypocenter data from the JMA hypocenter data portal.
- J-SHIS surface-ground information, including AVS30 and site amplification factors.
- A prefectural boundary GMT file used only for map rendering.

The repository includes only lightweight derived summaries. Raw JMA/J-SHIS data and large generated grids should be downloaded or regenerated separately.

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

## Authorship and AI Assistance

Codex/OpenAI was used as an AI-assisted coding, analysis, and drafting tool. Codex should not be listed as an author because it cannot take responsibility for the work, approve the final manuscript, or be accountable for data integrity and interpretation. The responsible author should be the human researcher who validates and submits the work.

Recommended repository citation is provided in [`CITATION.cff`](CITATION.cff). Update affiliations, ORCID, and corresponding-author information before journal submission.

## License

No reuse license has been selected yet. Until a license is added, all rights are reserved by default. Before public release, choose an explicit license for code and manuscript materials, for example MIT/BSD for code and CC BY 4.0 for text and figures if open reuse is intended.
