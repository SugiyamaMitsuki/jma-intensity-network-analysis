# 2018年大阪府北部地震に対する観測網密度の反実仮想分析

- 対象地震: `i2018_000842` (2018-06-18 07:58:34.14, M6.1, depth 13.0 km)
- 比較範囲: 震央から 100 km以内
- 反実仮想: 1995年末時点で稼働していた観測点位置に，2018年実観測の最近傍値を 10 km以内で割り当てた。

この分析は，1990年代相当の観測点密度・幾何が同じ地震の推計震度分布をどの程度変えるかを評価するための再サンプリングである。
したがって，当時の観測点固有の設置条件・機器特性・地盤応答を完全再現するものではない。

## 観測点支持

| scenario | n_station_used | density_per_10000km2 | station_count_within_10km | station_count_within_20km | station_count_within_50km | station_count_within_100km | nearest_epicentral_distance_km | station_nn_median_km | max_observed_intensity |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| current_2018_observed | 465 | 148.014 | 11 | 54 | 206 | 465 | 0.778 | 3.560 | 5.600 |
| counterfactual_1995_active | 59 | 18.780 | 2 | 5 | 27 | 59 | 0.778 | 9.320 | 5.300 |
| strict_1995_code_overlap | 21 | 6.685 | 1 | 1 | 7 | 21 | 0.778 | 19.939 | 5.400 |

## 推計震度分布の比較

![2018 current network versus 1995 counterfactual](../../paper/assets_ja/fig11_osaka_2018_counterfactual_maps.png)

震度は線形振幅ではないため，主比較には震度差ΔIを用いる。
右端の等価PGV比は `I = 2.68 + 1.72 log10(PGV)` に基づく参考換算であり，震度値そのものの比ではない。

![Counterfactual summary](../../paper/assets_ja/fig12_osaka_2018_counterfactual_summary.png)

## 2018年実観測点での検証

| scenario | validation_subset | n_validation | class_accuracy | class_within_1 | mae | rmse | bias | class_under_rate | class_over_rate |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| current_2018_observed | all_I>=1 | 465 | 0.877 | 1.000 | 0.103 | 0.150 | 0.002 | 0.082 | 0.041 |
| current_2018_observed | I>=5.0 | 28 | 0.643 | 1.000 | 0.132 | 0.181 | -0.061 | 0.286 | 0.071 |
| current_2018_observed | I>=5.5 | 5 | 0.400 | 1.000 | 0.152 | 0.209 | -0.152 | 0.600 | 0.000 |
| counterfactual_1995_active | all_I>=1 | 465 | 0.572 | 0.978 | 0.364 | 0.463 | -0.027 | 0.260 | 0.168 |
| counterfactual_1995_active | I>=5.0 | 28 | 0.357 | 0.786 | 0.350 | 0.434 | -0.316 | 0.643 | 0.000 |
| counterfactual_1995_active | I>=5.5 | 5 | 0.000 | 0.400 | 0.525 | 0.555 | -0.525 | 1.000 | 0.000 |
| strict_1995_code_overlap | all_I>=1 | 465 | 0.551 | 0.978 | 0.417 | 0.515 | -0.090 | 0.295 | 0.155 |
| strict_1995_code_overlap | I>=5.0 | 28 | 0.143 | 0.786 | 0.441 | 0.518 | -0.428 | 0.857 | 0.000 |
| strict_1995_code_overlap | I>=5.5 | 5 | 0.000 | 0.400 | 0.639 | 0.671 | -0.639 | 1.000 | 0.000 |

## 面積指標

| scenario | max_estimated_surface_intensity | area_ge_5p0_km2 | area_ge_5p5_km2 | area_ge_6p0_km2 | area_ge_6p5_km2 |
| --- | --- | --- | --- | --- | --- |
| current_2018_observed | 5.821 | 464.775 | 60.626 | 0.000 | 0.000 |
| counterfactual_1995_active | 5.505 | 283.026 | 4.041 | 0.000 | 0.000 |
| strict_1995_code_overlap | 5.547 | 133.384 | 8.083 | 0.000 | 0.000 |

## 原稿に入れるべき解釈

2018年の実観測網では，震央100 km以内に多数の観測点があり，震源近傍から大阪平野・京都盆地・奈良盆地にかけての局所的な震度勾配を直接拘束できる。
一方，1995年末の配置へ落とすと，近傍観測点数は大幅に減少し，推計震度分布は距離減衰式事前場と少数の残差観測に強く依存する。
この差は最大震度域そのものだけでなく，震度5弱以上・5強以上の面積評価にも現れる。
したがって，近年の地震被害域の面的把握を過去観測網へ外挿する場合，観測点密度の増加前後で同じ推計震度分布品質を仮定してはいけない。
