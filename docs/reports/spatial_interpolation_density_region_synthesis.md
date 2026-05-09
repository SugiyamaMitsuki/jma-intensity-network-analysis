# 観測点密度・距離減衰式比較・地域差の整理

## 結論

- 観測点0状態は距離減衰式のみ `gmpe_raw` として扱う。この基準の震度階級完全一致率は 45.5%、RMSEは 0.728。
- 空間補間が「十分」と言える水準を 70%の震度階級完全一致率 と置くと、必要密度は約 57.7点/10,000 km2、最近接観測点距離中央値は約 4.35 km。
- 65%程度を許容するなら約 6.6点/10,000 km2、最近接距離中央値約 13.0 kmで到達する。68%程度なら約 23.3点/10,000 km2、最近接距離中央値約 6.9 kmが目安。
- 70%水準では、距離減衰式のみと比べて階級一致率は約 +25ポイント、RMSEは約 52%低下する。
- 80%および90%の震度階級完全一致率は、試験した密度範囲では到達しない。個別試行の最大有効密度は 71.58点/10,000 km2 だが、保持率0.9でも主要補間法の中央値は約69-70%で頭打ちとなる。
- 地域差は残る。高密度側では関東、近畿・中国・四国、中部・北陸が相対的に高く、東北、北海道、九州は70%前後またはそれ未満にとどまる。

## 密度しきい値

| method | criterion | density | NN km | exact | status |
| --- | --- | --- | --- | --- | --- |
| gmpe_raw | exact>=0.65 |  |  |  | not_reached |
| gmpe_raw | exact>=0.68 |  |  |  | not_reached |
| gmpe_raw | exact>=0.70 |  |  |  | not_reached |
| gmpe_raw | exact>=0.80 |  |  |  | not_reached |
| gmpe_raw | exact>=0.90 |  |  |  | not_reached |
| gmpe_calibrated | exact>=0.65 |  |  |  | not_reached |
| gmpe_calibrated | exact>=0.68 |  |  |  | not_reached |
| gmpe_calibrated | exact>=0.70 |  |  |  | not_reached |
| gmpe_calibrated | exact>=0.80 |  |  |  | not_reached |
| gmpe_calibrated | exact>=0.90 |  |  |  | not_reached |
| idw | exact>=0.65 | 6.64 | 12.99 | 0.656 | reached |
| idw | exact>=0.68 | 23.32 | 6.89 | 0.690 | reached |
| idw | exact>=0.70 | 57.68 | 4.35 | 0.703 | reached |
| idw | exact>=0.80 |  |  |  | not_reached |
| idw | exact>=0.90 |  |  |  | not_reached |
| kriging | exact>=0.65 | 6.64 | 12.99 | 0.654 | reached |
| kriging | exact>=0.68 | 23.32 | 6.89 | 0.683 | reached |
| kriging | exact>=0.70 | 57.68 | 4.35 | 0.704 | reached |
| kriging | exact>=0.80 |  |  |  | not_reached |
| kriging | exact>=0.90 |  |  |  | not_reached |
| gmpe_kriging | exact>=0.65 | 6.64 | 12.99 | 0.654 | reached |
| gmpe_kriging | exact>=0.68 | 23.32 | 6.89 | 0.681 | reached |
| gmpe_kriging | exact>=0.70 | 57.68 | 4.35 | 0.703 | reached |
| gmpe_kriging | exact>=0.80 |  |  |  | not_reached |
| gmpe_kriging | exact>=0.90 |  |  |  | not_reached |

80%・90%基準は、密度ビン中央値では未到達である。したがって現在の検証データから「80%または90%に必要な密度」を経験的に推定することはできず、少なくとも試験済み最大有効密度 71.58点/10,000 km2 を超える。現行手法では高密度側で約70%に飽和しているため、80%以上を目標にする場合は、観測点密度だけでなく震源モデル、距離減衰参照場、地盤補正、不確実性重み付けの改善が必要である。

![density needed](/Users/mitsuki/02_作業/001_地震/2605_震度研究/outputs/png/gmpe_vs_spatial_interpolation_effective0/synthesis/density_needed_class_accuracy.png)

## 距離減衰式のみからの改善

| method | density bin | density | NN km | RMSE | exact | gain | RMSE red. | over | under |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| gmpe_raw | [0.0, 5.0) | 0.00 |  | 0.728 | 0.455 | 0.0 pt | 0.0% | 0.477 | 0.062 |
| idw | [0.0, 5.0) | 2.46 | 37.62 | 0.561 | 0.608 | 15.3 pt | 23.0% | 0.176 | 0.203 |
| idw | [5.0, 10.0) | 6.64 | 12.99 | 0.411 | 0.656 | 20.1 pt | 43.6% | 0.150 | 0.190 |
| idw | [10.0, 20.0) | 13.91 | 8.59 | 0.387 | 0.676 | 22.1 pt | 46.9% | 0.135 | 0.188 |
| idw | [20.0, 30.0) | 23.32 | 6.89 | 0.373 | 0.690 | 23.5 pt | 48.8% | 0.132 | 0.175 |
| idw | [50.0, 75.0) | 57.68 | 4.35 | 0.350 | 0.703 | 24.8 pt | 52.0% | 0.122 | 0.174 |
| kriging | [0.0, 5.0) | 2.46 | 37.62 | 0.588 | 0.602 | 14.7 pt | 19.3% | 0.170 | 0.214 |
| kriging | [5.0, 10.0) | 6.64 | 12.99 | 0.416 | 0.654 | 19.9 pt | 42.8% | 0.150 | 0.193 |
| kriging | [10.0, 20.0) | 13.91 | 8.59 | 0.396 | 0.669 | 21.3 pt | 45.7% | 0.139 | 0.194 |
| kriging | [20.0, 30.0) | 23.32 | 6.89 | 0.382 | 0.683 | 22.8 pt | 47.5% | 0.134 | 0.179 |
| kriging | [50.0, 75.0) | 57.68 | 4.35 | 0.354 | 0.704 | 24.8 pt | 51.3% | 0.122 | 0.170 |
| gmpe_kriging | [0.0, 5.0) | 2.46 | 37.62 | 0.591 | 0.608 | 15.3 pt | 18.8% | 0.173 | 0.219 |
| gmpe_kriging | [5.0, 10.0) | 6.64 | 12.99 | 0.417 | 0.654 | 19.9 pt | 42.8% | 0.149 | 0.190 |
| gmpe_kriging | [10.0, 20.0) | 13.91 | 8.59 | 0.396 | 0.668 | 21.3 pt | 45.6% | 0.139 | 0.196 |
| gmpe_kriging | [20.0, 30.0) | 23.32 | 6.89 | 0.384 | 0.681 | 22.5 pt | 47.3% | 0.133 | 0.181 |
| gmpe_kriging | [50.0, 75.0) | 57.68 | 4.35 | 0.354 | 0.703 | 24.8 pt | 51.4% | 0.122 | 0.171 |

![improvement](/Users/mitsuki/02_作業/001_地震/2605_震度研究/outputs/png/gmpe_vs_spatial_interpolation_effective0/synthesis/improvement_over_attenuation_by_density.png)

## 地域差

地域区分は震央/観測点の緯度経度による広域区分であり、行政界そのものではない。表は高密度側の空間補間 `[50,75)` と観測点0の `gmpe_raw` を比較した。

| station region | method | trial-regions | predictions | exact | RMSE | gain | RMSE red. |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Hokkaido | gmpe_raw | 282 | 27350 | 0.533 | 0.629 | 0.0 pt | 0.0% |
| Hokkaido | idw | 35 | 657 | 0.667 | 0.360 | 13.4 pt | 42.7% |
| Hokkaido | kriging | 35 | 657 | 0.696 | 0.372 | 16.3 pt | 40.8% |
| Hokkaido | gmpe_kriging | 35 | 657 | 0.694 | 0.373 | 16.2 pt | 40.6% |
| Tohoku | gmpe_raw | 496 | 101118 | 0.436 | 0.705 | 0.0 pt | 0.0% |
| Tohoku | idw | 75 | 3751 | 0.676 | 0.377 | 24.0 pt | 46.5% |
| Tohoku | kriging | 75 | 3751 | 0.691 | 0.382 | 25.5 pt | 45.8% |
| Tohoku | gmpe_kriging | 75 | 3751 | 0.686 | 0.383 | 25.0 pt | 45.6% |
| Kanto | gmpe_raw | 557 | 117003 | 0.446 | 0.633 | 0.0 pt | 0.0% |
| Kanto | idw | 80 | 4773 | 0.749 | 0.292 | 30.3 pt | 53.8% |
| Kanto | kriging | 80 | 4773 | 0.738 | 0.302 | 29.2 pt | 52.3% |
| Kanto | gmpe_kriging | 80 | 4773 | 0.738 | 0.302 | 29.2 pt | 52.3% |
| Chubu-Hokuriku | gmpe_raw | 560 | 137759 | 0.432 | 0.725 | 0.0 pt | 0.0% |
| Chubu-Hokuriku | idw | 80 | 5860 | 0.702 | 0.361 | 27.1 pt | 50.2% |
| Chubu-Hokuriku | kriging | 80 | 5860 | 0.707 | 0.363 | 27.6 pt | 50.0% |
| Chubu-Hokuriku | gmpe_kriging | 80 | 5860 | 0.707 | 0.364 | 27.6 pt | 49.8% |
| Kinki-Chugoku-Shikoku | gmpe_raw | 507 | 62352 | 0.464 | 0.704 | 0.0 pt | 0.0% |
| Kinki-Chugoku-Shikoku | idw | 77 | 1542 | 0.773 | 0.276 | 30.9 pt | 60.9% |
| Kinki-Chugoku-Shikoku | kriging | 77 | 1542 | 0.783 | 0.281 | 31.9 pt | 60.1% |
| Kinki-Chugoku-Shikoku | gmpe_kriging | 77 | 1542 | 0.778 | 0.283 | 31.4 pt | 59.8% |
| Kyushu | gmpe_raw | 208 | 44188 | 0.400 | 0.811 | 0.0 pt | 0.0% |
| Kyushu | idw | 19 | 932 | 0.641 | 0.404 | 24.1 pt | 50.3% |
| Kyushu | kriging | 19 | 932 | 0.645 | 0.412 | 24.5 pt | 49.2% |
| Kyushu | gmpe_kriging | 19 | 932 | 0.645 | 0.412 | 24.5 pt | 49.2% |

![region high density](/Users/mitsuki/02_作業/001_地震/2605_震度研究/outputs/png/gmpe_vs_spatial_interpolation_effective0/synthesis/high_density_accuracy_by_station_region.png)

![region density heatmap](/Users/mitsuki/02_作業/001_地震/2605_震度研究/outputs/png/gmpe_vs_spatial_interpolation_effective0/synthesis/kriging_accuracy_by_region_and_density.png)

## 読み取り

距離減衰式のみは、震源距離に沿った平均トレンドを与えるが、階級一致率は低く、過大評価が多い。観測点を数十点入れるだけでも空間補間は大きく改善し、最初の 0 -> 6.6点/10,000 km2 の増分で改善の大部分が得られる。一方、70%基準に乗せるには最近接距離4-5 km程度まで密にする必要がある。ただし80%・90%基準は試験密度範囲では未到達であり、密度増加だけでなくモデル側の誤差床を下げる必要がある。

地域差は、地盤条件のばらつきだけでなく、対象イベントが内陸直下か沖合か、強震域が観測網の内側にあるか外側にあるかにも依存する。したがって、全国一律の密度しきい値は第一近似であり、複雑地形・沿岸/島嶼・沖合イベントでは局所的により高密度な観測点、または断層面距離を使ったGMPE参照場の改善が必要になる。

## 出力

- CSV: `/Users/mitsuki/02_作業/001_地震/2605_震度研究/outputs/csv/gmpe_vs_spatial_interpolation_effective0/synthesis`
- PNG: `/Users/mitsuki/02_作業/001_地震/2605_震度研究/outputs/png/gmpe_vs_spatial_interpolation_effective0/synthesis`
