# 観測点0を距離減衰式のみとした推計震度補間比較

## 要点

- 観測点0状態は `gmpe_raw` と定義し、実効観測点数 `effective_n_train=0`、実効観測点密度 `effective_train_density_per_10000km2=0` に固定した。
- `gmpe_raw` は距離減衰式からPGVを予測し震度換算するだけで、観測点によるイベント校正や残差補間を行わない。
- 主評価指標は震度階級の完全一致率。距離減衰式のみは約45.5%、観測点で校正したGMPEは約58.8%、観測点残差を使うIDW/kriging/GMPE+krigingは高密度側で約70.3%まで改善した。
- 70%の階級一致率を満たすには、今回の設定では約57.7点/10,000 km2、最近接観測点距離中央値で約4.35 kmが必要だった。
- 80%・90%の階級一致率は、試験した密度範囲では到達しない。個別試行最大密度は71.58点/10,000 km2だが、保持率0.9でも主要補間法の中央値は約69-70%である。

## 手法の定義

| method | 観測点の使い方 | 解釈 |
| --- | --- | --- |
| `gmpe_raw` | 0点 | 距離減衰式のみ。観測点0状態の基準。 |
| `gmpe_calibrated` | 保持観測点 | 距離減衰式をイベントごとに線形校正。空間残差は補間しない。 |
| `idw` | 保持観測点 | 地盤増幅をキャンセルした基盤面震度をIDW補間し、評価点の増幅で地表に戻す。 |
| `kriging` | 保持観測点 | 平滑化ordinary krigingで基盤面震度を補間する。 |
| `gmpe_kriging` | 保持観測点 | GMPE参照場を校正し、観測点残差をkrigingする。 |

## 全体精度

| method | RMSE | MAE | within ±0.5 | exact class | within 1 class | over | under | class MAE |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| gmpe_raw | 0.728 | 0.610 | 0.450 | 0.455 | 0.942 | 0.477 | 0.062 | 0.601 |
| gmpe_calibrated | 0.513 | 0.406 | 0.683 | 0.588 | 0.976 | 0.176 | 0.238 | 0.442 |
| idw | 0.384 | 0.292 | 0.826 | 0.680 | 0.991 | 0.135 | 0.180 | 0.329 |
| kriging | 0.390 | 0.297 | 0.819 | 0.675 | 0.991 | 0.137 | 0.184 | 0.335 |
| gmpe_kriging | 0.391 | 0.297 | 0.819 | 0.674 | 0.991 | 0.138 | 0.185 | 0.336 |

## 実効観測点密度との関係

![震度階級一致率](/Users/mitsuki/02_作業/001_地震/2605_震度研究/outputs/png/gmpe_vs_spatial_interpolation_effective0/thinning_class_accuracy_vs_station_density.png)

| method | density bin | trials | density | n train | RMSE | exact | within 1 | over | under | NN km |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| gmpe_raw | [0.0, 5.0) | 660 | 0.00 | 0 | 0.728 | 0.455 | 0.942 | 0.477 | 0.062 |  |
| gmpe_calibrated | [0.0, 5.0) | 45 | 2.46 | 62 | 0.583 | 0.556 | 0.940 | 0.214 | 0.223 | 37.62 |
| gmpe_calibrated | [5.0, 10.0) | 110 | 6.64 | 146 | 0.517 | 0.585 | 0.976 | 0.180 | 0.237 | 12.99 |
| gmpe_calibrated | [10.0, 20.0) | 160 | 13.91 | 311 | 0.514 | 0.589 | 0.977 | 0.170 | 0.242 | 8.59 |
| gmpe_calibrated | [20.0, 30.0) | 85 | 23.32 | 539 | 0.506 | 0.591 | 0.974 | 0.176 | 0.241 | 6.89 |
| gmpe_calibrated | [30.0, 50.0) | 170 | 39.68 | 884 | 0.510 | 0.587 | 0.976 | 0.177 | 0.236 | 5.53 |
| gmpe_calibrated | [50.0, 75.0) | 90 | 57.68 | 1314 | 0.506 | 0.591 | 0.977 | 0.173 | 0.234 | 4.35 |
| idw | [0.0, 5.0) | 45 | 2.46 | 62 | 0.561 | 0.608 | 0.938 | 0.176 | 0.203 | 37.62 |
| idw | [5.0, 10.0) | 110 | 6.64 | 146 | 0.411 | 0.656 | 0.988 | 0.150 | 0.190 | 12.99 |
| idw | [10.0, 20.0) | 160 | 13.91 | 311 | 0.387 | 0.676 | 0.991 | 0.135 | 0.188 | 8.59 |
| idw | [50.0, 75.0) | 90 | 57.68 | 1314 | 0.350 | 0.703 | 0.994 | 0.122 | 0.174 | 4.35 |
| kriging | [0.0, 5.0) | 45 | 2.46 | 62 | 0.588 | 0.602 | 0.936 | 0.170 | 0.214 | 37.62 |
| kriging | [5.0, 10.0) | 110 | 6.64 | 146 | 0.416 | 0.654 | 0.987 | 0.150 | 0.193 | 12.99 |
| kriging | [10.0, 20.0) | 160 | 13.91 | 311 | 0.396 | 0.669 | 0.991 | 0.139 | 0.194 | 8.59 |
| kriging | [50.0, 75.0) | 90 | 57.68 | 1314 | 0.354 | 0.704 | 0.993 | 0.122 | 0.170 | 4.35 |
| gmpe_kriging | [0.0, 5.0) | 45 | 2.46 | 62 | 0.591 | 0.608 | 0.932 | 0.173 | 0.219 | 37.62 |
| gmpe_kriging | [5.0, 10.0) | 110 | 6.64 | 146 | 0.417 | 0.654 | 0.987 | 0.149 | 0.190 | 12.99 |
| gmpe_kriging | [10.0, 20.0) | 160 | 13.91 | 311 | 0.396 | 0.668 | 0.990 | 0.139 | 0.196 | 8.59 |
| gmpe_kriging | [50.0, 75.0) | 90 | 57.68 | 1314 | 0.354 | 0.703 | 0.993 | 0.122 | 0.171 | 4.35 |

## 70%・80%・90%階級一致率の到達密度

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

80%・90%基準は密度ビン中央値では未到達である。今回の設計で試験した個別試行の最大有効密度は71.58点/10,000 km2であり、これを超える密度が少なくとも必要である。ただし高密度側の完全一致率は約70%で飽和しており、密度だけを上げても80%以上に届くとは解釈しない。

## 解釈

観測点0の `gmpe_raw` は、震源距離に伴う平均的な減衰トレンドを与えるが、震度階級の完全一致率は45%程度にとどまる。過大評価率が約48%と高く、今回の震度換算・震源距離近似・断層距離近似だけでは、地震ごとの放射特性や局所的な地盤・伝播差を吸収できない。

保持観測点で距離減衰式を校正する `gmpe_calibrated` は、イベント全体の振幅バイアスを補正するため一致率が約59%に上がる。しかし密度を増やしてもほぼ頭打ちで、空間残差を使わない限り70%基準には届かない。

IDW、平滑化kriging、GMPE残差krigingはいずれも観測点密度とともに改善し、約57.7点/10,000 km2で70%基準に達する。ただし80%・90%基準は未到達であり、現行手法の誤差床を下げる改善が必要である。今回のクロスバリデーションでは `gmpe_kriging` が通常krigingを明確に上回るところまでは出ていない。これは震源モデルを点震源近似、断層距離を簡略化していることも影響している可能性がある。断層面距離・Mw・断層タイプをより厳密に入れると、疎な外挿領域でGMPE参照場の利点が出やすい。

## 出力

- 試行集計CSV: `/Users/mitsuki/02_作業/001_地震/2605_震度研究/outputs/csv/gmpe_vs_spatial_interpolation_effective0/thinning_trial_summary.csv`
- 密度別集計CSV: `/Users/mitsuki/02_作業/001_地震/2605_震度研究/outputs/csv/gmpe_vs_spatial_interpolation_effective0/thinning_density_bin_summary.csv`
- 予測誤差CSV: `/Users/mitsuki/02_作業/001_地震/2605_震度研究/outputs/csv/gmpe_vs_spatial_interpolation_effective0/thinning_prediction_errors.csv.gz`
- 手法別要約CSV: `/Users/mitsuki/02_作業/001_地震/2605_震度研究/outputs/csv/gmpe_vs_spatial_interpolation_effective0/gmpe_vs_spatial_method_summary.csv`
- 図PNG: `/Users/mitsuki/02_作業/001_地震/2605_震度研究/outputs/png/gmpe_vs_spatial_interpolation_effective0`
