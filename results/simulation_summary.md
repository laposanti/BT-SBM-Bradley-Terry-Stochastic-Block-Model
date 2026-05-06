# Simulation summary

Source CSV: `simulation_comparison_all_models.csv`

## Fit status

|model  | OK| SKIPPED| n_total|   ok_rate|
|:------|--:|-------:|-------:|---------:|
|BT     |  9|       0|       9| 1.0000000|
|BT-SBM |  9|       0|       9| 1.0000000|
|RCBTL  |  1|       8|       9| 0.1111111|

## Performance (OK fits only)

|model  | n_ok|  rmse_mean| rmse_median|  mae_mean| rmse_rel_mean| mae_rel_mean| mean_abs_log_ratio_mean| rms_log_ratio_mean| spearman_mean| kendall_mean| pearson_mean|
|:------|----:|----------:|-----------:|---------:|-------------:|------------:|-----------------------:|------------------:|-------------:|------------:|------------:|
|RCBTL  |    1|   1.057674|    1.057674|  0.392741|            NA|           NA|                      NA|                 NA|     0.8414355|    0.7184239|    0.9471209|
|BT     |    9|  97.705455|   51.058509| 20.088918|            NA|           NA|                      NA|                 NA|     0.9015512|    0.7931125|    0.8917631|
|BT-SBM |    9| 119.500619|   16.153057| 24.344359|            NA|           NA|                      NA|                 NA|     0.8935554|    0.7853884|    0.8155540|

## Timing

|model  |  n| elapsed_mean_s| elapsed_median_s| elapsed_p90_s| elapsed_mean_min|
|:------|--:|--------------:|----------------:|-------------:|----------------:|
|BT     |  9|       22.30944|           22.346|       22.5854|        0.3718241|
|BT-SBM |  9|       23.91744|           23.959|       24.1436|        0.3986241|
|RCBTL  |  1|    40003.45200|        40003.452|    40003.4520|      666.7242000|

## Clustering recovery (where available)

|model  | n_ok| K_median_mean| pr_K_true_mean| ari_minVI_mean| vi_minVI_mean| ari_binder_mean| vi_binder_mean|
|:------|----:|-------------:|--------------:|--------------:|-------------:|---------------:|--------------:|
|BT-SBM |    9|      4.555556|      0.4125556|      0.7673566|     0.7287145|       0.7915687|      0.7650334|
|RCBTL  |    1|     18.000000|      0.0000000|      0.1776790|     2.7359335|       0.0893360|      3.4954472|

