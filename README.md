# RF-GLS_article_code


Please download the folder https://github.com/ArkajyotiSaha/RF-GLS_article_code/tree/main/RF-GLS_article_code. TIn order to run the .R codes please set the working directory of R to this folder.

## Function description:

* cost_function.R: Cost function, cutoff selection, node retresentative in partitoin of single node with CART and DART criteria.
* m1_250.R, m2_250.R, m3_250.R: Estimation & Prediction with correctly specified model, with n = 250, parameters unknown for m = m_1, m_2 and m_3 respectively.
* m1_250_known.R, m2_250_known.R: Estimation & Prediction with correctly specified model, with n = 250, parameters known for m = m_1, and m_2 respectively.
* m1_1000.R, m2_1000.R: Estimation & Prediction with correctly specified model with NNGP approximation, with n = 1000, parameters unknown for m = m_1, and m_2 respectively.
* m1_matern_misspec.R, m2_matern_misspec.R, m3_matern_misspec.R: Estimation & Prediction with misspecified model, errors generated from matern covariance function with n = 250, for m = m_1, m_2 and m_3 respectively.
* m1_smooth_misspec.R, m2_smooth_misspec.R, m3_smooth_misspec.R: Estimation & Prediction with misspecified model, errors generated from smooth function with n = 250, for m = m_1, m_2 and m_3 respectively.

## Runtime:

* cost_function.R: 5 mins.

T = Average time required for running the following functions in cluster (24-core Intel(R) Xeon(R) CPU E5-2650 v4 processor with 2.20GHz) for one seed. To produce the results in the manuscript they were run for 100 seeds for each setup.

* m = m_1 & m_2 (M_{try} = 1), n = 250; T = 25 mins
* m = m_1 & m_2 (M_{try} = 1), n = 1000; T = 2 hrs.
* m = m_3, n = 250; T = 2 hrs.

## Plot vs. Function

The following functions were run to produce the figures in the manuscript and the supplementary material.

* Figure 1 : cost_function.R
* Figure 2 : m2_250.R
* Figure 3 : m2_250.R
* Figure 4 : m2_250.R
* Figure 5 : m2_smooth_misspec.R
* Figure 6a: m1_250_known.R
* Figure 6b: m2_250_known.R
* Figure 7 : m1_250.R
* Figure 8 : m1_1000.R
* Figure 9 : m2_1000.R
* Figure 10: m3_250.R
* Figure 11: m1_matern_misspec.R
* Figure 12: m2_matern_misspec.R
* Figure 13: m3_matern_misspec.R
* Figure 14: m1_smooth_misspec.R
* Figure 15: m1_smooth_misspec.R
* Figure 16: m2_smooth_misspec.R
* Figure 17: m3_smooth_misspec.R

The timings are approximate and can vary significantly based on the machine configurations. There might be some issues with compiling the .Rcpp code in macOS with clang. Please look for solutions to bypass the problem as these can be very machine/OS version specific issues. 
