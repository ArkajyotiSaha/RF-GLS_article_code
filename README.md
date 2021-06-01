# RF-GLS_article_code


Please download the folder https://github.com/ArkajyotiSaha/RF-GLS_article_code/tree/main/RF-GLS_article_code. TIn order to run the .R codes please set the working directory of R to this folder.

We briefly describe the codes corresponding to each of the plots in the anuscript; in parenthesis, we mention the average time required for running the function in cluster (24-core Intel(R) Xeon(R) CPU E5-2650 v4 processor with 2.20GHz).

* Figure 1 : cost_function.R (<5 mins)
* Figure 2 : m2_250.R (~25 mins)
* Figure 3 : m2_250.R (~25 mins)
* Figure 4 : m2_250.R (~25 mins)
* Figure 5 : m2_smooth_misspec.R (~25 mins)
* Figure 6a: m1_250_known.R (~25 mins)
* Figure 6b: m2_250_known.R (~25 mins)
* Figure 7 : m1_250.R (~25 mins)
* Figure 8 : m1_1000.R (~2 hrs.)
* Figure 9 : m2_1000.R (~2 hrs.)
* Figure 10: m3_250.R (~2 hrs. 5 mins)
* Figure 11: m1_matern_misspec.R (~25 mins)
* Figure 12: m2_matern_misspec.R (~25 mins)
* Figure 13: m3_matern_misspec.R (~2 hrs. 5 mins)
* Figure 14: m1_smooth_misspec.R (~25 mins)
* Figure 15: m1_smooth_misspec.R (~25 mins)
* Figure 16: m2_smooth_misspec.R (~25 mins)
* Figure 17: m3_smooth_misspec.R (~2 hrs. 5 mins)

The timings are approximate and can vary significantly based on the machine configurations. There might be some issues with compiling the .Rcpp code in macOS with clang. Please look for solutions to bypass the problem as these can be very machine/OS version specific issues. 
