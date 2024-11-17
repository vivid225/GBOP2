`.sourceCpp_1_DLLInfo` <- dyn.load('/Users/wanni/Desktop/BOP2_PSO_Shiny_App/code/latest code/cache/sourceCpp-aarch64-apple-darwin20-1.0.13/sourcecpp_f3fb67588637/sourceCpp_4.so')

set_seed <- Rcpp:::sourceCppFunction(function(seed) {}, TRUE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_set_seed')
Exacterror <- Rcpp:::sourceCppFunction(function(nobs, ncut, pnull, palter) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_Exacterror')
exact_error_recursive2_Rcpp <- Rcpp:::sourceCppFunction(function(nobs, ncut, pnull, palter, ntotal) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_exact_error_recursive2_Rcpp')
my_dbinom <- Rcpp:::sourceCppFunction(function(x, size, prob) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_my_dbinom')
my_pbinom <- Rcpp:::sourceCppFunction(function(q, size, prob) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_my_pbinom')
dbinom_product <- Rcpp:::sourceCppFunction(function(vXs, nobs, prob) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_dbinom_product')
GetocBiRcpp <- Rcpp:::sourceCppFunction(function(seed, nsim, contrast, nobs, b, b2, pow2, dprior, ptrue, phi, fff) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_GetocBiRcpp')
GridSearchBiRcpp <- Rcpp:::sourceCppFunction(function(seed, contrast, nobs, dprior, b1, b2, pow, pn, pa, cutstart, nsim, err1, ffff) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_GridSearchBiRcpp')

rm(`.sourceCpp_1_DLLInfo`)
