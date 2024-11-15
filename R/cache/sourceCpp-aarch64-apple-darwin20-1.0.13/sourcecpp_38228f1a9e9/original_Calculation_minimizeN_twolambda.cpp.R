`.sourceCpp_7_DLLInfo` <- dyn.load('/Users/wanni/Desktop/BOP2_PSO_Shiny_App/code/latest code/cache/sourceCpp-aarch64-apple-darwin20-1.0.13/sourcecpp_38228f1a9e9/sourceCpp_8.so')

set_seed <- Rcpp:::sourceCppFunction(function(seed) {}, TRUE, `.sourceCpp_7_DLLInfo`, 'sourceCpp_7_set_seed')
Exacterror <- Rcpp:::sourceCppFunction(function(nobs, ncut, pnull, palter) {}, FALSE, `.sourceCpp_7_DLLInfo`, 'sourceCpp_7_Exacterror')
exact_error_recursive2_Rcpp <- Rcpp:::sourceCppFunction(function(nobs, ncut, pnull, palter, ntotal) {}, FALSE, `.sourceCpp_7_DLLInfo`, 'sourceCpp_7_exact_error_recursive2_Rcpp')
my_dbinom <- Rcpp:::sourceCppFunction(function(x, size, prob) {}, FALSE, `.sourceCpp_7_DLLInfo`, 'sourceCpp_7_my_dbinom')
my_pbinom <- Rcpp:::sourceCppFunction(function(q, size, prob) {}, FALSE, `.sourceCpp_7_DLLInfo`, 'sourceCpp_7_my_pbinom')
dbinom_product <- Rcpp:::sourceCppFunction(function(vXs, nobs, prob) {}, FALSE, `.sourceCpp_7_DLLInfo`, 'sourceCpp_7_dbinom_product')
GetocBiRcpp <- Rcpp:::sourceCppFunction(function(seed, contrast, nobs, b, b2, pow2, dprior, ptrue, phi, fff) {}, FALSE, `.sourceCpp_7_DLLInfo`, 'sourceCpp_7_GetocBiRcpp')

rm(`.sourceCpp_7_DLLInfo`)
