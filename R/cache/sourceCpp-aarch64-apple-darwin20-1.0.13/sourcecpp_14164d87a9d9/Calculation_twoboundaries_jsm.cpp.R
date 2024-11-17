`.sourceCpp_11_DLLInfo` <- dyn.load('/Users/wanni/Desktop/BOP2_PSO_Shiny_App/GBOP2/R/cache/sourceCpp-aarch64-apple-darwin20-1.0.13/sourcecpp_14164d87a9d9/sourceCpp_12.so')

set_seed <- Rcpp:::sourceCppFunction(function(seed) {}, TRUE, `.sourceCpp_11_DLLInfo`, 'sourceCpp_11_set_seed')
Exacterror <- Rcpp:::sourceCppFunction(function(nobs, ncut, pnull, palter) {}, FALSE, `.sourceCpp_11_DLLInfo`, 'sourceCpp_11_Exacterror')
exact_error_recursive_Rcpp <- Rcpp:::sourceCppFunction(function(nobs, ncut1, ncut2, pnull, palter, ntotal) {}, FALSE, `.sourceCpp_11_DLLInfo`, 'sourceCpp_11_exact_error_recursive_Rcpp')
my_dbinom <- Rcpp:::sourceCppFunction(function(x, size, prob) {}, FALSE, `.sourceCpp_11_DLLInfo`, 'sourceCpp_11_my_dbinom')
my_pbinom <- Rcpp:::sourceCppFunction(function(q, size, prob) {}, FALSE, `.sourceCpp_11_DLLInfo`, 'sourceCpp_11_my_pbinom')
dbinom_product <- Rcpp:::sourceCppFunction(function(vXs, nobs, prob) {}, FALSE, `.sourceCpp_11_DLLInfo`, 'sourceCpp_11_dbinom_product')
GetocBiRcpp <- Rcpp:::sourceCppFunction(function(seed, nsim, contrast, nobs, b, b_grad1, b_grad2, pow1, pow2, pow3, dprior, ptrue, phi, delta0, delta1, fff) {}, FALSE, `.sourceCpp_11_DLLInfo`, 'sourceCpp_11_GetocBiRcpp')

rm(`.sourceCpp_11_DLLInfo`)
