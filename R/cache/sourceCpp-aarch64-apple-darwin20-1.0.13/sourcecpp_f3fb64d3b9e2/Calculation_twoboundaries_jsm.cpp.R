`.sourceCpp_5_DLLInfo` <- dyn.load('/Users/wanni/Desktop/BOP2_PSO_Shiny_App/code/latest code/cache/sourceCpp-aarch64-apple-darwin20-1.0.13/sourcecpp_f3fb64d3b9e2/sourceCpp_6.so')

set_seed <- Rcpp:::sourceCppFunction(function(seed) {}, TRUE, `.sourceCpp_5_DLLInfo`, 'sourceCpp_5_set_seed')
Exacterror <- Rcpp:::sourceCppFunction(function(nobs, ncut, pnull, palter) {}, FALSE, `.sourceCpp_5_DLLInfo`, 'sourceCpp_5_Exacterror')
exact_error_recursive_Rcpp <- Rcpp:::sourceCppFunction(function(nobs, ncut1, ncut2, pnull, palter, ntotal) {}, FALSE, `.sourceCpp_5_DLLInfo`, 'sourceCpp_5_exact_error_recursive_Rcpp')
my_dbinom <- Rcpp:::sourceCppFunction(function(x, size, prob) {}, FALSE, `.sourceCpp_5_DLLInfo`, 'sourceCpp_5_my_dbinom')
my_pbinom <- Rcpp:::sourceCppFunction(function(q, size, prob) {}, FALSE, `.sourceCpp_5_DLLInfo`, 'sourceCpp_5_my_pbinom')
dbinom_product <- Rcpp:::sourceCppFunction(function(vXs, nobs, prob) {}, FALSE, `.sourceCpp_5_DLLInfo`, 'sourceCpp_5_dbinom_product')
GetocBiRcpp <- Rcpp:::sourceCppFunction(function(seed, nsim, contrast, nobs, b, b_grad1, b_grad2, pow1, pow2, pow3, dprior, ptrue, phi, delta0, delta1, fff) {}, FALSE, `.sourceCpp_5_DLLInfo`, 'sourceCpp_5_GetocBiRcpp')

rm(`.sourceCpp_5_DLLInfo`)
