#include <cmath>
#include <algorithm>
#include <vector>
#include <limits>
#include <Rcpp.h>
#include <math.h>
#include <Rmath.h>

using namespace Rcpp;


//#include <RcppEigen.h>
//using namespace Eigen;
using namespace std;




// set seed
// [[Rcpp::export]]
void set_seed(unsigned int seed) {
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(seed);
}

double dbinom_product( const NumericVector & vXs,  NumericVector & nobs, double prob);
double my_dbinom(double x, double size, double prob);
double my_pbinom(double q, double size, double prob);


// [[Rcpp::export]]
List Exacterror(NumericVector nobs, NumericVector ncut, double pnull, double palter)
{
  double dropn(0);
  double dropa(0);
  double prea(0);
  double ptsa(0);
  
  if (nobs.size() == 1){
    dropn = my_pbinom(ncut[0],nobs[0],pnull);
    dropa = my_pbinom(ncut[0],nobs[0],palter);
    ptsa = nobs[0];
  }else if(nobs.size() == 2){
    for(int i=0;i<=nobs[0];i++){
      if(i<=ncut[0]){
        dropn = dropn + my_dbinom(i,nobs[0],pnull);
        dropa = dropa + my_dbinom(i,nobs[0],palter);
        prea = prea + my_dbinom(i,nobs[0],palter);
        ptsa = ptsa + my_dbinom(i,nobs[0],palter)*nobs[0];
      }else{
        for(int j=0;j<=nobs[1]-nobs[0];j++){
          if(i+j<=ncut[1]){
            dropn = dropn + my_dbinom(i,nobs[0],pnull)*my_dbinom(j,nobs[1]-nobs[0],pnull);
            dropa = dropa + my_dbinom(i,nobs[0],palter)*my_dbinom(j,nobs[1]-nobs[0],palter);
            ptsa = ptsa + my_dbinom(i,nobs[0],palter)*my_dbinom(j,nobs[1]-nobs[0],palter)*nobs[1];
          }else{
            ptsa = ptsa + my_dbinom(i,nobs[0],palter)*my_dbinom(j,nobs[1]-nobs[0],palter)*nobs[1];
          }
        }
      }
    }
  }else if(nobs.size() == 3){
    for(int i=0;i<=nobs[0];i++){
      if(i<=ncut[0]){
        dropn = dropn + my_dbinom(i,nobs[0],pnull);
        dropa = dropa + my_dbinom(i,nobs[0],palter);
        prea = prea + my_dbinom(i,nobs[0],palter);
        ptsa = ptsa + my_dbinom(i,nobs[0],palter)*nobs[0];
      }else{
        for(int j=0;j<=nobs[1]-nobs[0];j++){
          if(i+j<=ncut[1]){
            dropn = dropn + my_dbinom(i,nobs[0],pnull)*my_dbinom(j,nobs[1]-nobs[0],pnull);
            dropa = dropa + my_dbinom(i,nobs[0],palter)*my_dbinom(j,nobs[1]-nobs[0],palter);
            prea = prea + my_dbinom(i,nobs[0],palter)*my_dbinom(j,nobs[1]-nobs[0],palter);
            ptsa = ptsa + my_dbinom(i,nobs[0],palter)*my_dbinom(j,nobs[1]-nobs[0],palter)*nobs[1];
          }else{
            for(int k=0;k<=nobs[2]-nobs[1];k++){
              if(i+j+k<=ncut[2]){
                dropn = dropn + my_dbinom(i,nobs[0],pnull)*my_dbinom(j,nobs[1]-nobs[0],pnull)*my_dbinom(k,nobs[2]-nobs[1],pnull);
                dropa = dropa + my_dbinom(i,nobs[0],palter)*my_dbinom(j,nobs[1]-nobs[0],palter)*my_dbinom(k,nobs[2]-nobs[1],palter);
                ptsa = ptsa + my_dbinom(i,nobs[0],palter)*my_dbinom(j,nobs[1]-nobs[0],palter)*my_dbinom(k,nobs[2]-nobs[1],palter)*nobs[2];
              }else{
                ptsa = ptsa + my_dbinom(i,nobs[0],palter)*my_dbinom(j,nobs[1]-nobs[0],palter)*my_dbinom(k,nobs[2]-nobs[1],palter)*nobs[2];
              }
            }
          }
        }
      }
    }
  }else if(nobs.size() == 4){
    for(int i=0;i<=nobs[0];i++){
      if(i<=ncut[0]){
        dropn = dropn + my_dbinom(i,nobs[0],pnull);
        dropa = dropa + my_dbinom(i,nobs[0],palter);
        prea = prea + my_dbinom(i,nobs[0],palter);
        ptsa = ptsa + my_dbinom(i,nobs[0],palter)*nobs[0];
      }else{
        for(int j=0;j<=nobs[1]-nobs[0];j++){
          if(i+j<=ncut[1]){
            dropn = dropn + my_dbinom(i,nobs[0],pnull)*my_dbinom(j,nobs[1]-nobs[0],pnull);
            dropa = dropa + my_dbinom(i,nobs[0],palter)*my_dbinom(j,nobs[1]-nobs[0],palter);
            prea = prea + my_dbinom(i,nobs[0],palter)*my_dbinom(j,nobs[1]-nobs[0],palter);
            ptsa = ptsa + my_dbinom(i,nobs[0],palter)*my_dbinom(j,nobs[1]-nobs[0],palter)*nobs[1];
          }else{
            for(int k=0;k<=nobs[2]-nobs[1];k++){
              if(i+j+k<=ncut[2]){
                dropn = dropn + my_dbinom(i,nobs[0],pnull)*my_dbinom(j,nobs[1]-nobs[0],pnull)*my_dbinom(k,nobs[2]-nobs[1],pnull);
                dropa = dropa + my_dbinom(i,nobs[0],palter)*my_dbinom(j,nobs[1]-nobs[0],palter)*my_dbinom(k,nobs[2]-nobs[1],palter);
                prea = prea + my_dbinom(i,nobs[0],palter)*my_dbinom(j,nobs[1]-nobs[0],palter)*my_dbinom(k,nobs[2]-nobs[1],palter);
                ptsa = ptsa + my_dbinom(i,nobs[0],palter)*my_dbinom(j,nobs[1]-nobs[0],palter)*my_dbinom(k,nobs[2]-nobs[1],palter)*nobs[2];
              }else{
                for(int m=0;m<=nobs[3]-nobs[2];m++){
                  if(i+j+k+m<=ncut[3]) {
                    dropn = dropn + my_dbinom(i,nobs[0],pnull)*my_dbinom(j,nobs[1]-nobs[0],pnull)*my_dbinom(k,nobs[2]-nobs[1],pnull)*my_dbinom(m,nobs[3]-nobs[2],pnull);
                    dropa = dropa + my_dbinom(i,nobs[0],palter)*my_dbinom(j,nobs[1]-nobs[0],palter)*my_dbinom(k,nobs[2]-nobs[1],palter)*my_dbinom(m,nobs[3]-nobs[2],palter);
                    ptsa = ptsa + my_dbinom(i,nobs[0],palter)*my_dbinom(j,nobs[1]-nobs[0],palter)*my_dbinom(k,nobs[2]-nobs[1],palter)*my_dbinom(m,nobs[3]-nobs[2],palter)*nobs[3];
                  }else{
                    ptsa = ptsa + my_dbinom(i,nobs[0],palter)*my_dbinom(j,nobs[1]-nobs[0],palter)*my_dbinom(k,nobs[2]-nobs[1],palter)*my_dbinom(m,nobs[3]-nobs[2],palter)*nobs[3];
                  }
                }
              }
            }
          }
        }
      }
    }
  }else if(nobs.size() == 5){
    for(int i=0;i<=nobs[0];i++){
      if(i<=ncut[0]){
        dropn = dropn + my_dbinom(i,nobs[0],pnull);
        dropa = dropa + my_dbinom(i,nobs[0],palter);
        prea  = prea  + my_dbinom(i,nobs[0],palter);
        ptsa  = ptsa  + my_dbinom(i,nobs[0],palter)*nobs[0];
      }else{
        for(int j=0;j<=nobs[1]-nobs[0];j++){
          if(i+j<=ncut[1]){
            //dropn = dropn + my_dbinom(i,nobs[0],pnull, false) * my_dbinom(j,nobs[1]-nobs[0],pnull, false);
            dropn = dropn + dbinom_product( NumericVector::create(i,j), nobs, pnull );
            dropa = dropa + dbinom_product( NumericVector::create(i,j), nobs, palter);
            prea  = prea  + dbinom_product( NumericVector::create(i,j), nobs, palter);
            ptsa  = ptsa  + dbinom_product( NumericVector::create(i,j), nobs, palter) * nobs[1];
          }else{
            for(int k=0;k<=nobs[2]-nobs[1];k++){
              if(i+j+k<=ncut[2]){
                dropn = dropn + dbinom_product( NumericVector::create(i,j,k), nobs, pnull );
                dropa = dropa + dbinom_product( NumericVector::create(i,j,k), nobs, palter);
                prea  = prea  + dbinom_product( NumericVector::create(i,j,k), nobs, palter);
                ptsa  = ptsa  + dbinom_product( NumericVector::create(i,j,k), nobs, palter) * nobs[2];
              }else{
                for(int m=0;m<=nobs[3]-nobs[2];m++){
                  if(i+j+k+m<=ncut[3]){
                    dropn = dropn + dbinom_product( NumericVector::create(i,j,k,m), nobs, pnull );
                    dropa = dropa + dbinom_product( NumericVector::create(i,j,k,m), nobs, palter);
                    prea  = prea  + dbinom_product( NumericVector::create(i,j,k,m), nobs, palter);
                    ptsa  = ptsa  + dbinom_product( NumericVector::create(i,j,k,m), nobs, palter) * nobs[3];
                  }else{
                    for(int n=0;n<=nobs[4]-nobs[3];n++){
                      if(i+j+k+m+n<=ncut[4]){
                        dropn = dropn + my_dbinom(i,nobs[0],pnull )*my_dbinom(j,nobs[1]-nobs[0],pnull )*my_dbinom(k,nobs[2]-nobs[1],pnull )*my_dbinom(m,nobs[3]-nobs[2],pnull )*my_dbinom(n,nobs[4]-nobs[3],pnull );
                        dropa = dropa + my_dbinom(i,nobs[0],palter)*my_dbinom(j,nobs[1]-nobs[0],palter)*my_dbinom(k,nobs[2]-nobs[1],palter)*my_dbinom(m,nobs[3]-nobs[2],palter)*my_dbinom(n,nobs[4]-nobs[3],palter);
                        ptsa  = ptsa  + my_dbinom(i,nobs[0],palter)*my_dbinom(j,nobs[1]-nobs[0],palter)*my_dbinom(k,nobs[2]-nobs[1],palter)*my_dbinom(m,nobs[3]-nobs[2],palter)*my_dbinom(n,nobs[4]-nobs[3],palter)*nobs[4];
                      }else{
                        ptsa = ptsa + my_dbinom(i,nobs[0],palter)*my_dbinom(j,nobs[1]-nobs[0],palter)*my_dbinom(k,nobs[2]-nobs[1],palter)*my_dbinom(m,nobs[3]-nobs[2],palter)*my_dbinom(n,nobs[4]-nobs[3],palter)*nobs[4];
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  List result = List::create(Named("edrop")=prea, Named("pts")=ptsa, Named("t1err")=1-dropn, Named("power")=1-dropa);
  return result;
}




// Define a struct to hold data
struct Data {
  std::vector<int> num_event;
  std::vector<double> prob;
};

// Define a struct to hold result
struct Result {
  Data effective_df;
  double effective_prob;
  double cummu_effective;
  Data nonstop_df;
  double expected_size;
  double stopping_prob;
  double nonstopping_prob;
  NumericVector prob_size_futile;
  NumericVector prob_size_effective;
  NumericVector prob_size_nonstop;
};

// Function to calculate exact error recursively
Result exacterror_recursive(Rcpp::NumericVector nobs, Rcpp::NumericVector ncut1, Rcpp::NumericVector ncut2, double test_prob, int ntotal) {
  int n_c = nobs.size();
  double cummu_effective = 0;
  double ptsa = 0;
  double stopping_prob = 0;
  double nonstopping_prob = 0;
  std::vector<double> prob_sizes_futile;
  std::vector<double> prob_sizes_effective;
  std::vector<double> prob_sizes_nonstop;
  
  if (n_c == 1) {
    Data effective_df;
    std::vector<double> effective_prob_vec;
    for (int i = ncut2[0]; i < nobs[0] + 1; i++) {
      effective_df.num_event.push_back(i);
      effective_prob_vec.push_back(R::dbinom(i, nobs[0], test_prob, false));
    }
    effective_df.prob = effective_prob_vec;
    
    Data nonstop_df;
    std::vector<double> nonstop_prob_vec;
    for (int i = ncut1[0] + 1; i < ncut2[0]; i++) {
      nonstop_df.num_event.push_back(i);
      nonstop_prob_vec.push_back(R::dbinom(i, nobs[0], test_prob, false));
    }
    nonstop_df.prob = nonstop_prob_vec;
    
    Data futile_df;
    std::vector<double> futile_prob_vec;
    for (int i = 0; i < ncut1[0] + 1; i++) {
      futile_df.num_event.push_back(i);
      futile_prob_vec.push_back(R::dbinom(i, nobs[0], test_prob, false));
    }
    futile_df.prob = futile_prob_vec;
    
    double futile_prob = std::accumulate(futile_df.prob.begin(), futile_df.prob.end(), 0.0);
    double effective_prob = std::accumulate(effective_df.prob.begin(), effective_df.prob.end(), 0.0);
    double nonstop_prob = std::accumulate(nonstop_df.prob.begin(), nonstop_df.prob.end(), 0.0);
      
    if (n_c == ntotal) {
      ptsa = ptsa + (futile_prob + effective_prob + nonstop_prob) * nobs[0];
    } else {
      ptsa = ptsa + futile_prob * nobs[0] + effective_prob * nobs[0];
      stopping_prob = stopping_prob + futile_prob + effective_prob;
      nonstopping_prob = nonstopping_prob + nonstop_prob;
    }
    
    // Rcout << "Round" << n_c << ptsa << std::endl;
    // Rcout << "Effprob, fut prob" << futile_prob << effective_prob << std::endl;
    // Rcout << "Stopping prob: " << stopping_prob << std::endl;
    // Rcout << "1-Nonstopping prob " << 1-nonstopping_prob << std::endl;
    prob_sizes_futile.insert(prob_sizes_futile.end(), futile_prob); 
    prob_sizes_effective.insert(prob_sizes_effective.end(), effective_prob); 
    prob_sizes_nonstop.insert(prob_sizes_nonstop.end(), nonstop_prob); 
    
    Result res_list;
    res_list.effective_df = effective_df;
    res_list.effective_prob = effective_prob;
    res_list.cummu_effective = cummu_effective + effective_prob;
    res_list.nonstop_df = nonstop_df;
    res_list.expected_size = ptsa;
    res_list.stopping_prob = stopping_prob;
    res_list.nonstopping_prob = nonstopping_prob;
    res_list.prob_size_futile = wrap(prob_sizes_futile);
    res_list.prob_size_effective = wrap(prob_sizes_effective);
    res_list.prob_size_nonstop = wrap(prob_sizes_nonstop);
    
    return res_list;
  } else {
    Rcpp::NumericVector nobs_new(nobs.begin(), nobs.end() - 1);
    // nobs_new.erase(nobs_new.end() - 1);
    Result res_before = exacterror_recursive(nobs_new, ncut1, ncut2, test_prob, ntotal);
    Data effective_df = res_before.effective_df;
    Data nonstop_df = res_before.nonstop_df;
    cummu_effective = res_before.cummu_effective;
    ptsa = res_before.expected_size;
    stopping_prob = res_before.stopping_prob;
    nonstopping_prob = res_before.nonstopping_prob;
    NumericVector prob_sizes_futile = res_before.prob_size_futile;
    NumericVector prob_sizes_effective = res_before.prob_size_effective;
    NumericVector prob_sizes_nonstop = res_before.prob_size_nonstop;
    
    Data effective_df_new;
    for (int i = 0; i < nonstop_df.num_event.size(); i++) {
      std::vector<int> newevent;
      for (int j = ncut2[n_c - 1]; j <= nobs[n_c - 1]; j++) {
        int temp_newevent = j - nonstop_df.num_event[i];
        if (temp_newevent <= nobs[n_c - 1] - nobs[n_c - 2] && temp_newevent >= 0) {
          newevent.push_back(temp_newevent);
        }
      }
      // for (int j = ncut2[n_c - 1]; j <= nobs[n_c - 1]; j++) {
      //   int adjusted_x = j - nonstop_df.num_event[i];
      //   if (adjusted_x >= 0) {
      //     newevent.push_back(adjusted_x);
      //   } else {
      //     newevent.push_back(0);
      //   }
      // }
      
      std::sort(newevent.begin(), newevent.end());
      auto last = std::unique(newevent.begin(), newevent.end());
      newevent.erase(last, newevent.end());
      
      double prob = nonstop_df.prob[i];
      Data newdf;
      for (int k = 0; k < newevent.size(); k++) {
        int new_event_value = nonstop_df.num_event[i] + newevent[k];
        double dbinom_result = R::dbinom(newevent[k], nobs[n_c - 1] - nobs[n_c - 2], test_prob, false);
        
        newdf.num_event.push_back(new_event_value);
        newdf.prob.push_back(prob * dbinom_result);
      }
      
      effective_df_new.num_event.insert(effective_df_new.num_event.end(), newdf.num_event.begin(), newdf.num_event.end());
      effective_df_new.prob.insert(effective_df_new.prob.end(), newdf.prob.begin(), newdf.prob.end());
    }
    
    Data nonstop_df_new;
    for (int i = 0; i < nonstop_df.num_event.size(); i++) {
      std::vector<int> newevent;
      for (int j = ncut1[n_c-1] + 1; j < ncut2[n_c-1]; j++) {
        int temp_newevent = j - nonstop_df.num_event[i];
        if (temp_newevent <= nobs[n_c - 1] - nobs[n_c - 2] && temp_newevent >= 0) {
          newevent.push_back(temp_newevent);
        }
      }
      // for (int j = ncut1[n_c-1] + 1; j < ncut2[n_c-1]; j++) {
      //   int adjusted_x = j - nonstop_df.num_event[i];
      //   if (adjusted_x >= 0) {
      //     newevent.push_back(adjusted_x);
      //   } else {
      //     newevent.push_back(0);
      //   }
      // }

      std::sort(newevent.begin(), newevent.end());
      auto last = std::unique(newevent.begin(), newevent.end());
      newevent.erase(last, newevent.end());

      double prob = nonstop_df.prob[i];
      Data newdf;
      for (int k = 0; k < newevent.size(); k++) {
        int new_event_value = nonstop_df.num_event[i] + newevent[k];
        double dbinom_result = R::dbinom(newevent[k], nobs[n_c - 1] - nobs[n_c - 2], test_prob, false);

        newdf.num_event.push_back(new_event_value);
        newdf.prob.push_back(prob * dbinom_result);
      }

      nonstop_df_new.num_event.insert(nonstop_df_new.num_event.end(), newdf.num_event.begin(), newdf.num_event.end());
      nonstop_df_new.prob.insert(nonstop_df_new.prob.end(), newdf.prob.begin(), newdf.prob.end());
    }

    // aggregate non_df_new:
    // Create a hashmap to store num_event as keys and their corresponding index in nonstop_df_new_sum as values
    Data nonstop_df_new_sum;
    std::unordered_map<int, int> num_event_index_map_null;

    for (int i = 0; i < nonstop_df_new.num_event.size(); i++) {
      int current_num_event = nonstop_df_new.num_event[i];
      auto it = num_event_index_map_null.find(current_num_event);

      if (it != num_event_index_map_null.end()) {
        // If num_event is already in nonstop_df_new_sum, update the probability
        nonstop_df_new_sum.prob[it->second] += nonstop_df_new.prob[i];
      } else {
        // If num_event is not found, add it to nonstop_df_new_sum and update the hashmap
        nonstop_df_new_sum.num_event.push_back(current_num_event);
        nonstop_df_new_sum.prob.push_back(nonstop_df_new.prob[i]);
        num_event_index_map_null[current_num_event] = nonstop_df_new_sum.num_event.size() - 1;
      }
    }

    // futile probability:
    Data futile_df_new;
    for (int i = 0; i < nonstop_df.num_event.size(); i++) {
      std::vector<int> newevent;
      for (int j = 0; j <= ncut1[n_c-1]; j++) {
        int temp_newevent = j - nonstop_df.num_event[i];
        if (temp_newevent <= nobs[n_c - 1] - nobs[n_c - 2] && temp_newevent >= 0) {
          newevent.push_back(temp_newevent);
        }
      }
      // for (int j = 0; j <= ncut1[n_c-1]; j++) {
      //   int adjusted_x = j - nonstop_df.num_event[i];
      //   if (adjusted_x >= 0) {
      //     newevent.push_back(adjusted_x);
      //   } else {
      //     newevent.push_back(0);
      //   }
      // }

      std::sort(newevent.begin(), newevent.end());
      auto last = std::unique(newevent.begin(), newevent.end());
      newevent.erase(last, newevent.end());

      double prob = nonstop_df.prob[i];
      Data newdf;
      for (int k = 0; k < newevent.size(); k++) {
        int new_event_value = nonstop_df.num_event[i] + newevent[k];
        double dbinom_result = R::dbinom(newevent[k], nobs[n_c - 1] - nobs[n_c - 2], test_prob, false);

        newdf.num_event.push_back(new_event_value);
        newdf.prob.push_back(prob * dbinom_result);
      }

      futile_df_new.num_event.insert(futile_df_new.num_event.end(), newdf.num_event.begin(), newdf.num_event.end());
      futile_df_new.prob.insert(futile_df_new.prob.end(), newdf.prob.begin(), newdf.prob.end());
    }

    double futile_prob = std::accumulate(futile_df_new.prob.begin(), futile_df_new.prob.end(), 0.0);
    double effective_prob = std::accumulate(effective_df_new.prob.begin(), effective_df_new.prob.end(), 0.0);
    double nonstop_prob = std::accumulate(nonstop_df_new_sum.prob.begin(), nonstop_df_new_sum.prob.end(), 0.0);

    // double test_prob2 = std::accumulate(nonstop_df.prob.begin(), nonstop_df.prob.end(), 0.0);
    
    if (n_c == ntotal) {
      ptsa = ptsa + (futile_prob + effective_prob + nonstop_prob) * nobs[n_c - 1];
      // Rcout << "Final Round Prob. " << futile_prob + effective_prob + nonstop_prob << std::endl;
      // Rcout << "Final Round Prob. (plus) " << futile_prob + effective_prob + nonstop_prob << std::endl;
      // Rcout << "Final Round Prob. (nonstop) " << test_prob2 << std::endl;
      
    } else {
      ptsa = ptsa + (effective_prob + futile_prob) * nobs[n_c - 1];
      stopping_prob = stopping_prob + futile_prob + effective_prob;
      nonstopping_prob = nonstopping_prob + nonstop_prob;
    }
    
    // Rcout << "Round " << n_c << ptsa << std::endl;
    // Rcout << "Stopping prob: " << stopping_prob << std::endl;
    // Rcout << "1-Nonstopping prob: " << 1-nonstopping_prob << std::endl;

    cummu_effective = cummu_effective + effective_prob;
    prob_sizes_futile.insert(prob_sizes_futile.end(), futile_prob); 
    prob_sizes_effective.insert(prob_sizes_effective.end(), effective_prob); 
    prob_sizes_nonstop.insert(prob_sizes_nonstop.end(), nonstop_prob); 

    Result res_list;
    res_list.effective_df = effective_df_new;
    res_list.effective_prob = effective_prob;
    res_list.cummu_effective = cummu_effective;
    res_list.nonstop_df = nonstop_df_new_sum;
    res_list.expected_size = ptsa;
    res_list.stopping_prob = stopping_prob;
    res_list.nonstopping_prob = nonstopping_prob;
    res_list.prob_size_futile = wrap(prob_sizes_futile);
    res_list.prob_size_effective = wrap(prob_sizes_effective);
    res_list.prob_size_nonstop = wrap(prob_sizes_nonstop);
    
    return res_list;
  }
}


// [[Rcpp::export]]
List exact_error_recursive_Rcpp(Rcpp::NumericVector nobs, Rcpp::NumericVector ncut1, Rcpp::NumericVector ncut2, double pnull, double palter, int ntotal) {
  // Call the C++ function
  Result result_null = exacterror_recursive(nobs, ncut1, ncut2, pnull, ntotal);
  Result result_alter = exacterror_recursive(nobs, ncut1, ncut2, palter, ntotal);
  
  // Convert Result to Rcpp List
  List result_list = List::create(
    Named("t1err") = result_null.cummu_effective,
    Named("power") = result_alter.cummu_effective,
    Named("pts_H0") = result_null.expected_size,
    Named("pts_Ha") = result_alter.expected_size,
    Named("prob_futile_null") = result_null.prob_size_futile,
    Named("prob_effective_null") = result_null.prob_size_effective,
    Named("prob_nonstop_null") = result_null.prob_size_nonstop,
    Named("prob_futile_alter") = result_alter.prob_size_futile,
    Named("prob_effective_alter") = result_alter.prob_size_effective,
    Named("prob_nonstop_alter") = result_alter.prob_size_nonstop
  );
  
  return result_list;
}



// [[Rcpp::export]]
double my_dbinom(double x, double size, double prob) {
  if (x < 0 || x > size) {
    return 0;
  }
  return R::choose(size, x) * std::pow(prob, x) * std::pow(1 - prob, size - x);
}

// [[Rcpp::export]]
double my_pbinom(double q, double size, double prob) {
  if (q < 0) {
    return 0;
  }
  if (q >= size) {
    return 1;
  }
  double sum = 0;
  double step = 1; // Increment step size
  for (double i = 0; i <= q; i += step) {
    Rcpp::Rcout << i << std::endl;
    sum += my_dbinom(i, size, prob) * step; // Multiply by step size
  }
  return sum;
}

// [[Rcpp::export]]
double dbinom_product(  const NumericVector & vXs,  NumericVector & nobs, double prob)
{
  int iMax = vXs.size();   // nobs.size() must be at least iMax also
  double dResult = my_dbinom(vXs[0], nobs[0], prob);
  for (int i = 1; i < iMax; ++i)
  {
    dResult *= my_dbinom(vXs[i], nobs[i] - nobs[i-1], prob);
  }
  return dResult;
}



//Binary Efficacy
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
List GetocBiRcpp(int seed, double nsim, NumericMatrix contrast, NumericVector nobs, double b, 
                 double b_grad1, double b_grad2, double pow1, double pow2, double pow3,
                 NumericVector dprior, double ptrue, double phi, double delta0, double delta1, Function fff)
  {

  int n1 =nobs.size();
  NumericVector stopbound1(2);
  // std::map<std::string, int> stopbound1;
  NumericVector stopbound_fut(n1);
  NumericVector stopbound_grad(n1);
  double nmax=max(nobs);
  
  List resultslist(5);

  set_seed(seed);

  for (int i=0; i<n1; i++)
  {
    double n=nobs[i];
    double lambda2 = b_grad1*pow(n/nmax,pow3); //cutoff1,cutoff2,cutoff2_lambda1,cutoff2_lambda2
    stopbound1 = fff(b*pow(n/nmax,pow1),pow(n/nmax,pow2),lambda2,b_grad2,dprior,n,phi,delta0,delta1,contrast);
    
    if (stopbound1[0] < 0 || stopbound1[1] < 0 ) { // check if there is -Inf bound
      resultslist[0]= 0; // edrop
      resultslist[1]= 0; // t1err
      resultslist[2]= 0; // power
      resultslist[3]= 1; // pts
      resultslist[4]= 1;
      return resultslist;
    }
    
    stopbound_fut[i]=stopbound1[0];
    stopbound_grad[i]=stopbound1[1];
  }
  
  for (int i=0; i<n1; i++){
    nobs[i] = std::ceil(nobs[i]);
  }
  // Rcpp::Rcout << "nobs " << nobs << std::endl;
  // Rcpp::Rcout << "Fut Boundary " << stopbound_fut << std::endl;
  // Rcpp::Rcout << "Grad Boundary " << stopbound_grad << std::endl;
  
  List temp(4);
  temp = exact_error_recursive_Rcpp(nobs, stopbound_fut, stopbound_grad, phi, ptrue, n1);

  resultslist[0]= 0; // edrop
  resultslist[1]= temp[0]; // t1err
  resultslist[2]= temp[1]; // power
  resultslist[3]= temp[2]; // pts
  resultslist[4]= temp[3]; // pts_Ha

  return resultslist;
}


// NumericMatrix GridSearchBiRcpp(int seed, NumericMatrix contrast, NumericVector nobs, NumericVector dprior,
//                                NumericVector b, NumericVector b_grad1, NumericVector pow1, NumericVector pow2, 
//                                double delta0, double delta1 double pn, double pa,
//                                double cutstart, double nsim, double err1, Function ffff)
// {
//   List temp1;
//   List temp2;
//   int k1=cutstart-1;
//   NumericVector ocmatbinary_v(6);
//   NumericMatrix ocmatbinary(6,1);
//   //ListNumericVector ocmatbinary = c();
//   // if nobs.size()<5, the calculation is based on exact computation.
//   // Rcpp::Rcout << "exact method " << std::endl;
//   for (int j=0; j<b.size(); j++)
//   {
//     for (int l=0; l<b2.size(); l++){
//       
//       for (int k = 0; k<pow.size(); k++)
//       {
//         
//         k1=k;
//         double b_mid=b[j];
//         double b_grad1_mid = b_grad1[l];
//         double pow1_mid=pow1[k];
//         double pow2_mid=pow2[m];
//         
//         temp1=GetocBiRcpp(seed, nsim, contrast, nobs, b_mid, b_grad1_mid, pow1_mid, pow2_mid, dprior, pa, pn, delta0, delta1, ffff);
//         // power_default = GetocBiRcpp(seed=input$seed, nsim=input$numOfSimForTiralSetting,
//         //                             contrast=as.matrix(1), nobs=(nobs2),
//         //                             b=res$par[,1], b_grad1=res$par[,6],
//         //                                                           
//         //                             pow1=res$par[,2], pow2 = res$par[,5],dprior= c(p.n,1-p.n), 
//         //                             ptrue=p.a, phi=p.n,delta0=inputlist$delta0,delta1=inputlist$delta1, 
//         //                             maxresp)
//         
//         // Rcpp::Rcout << "temp1: " << temp1_pts << std::endl;
//         ocmatbinary_v[0]=b[j];
//         ocmatbinary_v[1]=pow[k];
//         ocmatbinary_v[2]=temp1[1]; // t1err
//         ocmatbinary_v[3]=temp1[2]; // power
//         ocmatbinary_v[4]=temp1[3]; // pts
//         ocmatbinary_v[5]=b2[l];
//         //ocmatbinary(_,0)=ocmatbinary_v;
//         ocmatbinary = cbind(ocmatbinary, ocmatbinary_v);
//         
//         
//       }
//     }
//     
//   }
//   
//   return(ocmatbinary);
// }




#include <Rcpp.h>
#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// set_seed
void set_seed(unsigned int seed);
RcppExport SEXP sourceCpp_11_set_seed(SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type seed(seedSEXP);
    set_seed(seed);
    return R_NilValue;
END_RCPP
}
// Exacterror
List Exacterror(NumericVector nobs, NumericVector ncut, double pnull, double palter);
RcppExport SEXP sourceCpp_11_Exacterror(SEXP nobsSEXP, SEXP ncutSEXP, SEXP pnullSEXP, SEXP palterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type nobs(nobsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ncut(ncutSEXP);
    Rcpp::traits::input_parameter< double >::type pnull(pnullSEXP);
    Rcpp::traits::input_parameter< double >::type palter(palterSEXP);
    rcpp_result_gen = Rcpp::wrap(Exacterror(nobs, ncut, pnull, palter));
    return rcpp_result_gen;
END_RCPP
}
// exact_error_recursive_Rcpp
List exact_error_recursive_Rcpp(Rcpp::NumericVector nobs, Rcpp::NumericVector ncut1, Rcpp::NumericVector ncut2, double pnull, double palter, int ntotal);
RcppExport SEXP sourceCpp_11_exact_error_recursive_Rcpp(SEXP nobsSEXP, SEXP ncut1SEXP, SEXP ncut2SEXP, SEXP pnullSEXP, SEXP palterSEXP, SEXP ntotalSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type nobs(nobsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type ncut1(ncut1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type ncut2(ncut2SEXP);
    Rcpp::traits::input_parameter< double >::type pnull(pnullSEXP);
    Rcpp::traits::input_parameter< double >::type palter(palterSEXP);
    Rcpp::traits::input_parameter< int >::type ntotal(ntotalSEXP);
    rcpp_result_gen = Rcpp::wrap(exact_error_recursive_Rcpp(nobs, ncut1, ncut2, pnull, palter, ntotal));
    return rcpp_result_gen;
END_RCPP
}
// my_dbinom
double my_dbinom(double x, double size, double prob);
RcppExport SEXP sourceCpp_11_my_dbinom(SEXP xSEXP, SEXP sizeSEXP, SEXP probSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< double >::type prob(probSEXP);
    rcpp_result_gen = Rcpp::wrap(my_dbinom(x, size, prob));
    return rcpp_result_gen;
END_RCPP
}
// my_pbinom
double my_pbinom(double q, double size, double prob);
RcppExport SEXP sourceCpp_11_my_pbinom(SEXP qSEXP, SEXP sizeSEXP, SEXP probSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    Rcpp::traits::input_parameter< double >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< double >::type prob(probSEXP);
    rcpp_result_gen = Rcpp::wrap(my_pbinom(q, size, prob));
    return rcpp_result_gen;
END_RCPP
}
// dbinom_product
double dbinom_product(const NumericVector& vXs, NumericVector& nobs, double prob);
RcppExport SEXP sourceCpp_11_dbinom_product(SEXP vXsSEXP, SEXP nobsSEXP, SEXP probSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type vXs(vXsSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type nobs(nobsSEXP);
    Rcpp::traits::input_parameter< double >::type prob(probSEXP);
    rcpp_result_gen = Rcpp::wrap(dbinom_product(vXs, nobs, prob));
    return rcpp_result_gen;
END_RCPP
}
// GetocBiRcpp
List GetocBiRcpp(int seed, double nsim, NumericMatrix contrast, NumericVector nobs, double b, double b_grad1, double b_grad2, double pow1, double pow2, double pow3, NumericVector dprior, double ptrue, double phi, double delta0, double delta1, Function fff);
RcppExport SEXP sourceCpp_11_GetocBiRcpp(SEXP seedSEXP, SEXP nsimSEXP, SEXP contrastSEXP, SEXP nobsSEXP, SEXP bSEXP, SEXP b_grad1SEXP, SEXP b_grad2SEXP, SEXP pow1SEXP, SEXP pow2SEXP, SEXP pow3SEXP, SEXP dpriorSEXP, SEXP ptrueSEXP, SEXP phiSEXP, SEXP delta0SEXP, SEXP delta1SEXP, SEXP fffSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< double >::type nsim(nsimSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type contrast(contrastSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type nobs(nobsSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type b_grad1(b_grad1SEXP);
    Rcpp::traits::input_parameter< double >::type b_grad2(b_grad2SEXP);
    Rcpp::traits::input_parameter< double >::type pow1(pow1SEXP);
    Rcpp::traits::input_parameter< double >::type pow2(pow2SEXP);
    Rcpp::traits::input_parameter< double >::type pow3(pow3SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type dprior(dpriorSEXP);
    Rcpp::traits::input_parameter< double >::type ptrue(ptrueSEXP);
    Rcpp::traits::input_parameter< double >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< double >::type delta0(delta0SEXP);
    Rcpp::traits::input_parameter< double >::type delta1(delta1SEXP);
    Rcpp::traits::input_parameter< Function >::type fff(fffSEXP);
    rcpp_result_gen = Rcpp::wrap(GetocBiRcpp(seed, nsim, contrast, nobs, b, b_grad1, b_grad2, pow1, pow2, pow3, dprior, ptrue, phi, delta0, delta1, fff));
    return rcpp_result_gen;
END_RCPP
}
