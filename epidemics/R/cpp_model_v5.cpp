#include <Rcpp.h>
#include <vector>
#include <unordered_map>
using namespace Rcpp;
using namespace std;

// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
List cpp_model(IntegerVector population, double dog_ratio, NumericVector beta,
               NumericMatrix distance, IntegerVector init_inf, 
               int time_max, NumericVector inf_period_dist, 
               double birth_rate, double death_rate, double nu = 0, double mu = 0, 
               double epsilon = 0, double iota = 0, double connection_strength = 1,
               bool mobility = false) {

  ////////////////////////////////////////////////
  // Define objects
  int n_patch = population.size();
  IntegerMatrix sus(n_patch, time_max);
  IntegerMatrix exp(n_patch, time_max);
  IntegerMatrix inf(n_patch, time_max);
  IntegerMatrix exp_s(n_patch, time_max); // New exposed
  IntegerMatrix exp_d(n_patch, time_max); // Deaths among exposed
  IntegerMatrix inf_s(n_patch, time_max); // Start of infectiousness
  IntegerVector num(n_patch);
  NumericMatrix beta_matrix(n_patch, n_patch);
  unordered_map<string, vector<int>> data; // [unique ID]
  vector<vector<vector<int>>> inf_ids; // [patch][time]
  inf_ids.resize(n_patch);
  for (int patch = 0; patch < n_patch; ++patch) {
    inf_ids[patch].resize(time_max);
  }
  
  // Compute the beta_matrix
  NumericVector density_dependent(n_patch);
    
  if (iota) { 
    for (int i = 0; i < n_patch; ++i) {
      for (int j = 0; j < n_patch; ++j) {
        if (i != j && distance(i,j) != 0) {
          density_dependent[i] += pow(population[j], nu) / pow(distance(i,j), epsilon);
        }
      }
    } 
  }
  
  
  if (iota && epsilon) {
  for (int i = 0; i < n_patch; ++i) {
    for (int j = 0;j < n_patch; ++j) {
      if (j != i && distance(i,j) != 0) {
        beta_matrix(i,j) = beta[j] * connection_strength * pow(population[i], nu) * pow(population[j], mu) / (pow(distance(i,j), epsilon) * pow(density_dependent[i], iota)); 
      } else if (j == i) {
        beta_matrix(i,j) = beta[i];
      }
    }
  }
  }
  
  if (iota && !epsilon) {
    for (int i = 0; i < n_patch; ++i) {
      for (int j = 0;j < n_patch; ++j) {
        if (j != i && distance(i,j) != 0) {
          beta_matrix(i,j) = beta[j] * connection_strength * pow(population[i], nu) * pow(population[j], mu) / pow(density_dependent[i], iota); 
        } else if (j == i){
          beta_matrix(i,j) = beta[i];
        }
      }
    }
  }
  
  if (!iota && epsilon) {
    for (int i = 0; i < n_patch; ++i) {
      for (int j = 0;j < n_patch; ++j) {
        if (j != i && distance(i,j) != 0) {
          if (mobility) {
            beta_matrix(i,j) = beta[j] * connection_strength * pow(population[i], nu) * pow(population[j], mu) * pow(distance(i,j), epsilon); 
          } else {
            beta_matrix(i,j) = beta[j] * connection_strength * pow(population[i], nu) * pow(population[j], mu) / pow(distance(i,j), epsilon);
          }
        } else if (j == i){
          beta_matrix(i,j) = beta[i];
        }
      }
    }
  }
  
  if (!iota && !epsilon) {
    for (int i = 0; i < n_patch; ++i) {
      for (int j = 0;j < n_patch; ++j) {
        if (j != i && distance(i,j) != 0) {
          beta_matrix(i,j) = beta[j] * connection_strength * pow(population[i], nu) * pow(population[j], mu); 
        } else if (j == i){
          beta_matrix(i,j) = beta[i];
        }
      }
    }
  }
  
  // Rcout << beta_matrix << endl;  // DEBUG

  ////////////////////////////////////////////////
  // Set initial conditions
  for (int i = 0; i < n_patch; ++i) {
    sus(i, 0) = population[i]*dog_ratio - init_inf[i];
  }
  
  inf(_, 0) = init_inf;
  int curr_time = 0;
  
  for (int patch = 0; patch < n_patch; ++patch) {
    num[patch] = init_inf[patch];
    
    for (int id = 0; id < init_inf[patch]; ++id) {
      
      int death_delay = max(int(round(rexp(1, death_rate), 0)[0]), 1);
      int inf_period = sample(inf_period_dist.size(), 1, false, inf_period_dist, true)[0];
      int natural_death = 0;
      int stop = curr_time + inf_period;
      int hard_stop = min(stop, time_max);
      
      if (death_delay < inf_period) {
        natural_death = 1;
        stop = curr_time + death_delay;
      }
      
      for (int step = curr_time; step < hard_stop; ++step) {
        inf_ids[patch][step].push_back(id);
      }
      
      string unique_id = to_string(patch) + "_" + to_string(id);
      data[unique_id] = {-1, -1, -1, curr_time + death_delay, curr_time,
                         inf_period, -1, natural_death, int(hard_stop != stop)};
    }
  }
  
  ////////////////////////////////////////////////
  // Update loop
  NumericMatrix lambda(n_patch, n_patch);
  NumericVector proba_patch(n_patch);
  IntegerVector pop(n_patch);
  NumericMatrix foi(n_patch, n_patch);
  
  for (int curr_time = 1; curr_time < time_max; ++curr_time) {
    
    // Total population vector
    pop = sus(_, curr_time - 1) + exp(_, curr_time - 1) + inf(_, curr_time - 1);
    
    // Compute lambda
    for (int i = 0; i < n_patch; ++i) {
      for (int j = 0; j < n_patch; ++j) {
        // if (pop[j] == 0) {
        //   lambda(i,j) = 0;
        // } else {
          lambda(i, j) = beta_matrix(i, j) * sus(j, curr_time - 1) * inf(i, curr_time - 1); // / pop[j];
        // }
      }
    }
    
    if (curr_time == 1) {
      foi = clone(lambda);
    }
    
    // Update
    for (int patch = 0; patch < n_patch; ++patch) {
      
      double sum_lambda = sum(lambda(_, patch));
      int newborn = round(rpois(1, birth_rate * pop[patch]), 0)[0];
      int s_to_e = min(int(round(rpois(1, sum_lambda), 0)[0]), sus(patch, curr_time - 1));
      // Rcout << "Time: " << curr_time << " || Patch: " << patch << " || s_to_e: " << s_to_e << endl; // DEBUG
      int death_s = min(int(round(rpois(1, death_rate * sus(patch, curr_time - 1)), 0)[0]),
                        sus(patch, curr_time - 1) - s_to_e);
      
      proba_patch = lambda(_, patch);
      if (sum_lambda > 0.0) {
        proba_patch = proba_patch / sum_lambda;
      }
      
      for (int l = 0; l < n_patch; ++l) {
        if (proba_patch[l] < 0) {
          IntegerVector S = sus(_, curr_time - 1);
          IntegerVector E = exp(_, curr_time - 1);
          IntegerVector I = inf(_, curr_time - 1);
          Rcout << "Time = " << curr_time << " || " << S << " || " << E << " || " << I << endl;
          stop("Negative lambda!");
        } 
      }

      // Loop on new infected in patch patch and at time curr_time
      for (int id = num[patch]; id < num[patch] + s_to_e; ++id) {
        
        // ID
        string unique_id = to_string(patch) + "_" + to_string(id);
        
        // Times of events
        int death_delay = max(int(round(rexp(1, death_rate), 0)[0]), 1);
        int inc_period = max(int(round(rgamma(1, 2, 11.055), 0)[0]), 1);
        int inf_period = sample(inf_period_dist.size(), 1, false, inf_period_dist, true)[0];
        
        // Source 
        int patch_source = sample(n_patch, 1, false, proba_patch, false)[0];
        vector<int> ids = inf_ids[patch_source][curr_time - 1];
        IntegerVector vec(ids.begin(), ids.end());
        int id_source = sample(vec, 1)[0];
        
        //  Generation time
        string unique_id_source = to_string(patch_source) + "_" + to_string(id_source);
        int generation_time = curr_time + inc_period - data[unique_id_source][4];
        
        // Parameters to manage mortality
        int natural_death_e = 0;
        int natural_death_i = 0;
        int stop = curr_time + inc_period + inf_period;
        if (curr_time + inc_period < time_max) inf_s(patch, curr_time + inc_period) += 1;
        
        // Modify parameters according to the situation 
        if (death_delay < inc_period) {
          
          natural_death_e = 1;
          stop = curr_time + death_delay;
          if (curr_time + inc_period < time_max) inf_s(patch, curr_time + inc_period) -= 1;
          
        } else if (death_delay >= inc_period && death_delay < (inc_period + inf_period)) {
          
          natural_death_i = 1;
          stop = curr_time + death_delay;
          
        }
        
        int hard_stop = min(stop, time_max);
        
        // Fill data
        if (natural_death_e) {
          
          if (stop < time_max) exp_d(patch, stop) += 1;
          data[unique_id] = {patch_source, id_source, curr_time, stop, -1,
                             -1, -1, natural_death_e, int(hard_stop != stop)};
          
        } else if (natural_death_i){
          
          for (int step = curr_time + inc_period; step < hard_stop; ++step) {
            inf_ids[patch][step].push_back(id);
          }

          data[unique_id] = {patch_source, id_source, curr_time, stop,
                             curr_time + inc_period, inf_period,
                             generation_time, natural_death_i,
                             int(hard_stop != stop)};
          
        } else {
          
          for (int step = curr_time + inc_period; step < hard_stop; ++step) {
            inf_ids[patch][step].push_back(id);
          }
          
          data[unique_id] = {patch_source, id_source, curr_time, stop,
                             curr_time + inc_period, inf_period,
                             generation_time, natural_death_i,
                             int(hard_stop != stop)};
          
        }
      }
      
      // Counter of infected 
      num[patch] += s_to_e;
      
      //  Update counts
      int e_to_i = inf_s(patch, curr_time);
      sus(patch, curr_time) = sus(patch, curr_time - 1) - s_to_e + newborn - death_s;
      exp(patch, curr_time) = exp(patch, curr_time - 1) + s_to_e - e_to_i - exp_d(patch, curr_time);
      exp_s(patch, curr_time) = s_to_e;
      inf(patch, curr_time) = inf_ids[patch][curr_time].size();
    }
  }
  
  ////////////////////////////////////////////////
  // Outputs
  List dynamics = List::create(Named("sus") = sus, Named("exp") = exp, Named("inf") = inf, 
                                     Named("inc") = inf_s, Named("exp_s") = exp_s); 
  List sourcing = List::create(Named("data") = data, Named("inf_ids") = inf_ids);
  int epi_size = sum(inf_s);
  
  return List::create(Named("dynamics") = dynamics,
                      Named("sourcing") = sourcing,
                      Named("size") = epi_size, 
                      Named("foi") = foi, 
                      Named("beta_mat") = beta_matrix);

}


