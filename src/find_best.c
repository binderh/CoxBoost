#include <math.h>
//#include <stdio.h>

void get_I_vec(double *x, int *n, int *p, int *uncens_n, 
               double *weighted_risk, double *weighted_risk_sum,double *I_vec) 
{
    int actual_x, i, j;
    double x_bar,V;
    
    for (actual_x = 0; actual_x < *p; actual_x++) {        
        for (i = 0; i < *uncens_n; i++) {
            x_bar = 0;
            for (j = 0; j < *n; j++) x_bar += weighted_risk[(i * *n) + j] * x[(actual_x * *n) + j];
            x_bar /= weighted_risk_sum[i];
            
            V = 0;
            for (j = 0; j < *n; j++) V += weighted_risk[(i * *n) + j] * (x[(actual_x * *n) + j] - x_bar) * (x[(actual_x * *n) + j] - x_bar);
            V /= weighted_risk_sum[i]; 
            I_vec[actual_x] += V;
        }
    }
}


void find_best_candidate(double *x, int *n, int *p, int *uncens, int *uncens_n, double *beta, double *risk, double *eta, 
               double *weight_at_event, int *max_nz, int *max_1, double *weighted_risk, double *weighted_risk_sum,
               int *weight_time_dependent, 
               double *penalty, int *penscore, int *candidates, int *candno, int *is_01,
               int *warncount, int *min_index, double *min_beta_delta, double *score_vec, double *beta_delta_vec,
               double *U_vec, double *I_vec)
{
    int actual_cand,actual_x, i, j;
    double x_bar,prev_x_bar_base,V;
    double U, I, beta_delta;

    double buffer1, prev_buffer1;
    
    double score, max_score, max_score_beta_delta;
    int max_score_index;

    const double stability_penalty = 1.0/0.1 - 1; 
    
    *warncount = 0;
    
    for (actual_cand = 0; actual_cand < *candno; actual_cand++) {
        actual_x = candidates[actual_cand];
        
        U = 0;
        I = 0;
        
        if (*is_01) {
            prev_buffer1 = 0;
        
            for (i = 0; i < *uncens_n; i++) {
                if (*weight_time_dependent) {
                    buffer1 = 0;
                    for (j = 0; j < max_nz[i]; j++) {
                        buffer1 += weighted_risk[(i * *n) + j]*x[(actual_x * *n) + j];
                    }                    
                } else {
                    buffer1 = prev_buffer1;
    
                    for (j = max_1[i]; j < max_1[i+1]; j++) {
                        buffer1 += weighted_risk[(i * *n) + j]*x[(actual_x * *n) + j];
                    }
    
                    prev_buffer1 = buffer1;            
    
                    for (j = max_1[i+1]; j < max_nz[i]; j++) {
                        buffer1 += weighted_risk[(i * *n) + j]*x[(actual_x * *n) + j];
                    }
                }
    
                x_bar = buffer1/weighted_risk_sum[i];
                U += weight_at_event[i]*(x[(actual_x * *n) + uncens[i]] - x_bar); 
                I += weight_at_event[i]*(x_bar*x_bar*(weighted_risk_sum[i] - buffer1) + (1-x_bar)*(1-x_bar)*buffer1) / weighted_risk_sum[i];
            }        
        } else {        
            prev_x_bar_base = 0;
            
            for (i = 0; i < *uncens_n; i++) {
                if (*weight_time_dependent) {
                    x_bar = 0;
                    for (j = 0; j < max_nz[i]; j++) x_bar += weighted_risk[(i * *n) + j] * x[(actual_x * *n) + j];
                } else {
                    x_bar = prev_x_bar_base;
                    for (j = max_1[i]; j < max_1[i+1]; j++) x_bar += weighted_risk[(i * *n) + j] * x[(actual_x * *n) + j];
                    prev_x_bar_base = x_bar;
                    for (j = max_1[i+1]; j < max_nz[i]; j++) x_bar += weighted_risk[(i * *n) + j] * x[(actual_x * *n) + j];
                }
    
                x_bar /= weighted_risk_sum[i];
                U += weight_at_event[i]*(x[(actual_x * *n) + uncens[i]] - x_bar); 
                
                V = 0;
                for (j = 0; j < max_nz[i]; j++) V += weighted_risk[(i * *n) + j] * (x[(actual_x * *n) + j] - x_bar) * (x[(actual_x * *n) + j] - x_bar);
                V /= weighted_risk_sum[i]; 
                I += weight_at_event[i]*V;
            }
        }
        
        U_vec[actual_x] = U;
        I_vec[actual_x] = I;
        
        beta_delta_vec[actual_x] = beta_delta = U/(I + penalty[actual_x]);
        
        if (*penscore) {
            score_vec[actual_x] = score = U*U/(I + penalty[actual_x]);
        } else {
            score_vec[actual_x] = score = U*U/(I + stability_penalty);
        }

        if (actual_cand == 0 || score > max_score) {
            max_score_index = actual_x + 1;
            max_score = score;
            max_score_beta_delta = beta_delta;
        }        
    }

    *min_index = max_score_index;
    *min_beta_delta = max_score_beta_delta;
    
    //*min_deviance = 0; 
    //
    //for (i = 0; i < *uncens_n; i++) {
    //    ldenom = 0;
    //    for (j = 0; j < *n; j++) ldenom += weights[(i * *n) + j] * risk[j] * exp(x[((*min_index - 1) * *n) + j]* *min_beta_delta);
    //    *min_deviance += eta[uncens[i]] + x[((*min_index - 1) * *n) + uncens[i]]* *min_beta_delta - log(ldenom); 
    //}

    //printf("min %d %d\n",max_score_index,*min_index);
}

