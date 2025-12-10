
// saved as ACS_sensitivity_early_late.stan

functions{
  real tprob(vector dtime_early, vector tpos_early, vector wgt_early, vector dtime_late, vector tpos_late, vector wgt_late, real sens_early, real sens_late, real lamb_sum, real lamb_end, real p_late){
    vector[num_elements(dtime_early)] pp1;
    vector[num_elements(dtime_early)] logl1;
  
    vector[num_elements(dtime_late)] pp2;
    vector[num_elements(dtime_late)] logl2;
  
    real ll;
    int n_early;
    int n_late;
    
    n_early = num_elements(dtime_early);
    n_late = num_elements(dtime_late);
    
    for (i in 1:num_elements(dtime_early)){
      pp1[i] =  sens_early*exp(-(dtime_early[i]*(lamb_sum)));
      logl1[i] =  tpos_early[i] == 1 ? tpos_early[i]*log(pp1[i])*wgt_early[i] 
                               : (1-tpos_early[i])*log(1-pp1[i])*wgt_early[i];
    } 
    
    for (i in 1:num_elements(dtime_late)){
      pp2[i] =  sens_late*exp(-(dtime_late[i]*lamb_end)) + sens_early*(lamb_end/(lamb_sum-lamb_end))*(exp(-(dtime_late[i]*lamb_end)) - exp(-(dtime_late[i]*lamb_sum))) ;
      logl2[i] =  tpos_late[i] == 1 ? tpos_late[i]*log(pp2[i])*wgt_late[i] 
                               : (1-tpos_late[i])*log(1-pp2[i])*wgt_late[i];
    } 
    
    ll = n_early*log(1-p_late) + n_late*log(p_late) + sum(logl1) + sum(logl2);
    return ll;
  }   
}

data{
  real sens_early_prior[2]; // beta prior
  real sens_late_prior[2]; //beta prior
  real lamb_sum_prior[2]; // gamma prior
  real lamb_end_prior[2]; // gamma prior
  
  int<lower=0> J1;                         // number of early stage cancers in ACS study 
  int<lower=0> J2;                         // number of late stage cancers in ACS study 
  
  vector[J1] tpos_early;            // indicator of detection or no detection by blood test for early stage
  vector[J1] dtime_early;           // time of the clinical diagnosis for early stage
  vector[J1] wgt_early;             // inverse probability weight for early stage
  vector[J2] tpos_late;             // indicator of detection or no detection by blood test for late stage
  vector[J2] dtime_late;            // time of the clinical diagnosis for late stage
  vector[J2] wgt_late;              // inverse probability weight for late stage
}

parameters{
  real<lower=0, upper=1> sens_early;                // population treatment effect
  real<lower=0, upper=1> sens_late;                 // population treatment effect
  
  real<lower=0.1> lamb_sum;                        // the hazard for exponential distribution 
  real<lower=0.1> lamb_end;
  real<lower=0, upper=1> p_late; 
}

model {
  sens_early ~ beta(sens_early_prior[1], sens_early_prior[2]); // CCGA stage I/II/III 40.7%; mean 0.33 se 0.24
  sens_late  ~ beta(sens_late_prior[1], sens_late_prior[2]); // mean 0.67 se 0.24
  #sens_early ~ beta(3,7); 
  #sens_late ~ beta(7,3); 
  lamb_sum ~ gamma(lamb_sum_prior[2],lamb_sum_prior[2]);
  lamb_end ~ gamma(lamb_end_prior[2],lamb_end_prior[2]);
  #lamb_sum ~ gamma(0.5,0.5); // mean 1 se 1.37
  #lamb_end ~ gamma(0.25,0.5);
  p_late ~ beta(1,4);
  target += tprob(dtime_early,tpos_early, wgt_early,dtime_late,tpos_late,wgt_late,sens_early, sens_late, lamb_sum, lamb_end, p_late); // log-likelihood
}

generated quantities {
  real<lower=0> mst1;
  real<lower=0> mst2;
  real<lower=0> mst3;
  real<lower=0> mst4;
  real<lower=0> mst;
  
  mst1 = 1/lamb_sum; // mean sojourn time in early stage
  mst2 = 1/lamb_end; // mean sojourn time in late stage
  mst3 = 1/(lamb_sum*(1-p_late)); // mean sojourn time in late stage
  mst4 = 1/(lamb_sum*p_late); // mean sojourn time in late stage
  mst = mst1 + p_late*mst2; // mean sojourn time
}

