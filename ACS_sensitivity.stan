
// saved as ACS_sensitivity.stan

functions{
  real tprob(vector dtime, vector tpos, vector wgt, real sens, real lamb){
    vector[num_elements(dtime)] pp;
    vector[num_elements(dtime)] logl;
    real ll;
    for (i in 1:num_elements(dtime)){
      pp[i] =  sens*exp(-(dtime[i]*lamb));
      logl[i] =  tpos[i] == 1 ? tpos[i]*log(pp[i])*wgt[i] 
                               : (1-tpos[i])*log(1-pp[i])*wgt[i];
    } 
    ll = sum(logl);
    return ll;
  }   
}

data{
  real sens_prior[2]; // beta prior
  real lamb_prior[2]; // gamma prior
  int<lower=0> J;            // number of cancers in ACS study 
  vector[J] tpos;            // indicator of detection or no detection by blood test
  vector[J] dtime;           // time of the clinical diagnosis 
  vector[J] wgt;             // inverse probability weight
}

parameters{
  real<lower=0, upper=1> sens;                // population treatment effect
  real<lower=0.05> lamb;                         // the hazard for exponential distribution 
}

model {
  lamb ~ gamma(lamb_prior[1],lamb_prior[2]); //mean 0.5 se 1.5
  sens ~ beta(sens_prior[1],sens_prior[2]); // mean 0.4, se 0.2
  #lamb ~ gamma(0.125,0.25); //mean 0.5 se 1.5
  #sens ~ beta(2,3); // mean 0.4, se 0.2
  target += tprob(dtime,tpos, wgt, sens, lamb); // log-likelihood
}
generated quantities {
  real<lower=0> mst;
  real<lower=0,upper=1> year1_rate;
  real<lower=0,upper=1> year2_rate;
  real<lower=0,upper=1> year3_rate;
  mst = 1/lamb; // mean sojourn time
  year1_rate = sens/lamb * (1- exp(-lamb));
  year2_rate = sens/lamb * (1- exp(-2*lamb))/2;
  year3_rate = sens/lamb * (1- exp(-3*lamb))/3;
}

