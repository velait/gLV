// function for the state equations
real[] dz_dt(real t, //time
             real[] z,  //state  
             real[] theta, //parameters
             real[] x_l, //data (unused) 
             int[] x_i //unused 
             ) {
             
             real u = z[1];
             real v = z[2];
             real w = z[3];
             real r1 = theta[1];
             real r2 = theta[2];
             real r3 = theta[3];
             
             real du_dt = (r1 + theta[1]*u + theta[2]*v + theta[3]*w)*u;
             real dv_dt = (r2 + theta[4]*u + theta[5]*v + theta[6]*w)*v;
             real dw_dt = (r3 + theta[7]*u + theta[8]*v + theta[9]*w)*w;
             
             return { du_dt, dv_dt, dw_dt };
             }

data {
  int<lower=0> N;      //number of measurements
  real ts[N];          //times
  real y0[3];          //initial measured populations
  real<lower=0> y[N,3]; // measured populations
}

parameters {
  real theta[9];
  real<lower=0> z0[3];    //initial populations
  real<lower=0> sigma[3]; //measurement errors
}

transformed parameters {
  //population of measured years
  real z[N,3] 
    = integrate_ode_rk45(dz_dt, dz0, 0, ts, theta,
                          rep_array(0.0, 0), rep_array(0,0), 
                          1e-6, 1e-5, 1e3);
}


model{
//priors

//likelihood
y0 ~ lognormal(log(z0), sigma);

for(k in 1:3) {
  y[ ,k] ~ lognormal(log(z[,k]), sigma[k]);
}
}

generated quantities {
  real y0_sim[3];
  real y_sim[N, 3];

  for (k in 1:3) {
    y0_sim[k] = lognormal_rng(log(z0[k]), sigma[k]);
    for (n in 1:N)
      y_sim[n, k] = lognormal_rng(log(z[n, k]), sigma[k]);
  }
}

