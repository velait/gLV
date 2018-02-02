
functions {
// function for the state equations
  real[] dz_dt(real t, 
               real[] z,  //state  
               real[] theta, //parameters
               real[] x_r, //data (unused) 
               int[] x_i) {
   real u = z[1];
   real v = z[2];
   real w = z[3];
   real r1 = theta[1];
   real r2 = theta[2];
   real r3 = theta[3];
   //
   real A11 = theta[4];
   real A12 = theta[5];
   real A13 = theta[6];
   real A21 = theta[7];
   real A22 = theta[8];
   real A23 = theta[9];
   real A31 = theta[10];
   real A32 = theta[11];
   real A33 = theta[12];
   
   real du_dt = (r1 + A11*u + A12*v + A13*w)*u;
   real dv_dt = (r2 + A21*u + A22*v + A23*w)*v;
   real dw_dt = (r3 + A31*u + A32*v + A33*w)*w;
   
   return { du_dt, dv_dt, dw_dt };
  }
}             

data {
  int<lower=0> N;      //number of measurements
  real ts[N];          //times
  real y0[3];          //initial measured populations
  real<lower=0> y[N, 3]; // measured populations
}

parameters {
  real theta[12];
  real<lower=0> z0[3];    //initial populations
  real<lower=0> sigma[3]; //measurement errors
}

transformed parameters {
  //population of measured years
  real z[N,3] 
    = integrate_ode_rk45(dz_dt, z0, 0, ts, theta,
                          rep_array(0.0, 0), rep_array(0,0), 
                          1e-6, 1e-5, 1e3);
}


model{
//priors
  sigma ~ lognormal(0, 0.5);
  theta[{1, 2, 3}] ~ normal(1, 0.5);
  theta[{4, 5, 6, 7, 8, 9, 10, 11, 12}] ~ normal(0.05, 0.05);

  z0[1] ~ lognormal(log(30), 5);
  z0[2] ~ lognormal(log(5), 5);
  z0[3] ~ lognormal(log(5), 5);
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
    for (n in 1:N) {
      y_sim[n, k] = lognormal_rng(log(z[n, k]), sigma[k]);
    } 
  }
}

