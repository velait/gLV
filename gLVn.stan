// function for the state equations
functions {

// gLV: z = z*(r + Az)
// theta holds first the vector r, then the rows of matrix A
// s, S is the number of species, |theta| = s^2 + s

// nSpecies gives the number of species from theta
int nSpecies(real[] theta) {
int a = 1;
real D = 1+4*size(theta);
real S = pow(D,0.5);
S = S/2;

while(a <= S-1) {
a = a+1;
}
return a;
}

// this calculates the state for solver:
real[] dz_dt(real t, real[] z, real[] theta, real[] x_r, int[] x_i) {

// number of species
int S = nSpecies(theta);

//state array
real dz[S];

// states according to gLV
for(i in 1:S) {

// dotP is ith element of Az:
real dotP = 0;
for(j in 1:S) {
dotP = dotP + z[j]*theta[i*S+j];
}

dz[i] = (theta[i] + dotP)*z[i];
}

//return array
return dz;
}

}

data {
  int<lower=0> N;      //number of measurements
  real ts[N];          //times
  int<lower=0> S;      //number of species
  real y0[S];          //initial measured populations
  real<lower=0> y[N,S]; // measured populations
}

parameters {
  real<lower=0> z0[S];    //initial populations
  real<lower=0> sigma[S]; //measurement errors
}

transformed parameters {
  // |theta| = S + S^2
  int length_theta;
  length_theta = S + S*S;
  real theta[length_theta];
  
  //population of measured years
  real z[N,S] 
    = integrate_ode_rk45(dz_dt, z0, 0, ts, theta,
                          rep_array(0.0, 0), rep_array(0,0), 
                          1e-6, 1e-5, 1e3);
}


model{
//priors
//likelihood
y0 ~ lognormal(log(z0), sigma);

for(i in 1:S) {
  y[ ,i] ~ lognormal(log(z[,i]), sigma[i]]);
}
}


