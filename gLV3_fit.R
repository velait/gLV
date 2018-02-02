### Simulate data for gLV3
library(deSolve)

LotVmod <- function (Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    ds1 = s1*(a - b*s2)
    ds2 = -s2*(c - d*s1)
    ds3 = s3*(e + f*s2 - g*s3)
    return(list(c(ds1, ds2, ds3)))
  })
}

#parameters:
Pars <- c(a=0.6, b=0.03, c=0.8, d=0.03, e=0.2, f=0.05, g=0.2)
State <- c(s1 = 10, s2 = 10, s3 = 10)
Time <- seq(0, 20, by = 1)

#make data frame
out <- as.data.frame(ode(func = LotVmod, y = State, parms = Pars, times = Time))

#plot
ggplot(data=out) + geom_line(aes(x=time, y=s1, colour="red")) + geom_line(aes(x=time, y=s2, colour="green")) + geom_line(aes(x=time, y=s3))

### Does STAN figure the parameters out?
# prepare data for stan
ts <- 1:(nrow(out)-1)
N <- nrow(out)-1
y0 <- c(out$s1[1], out$s2[1], out$s3[1])
y <- as.matrix(out[2:nrow(out), 2:4])
gLV3data <- list(N, ts, y0, y)

model <- stan_model("gLV3.stan")

fit <- sampling(model, data=gLV3data, chains=1, iter=1000, control=list(stepsize=5, adapt_delta=0.9), seed=123)

