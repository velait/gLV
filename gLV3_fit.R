### Simulate data for gLV3
library(deSolve)

LotVmod <- function (Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    dx = x*(a - b*y)
    dy = -y*(c - d*x)
    dz = z*(e + f*y - g*z)
    return(list(c(dx, dy, dz)))
  })
}

#parameters:
Pars <- c(a=0.6, b=0.03, c=0.8, d=0.03, e=0.2, f=0.05, g=0.2)
State <- c(x = 10, y = 10, z = 10)
Time <- seq(0, 20, by = 1)

#make data frame
out <- as.data.frame(ode(func = LotVmod, y = State, parms = Pars, times = Time))

#plot
ggplot(data=out) + geom_line(aes(x=time, y=x, colour="red")) + geom_line(aes(x=time, y=y, colour="green")) + geom_line(aes(x=time, y=z))

### Does STAN figure the parameters out?

model <- stan_model("gLV3.stan")

