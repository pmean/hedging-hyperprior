# These are all the preliminary libraries and the custom built functions that I need.

library(dplyr)
library(magrittr)
library(plot3D)
library(rjags)
library(rstan)
library(tidyr)

# start with a blank slate
rm(list=ls())

# these are some recommended setting for stan
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

cycle_fast <- function(g=1001, cmax=1, rep=g) {
  cyc <- seq(0, cmax, length=g)
  return(rep(cyc, times=rep))
}

cycle_slow <- function(g=1001, cmax=1, rep=g) {
  cyc <- seq(0, cmax, length=g)
  return(rep(cyc, each=rep))
}
cycle_fast(3)
cycle_slow(3)

matrix(cycle_fast(3), nrow=3, byrow=TRUE)

calculate_prior <- function(pi, tau, g=21, alpha=4, beta=16, d=1, hp=dens.uniform) {
  v <- dbeta(pi, d+(alpha-d)*tau, d+(beta-d)*tau)*hp(tau)
  # with byrow=TRUE, cycle_slow is the rows (x) and cycle_fast is the columns (y)
  matrix(v, nrow=g, byrow=TRUE)
}

draw_prior <- function(hp=dens.uniform, g=33, alpha=4, beta=16,
                       tau_max=1, d=1, zmax=4.7) {
  tau <- cycle_slow(g, tau_max)
  pi  <- cycle_fast(g)
  d <- calculate_prior(pi, tau, g, alpha, beta, d, hp)
  persp3D(z=d, xlab="tau", ylab="pi",zlab="Density", zlim=c(0, zmax),
          theta=th, phi=ph, colvar=NULL, col="white", border="black",
          bty="u", col.axis="gray90")
  return(d)
}

tau_slice <- function(tau_cut, hp=dens.uniform, g=33, alpha=4, beta=16,
                      tau_max=1, d=1, zmax=4.7) {
  delta <- 1E-2
  tau <- cycle_slow(g, tau_max)
  pi  <- cycle_fast(g)
  d <- calculate_prior(pi, tau, g, alpha, beta, d, hp)
  d[tau > tau_cut+delta] <- NA
  persp3D(z=d, xlab="tau", ylab="pi",zlab="Density", zlim=c(0, zmax),
          theta=th, phi=ph, colvar=NULL, col="white", border="gray90", 
          bty="u", col.axis="gray90")
  d[tau < tau_cut-delta] <- NA
  persp3D(z=d, xlab="tau", ylab="pi",zlab="Density", zlim=c(0, zmax),
          theta=th, phi=ph, colvar=NULL, col="white", border="black", add=TRUE)
  return(d)
}

pi_slice <- function(pi_cut, hp=dens.uniform, g=33, alpha=4, beta=16,
                     tau_max=1, d=1, zmax=4.7) {
  delta <- 1E-2
  tau <- cycle_slow(g, tau_max)
  pi  <- cycle_fast(g)
  d <- calculate_prior(pi, tau, g, alpha, beta, d, hp)
  d[pi > pi_cut+delta] <- NA
  persp3D(z=d, xlab="tau", ylab="pi",zlab="Density", zlim=c(0, zmax),
          theta=th, phi=ph, colvar=NULL, col="white", border="gray90",
          bty="u", col.axis="gray90")
  d[pi < pi_cut-delta] <- NA
  persp3D(z=d, xlab="tau", ylab="pi",zlab="Density", zlim=c(0, zmax),
          theta=th, phi=ph, colvar=NULL, col="white", border="black", add=TRUE)
  return(d)
}

ct <- function(la,sr=0) {
  plot(0:1,0:1,type="n",axes=FALSE,xlab=" ",ylab=" ")
  text(0.5,0.5,la,srt=sr,cex=2)
}

pl2 <- function(y) {
  plot(x,y,type="l",axes=FALSE,xlab=" ",ylab=" ")
  axis(side=1,at=(0:10)/10,cex=2)
}

dens.uniform <- function(tau) {dunif(tau,0.00,1.00)}
dens.beta12 <- function(tau) {dbeta(tau,1,2)}
dens.beta19 <- function(tau) {dbeta(tau,1,9)}
dens.beta1T <- function(tau) {dbeta(tau,1,10)}
dens.beta21 <- function(tau) {dbeta(tau,2,1)}
dens.beta91 <- function(tau) {dbeta(tau,9,1)}
dens.betaT1 <- function(tau) {dbeta(tau,10,1)}

calc_likelihood <- function(pi, x=54, n=60, g=33) {
  v <- dbinom(x, n, pi)
  matrix(v, nrow=g, byrow=TRUE)
}

draw_likelihood <- function(x=54, n=60,
                            g=33, xmax=1) {
  pi <- cycle_fast(g)
  d <- calc_likelihood(pi, x, n, g)
  persp(z=d, theta=th, phi=ph,
        xlab="tau", ylab="pi", zlab="Likelihood",
        zlim=c(0, max(d)))
  return(d)
}

calc_hedge_trace <- function(g=99, alpha=4, beta=16, d=1, hp=dens.uniform) {
  tau <- cycle_slow(g)
  pi <- cycle_fast(g)
  v <- dbeta(pi, d+(alpha-d)*tau, d+(beta-d)*tau)*hp(tau)
  data.frame(pi=pi, tau=tau, v=v) %>% 
    mutate(prod=tau*v) %>%
    group_by(pi) %>%
    summarize(num=sum(prod, na.rm=TRUE), den=sum(v, na.rm=TRUE)) %>%
    ungroup %>%
    mutate(mn=num/den) %>%
    return
}


plot_hedge_trace <- function(g=99, alpha=4, beta=16, d=1, hp=dens.uniform) {
  hedge_trace <- calc_hedge_trace(g, alpha, beta, d, hp)
  plot(hedge_trace$pi, hedge_trace$mn, type="l", ylim=c(0,1),
       xlab="pi", ylab="Expected value of tau")
  return(hedge_trace)
}