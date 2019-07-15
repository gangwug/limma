/*
normexp fitting
Jeremy Silver, Dec 2007
Minor modifications by Gordon Smyth, Sep 2008
Minor change to memory allocation by Jeremy Silver, Sep 2008
*/

#include <stdio.h>
#include <stdlib.h>
#include <Rmath.h>
#include <R.h>
#include <R_ext/Applic.h>

double *x;
int *n;

void ex(){}

double normexp_m2loglik_saddle(int m, double *par, void *ex){
// normexp minus-twice log-likelihood
// Function of mu, log(sigma) and log(alpha)

  extern int *n;
  extern double *x;

  double mu;
  mu = par[0];
  double sigma, sigma2;
  sigma = exp(par[1]);
  sigma2 = sigma*sigma;
  double alpha;
  alpha = exp(par[2]);

  double upperbound1;
  double upperbound2;
  double *upperbound;
  double *theta;
  double k1,k2,k3,k4;
  double err;
  double c0;
  double c1;
  double c2;
  double logf;
  double omat, omat2;
  double loglik;
  double thetaQuadratic;
  int keepRepeating = 1;
  int j,i;
  int *hasConverged; // vector of 0/1 indicators,
  // hasConverged[i] = 0 if theta[i] has not yet converged,
  // hasConverged[i] = 1 if theta[i] has converged,
  int nConverged = 0;   // the sum of hasConverged
  double alpha2 = alpha * alpha;
  double alpha3 = alpha * alpha2;
  double alpha4 = alpha2 * alpha2;
  double dK, ddK, delta;

  upperbound = (double *) Calloc(*n, double);
  theta = (double *) Calloc(*n, double);
  hasConverged = (int *) Calloc(*n, int);

  c2 = sigma2 * alpha;

  for(i = 0; i < *n; i++){

    err = x[i] - mu;
    //	Sigma small approximation
    upperbound1 = fmax(0.0, ((err - alpha)/( alpha * fabs(err))));
    //	alpha small approximation
    upperbound2 = err/sigma2;
    upperbound[i] = fmin(upperbound1, upperbound2);
    c1 = -sigma2 - err * alpha;
    c0 = -alpha + err;
    //	Solve quadratic approximation
    //	Theoretically exact, but subject to subtractive cancellation
    thetaQuadratic = (-c1 - sqrt(c1*c1 - 4 * c0 * c2) ) / (2*c2);
    theta[i] = fmin(thetaQuadratic,upperbound[i]);
    hasConverged[i] = 0;
  }


  //	Globally convergence Newton iteration
  j = 0;
  while(keepRepeating == 1){
    j++;
    for(i = 0; i < *n; i++){
      // Only loop over entries of theta[.] that haven't yet converged
      if(hasConverged[i] == 0){
	omat = 1 - alpha * theta[i];
        dK = mu + sigma2 * theta[i] + alpha / omat;
        ddK = sigma2 + alpha2 / (omat*omat);
        delta = (x[i] - dK)/ddK;
        theta[i] += delta;
        if (j == 1){
          theta[i] = fmin(theta[i], upperbound[i]);
        }
        if(fabs(delta) < 1e-10){
          hasConverged[i] = 1;
          nConverged++;
        }
      }
    }

    if(nConverged == *n || j > 50){
      keepRepeating = 0;
    }
  }

  R_CheckUserInterrupt();

  loglik = 0.0;
  for(i = 0; i < *n; i++){
    omat = 1 - alpha * theta[i];
    omat2 = omat*omat;
    k1 = mu * theta[i] + 0.5 * sigma2 * theta[i] * theta[i] - log(omat);
    k2 = sigma2 + alpha2/omat2;
    logf = -0.5 * log(2.0 * M_PI * k2) - x[i] * theta[i] + k1;
    k3 = 2.0 * alpha3/(omat * omat2);
    k4 = 6.0 * alpha4/(omat2 * omat2);
    logf += k4/(8.0 * k2 * k2) - (5.0 * k3 * k3)/(24.0 * k2 * k2 * k2);
    loglik += logf;
  }

  Free(upperbound);
  Free(theta);
  Free(hasConverged);

  return -2.0 * loglik;

}

/*
    k <- mu * theta + 0.5 * sigma^2 * theta^2 - log(1 - alpha * theta)
    k2 <- sigma^2 + alpha^2/(1 - alpha * theta)^2
    logf <- -0.5 * log(2 * pi * k2) - x * theta + k
    if (secondorder) {
        k3 <- 2 * alpha^3/(1 - alpha * theta)^3
        k4 <- 6 * alpha^4/(1 - alpha * theta)^4
        logf <- logf + 1/8 * k4/k2^2 - 5/24 * k3^2/k2^3
    }
*/

void fit_saddle_nelder_mead(double *par, double *X, int *N, int *fail, int *fncount, double *Fmin){
// Minimize normexp m2loglik by Nelder-Mead
// as a function of mu, log(sigma) and log(alpha)

  double parsOut[3];
  parsOut[0] = par[0];
  parsOut[1] = par[1];
  parsOut[2] = par[2];
  double abstol = -1e308; // infinity
  double intol = 1.490116e-08; // square root of machine precision
  void ex();
  double alpha = 1.0;
  double beta = 0.5;
  double gamma = 2.0;
  int trace = 0;
  int maxit = 500;

  extern int *n;
  extern double *x;
  n = N;
  x = X;

  nmmin(3, par, parsOut, Fmin, normexp_m2loglik_saddle, fail, abstol, intol, &ex, alpha, beta, gamma, trace, fncount, maxit);

  par[0] = parsOut[0];
  par[1] = parsOut[1];
  par[2] = parsOut[2];

}

void normexp_m2loglik(double *mu, double *s2, double *al, int *n, double *f, double *m2LL){
// normexp minus-twice log-likelihood
// as a function of mu, sigma^2 and alpha

  double e;
  double mu_sf;
  double s2onal = *s2/ *al;
  double logal = log(*al);
  double s2on2al2 = 0.5 * *s2/(*al * *al);
  double s = sqrt(*s2);
  int i;
  *m2LL = 0.0;

  for(i = 0; i < *n; i++){
    e = f[i] - *mu;
    mu_sf = e - s2onal;
    *m2LL += -logal - e/ *al + s2on2al2 + pnorm(0.0,mu_sf,s,0,1);
  }
  //  -2*sum(-log(al) - e/al + 0.5*s2/(al^2) + pnorm(0,mu.sf,sqrt(s2),lower.tail = FALSE,log.p = TRUE))

  *m2LL *= -2.0;

}

void normexp_gm2loglik(double *mu, double *s2, double *al, int *n, double *f, double *dm2LL){
// gradient of normexp m2loglik
// with respect to mu, log(sigma^2) and log(alpha)

  double e;
  double mu_sf;
  double s2onal = *s2/ *al;
  double s = sqrt(*s2);
  double v1onal = 1/ *al;
  double psionPsi;
  double al2 = *al * *al;
  double s2onal3 = *s2/(al2 * *al);
  double v1on2al2 = 0.5/al2;
  double s2onal2 = *s2/al2;
  double v1on2s2 = 0.5/ *s2;
  int i;

  for(i = 0; i < 3; i++){
    dm2LL[i] = 0.0;
  }
  // Calculate derivatives
  for(i = 0; i < *n; i++){
    e = f[i] - *mu;
    mu_sf = e - s2onal;
    psionPsi = exp(dnorm(0.0,mu_sf,s,1) - pnorm(0.0,mu_sf,s,0,1));
    dm2LL[0] += v1onal - psionPsi;
    dm2LL[1] += v1on2al2 - (v1onal + v1on2s2 * mu_sf) * psionPsi;
    dm2LL[2] += e/al2 - v1onal - s2onal3 + psionPsi * s2onal2;
  }
  // correct for taking derivatives wrt -2 times the log(likelihood)
  for(i = 0; i < 3; i++){
    dm2LL[i] *= -2.0;
  }
  // correct for differentiation wrt log(alpha) and log(sigma^2)
  dm2LL[1] *= *s2;
  dm2LL[2] *= *al;

}

void normexp_hm2loglik(double *mu, double *s2, double *al, int *n, double *f, double *d2m2LL){
// Hessian of normexp m2loglik
// with respect to mu, log(sigma^2) and log(alpha)

  double e;
  double mu_sf;
  double s2onal = *s2/ *al;
  double s2onalsq = s2onal*s2onal;
  double s2onalcu = s2onalsq*s2onal;
  double s = sqrt(*s2);
  double v1onal = 1/ *al;
  double v1onal2 = v1onal*v1onal;
  double v1onal3 = v1onal2*v1onal;
  double v1onal4 = v1onal3*v1onal;
  double v3s2onal4 = 3.0 * *s2*v1onal4;
  double psionPsi;
  double psionPsi2;
  double al2 = *al * *al;
  double s2onal3 = *s2/(al2 * *al);
  double v1on2al2 = 0.5/al2;
  double s2onal2 = *s2/al2;
  double s4onal4 = s2onal2*s2onal2;
  double s2onal4 = *s2*v1onal4;
  double v1on2s2 = 0.5/ *s2;
  double v1on2s4 = v1on2s2/ *s2;
  double v1on2s2sq = v1on2s2*v1on2s2;
  double v2onal3 = 2.0*v1onal3;
  double v3al = 3.0 * *al;
  double v2al = 2.0 * *al;
  double v1on4s6 = v1on2s2sq/ *s2;
  double eps2onal;
  double e2;
  double eps2onalsq;
  int i;
  double dL_dal = 0.0;
  double dL_ds2 = 0.0;
  double d2L_dbtdbt = 0.0;
  double d2L_dbtds2 = 0.0;
  double d2L_dbtdal = 0.0;
  double d2L_dalds2 = 0.0;
  double d2L_ds2ds2 = 0.0;
  double d2L_daldal = 0.0;

  for(i = 0; i < *n; i++){
    e = f[i] - *mu;
    mu_sf = e - s2onal;
    eps2onal = e + s2onal;
    e2 = e*e;
    eps2onalsq = eps2onal*eps2onal;
    psionPsi = dnorm(0.0,mu_sf,s,1) - pnorm(0.0,mu_sf,s,0,1);
    psionPsi2 = 2.0 * psionPsi;
    psionPsi = exp(psionPsi);
    psionPsi2 = exp(psionPsi2);

    dL_dal += v1on2al2 - (v1onal + v1on2s2 * mu_sf) * psionPsi;
    dL_ds2 += e/al2 - v1onal - s2onal3 + psionPsi * s2onal2;
    d2L_dbtdbt += - psionPsi2 - psionPsi*mu_sf/ *s2;
    d2L_dbtds2 +=  -0.5*eps2onal*psionPsi2/ *s2 + (-eps2onalsq + 2.0*s2onal*eps2onal + *s2)*psionPsi*v1on2s4;
    d2L_dbtdal +=  -v1onal2 + s2onal2*psionPsi2 + mu_sf*psionPsi*v1onal2;
    d2L_ds2ds2 +=  -v1on2s2sq*eps2onalsq*psionPsi2 + psionPsi*(-e2*e + e*(v3al - e)*s2onal + (e + *al)*s2onalsq + s2onalcu)*v1on4s6;
    d2L_dalds2 +=  -v1onal3 + v1on2al2*(psionPsi2*eps2onal + (e2 + *s2 - s2onalsq)*psionPsi/ *s2);
    d2L_daldal +=  v1onal2 - v2onal3*e + v3s2onal4 - psionPsi2*s4onal4 - psionPsi*(mu_sf + v2al)*s2onal4;
  }

  /*
    d2L.dbtdbt <- sum(-psionPsi2 - psionPsi*mu.sf/s2)
    d2L.dbtds2 <- sum( -0.5*(e + s2onal)*psionPsi2/s2 + 0.5*(-((e + s2onal)^2) + 2*s2onal*(e + s2onal) + s2)*psionPsi/(s2^2)) # OK TO 3 DP
    d2L.dbtdal <- sum( -al^-2 + s2onal*psionPsi2/al + mu.sf*psionPsi/(al^2))
    d2L.ds2ds2 <- sum( -(0.25/(s2^2))*((e + s2onal)^2)*psionPsi2 + psionPsi*(-e^3 + e*(3*al - e)*s2onal + (e + al)*(s2onal^2) + (s2onal^3))/(4*(s2^3)) )
    d2L.dalds2 <- sum( -1/(al^3) + (al^-2)*0.5*(psionPsi2*(e + s2onal) + (e^2 + s2 - (s2onal^2))*psionPsi/s2))
    d2L.daldal <- sum( (al^-2) - 2*e/(al^3) + 3*s2/(al^4) - psionPsi2*((s2^2)/(al^4)) - psionPsi*(mu.sf + 2*al)*((s2)/(al^4)))
  */
  d2m2LL[0] = -2.0 * d2L_dbtdbt;
  d2m2LL[1] = -2.0 * *s2 * d2L_dbtds2;
  d2m2LL[2] = -2.0 * *al * d2L_dbtdal;
  d2m2LL[3] = -2.0 * *s2 * d2L_dbtds2;
  d2m2LL[4] = -2.0 * ( *s2 * *s2 * d2L_ds2ds2 + *s2 * dL_ds2);
  d2m2LL[5] = -2.0 * *al * *s2 * d2L_dalds2;
  d2m2LL[6] = -2.0 * *al * d2L_dbtdal;
  d2m2LL[7] = -2.0 * *al* *s2 * d2L_dalds2;
  d2m2LL[8] = -2.0 * ( *al * *al * d2L_daldal + *al * dL_dal);

  /*
    -2*rbind(
    c(d2L.dbtdbt,    s2*d2L.dbtds2,                 al*d2L.dbtdal),
    c(s2*d2L.dbtds2, (s2^2)*d2L.ds2ds2 + s2*dL.ds2, al*s2*d2L.dalds2),
    c(al*d2L.dbtdal, al*s2*d2L.dalds2,              (al^2)*d2L.daldal + al*dL.dal))
  */

}

