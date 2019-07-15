#include "R.h"
#include "Rdefines.h"
#include <math.h>

#define THRESHOLD 0.0000001

/* This function determines the number of points to be used in the analysis,
 * based on the chosen delta. It returns that number as well as the index
 * values of those points in a pointer reference. Note that the first
 * and last point (i.e. smallest and largest) are always considered.
 */

void find_seeds (int ** indices, int * number, const double* xptr, const int npts, const double delta) {
	int pt, last_pt=0;
	int total=2;
	for (pt=1; pt<npts-1; ++pt) {
		if (xptr[pt] - xptr[last_pt] > delta) {
			++total;
			last_pt=pt;
		}
	}
	(*number)=total;

	/* Second pass to actually record the indices at which these events happen. */
	int* idptr=(int*)R_alloc(*number, sizeof(int));
	idptr[0]=0;
	total=1;
	last_pt=0;
	for (pt=1; pt<npts-1; ++pt) {
        if (xptr[pt] - xptr[last_pt] > delta) {
			idptr[total]=pt;
			last_pt=pt;
			++total;
		}
	}
	idptr[total]=npts-1;
	++total;
	(*indices)=idptr;
	return;
}

/* This function identifies the start and end index in the span for each chosen sampling
 * point. It returns two arrays via reference containing said indices. It also returns
 * an array containing the maximum distance between points at each span.
 *
 * We don't use the update-based algorithm in Cleveland's paper, as it ceases to be
 * numerically stable once you throw in double-precision weights. It's not particularly
 * amenable to updating through cycles of addition and subtraction. At any rate, the
 * algorithm as a whole remains quadratic (as weights must be recomputed) so there's no
 * damage to scalability.
 */

void find_limits (const int* indices, const int num, const double* xptr, const double* wptr,
		const int npts, const double spanweight, int** start, int** end, double** dist) {
	int* spbegin=(int*)R_alloc(num, sizeof(int));
	int* spend=(int*)R_alloc(num, sizeof(int));
	double* spdist=(double*)R_alloc(num, sizeof(double));

	int curx;
	for (curx=0; curx<num; ++curx) {
		const int curpt=indices[curx];
		int left=curpt, right=curpt;
		double curw=wptr[curpt];
		int ende=(curpt==npts-1), ends=(curpt==0);
		double mdist=0, ldist, rdist;

		while (curw < spanweight && (!ende || !ends)) {
			if (ende) {
				/* Can only extend backwards. */
				--left;
				curw+=wptr[left];
				if (left==0) { ends=1; }
				ldist=xptr[curpt]-xptr[left];
				if (mdist < ldist) { mdist=ldist; }
			} else if (ends) {
				/* Can only extend forwards. */
				++right;
				curw+=wptr[right];
				if (right==npts-1) { ende=1; }
				rdist=xptr[right]-xptr[curpt];
				if (mdist < rdist) { mdist=rdist; }
			} else {
				/* Can do either; extending by the one that minimizes the curpt mdist. */
				ldist=xptr[curpt]-xptr[left-1];
				rdist=xptr[right+1]-xptr[curpt];
				if (ldist < rdist) {
					--left;
					curw+=wptr[left];
					if (left==0) { ends=1; }
					if (mdist < ldist) { mdist=ldist; }
				} else {
					++right;
					curw+=wptr[right];
					if (right==npts-1) { ende=1; }
					if (mdist < rdist) { mdist=rdist; }
				}
			}
		}

		/* Extending to ties. */
		while (left>0 && xptr[left]==xptr[left-1]) { --left; }
		while (right<npts-1 && xptr[right]==xptr[right+1]) { ++right; }

		/* Recording */
		spbegin[curx]=left;
		spend[curx]=right;
		spdist[curx]=mdist;
	}

	(*start)=spbegin;
	(*end)=spend;
	(*dist)=spdist;
	return;
}

/* Computes the lowess fit at a given point using linear regression with a combination of tricube,
 * prior and robustness weighting. Some additional effort is put in to avoid numerical instability
 * and undefined values when divisors are near zero.
 */

double lowess_fit (const double* xptr, const double* yptr, const double* wptr, const double* rwptr,
		const int npts, const int curpt, const int left, const int right, const double dist, double* work) {
	double ymean=0, allweight=0;
	int pt;
	if (dist < THRESHOLD) {
		for (pt=left; pt<=right; ++pt) {
			work[pt]=wptr[pt]*rwptr[pt];
			ymean+=yptr[pt]*work[pt];
			allweight+=work[pt];
		}
		ymean/=allweight;
		return ymean;
	}
	double xmean=0;
	for (pt=left; pt<=right; ++pt) {
		work[pt]=pow(1-pow(fabs(xptr[curpt]-xptr[pt])/dist, 3.0), 3.0)*wptr[pt]*rwptr[pt];
		xmean+=work[pt]*xptr[pt];
		ymean+=work[pt]*yptr[pt];
		allweight+=work[pt];
	}
	xmean/=allweight;
	ymean/=allweight;

	double var=0, covar=0, temp;
	for (pt=left; pt<=right; ++pt) {
		temp=xptr[pt]-xmean;
		var+=temp*temp*work[pt];
		covar+=temp*(yptr[pt]-ymean)*work[pt];
	}
	if (var < THRESHOLD) { return ymean; }

	const double slope=covar/var;
	const double intercept=ymean-slope*xmean;
	return slope*xptr[curpt]+intercept;
}

/* This is a C version of the local weighted regression (lowess) trend fitting algorithm,
 * based on the Fortran code in lowess.f from http://www.netlib.org/go written by Cleveland.
 * Consideration of non-equal prior weights is added to the span calculations and linear
 * regression. These weights are intended to have the equivalent effect of frequency weights
 * (at least, in the integer case; extended by analogy to all non-negative values).
 */

SEXP weighted_lowess(SEXP covariate, SEXP response, SEXP weight, SEXP span, SEXP iter, SEXP delta) {
    if (!IS_NUMERIC(covariate)) { error("covariates must be double precision"); }
    if (!IS_NUMERIC(response)) { error("responses must be double precision"); }
    if (!IS_NUMERIC(weight)) { error("weights must be double precision"); }

	const int npts=LENGTH(covariate);
	if (npts!=LENGTH(response) || npts!=LENGTH(weight)) { error("weight, covariate and response vectors have unequal lengths"); }
	if (npts<2) { error("need at least two points"); }
	const double* covptr=NUMERIC_POINTER(covariate);
	const double* resptr=NUMERIC_POINTER(response);
	const double* weiptr=NUMERIC_POINTER(weight);

	if (!IS_NUMERIC(span) || LENGTH(span)!=1) { error("span should be a double-precision scalar"); }
	const double spv=NUMERIC_VALUE(span);
	if (!IS_INTEGER(iter) || LENGTH(iter)!=1) { error("number of robustness iterations should be an integer scalar"); }
	const int niter=INTEGER_VALUE(iter);
	if (niter<=0) { error("number of robustness iterations should be positive"); }
	if (!IS_NUMERIC(delta) || LENGTH(delta)!=1) { error("delta should be a double-precision scalar"); }
	const double dv=NUMERIC_VALUE(delta);

	/*** NO MORE ERRORS AT THIS POINT, MEMORY ASSIGNMENTS ARE ACTIVE. ***/

	/* Computing the span weight that each span must achieve. */
	double totalweight=0;
	int pt;
	for (pt=0; pt<npts; ++pt) { totalweight+=weiptr[pt]; }
	double spanweight=totalweight*spv;
	const double subrange=(covptr[npts-1]-covptr[0])/npts;

	/* Setting up the indices of points for sampling; the frame start and end for those indices, and the max dist. */
	int *seed_index;
	int nseeds;
	find_seeds(&seed_index, &nseeds, covptr, npts, dv);
   	int *frame_start, *frame_end;
	double* max_dist;
	find_limits (seed_index, nseeds, covptr, weiptr, npts, spanweight, &frame_start, &frame_end, &max_dist);

	/* Setting up arrays to hold the fitted values, residuals and robustness weights. */
	SEXP output=PROTECT(NEW_LIST(2));
	SET_VECTOR_ELT(output, 0, NEW_NUMERIC(npts));
	double* fitptr=NUMERIC_POINTER(VECTOR_ELT(output, 0));
	double* rsdptr=(double*)R_alloc(npts, sizeof(double));
	SET_VECTOR_ELT(output, 1, NEW_NUMERIC(npts));
	double* robptr=NUMERIC_POINTER(VECTOR_ELT(output, 1));
	int* rorptr=(int*)R_alloc(npts, sizeof(int));
	for (pt=0; pt<npts; ++pt) { robptr[pt]=1; }

	/* Robustness iterations. */
	int it=0;
	for (it=0; it<niter; ++it) {
		int cur_seed, last_pt=0, subpt;
		double current;

		/* Computing fitted values for seed points, and interpolating to the intervening points. */
		fitptr[0]=lowess_fit(covptr, resptr, weiptr, robptr, npts, 0, frame_start[0], frame_end[0], max_dist[0], rsdptr);
		for (cur_seed=1; cur_seed<nseeds; ++cur_seed) {
			pt=seed_index[cur_seed];
			fitptr[pt]=lowess_fit(covptr, resptr, weiptr, robptr, npts, pt, frame_start[cur_seed],
				frame_end[cur_seed], max_dist[cur_seed], rsdptr); /* using rsdptr as a holding cell. */

			if (pt-last_pt > 1) {
	 			/* Some protection is provided against infinite slopes. This shouldn't be
 				 * a problem for non-zero delta; the only concern is at the final point
 				 * where the covariate distance may be zero. Besides, if delta is not
 				 * positive, pt-last_pt could never be 1 so we'd never reach this point.
 				 */
				current = covptr[pt]-covptr[last_pt];
				if (current > THRESHOLD*subrange) {
					const double slope=(fitptr[pt]-fitptr[last_pt])/current;
					const double intercept=fitptr[pt] - slope*covptr[pt];
					for (subpt=last_pt+1; subpt<pt; ++subpt) { fitptr[subpt]=slope*covptr[subpt]+intercept; }
				} else {
					const double endave=0.5*(fitptr[pt]+fitptr[last_pt]);
					for (subpt=last_pt+1; subpt<pt; ++subpt) { fitptr[subpt]=endave; }
				}
			}
			last_pt=pt;
		}

		/* Computing the weighted MAD of the absolute values of the residuals. */
		double resid_scale=0;
		for (pt=0; pt<npts; ++pt) {
			rsdptr[pt]=fabs(resptr[pt]-fitptr[pt]);
			resid_scale+=rsdptr[pt];
			rorptr[pt]=pt;
		}
		resid_scale/=npts;
		rsort_with_index(rsdptr, rorptr, npts);

		current=0;
		double cmad=0;
		const double halfweight=totalweight/2;
		for (pt=0; pt<npts; ++pt) {
			current+=weiptr[rorptr[pt]];
			if (current==halfweight) {  /* In the unlikely event of an exact match. */
				cmad=3*(rsdptr[pt]+rsdptr[pt+1]);
				break;
			} else if (current>halfweight) {
				cmad=6*rsdptr[pt];
				break;
			}
		}

		/* If it's too small, then robustness weighting will have no further effect.
		 * Any points with large residuals would already be pretty lowly weighted.
		 * This is based on a similar step in lowess.c in the core R code.
		 */
		if (cmad <= THRESHOLD * resid_scale) { break; }

		/* Computing the robustness weights. */
		for (pt=0; pt<npts; ++pt) {
			if (rsdptr[pt]<cmad) {
				robptr[rorptr[pt]]=pow(1-pow(rsdptr[pt]/cmad, 2.0), 2.0);
			} else { robptr[rorptr[pt]]=0; }
		}
	}

	UNPROTECT(1);
	return output;
}
