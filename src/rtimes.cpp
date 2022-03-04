#include "rts.h"
// authors: Christoph Klauer
//#include "norminc.h"

// namespace rtsNS {

	#define NNODES(I,J) nnodes[I*kernpar+J]
	#define NPPR(T,R) nppr[(T)*respno + R]
	#define FACTOR(T,R) factor[T*respno+R]

	//restpars: igroup*respno; 1 = rsig (vorher + 2); respno*(respno+1)/2 = tau; indi*respno= alphas; indi = rsigalphas;
	//n_all_parameters n_value_store


	double malpha(int t, int r, double* restpars, double* slams) {
		return slams[r] * restpars[alphaoff + t * respno + r];
	}


	void make_rtau(double *restpars, double *taui, double* slams, gsl_rng *rst) {
	#define XY(I,J) xy[(I)*(respno) + J]


		double *xy = 0;	xy = (double *)R_Calloc((indi + respno + 1 + pr_df_add_inv_wish)*(respno), double);
		double *tau = 0;	tau = (double *)R_Calloc(respno*(respno), double);

		for (int i = 0; i != indi; i++) {
			for (int j = 0; j != respno; j++) XY(i, j) = restpars[alphaoff + i * respno + j];
		}
		invwis(indi, respno, xy, tau, taui, pr_sf_scale_matrix_TAU, rst);
		int iz = -1;
		for (int i = 0; i != respno; i++) for (int j = i; j != respno; j++) restpars[igroup*respno + 1 + (++iz)] = tau[i*respno + j];
		if (xy) R_Free(xy);
		if (tau) R_Free(tau);
	}

	void make_rsigalpha(vector<trial> daten, double* factor, double *rest, double *restpar, double* slams, bool xflag, gsl_rng *rst) {

		double *u = 0; u = (double *)R_Calloc(indi, double);
		double *fn = 0; fn = (double *)R_Calloc(respno, double);
		double *n_per_person = 0; n_per_person = (double *)R_Calloc(indi, double);

		double sigsquar = restpar[1 + respno * igroup - 1];
		int no_trials = static_cast<int>(daten.size());
		for (int it = 0; it != indi; it++) {
			u[it] = pr_df_sigma_sqr*sigsquar; n_per_person[it] = pr_df_sigma_sqr;
			for (int r = 0; r != respno; r++) n_per_person[it] += (double) NPPR(it, r);
		}

		for (int i = 0; i != no_trials; i++) {
			int t = daten[i].person; int r = cat2resp[daten[i].category];
			u[t] += gsl_pow_2(rest[i] - malpha(t, r, restpar, slams) - restpar[t2group[t] * respno + r]);
		}

		for (int t = 0; t != indi; t++) {
			int ihelp = 0;
		STEP0:		double factornew = 0, factorold = 0, store = restpar[sigalphaoff + t];
			for (int r = 0; r != respno; r++) {
				factorold += FACTOR(t, r);
			}

			double x[1];

			x[0] = gsl_ran_chisq(rst, n_per_person[t]);

			restpar[sigalphaoff + t] = 1 / x[0] * u[t];

			if ((xflag) && (ihelp == 0)) {// std::cout<< "Step 0";
				ihelp++;
					for (int r = 0; r != respno; r++) {
						double mu = malpha(t, r, restpar, slams) + restpar[t2group[t] * respno + r];
						double rsig = sqrt(restpar[sigalphaoff + t]);
						FACTOR(t, r) = lnnorm(mu / rsig)*NPPR(t, r);
					}
				goto STEP0;
			}

			for (int r = 0; r != respno; r++) {
				double mu = malpha(t, r, restpar, slams) + restpar[t2group[t] * respno + r];
				double rsig = sqrt(restpar[sigalphaoff + t]);
				fn[r] = lnnorm(mu / rsig)*NPPR(t, r);
				factornew += fn[r];
			}

			double ux = oneuni(rst);
			if (log(ux) > factorold - factornew) {
				restpar[sigalphaoff + t] = store; //std::cout << "sigalph";
			}
			else for (int r = 0; r != respno; r++) FACTOR(t, r) = fn[r];


		}
		if (u) R_Free(u);
		if (fn) R_Free(fn);
		if (n_per_person) R_Free(n_per_person);
	}



	void make_ralpha(vector<trial> daten, double* factor, double *rest, double *restpars, double* slams, double *taui, gsl_rng *rst) {



		double *w = 0;	w = (double *)R_Calloc(respno, double);
		double *hba = 0;	hba = (double *)R_Calloc(respno, double);
		double *fig = 0;	fig = (double *)R_Calloc(indi*respno, double);
		double *xfig = 0;	xfig = (double *)R_Calloc(respno*respno, double);
		double *ba = 0;	ba = (double *)R_Calloc(indi*respno, double);
		double *fn = 0; fn = (double *)R_Calloc(respno, double);
	#define FIG(T,I) fig[T*respno+I]
	#define XFIG(I,J) xfig[I*respno+J]
	#define BA(T,I) ba[T*respno+I]
	#define TAUI(I,J) taui[I*respno+J]

		// int no_trials = int(daten.size());

		for (int t = 0; t != indi; t++) for (int iz = 0; iz != respno; iz++) { BA(t, iz) = 0.0; FIG(t, iz) = 0.0; }



		for (int i = 0; i != static_cast<int>(daten.size()); i++) {
			trial one = daten[i]; int t = one.person, r = cat2resp[one.category];
			double be = restpars[t2group[t] * respno + r];
			BA(t, r) += (rest[i] - be);
		}
		for (int r = 0; r != respno; r++) {
			double slam = slams[r];
			for (int t = 0; t != indi; t++) {
				FIG(t, r) = NPPR(t, r)*gsl_pow_2(slam);
				BA(t, r) *= slam;
			}
		}

		for (int iz = 0; iz != respno; iz++) for (int jz = 0; jz != respno; jz++) if (iz != jz) XFIG(iz, jz) = TAUI(iz, jz);

		for (int t = 0; t != indi; t++) {

			double factornew = 0.0, factorold = 0;
			double *store = 0;	store = (double *)R_Calloc(respno, double);
			for (int ir = 0; ir != respno; ir++) store[ir] = restpars[alphaoff + t * respno + ir];

			for (int r = 0; r != respno; r++) {
				factorold += FACTOR(t, r);
			}

			for (int iz = 0; iz != respno; iz++) {
				w[iz] = BA(t, iz) / restpars[sigalphaoff + t]; XFIG(iz, iz) = FIG(t, iz) / restpars[sigalphaoff + t] + TAUI(iz, iz);
			}
			bayesreg(respno, w, xfig, hba, rst);
			for (int ip = 0; ip != respno; ip++) restpars[alphaoff + t * respno + ip] = hba[ip];

			for (int r = 0; r != respno; r++) {
				double mu = malpha(t, r, restpars, slams) + restpars[t2group[t] * respno + r];
				double rsig = sqrt(restpars[sigalphaoff + t]);
				fn[r] = lnnorm(mu / rsig)*NPPR(t, r);
				factornew += fn[r];
			}

			double u = oneuni(rst);
			if (log(u) > factorold - factornew)
			{
				for (int ir = 0; ir != respno; ir++) restpars[alphaoff + t * respno + ir] = store[ir]; // std::cout << " alpha";
			}
			else for (int ir = 0; ir != respno; ir++) FACTOR(t, ir) = fn[ir];
			if (store) R_Free(store);

		}


		if (w) R_Free(w);
		if (hba) R_Free(hba);
		if (fig) R_Free(fig);
		if (xfig) R_Free(xfig);
		if (ba) R_Free(ba);
		if (fn) R_Free(fn);

	}


	void make_slams(vector<trial> daten, double* factor, double *rest, double *restpars, double *slams, gsl_rng *rst) {

		double *fig = 0;	fig = (double *)R_Calloc(indi*respno, double);
		double *ba = 0;	ba = (double *)R_Calloc(indi*respno, double);
		double *fn = 0; fn = (double *)R_Calloc(indi, double);

		for (int t = 0; t != indi; t++) for (int iz = 0; iz != respno; iz++) { BA(t, iz) = 0.0; FIG(t, iz) = 0.0; }

		for (int i = 0; i != static_cast<int>(daten.size()); i++) {
			trial one = daten[i]; int t = one.person, r = cat2resp[one.category];
			double be = restpars[t2group[t] * respno + r];
			BA(t, r) += (rest[i] - be);
		}
		for (int r = 0; r != respno; r++) {
			for (int t = 0; t != indi; t++) {
				double alpha = restpars[alphaoff + t * respno + r];
				FIG(t, r) = NPPR(t, r)*gsl_pow_2(alpha) / restpars[sigalphaoff + t];
				BA(t, r) *= alpha / restpars[sigalphaoff + t];
			}


			double factornew = 0.0, factorold = 0.0, store = slams[r];
			double u = PRIOR; double w = 0;
			for (int t = 0; t != indi; t++) {
				w += BA(t, r);
				u += FIG(t, r);

				factorold += FACTOR(t, r);
			}

			if (u <= 0) u = DBL_MIN;
			slams[r] = (PRIOR + w) / u + onenorm(rst) / sqrt(u);

			for (int t = 0; t != indi; t++) {
				double mu = malpha(t, r, restpars, slams) + restpars[t2group[t] * respno + r];
				double rsig = sqrt(restpars[sigalphaoff + t]);
				fn[t] = lnnorm(mu / rsig)*NPPR(t, r);
				factornew += fn[t];
			}
			double ux = oneuni(rst);
			if (log(ux) > factorold - factornew) {
				slams[r] = store;
			}
			else for (int t = 0; t != indi; t++) FACTOR(t, r) = fn[t];
		}
		if (fig) R_Free(fig);
		if (ba) R_Free(ba);
		if (fn) R_Free(fn);
	}







	void make_rmu(vector<trial> daten, double* factor, double *rest, double *restpar, double *slams, gsl_rng *rst) {

		int no_trials = static_cast<int>(daten.size()); double sig_prior = pr_var_mu_gamma;//1.0 / 0.1;
		double *u = 0; u = (double *)R_Calloc(igroup*respno, double);
		double *rsig = 0; rsig = (double *)R_Calloc(igroup*respno, double);
		double *spostg = 0; spostg = (double *)R_Calloc(igroup*respno, double);
		double *factorold = 0; factorold = (double *)R_Calloc(igroup*respno, double);
		double *fn = 0; fn = (double *)R_Calloc(indi * respno, double);

		double *store = 0; store = (double *)R_Calloc(igroup*respno, double);
		double *factornew = 0; factornew = (double *)R_Calloc(igroup*respno, double);
		bool *keepold = 0; keepold = (bool *)R_Calloc(igroup*respno, bool);

		for (int ig = 0; ig != igroup * respno; ig++) { u[ig] = pr_mean_mu_gamma/sig_prior; rsig[ig] = 0.0; spostg[ig] = 1.0 / sig_prior; }



		for (int i = 0; i != no_trials; i++) {
			int t = daten[i].person; int r = cat2resp[daten[i].category]; int ig = t2group[t];
			u[ig*respno + r] += (rest[i] - malpha(t, r, restpar, slams)) / restpar[sigalphaoff + t];
		}
		for (int it = 0; it != indi; it++) for (int ir = 0; ir != respno; ir++)
			spostg[t2group[it] * respno + ir] += NPPR(it, ir)*1.0 / restpar[sigalphaoff + it];


		for (int ig = 0; ig != igroup; ig++) for (int r = 0; r != respno; r++)
			factorold[ig*respno + r] = 0.0;
		for (int t = 0; t != indi; t++) for (int r = 0; r != respno; r++)
			factorold[t2group[t] * respno + r] += FACTOR(t, r);

		for (int ig = 0; ig != igroup; ig++)
			for (int ir = 0; ir != respno; ir++) {

				factornew[ig*respno + ir] = 0.0;
				store[ig*respno + ir] = restpar[ig*respno + ir];



				u[ig*respno + ir] = u[ig*respno + ir] / (spostg[ig*respno + ir]);
				rsig[ig*respno + ir] = sqrt(1 / spostg[ig*respno + ir]);

				restpar[ig*respno + ir] = u[ig*respno + ir] + onenorm(rst)*rsig[ig*respno + ir];
			}
		//proposals generated; now Metropolis-Hastings

		for (int t = 0; t != indi; t++)
			for (int ir = 0; ir != respno; ir++) {
				int ig = t2group[t];
				double mu = malpha(t, ir, restpar, slams) + restpar[ig * respno + ir];
				double xsig = sqrt(restpar[sigalphaoff + t]);
				fn[t*respno+ir] = lnnorm(mu / xsig)*NPPR(t, ir);
				factornew[ig*respno + ir] += fn[t*respno+ir];
			}
		for (int ig = 0; ig != igroup; ig++)
			for (int ir = 0; ir != respno; ir++) {
				double ux = oneuni(rst);
				if (log(ux) > factorold[ig*respno + ir] - factornew[ig*respno + ir]) {
					restpar[ig*respno + ir] = store[ig*respno + ir];// std::cout << " rmu";
					keepold[ig*respno + ir] = true;
				}
				else keepold[ig*respno + ir] = false;
			}
		for (int t = 0; t != indi; t++) for (int ir = 0; ir != respno; ir++)
			if (!(keepold[t2group[t] * respno + ir])) FACTOR(t, ir) = fn[t*respno+ir];


		if (u) R_Free(u);
		if (rsig) R_Free(rsig);
		if (spostg) R_Free(spostg);
		if (factorold) R_Free(factorold);
		if (fn) R_Free(fn);
		R_Free(store);
		R_Free(factornew);
		R_Free(keepold);
	}

	void make_rsig(vector<trial> daten, double *rest, double *restpar, gsl_rng *rst) {

		// double prior_gamma = 0.0;
		double s = 0.0;
		for (int t = 0; t != indi; t++) {
			s += 1 / restpar[t + sigalphaoff];
		}
		s = s * pr_df_sigma_sqr; //	s=s*adf;
		double x[1];

		{
			double alpha = (indi*pr_df_sigma_sqr) / 2.0 + pr_shape_omega_sqr, beta = s / 2.0 + pr_rate_omega_sqr;

			x[0] = gsl_ran_gamma(rst, alpha, 1.0 / beta);
			restpar[1 + igroup * respno - 1] = x[0];
		}
	}

// }
