//implements summary statistics on posterior distribution
//goodness-of-fit tests
//specific contrast
// authors: Christoph Klauer and Raphael Hartmann

#include "rts.h"

namespace ertmpt {
  
  std::ofstream tests_out;
  
  #define SAMPLE(I,IP) sample[(I)*(n_all_parameters+1)+IP]
  #define TREE_AND_NODE2PAR(ITREE,R) tree_and_node2par[ITREE*nodemax+R]
  #define DRIN(J,K,X) drin[J*zweig*nodemax+K*nodemax + X]
  #define NDRIN(J,K) ndrin[J*zweig+K]
  #define AR(I,J,R) ar[I*zweig*nodemax + J*nodemax + R]
  #define LAMBDAS(T,IZ) lambdas[T*ilamfree+IZ] //PM=0 negativ; PM=1 positiv
  #define RHOS(IG,IZ) rhos[(IG)*ilamfree+IZ]
  
  void lies(int n_all_parameters, double *sample) {
  
  	std::ifstream rein(RAUS);
  	int is, in;
  	rein >> is >> in;
  	if (is != SAMPLE_SIZE) Rprintf("HM\n");
  	if (in != (n_all_parameters+1)) Rprintf("HO\n");
  	for (int i = 0; i != is; i++)
  		for (int j = 0; j != in; j++) rein >> SAMPLE(i, j);
  	rein.close();
  }
  
  void hdi(int length, double * parameter, double p, double iv[2]) {
  	int inc = static_cast<int>(p*length) + 1;
  	int nci = length - inc;
  	int imin = -1;
  	double tempold = parameter[length - 1] - parameter[0];
  	for (int i = 0; i != nci; i++) {
  		double temp = parameter[i + inc] - parameter[i];
  		if (temp < tempold) {
  			tempold = temp;
  			imin = i;
  		}
  	}
  	iv[0] = parameter[imin];
  	iv[1] = parameter[imin + inc];
  
  }
  
  #define BETA(T,IP) beta[T*ifree+IP]
  
  void belege_beta(double *sample, int is, double *beta) {
  	for (int t = 0; t != indi; t++) for (int iz = 0; iz != ifree; iz++) BETA(t, iz) = SAMPLE(is, iz + ifree * t2group[t]) +
  		SAMPLE(is, (iz + ilamfree * igroup + ((ifree + ilamfree)*(ifree + ilamfree + 1)) / 2 + (t + igroup)*ifree));
  }
  
  
  void belege_lambdas_rhos(double *sample, int is, double *rhos, double *lambdas) {
  	for (int iz = 0; iz != ilamfree * igroup; iz++) rhos[iz] = SAMPLE(is, ifree*igroup + iz);
  	for (int t = 0; t != indi; t++) for (int iz = 0; iz != ilamfree; iz++)
  		lambdas[t*ilamfree + iz] = SAMPLE(is, ifree*igroup + ilamfree * igroup + ((ifree + ilamfree)*(ifree + ilamfree + 1)) / 2 + indi * ifree + t * ilamfree + iz);
  }
  
  void belege_nur_lambdas(double *sample, int is, double *lambdas) {
  	for (int t = 0; t != indi; t++)
  		for (int iz = 0; iz != ilamfree; iz++)
  			lambdas[t*ilamfree + iz] = SAMPLE(is, ifree*igroup + iz + t2group[t] * ilamfree)*
  			SAMPLE(is, ifree*igroup + ilamfree * igroup + ((ifree + ilamfree)*(ifree + ilamfree + 1)) / 2 + indi * ifree + t * ilamfree + iz);
  }
  
  void quantiles(std::vector<trial> daten, int n_all_parameters, double *sample) {
  	double qv[5];
	std::vector<double> temp(SAMPLE_SIZE);
	//std::streamsize prec = cout.precision(); std::cout << std::setprecision(4);
	if (save_diagnose) tests_out << std::setprecision(4);
	Rprintf("theta per group [median, 95 and 99%% HDI]\n"); if (save_diagnose) tests_out << "MUs per group" << std::endl;
	for (int ig = 0; ig != igroup; ig++)
		for (int ip = 0; ip != kernpar; ip++) if (comp[ip]) if(free2kern[kern2free[ip]] == ip) {
			for (int j = 0; j != SAMPLE_SIZE; j++) temp[j] = gsl_cdf_ugaussian_P(SAMPLE(j, kern2free[ip] + ig * ifree));
  			gsl_sort(temp.data(), 1, SAMPLE_SIZE);
  			qv[2] = gsl_stats_median_from_sorted_data(temp.data(), 1, SAMPLE_SIZE);
  			double iv[2]; hdi(SAMPLE_SIZE, temp.data(), 0.95, iv); qv[1] = iv[0]; qv[3] = iv[1]; hdi(SAMPLE_SIZE, temp.data(), 0.99, iv); qv[0] = iv[0]; qv[4] = iv[1];
  			Rprintf("%3d", ip + 1); for (int iq = 0; iq != 5; iq++) Rprintf("%12.4g", qv[iq]); Rprintf("\n");
  			if (save_diagnose) { tests_out << std::setw(3) << ip + 1; for (int iq = 0; iq != 5; iq++) tests_out << std::setw(12) << qv[iq]; tests_out << std::endl; }
  		}
  
  	Rprintf("tau- in ms per group [median, 95 and 99%% HDI]\n"); if (save_diagnose) tests_out << "rho minus per group" << std::endl;
  	for (int ig = 0; ig != igroup; ig++)
  		for (int ip = 0; ip != kernpar; ip++) if (comp[kernpar + ip]) if(free2kern[kern2free[ip+kernpar]] == ip+kernpar) {
  			for (int j = 0; j != SAMPLE_SIZE; j++) temp[j] = 1000.0 / SAMPLE(j, igroup*ifree + ig * ilamfree + kern2free[kernpar + ip] - ifree);
  			gsl_sort(temp.data(), 1, SAMPLE_SIZE);
  			qv[2] = gsl_stats_median_from_sorted_data(temp.data(), 1, SAMPLE_SIZE);
  			double iv[2]; hdi(SAMPLE_SIZE, temp.data(), 0.95, iv); qv[1] = iv[0]; qv[3] = iv[1]; hdi(SAMPLE_SIZE, temp.data(), 0.99, iv); qv[0] = iv[0]; qv[4] = iv[1];
  			Rprintf("%3d", ip + 1); for (int iq = 0; iq != 5; iq++) Rprintf("%12.4g", qv[iq]); Rprintf("\n");
  			if (save_diagnose) { tests_out << std::setw(3) << ip + 1; for (int iq = 0; iq != 5; iq++) tests_out << std::setw(12) << qv[iq]; tests_out << std::endl; }
  		}
  	Rprintf("tau+ in ms per group [median, 95 and 99%% HDI]\n"); if (save_diagnose) tests_out << "rho plus per group" << std::endl;
  	for (int ig = 0; ig != igroup; ig++)
  		for (int ip = 0; ip != kernpar; ip++) if (comp[ip + 2 * kernpar]) if(free2kern[kern2free[ip+2*kernpar]] == ip+2*kernpar) {
  			for (int j = 0; j != SAMPLE_SIZE; j++) temp[j] = 1000.0 / SAMPLE(j, igroup*ifree + ig * ilamfree + kern2free[2 * kernpar + ip] - ifree);
  			gsl_sort(temp.data(), 1, SAMPLE_SIZE);
  			qv[2] = gsl_stats_median_from_sorted_data(temp.data(), 1, SAMPLE_SIZE);
  			double iv[2]; hdi(SAMPLE_SIZE, temp.data(), 0.95, iv); qv[1] = iv[0]; qv[3] = iv[1]; hdi(SAMPLE_SIZE, temp.data(), 0.99, iv); qv[0] = iv[0]; qv[4] = iv[1];
  			Rprintf("%3d", ip + 1); for (int iq = 0; iq != 5; iq++) Rprintf("%12.4g", qv[iq]); Rprintf("\n");
  			if (save_diagnose) { tests_out << std::setw(3) << ip + 1; for (int iq = 0; iq != 5; iq++) tests_out << std::setw(12) << qv[iq]; tests_out << std::endl; }
  		}
  
  	Rprintf("SD and CORR of process params [median, 95 and 99%% HDI]\n"); if (save_diagnose) tests_out << "SIG" << std::endl;
  	int iz = igroup * (ifree + ilamfree) - 1;
  	for (int ix = 0; ix != ifree + ilamfree; ix++)
  		for (int jz = ix; jz != ifree + ilamfree; jz++) {
  			iz++;
  			for (int j = 0; j != SAMPLE_SIZE; j++) temp[j] = SAMPLE(j, iz);
  			gsl_sort(temp.data(), 1, SAMPLE_SIZE);
  			qv[2] = gsl_stats_median_from_sorted_data(temp.data(), 1, SAMPLE_SIZE);
  			double iv[2]; hdi(SAMPLE_SIZE, temp.data(), 0.95, iv); qv[1] = iv[0]; qv[3] = iv[1]; hdi(SAMPLE_SIZE, temp.data(), 0.99, iv); qv[0] = iv[0]; qv[4] = iv[1];
  			Rprintf("%3d%3d", free2kern[ix]+1, free2kern[jz] + 1); for (int iq = 0; iq != 5; iq++) Rprintf("%12.4g", qv[iq]); Rprintf("\n");
  			if (save_diagnose) { tests_out << std::setw(3) << free2kern[ix] + 1 << std::setw(3) << free2kern[jz] + 1; for (int iq = 0; iq != 5; iq++) tests_out << std::setw(12) << qv[iq]; tests_out << std::endl; }
  		}
  
  	//Rprintf("Restpars: Mean per group, Residual Variance in motor/encoding, Due to individual differences in motor/encoding\n");
  	Rprintf("mu_gamma in ms per group [median, 95 and 99%% HDI]\n");
  	if (save_diagnose) tests_out << "Restpars: Mean per group, Residual Variance in motor/encoding, Due to individual differences in motor/encoding" << std::endl;
  	iz = ifree * igroup + ilamfree * igroup + ((ifree + ilamfree)*(ifree + ilamfree + 1)) / 2 + indi * ifree + indi * ilamfree;
  
  	for (int ir = 0; ir != 1 + igroup * respno; ir++) {
  		if (ir != igroup*respno) for (int j = 0; j != SAMPLE_SIZE; j++) temp[j] = 1000.0*SAMPLE(j, (iz + ir));
  		if (ir == igroup*respno) for (int j = 0; j != SAMPLE_SIZE; j++) temp[j] = SAMPLE(j, (iz + ir));
  		gsl_sort(temp.data(), 1, SAMPLE_SIZE);
  		qv[2] = gsl_stats_median_from_sorted_data(temp.data(), 1, SAMPLE_SIZE);
  		double iv[2]; hdi(SAMPLE_SIZE, temp.data(), 0.95, iv); qv[1] = iv[0]; qv[3] = iv[1]; hdi(SAMPLE_SIZE, temp.data(), 0.99, iv); qv[0] = iv[0]; qv[4] = iv[1];
  		if (ir == igroup*respno) {
  			Rprintf("omega^2 [median, 95 and 99%% HDI]\n");
  		}
  		for (int iq = 0; iq != 5; iq++) Rprintf("%12.4g", qv[iq]);
  		Rprintf("\n");
  		if (save_diagnose) { for (int iq = 0; iq != 5; iq++) tests_out << std::setw(12) << qv[iq]; tests_out << std::endl; }
  	}
  
  	Rprintf("SD and CORR of motor/encoding params [median, 95 and 99%% HDI]\n"); if (save_diagnose) tests_out << "RSIG" << std::endl;
  	iz = n_all_parameters - restparsno + igroup * respno + 1 - 1;
  	for (int ip = 0; ip != respno; ip++)
  		for (int jp = ip; jp != respno; jp++) {
  			iz++;
  			for (int j = 0; j != SAMPLE_SIZE; j++) temp[j] = SAMPLE(j, iz);
  			gsl_sort(temp.data(), 1, SAMPLE_SIZE);
  			qv[2] = gsl_stats_median_from_sorted_data(temp.data(), 1, SAMPLE_SIZE);
  			double iv[2]; hdi(SAMPLE_SIZE, temp.data(), 0.95, iv); qv[1] = iv[0]; qv[3] = iv[1]; hdi(SAMPLE_SIZE, temp.data(), 0.99, iv); qv[0] = iv[0]; qv[4] = iv[1];
  			Rprintf("%3d%3d", ip + 1, jp + 1); for (int iq = 0; iq != 5; iq++) Rprintf("%12.4g", qv[iq]); Rprintf("\n");
  			if (save_diagnose) { tests_out << std::setw(3) << ip + 1 << std::setw(3) << jp + 1; for (int iq = 0; iq != 5; iq++) tests_out << std::setw(12) << qv[iq]; tests_out << std::endl; }
  		}
  
  	// Daten zum Vergleich:
  	double s = 0.0; int no_trials = static_cast<int>(daten.size());
	std::vector<double> u(indi);
	std::vector<int> nj(indi);
  	for (int t = 0; t != indi; t++) { u[t] = 0.0; nj[t] = 0; }
  
  	for (int i = 0; i != no_trials; i++) { u[daten[i].person] += daten[i].rt / 1000.0; nj[daten[i].person]++; }
  	for (int t = 0; t != indi; t++) { u[t] /= nj[t]; }
  
  	for (int i = 0; i != no_trials; i++) s += gsl_pow_2(daten[i].rt / 1000.0 - u[daten[i].person]);
  	s = s / (no_trials - 1);
  	double grand = 0.0; for (int t = 0; t != indi; t++) grand += (u[t] * nj[t]) / no_trials;
  	double salph = 0.0; for (int t = 0; t != indi; t++) salph += gsl_pow_2(u[t] - grand); salph = salph / (indi - 1);
  	//Rprintf("Daten: Mean, Residual Variance, Due to Individual differences\n");
  	Rprintf("RT: mean, variance, residual variance\n");
  	Rprintf("%12.4g%12.4g%12.4g\n", grand, s, salph);
  	if (save_diagnose) {
  		tests_out << "Daten: Mean, Residual Variance, Due to Individual differences " << std::endl;
   		tests_out << std::setw(12) << grand << std::setw(12) << s << std::setw(12) << salph << std::endl;
  	}

  	R_CheckUserInterrupt();

  }
   
  void make_pij_for_one_trial_new(trial one, double *x_for_all, double *pij, double &pj) {
  	// berechnet  p
  
  
  #define X_FOR_ALL(T,IP) x_for_all[T*kernpar+IP]
  
  
  	double d0ha;
  	int j = one.category, t = one.person, itree = one.tree;
  
  
  	// pj = 0.0;
  	for (int k = 0; k != branch[j]; k++) {
  //			 pij[k]=1.0;
  		for (int ir = 0; ir != NDRIN(j, k); ir++) {
  			int r = DRIN(j, k, ir); int ix = AR(j, k, r); int ip = TREE_AND_NODE2PAR(itree, r);
  			d0ha = (ix > 0) ? lnnorm(X_FOR_ALL(t, ip)) : lnnorm(- X_FOR_ALL(t, ip));
  			pij[k] = pij[k] + d0ha;
  		}
  		if (k == 0) pj = pij[k]; else pj = logsum(pj,pij[k]);
  	}
  	if (DEBUG) if (!(std::isfinite(pj))) {
  		Rprintf("pj in diagnosis is not finite %d\n", branch[j]);
  	}
  }
  
  #define PFAD_INDEX(J,K) pfad_index[J*zweig+K]
  
  void make_tij_for_one_trial_new(trial one, double *rhos, double *lambdas, double *restpars, double *pij) {
  
  	int t = one.person; /*int itree = one.tree;*/ int j = one.category; double rt = one.rt / 1000.0; int resp = cat2resp[j];
  
  	double rmu = restpars[t2group[t] * respno + resp] + restpars[alphaoff + t * respno + resp]; double rsig = sqrt(restpars[t + sigalphaoff]);
  
  	for (int k = 0; k != branch[j]; k++) {
  		int pfadlength = NDRIN(j, k);
		std::vector<double> lams(pfadlength);
		int complength = 0;

  		int ipfad = PFAD_INDEX(j, k);
  		pfadinfo akt_pfad = path_info[ipfad];
  		complength = akt_pfad.a;
  		for (int ir = 0; ir != akt_pfad.a; ir++) {
  			int ip = akt_pfad.pfad_par[ir];
  			int pm = akt_pfad.pm[ir];
  			int iz = kern2free[(1 + pm)*kernpar + ip] - ifree; lams[ir] = LAMBDAS(t,iz)*RHOS(t2group[t], iz);
  		}
  		// }
  		if (complength == 0) { pij[k] = (-0.5*gsl_pow_2((rt - rmu) / rsig)) / sqr2pi / rsig; }

  		if (complength == 1 /*&& (PFAD_INDEX(j, k) > -1)*/) {

  			double lam = lams[0];
  			double hplus, hminus;
  			loggammagaussian(akt_pfad.r[0] - 1, lam, rmu, rsig, rt, hplus, hminus);
  			double temp = logdiff(hplus, hminus) + log(lam) * akt_pfad.r[0];
  			if (temp == GSL_NEGINF) {
  				pij[k] = temp;
  			}
  			else pij[k] = (temp);
  		}

		if ((complength >= 2) /*&& (PFAD_INDEX(j, k) > -1)*/) {
			std::vector<double> loglams(complength);
  			for (int ir = 0; ir != complength; ir++) loglams[ir] = log(lams[ir]);
  			// int ipfad = PFAD_INDEX(j, k);
  			// pfadinfo akt_pfad = path_info[ipfad];
  			double temp = logf_tij(akt_pfad.a, akt_pfad.r, lams.data(), loglams.data(), rmu, rsig, rt);
  			if (temp ==GSL_NEGINF) {
  				pij[k] = temp;
  
  			}
  			else pij[k] = temp;
  
  		}
  
  
  		//	 if ((pfadlength>5) || (pfadlength==0)) cout<< "pfadlength" << std::endl;
  	}
  }
  
  
  void dic(int N_all_parameters, std::vector <trial> daten, double *beta, double *sample) {
  
  	// std::ofstream log_lik(LOGLIK);
  	// loglik_vec = (double *)malloc(SAMPLE_SIZE * datenzahl * sizeof(double));
    
  	if (log_lik_flag) {
  	  // log_lik << std::setprecision(12);
  	}
  	double dbar = 0.0, pd = 0.0, pv = 0.0;
  
	std::vector<double> x_for_all(indi * kernpar);
	std::vector<double> xbar(n_all_parameters);
	std::vector<double> pij(zweig);
	std::vector<double> lambdas(ilamfree * indi);
	std::vector<double> rhos(ilamfree * igroup);
	std::vector<double> restpars(restparsno);
  	//		double *temp=0; if (!(temp=NAG_ALLOC(SAMPLE_SIZE*sizeof(double)))){printf("Allocation failure\n");exit_status = -1;}
  
  	int trialno = static_cast<int>(daten.size());
  	for (int i = 0; i != n_all_parameters; i++) xbar[i] = 0.0;
  	for (int is = 0; is != SAMPLE_SIZE; is++) {
  
  		for (int i = 0; i != n_all_parameters; i++) xbar[i] += SAMPLE(is, i) / (SAMPLE_SIZE);
  
  		// belege x_for_all:
  		belege_beta(sample, is, beta);
  		for (int t = 0; t != indi; t++)  for (int ip = 0; ip != kernpar; ip++) x_for_all[ip + t * kernpar] = comp[ip] ? /*gsl_cdf_ugaussian_P*/(BETA(t, kern2free[ip])) : /*gsl_cdf_ugaussian_P*/(consts[ip]);
  		// belege rhos, lambdas, restpars
  		belege_lambdas_rhos(sample, is, rhos.data(), lambdas.data());
  		for (int ix = 0; ix != restparsno; ix++) restpars[ix] = SAMPLE(is, n_all_parameters - restparsno + ix);
  		double persample = 0.0;
  		for (int x = 0; x != trialno; x++) {
  			trial one = daten[x];  int t = one.person; int r = cat2resp[one.category]; double rmu = restpars[t2group[t] * respno + r] + restpars[t*respno + r + alphaoff]; double rsig = sqrt(restpars[t + sigalphaoff]);
  			double xsi = /*gsl_cdf_ugaussian_P*/lnnorm(rmu / rsig);
			make_tij_for_one_trial_new(one, rhos.data(), lambdas.data(), restpars.data(), pij.data());
			double p;
			make_pij_for_one_trial_new(one, x_for_all.data(), pij.data(), p);
			p -= xsi;
			if (p==GSL_NEGINF)
				Rprintf("DIC loglik Problem\n");
  			if (log_lik_flag) {
  			  // log_lik << std::setw(20) << p;
  			  loglik_vec[is + SAMPLE_SIZE * x] = p;
  			}
  			persample += -2 * (p);
  		}
  		// if (log_lik_flag) log_lik << std::endl;
  		dbar += persample / (SAMPLE_SIZE);
  		pv += gsl_pow_2(persample) / (SAMPLE_SIZE);
		R_CheckUserInterrupt();
	}

  
  	double xn = SAMPLE_SIZE * 1.0;
  	pv = pv - gsl_pow_2(dbar); pv = xn / (xn - 1.0)*pv; pv = 0.5*pv;
  
  	//belege x_for_all und so fort mit xbar
  	for (int t = 0; t != indi; t++) for (int iz = 0; iz != ifree; iz++)  BETA(t, iz) = xbar[iz + t2group[t] * ifree] +
  		xbar[(iz + ilamfree * igroup + ((ifree + ilamfree)*(ifree + ilamfree + 1)) / 2 + (t + igroup)*ifree)];
  	for (int t = 0; t != indi; t++) 	for (int ip = 0; ip != kernpar; ip++) if (comp[ip])
  		x_for_all[ip + t * kernpar] = /*gsl_cdf_ugaussian_P*/(BETA(t, kern2free[ip])); else
  		x_for_all[ip + t * kernpar] = /*gsl_cdf_ugaussian_P*/(consts[ip]);
  	// belege rhos, lambdas, restpars mit xbar
  	for (int iz = 0; iz != ilamfree * igroup; iz++) rhos[iz] = xbar[ifree*igroup + iz];
  	for (int t = 0; t != indi; t++) for (int iz = 0; iz != ilamfree; iz++)
  		lambdas[t*ilamfree + iz] = xbar[ifree*igroup + ilamfree * igroup + ((ifree + ilamfree)*(ifree + ilamfree + 1)) / 2 + indi * ifree + t * ilamfree + iz];
  	for (int ix = 0; ix != restparsno; ix++) restpars[ix] = xbar[n_all_parameters - restparsno + ix];
  	//
  	for (int x = 0; x != trialno; x++) {
  		trial one = daten[x];  int t = one.person;  int r = cat2resp[one.category]; double rmu = restpars[t2group[t] * respno + r] + restpars[t*respno + r + alphaoff]; double rsig = sqrt(restpars[t + sigalphaoff]);
  		double xsi = /*gsl_cdf_ugaussian_P*/lnnorm(rmu / rsig);
		make_tij_for_one_trial_new(one, rhos.data(), lambdas.data(), restpars.data(), pij.data());
		double p;
		make_pij_for_one_trial_new(one, x_for_all.data(), pij.data(), p);
		p -= xsi;
		if (p==GSL_NEGINF) Rprintf("DIC loglik Problem in pd\n");
  		pd += -2 * p;
  	}
  
  
  	pd = dbar - pd;
  
  	//std::cout << std::setprecision(8);
  	if (save_diagnose) tests_out << std::setprecision(8);
  
  	Rprintf("DIC1, DIC2, pd, pv:\n");
  	Rprintf("%15.8g%15.8g%15.8g%15.8g\n", pd + dbar, pv + dbar, pd, pv);
  	if (save_diagnose) {
  		tests_out << "DIC1, DIC2, pd, pv: " << std::endl;
  		tests_out << std::setw(15) << pd + dbar << std::setw(15) << pv + dbar << std::setw(15) << pd << std::endl << std::setw(15) << pv << std::endl;
  	}
  	// log_lik.close();

  }
  
  
  
  void test(double *t1, double *t2, std::string what) {
  	double p = 0.0; double out[2] = { 0.0,0.0 };
  	double r, count;
  	for (int i = 0; i != SAMPLE_SIZE; i++) {
  		r = 1.0 / (i + 1);
  		count = (t1[i] < t2[i]) ? 1.0 : 0.0;
  		out[0] += (t1[i] - out[0])*r; out[1] += (t2[i] - out[1])*r; p += (count - p)*r;
  	}
  	Rprintf("\n");
  	Rprintf("%s\n", what.c_str());
  	//std::cout << std::setprecision(4);
  	Rprintf("%12.4g%12.4g%12.4g\n", out[0], out[1], p);
  	if (save_diagnose) {
  		tests_out << std::endl;
  		tests_out << what << std::endl;
  		tests_out << std::setprecision(4);
  		tests_out << std::setw(12) << out[0] << std::setw(12) << out[1] << std::setw(12) << p << std::endl;
  	}
  	double iv[2];
  	for (int is = 0; is != SAMPLE_SIZE; is++) t1[is] = (t1[is] - t2[is]);
  	gsl_sort(t1, 1, SAMPLE_SIZE);
  	hdi(SAMPLE_SIZE, t1, .95, iv);
  	Rprintf("95%% HDI\n"); if (save_diagnose) tests_out << "95% HDI" << std::endl;
  	for (int iq = 0; iq != 2; iq++) Rprintf("%12.4g", iv[iq]);
  	Rprintf("\n");
  	if (save_diagnose) {
  		for (int iq = 0; iq != 2; iq++) tests_out << std::setw(12) << iv[iq];
  		tests_out << std::endl;
  	}
  }
  
  void aggregate(int n_all_parameters, int kerntree, int *idaten, std::vector<trial> daten, int *nks, int *jks, int* tree2cat, double *beta, double *sample, gsl_rng *rst) {
  
  #define IDATEN(T,I) idaten[T*kerncat + I]
  
  #define NKS(T,IT) nks[T*kerntree+IT]
  #define DRIN(J,K,X) drin[J*zweig*nodemax+K*nodemax + X]
  #define NDRIN(J,K) ndrin[J*zweig+K]
  
  //#define AR(I,J,K) ar[I*zweig*nodemax + J*nodemax + K]
  //#define TREE_AND_NODE2PAR(I,J) tree_and_node2par[I*nodemax+J]
  #define PIJ(I,J) pij[I*zweig+J]
  #define TREE2CAT(IT,J) tree2cat[IT*kerncat+J]
  
  	//	Nag_ModeRNG mode = Nag_GenerateWithoutReference;Nag_OrderType order = Nag_RowMajor;NagError fail;	INIT_FAIL(fail);
  
	std::vector<double> t1(SAMPLE_SIZE);
	std::vector<double> t2(SAMPLE_SIZE);
	std::vector<int> obs(kerncat);
	std::vector<double> expe(kerncat);
	std::vector<int> rep(kerncat);
	std::vector<int> sobs(kerncat * igroup);
	std::vector<double> sexp(kerncat * igroup);
	std::vector<int> srep(kerncat * igroup);
	std::vector<double> tobs(kerncat);
	std::vector<double> texp(kerncat);
	std::vector<double> trep(kerncat);
	std::vector<double> stobs(kerncat * igroup);
	std::vector<double> stexp(kerncat * igroup);
	std::vector<double> strep(kerncat * igroup);
	std::vector<double> pij(zweig * kerncat);
	std::vector<double> onepij(zweig);
	std::vector<double> x(kernpar);
	std::vector<double> lambdas(ilamfree * indi);
	std::vector<double> tdaten(indi * kerncat);
	std::vector<int> nobs(kerncat * igroup);
	std::vector<int> nrep(kerncat * igroup);
	std::vector<double> d(kerncat);
	std::vector<double> x1(kerncat * igroup * SAMPLE_SIZE);
	std::vector<double> x2(kerncat * igroup * SAMPLE_SIZE);
	std::vector<unsigned int> drep(kerncat);
	std::vector<int> ng(igroup);
  
  #define X1(IS,IG,J) x1[IS*igroup*kerncat+ IG*kerncat+J]
  #define X2(IS,IG,J) x2[IS*igroup*kerncat+IG*kerncat+J]
  #define SOBS(IG,J) sobs[IG*kerncat+J]
  #define SEXP(IG,J) sexp[IG*kerncat+J]
  #define SREP(IG,J) srep[IG*kerncat+J]
  #define STOBS(IG,J) stobs[IG*kerncat+J]
  #define STEXP(IG,J) stexp[IG*kerncat+J]
  #define STREP(IG,J) strep[IG*kerncat+J]
  #define NOBS(IG,J) nobs[IG*kerncat+J]
  #define NREP(IG,J) nrep[IG*kerncat+J]
  
  
  	for (int is = 0; is != SAMPLE_SIZE; is++) {
  		for (int j = 0; j != kerncat * igroup; j++) { sexp[j] = 0.0; sobs[j] = srep[j] = 0; x1[j + is * kerncat*igroup] = 0.0; x2[j + is * kerncat*igroup] = 0.0; }
  		// compute exp pro Person und rep pro Person; aggregate, compute chi-square
  		belege_beta(sample, is, beta);
  		for (int t = 0; t != indi; t++) {
  			for (int ip = 0; ip != kernpar; ip++) x[ip] = comp[ip] ? gsl_cdf_ugaussian_P(BETA(t, kern2free[ip])) : gsl_cdf_ugaussian_P(consts[ip]);
  			make_pij_for_individual(x.data(), pij.data(), expe.data()); for (int j = 0; j != kerncat; j++) drep[j] = 0;
  			for (int it = 0; it != kerntree; it++) {
  				for (int j = 0; j != jks[it]; j++) d[j] = expe[TREE2CAT(it, j)];
  				gsl_ran_multinomial(rst, jks[it], NKS(t, it), d.data(), drep.data());
  				for (int j = 0; j != jks[it]; j++) rep[TREE2CAT(it, j)] = drep[j];
  			}
  			int ig = t2group[t];
  			for (int j = 0; j != kerncat; j++) { expe[j] = NKS(t, cat2tree[j])*expe[j]; obs[j] = IDATEN(t, j); }
  			for (int j = 0; j != kerncat; j++) { SEXP(ig, j) += expe[j]; SOBS(ig, j) += obs[j]; SREP(ig, j) += rep[j]; X1(is, ig, j) += (obs[j] * 1.0) / NKS(t, cat2tree[j]); X2(is, ig, j) += (rep[j] * 1.0) / NKS(t, cat2tree[j]); }
  			//			for (int it=0;it!=kerntree;it++){int no=0,ne=0; for (int j=0;j!=jks[it];j++) {no+=obs[TREE2CAT(it,j)];ne+=rep[TREE2CAT(it,j)];}; if (no!=ne) std::cout<< "hello" << std::endl;}
  
  		}
  		t1[is] = 0.0; t2[is] = 0.0; for (int ig = 0; ig != igroup; ig++) for (int j = 0; j != kerncat; j++) { t1[is] += gsl_pow_2(SOBS(ig, j) - SEXP(ig, j)) / SEXP(ig, j); t2[is] += gsl_pow_2(SREP(ig, j) - SEXP(ig, j)) / SEXP(ig, j); }
  	}
  	test(t1.data(), t2.data(), "Posterior predictive check: frequencies");
  
  	for (int ig = 0; ig != igroup; ig++) { ng[ig] = 0; for (int t = 0; t != indi; t++) ng[ig] += (t2group[t] == ig); }
  
  	double qv[5];
  	// int *correct = 0; correct = (int *)malloc(kerncat * sizeof(int));
  	// for (int j = 0; j != kerncat; j++) correct[j] = ((j % 6 == 1) || (j % 6 == 2) || (j % 6 == 4)) ? 1 : 0;
  	//std::cout << std::setprecision(4);
  	for (int ig = 0; ig != igroup; ig++)
  		for (int j = 0; j != kerncat; j++) {
for (int is = 0; is != SAMPLE_SIZE; is++) t1[is] = X1(is, ig, j) / ng[ig];
			gsl_sort(t1.data(), 1, SAMPLE_SIZE);
  			qv[2] = gsl_stats_median_from_sorted_data(t1.data(), 1, SAMPLE_SIZE);
  			Rprintf("%3d", j); Rprintf("%12.4g", qv[2]);
  			for (int is = 0; is != SAMPLE_SIZE; is++) t2[is] = X2(is, ig, j) / ng[ig];
  			gsl_sort(t2.data(), 1, SAMPLE_SIZE);
  			qv[2] = gsl_stats_median_from_sorted_data(t2.data(), 1, SAMPLE_SIZE);
  			double iv[2]; hdi(SAMPLE_SIZE, t2.data(), 0.95, iv); qv[1] = iv[0]; qv[3] = iv[1]; hdi(SAMPLE_SIZE, t2.data(), 0.99, iv); qv[0] = iv[0]; qv[4] = iv[1];
  			Rprintf("%12.4g%12.4g%12.4g\n", qv[1], qv[2], qv[3]);
  		}
  
  
  	// Das Ganze f�r die Zeiten
  
  	for (int j = 0; j != kerncat * indi; j++) tdaten[j] = 0.0;
  #define TDATEN(T,J) tdaten[T*kerncat+J]
  	for (int it = 0; it != static_cast<int>(daten.size()); it++) {
  		trial one = daten[it]; TDATEN(one.person, one.category) += one.rt / 1000.0;
  	}
  	for (int t = 0; t != indi; t++) for (int j = 0; j != kerncat; j++) if (IDATEN(t, j) > 0) TDATEN(t, j) /= IDATEN(t, j);
  
  
  
  	for (int is = 0; is != SAMPLE_SIZE; is++) {
  		for (int j = 0; j != kerncat * igroup; j++) { stobs[j] = stexp[j] = strep[j] = 0.0; nobs[j] = nrep[j] = 0; }
  		// compute exp pro Person und rep pro Person; aggregate, compute chi-square
  		belege_beta(sample, is, beta);
  		belege_nur_lambdas(sample, is, lambdas.data());
  
  		for (int j = 0; j != kerncat; j++) nobs[j] = nrep[j] = 0;
  
  		for (int t = 0; t != indi; t++) {
for (int ip = 0; ip != kernpar; ip++) x[ip] = comp[ip] ? gsl_cdf_ugaussian_P(BETA(t, kern2free[ip])) : gsl_cdf_ugaussian_P(consts[ip]);
			make_pij_for_individual(x.data(), pij.data(), expe.data());
  
  			for (int j = 0; j != kerncat; j++) {
  				int r = cat2resp[j];
  				double mu = SAMPLE(is, (n_all_parameters - restparsno + t2group[t] * respno + r)) + SAMPLE(is, (n_all_parameters - restparsno + alphaoff + t * respno + r)), sig = sqrt(SAMPLE(is, (n_all_parameters - indi + t))); double z = -mu / sig;
  				z = 1.0 / sqrt(2.0*M_PI)*exp(-0.5*gsl_pow_2(z)) / (1.0 - gsl_cdf_ugaussian_P(z));
  				texp[j] = 0.0;
  				for (int k = 0; k != branch[j]; k++)
  					for (int xr = 0; xr != NDRIN(j, k); xr++) {
  						int r = DRIN(j, k, xr); int ip = TREE_AND_NODE2PAR(cat2tree[j], r);  int pm = (AR(j, k, r) > 0) ? 1 : 0;
  						if (comp[(1 + pm)*kernpar + ip]) {
  							double lambda = LAMBDAS(t, (kern2free[(1 + pm)*kernpar + ip] - ifree)); texp[j] += PIJ(j, k)*1.0 / lambda;
  						}
  					}
  				texp[j] += mu + sig * z;
  
  				// compute exp und variance
  			}
  			for (int j = 0; j != kerncat; j++) drep[j] = 0;
  			for (int it = 0; it != kerntree; it++) {
  				for (int j = 0; j != jks[it]; j++) d[j] = expe[TREE2CAT(it, j)];
  				gsl_ran_multinomial(rst, jks[it], NKS(t, it), d.data(), drep.data());
  				//				nag_rand_gen_multinomial(order, mode, 1,NKS(t,it),jks[it],d, r, lr,rst,drep,jks[it],&fail); if (fail.code != NE_NOERROR){printf("Error from nag_rand_gen_multinomial (g05tgc).\n%s\n",fail.message);exit_status = 1;}
  				for (int j = 0; j != jks[it]; j++) rep[TREE2CAT(it, j)] = drep[j];
  			}
  			for (int j = 0; j != kerncat; j++) {
  				int r = cat2resp[j];
  				double mu = SAMPLE(is, (n_all_parameters - restparsno + t2group[t] * respno + r)) + SAMPLE(is, (n_all_parameters - restparsno + alphaoff + t * respno + r)), sig = sqrt(SAMPLE(is, (n_all_parameters - indi + t)));
  				trep[j] = 0.0;
  				tobs[j] = TDATEN(t, j);
  				for (int ir = 0; ir != rep[j]; ir++) {
  					double temp = 0.0;
  					for (int k = 0; k != branch[j]; k++) onepij[k] = log(PIJ(j, k));
  					int ipath = make_path_for_one_trial(branch[j], onepij.data(), 0.0, rst);
  					for (int xr = 0; xr != NDRIN(j, ipath); xr++) {
  						int r = DRIN(j, ipath, xr); int ip = TREE_AND_NODE2PAR(cat2tree[j], r);  int pm = (AR(j, ipath, r) > 0) ? 1 : 0;
  						if (comp[(1 + pm)*kernpar + ip]) {
  							double lambda = LAMBDAS(t, (kern2free[(1 + pm)*kernpar + ip] - ifree));
  							temp += oneexp(lambda, rst);
  						}
  					}
  					temp += truncnorm(mu / sig, rst)*sig;
  					trep[j] += temp;
  				}
  				if (rep[j] > 0) trep[j] /= rep[j]; else trep[j] = 0.0;
  			}
  			int ig = t2group[t];
  			for (int j = 0; j != kerncat; j++) {
  				STEXP(ig, j) += texp[j]; if (tobs[j] > 0) { STOBS(ig, j) += tobs[j]; NOBS(ig, j)++; } if (trep[j] > 0) { STREP(ig, j) += trep[j]; NREP(ig, j)++; }
  			}
  		}
  
  		for (int ig = 0; ig != igroup; ig++) for (int j = 0; j != kerncat; j++) { STEXP(ig, j) /= ng[ig]; STOBS(ig, j) /= NOBS(ig, j); if (NREP(ig, j) > 0) STREP(ig, j) /= NREP(ig, j); X1(is, ig, j) = STOBS(ig, j); X2(is, ig, j) = STREP(ig, j); }
  		t1[is] = 0.0; t2[is] = 0.0;
  
  		for (int j = 0; j != kerncat * igroup; j++) {
  			t1[is] += gsl_pow_2(stobs[j] - stexp[j]) / stexp[j]; t2[is] += gsl_pow_2(strep[j] - stexp[j]) / stexp[j];
  		}
  	}
  	test(t1.data(), t2.data(), "Posterior predictive check: latencies");
  
  
  
  	//std::cout << std::setprecision(4);
  	for (int ig = 0; ig != igroup; ig++)
  		for (int j = 0; j != kerncat; j++) {
for (int is = 0; is != SAMPLE_SIZE; is++) t1[is] = X1(is, ig, j);
			gsl_sort(t1.data(), 1, SAMPLE_SIZE);
  			qv[2] = gsl_stats_median_from_sorted_data(t1.data(), 1, SAMPLE_SIZE);
  			Rprintf("%3d", j); Rprintf("%12.4g", qv[2]);
  			for (int is = 0; is != SAMPLE_SIZE; is++) t2[is] = X2(is, ig, j);
  			gsl_sort(t2.data(), 1, SAMPLE_SIZE);
  			qv[2] = gsl_stats_median_from_sorted_data(t2.data(), 1, SAMPLE_SIZE);
  			double iv[2]; hdi(SAMPLE_SIZE, t2.data(), 0.95, iv); qv[1] = iv[0]; qv[3] = iv[1]; hdi(SAMPLE_SIZE, t2.data(), 0.99, iv); qv[0] = iv[0]; qv[4] = iv[1];
  			Rprintf("%12.4g%12.4g%12.4g\n", qv[1], qv[2], qv[3]);
  		}
  	// if (correct) free(correct);
  

  }
  

  void correlation(double *sample, int *index1, int *index2)
  {
  #define SIGMA(I,J) sigma[I*(ilamfree+ifree)+J]
	double qv[5];
	std::vector<double> temp(SAMPLE_SIZE);
	std::vector<double> sigma((ilamfree + ifree)*(ilamfree + ifree));
  	for (int is = 0; is != SAMPLE_SIZE; is++) {
  		int iz = (igroup)*(ifree + ilamfree) - 1;
  		for (int ip = 0; ip != ilamfree + ifree; ip++)
  			for (int jp = ip; jp != ilamfree + ifree; jp++) {
  				iz++;
  				SIGMA(ip, jp) = SAMPLE(is, iz); SIGMA(jp, ip) = SIGMA(ip, jp);
  			}
  		for (int ip = 0; ip != ilamfree + ifree; ip++) for (int jp = 0; jp != ilamfree + ifree; jp++) if (ip != jp)
  			SIGMA(ip, jp) = sqrt(SIGMA(ip, ip))*SIGMA(ip, jp)*sqrt(SIGMA(jp, jp));
  		double sig1 = 0.0, sig2 = 0.0, cov = 0.0;
  		for (int ip = 0; ip != ilamfree + ifree; ip++) for (int jp = 0; jp != ilamfree + ifree; jp++) {
  			sig1 += index1[ip] * index1[jp] > 0 ? SIGMA(ip, jp) : 0;
  			sig2 += index2[ip] * index2[jp] > 0 ? SIGMA(ip, jp) : 0;
  			cov += index1[ip] * index2[jp] > 0 ? SIGMA(ip, jp) : 0;
  		}
  		temp[is] = cov / sqrt(sig1*sig2);
  	}
  	gsl_sort(temp.data(), 1, SAMPLE_SIZE);
  	qv[2] = gsl_stats_median_from_sorted_data(temp.data(), 1, SAMPLE_SIZE);
  	double iv[2]; hdi(SAMPLE_SIZE, temp.data(), 0.95, iv); qv[1] = iv[0]; qv[3] = iv[1]; hdi(SAMPLE_SIZE, temp.data(), 0.99, iv); qv[0] = iv[0]; qv[4] = iv[1];
  	Rprintf("Corr"); for (int iq = 0; iq != 5; iq++) Rprintf("%12.4g", qv[iq]); Rprintf("\n");
  	if (save_diagnose) { tests_out << "Corr "; for (int iq = 0; iq != 5; iq++) tests_out << std::setw(12) << qv[iq]; tests_out << std::endl; }
  
  }
  
	void groupwise(double *sample) {

	std::vector<double> t1(SAMPLE_SIZE);
	std::vector<double> t2(SAMPLE_SIZE);
  
  	for (int ip = 0; ip != ifree; ip++)
  	{
  		for (int is = 0; is != SAMPLE_SIZE; is++) {
  			//				par.push_back(gsl_cdf_ugaussian_P(SAMPLE(is,ip))-gsl_cdf_ugaussian_P(SAMPLE(is,ip+ifree)));
  			t1[is] = gsl_cdf_ugaussian_P(SAMPLE(is, ip));
		t2[is] = gsl_cdf_ugaussian_P(SAMPLE(is, ip + ifree));
		}
		test(t1.data(), t2.data(), "group-tests mu");
  	}
  
  	for (int ip = 0; ip != ilamfree; ip++)
  	{
  
  		for (int is = 0; is != SAMPLE_SIZE; is++) {
  			//				par.push_back(1000.0/SAMPLE(is,ip+ifree*igroup)-1000.0/SAMPLE(is,ip+ifree*igroup+2*ilamfree));
  			t1[is] = 1000.0 / SAMPLE(is, ip + ifree * igroup);
		t2[is] = 1000.0 / SAMPLE(is, ip + ifree * igroup + ilamfree);
		}
		test(t1.data(), t2.data(), "group-tests pho");
  	}
  	int    iz = ifree * igroup + ilamfree * igroup + ((ifree + ilamfree)*(ifree + ilamfree + 1)) / 2 + indi * ifree + indi * ilamfree;
  	for (int ip = 0; ip != respno; ip++)
  	{
  
  		for (int is = 0; is != SAMPLE_SIZE; is++) {
  			//				par.push_back(SAMPLE(is,iz+ip)-SAMPLE(is,ip+iz+ respno));
  			t1[is] = SAMPLE(is, iz + ip);
		t2[is] = SAMPLE(is, iz + ip + respno);
		}
		test(t1.data(), t2.data(), "group-tests residual");
  	}
  	iz = ifree * igroup + ilamfree * igroup + ((ifree + ilamfree)*(ifree + ilamfree + 1)) / 2 + indi * ifree + indi * ilamfree;

  }
  
  void diagnosis(std::vector<trial> daten, int *idaten, int kerntree, gsl_rng *rst) {
	std::vector<int> nks(indi * kerntree);
	std::vector<int> jks(kerntree);
	std::vector<int> tree2cat(kerntree * kerncat);
	std::vector<double> beta(indi * ifree);

	n_all_parameters = ifree * igroup + ilamfree * igroup + ((ifree + ilamfree)*(ifree + ilamfree + 1)) / 2 + indi * ifree + indi * ilamfree + restparsno;

	std::vector<double> sample(SAMPLE_SIZE * (n_all_parameters + 1));
	lies(n_all_parameters, sample.data());
	if (save_diagnose) tests_out.open(diagn_tests);
	quantiles(daten, n_all_parameters, sample.data());
  	// make nks
  
  	for (int t = 0; t != indi; t++) for (int it = 0; it != kerntree; it++) NKS(t, it) = 0;
  	for (int t = 0; t != indi; t++) for (int j = 0; j != kerncat; j++) NKS(t, cat2tree[j]) += IDATEN(t, j);
  	for (int it = 0; it != kerntree; it++) jks[it] = 0;
  	for (int j = 0; j != kerncat; j++) { TREE2CAT(cat2tree[j], jks[cat2tree[j]]) = j; jks[cat2tree[j]]++; }

	dic(n_all_parameters, daten, beta.data(), sample.data());
	aggregate(n_all_parameters, kerntree, idaten, daten, nks.data(), jks.data(), tree2cat.data(), beta.data(), sample.data(), rst);
  	//    SSE_by_individual(n_all_parameters,kerntree, idaten, daten,nks,jks,tree2cat, beta,index,sample,rst);
  
  	// if (igroup > 1) groupwise(sample);
  
  
  	if (save_diagnose) tests_out.close();
  

  }

}



namespace drtmpt {
  
  //implements summary statistics on posterior distribution
  //goodness-of-fit tests
  //specific contrast
  
  
  //save results in tests.out
  std::ofstream tests_out;
  
  //read file with posterior samples
  void lies_sample(int n_all_parameters, std::vector<double>& sample) {
    
    std::ifstream rein(RAUS);
    int is, in;
    rein >> is >> in;
    
    sample.resize(is * n_all_parameters);
    Rprintf("\nTotal sample size is %d x %d = %d\n", NOTHREADS, is / NOTHREADS, is);
    sample_size = is;
    if (in != (n_all_parameters)) Rprintf("HO\n");
    for (int i = 0; i != is; i++)
      for (int j = 0; j != in; j++) rein >> dSAMPLE(i, j);
    rein.close();
  }
  
  //compute highest density intervals
  void hdi(int length, double * parameter, double p, double iv[2]) {
    //	int inc = static_cast<int>(p*length) + 1;
    // corrected:
    int inc = static_cast<int>(p * length);
    int nci = length - inc;
    int imin = -1;
    double tempold = parameter[length - 1] - parameter[0];
    for (int i = 0; i != nci; i++) {
      double temp = parameter[i + inc] - parameter[i];
      if (temp < tempold) {
        tempold = temp;
        imin = i;
      }
    }
    iv[0] = parameter[imin];
    iv[1] = parameter[imin + inc];
  }
  
  
  //transform parameters to diffusion-model parameters
  void belege_ts(double* sample, int is, double* tavw) {
    for (int t = 0; t != indi; t++) {
      int jj = 0;
      for (int type = 0; type != 3; type++) {
        int ift = ifree[type];
        for (int iz = 0; iz != ift; iz++) {
          if (dCOMP(type, iz)) {
            dTAVW(t, type, iz) = logit(avwtrans[type], dSAMPLE(is, jj + icompg * t2group[t]) + dSAMPLE(is, jj + icompg * igroup + t * icompg));
            jj++;
          }
          else dTAVW(t, type, iz) = dCONSTS(type, iz);
        }
      }
    }
  }
  
  //transform parameters to motor time means per person and personwise standard deviations of motor times
  void belege_lambdas_mus(double* sample, int is, double* lambdas) {
    
    for (int t = 0; t != indi; t++) {
      for (int r = 0; r != respno; r++) lambdas[t * respno + r] = dSAMPLE(is, irmuoff + t2group[t] * respno + r) +
        dSAMPLE(is, ilamoff + t * respno + r);
      lambdas[indi * respno + t] = exp(dSAMPLE(is, isigoff + t));
    }
  }
  
  #define dMTAVW(T, Type, IP) mtavw[T * 3 * ifreemax + Type * ifreemax + IP]
  
  void write_ind_estimates(double* sample) {
    std::vector<double> tavw(ifreemax * 3 * indi, 0.0);
    std::vector<double> lambdas((respno + 1) * indi, 0.0);
    std::vector<double> mtavw(ifreemax * 3 * indi, 0.0);
    std::vector<double> mlambdas((respno + 1) * indi, 0.0);
    
    std::ofstream person("persons");
    for (int is = 0; is != sample_size; is++) {
      double r = 1.0 / (is + 1);
      belege_ts(sample, is, tavw.data());
      for (int t = 0; t != indi; t++)
        for (int type = 0; type != 3; type++) {
          int ift = ifree[type];
          for (int iz = 0; iz != ift; iz++) {
            if (dCOMP(type, iz)) {
              dMTAVW(t, type, iz) += (dTAVW(t, type, iz) - dMTAVW(t, type, iz)) * r;
            }
          }
        }
        
      belege_lambdas_mus(sample, is, lambdas.data());
      int resin = (respno + 1) * indi;
      for (int ir = 0; ir != resin; ir++) mlambdas[ir] += (lambdas[ir] - mlambdas[ir]) * r;
    }
    
    for (int t = 0; t != indi; t++) {
      person << std::setw(5) << t;
      for (int type = 0; type != 3; type++) {
        int ift = ifree[type];
        for (int iz = 0; iz != ift; iz++) {
          if (dCOMP(type, iz)) person << std::setw(20) << dMTAVW(t, type, iz);
        }
        for (int ir = 0; ir != respno; ir++) person << std::setw(20) << mlambdas[t * respno + ir];
        person << std::setw(20) << mlambdas[indi * respno + t] << std::endl;
      }
    }
    person.close();

  }
  
  
  //summary statistics for population means
  void quantiles(const std::vector<trial> & daten, int n_all_parameters, double* sample) {
    
    double  qv[5];
    std::vector<double> temp(sample_size);
    //std::streamsize prec = std::cout.precision(); std::cout << std::setprecision(4);
    if (save_diagnose) tests_out << std::setprecision(4);
    Rprintf("mean thresholds per group [median, 96 and 99%% HDI]\n");
    if (save_diagnose) tests_out << "Thresholds per group" << std::endl;
    int if0 = ifree[0], if1 = ifree[1], if2 = ifree[2];
    for (int ig = 0; ig != igroup; ig++) {
      int jz = 0;
      for (int iz = 0; iz != if0; iz++) if (dCOMP(0, iz)) {
        for (int j = 0; j != sample_size; j++) temp[j] = logit(avwtrans[0], dSAMPLE(j, jz + ig * icompg));
        gsl_sort(temp.data(), 1, sample_size);
        qv[2] = gsl_stats_median_from_sorted_data(temp.data(), 1, sample_size);
        double iv[2]; hdi(sample_size, temp.data(), 0.95, iv); qv[1] = iv[0]; qv[3] = iv[1]; hdi(sample_size, temp.data(), 0.99, iv); qv[0] = iv[0]; qv[4] = iv[1];
        Rprintf("%3d%3d", ig + 1, jz + 1);
        for (int iq = 0; iq != 5; iq++) Rprintf("%12.4g", qv[iq]);
        Rprintf("\n");
        if (save_diagnose) { tests_out << std::setw(3) << ig + 1 << std::setw(3) << jz + 1; for (int iq = 0; iq != 5; iq++) tests_out << std::setw(12) << qv[iq]; tests_out << std::endl; }
        jz++;
      }
    }
    Rprintf("mean drift rates per group [median, 96 and 99%% HDI]\n");
    if (save_diagnose) tests_out << "Drift rates per group" << std::endl;
    for (int ig = 0; ig != igroup; ig++) {
      int jz = 0;
      for (int iz = 0; iz != if1; iz++) if (dCOMP(1, iz)) {
        for (int j = 0; j != sample_size; j++) temp[j] = logit(avwtrans[1], dSAMPLE(j, jz + icomp[0] + ig * icompg));
        gsl_sort(temp.data(), 1, sample_size);
        qv[2] = gsl_stats_median_from_sorted_data(temp.data(), 1, sample_size);
        double iv[2]; hdi(sample_size, temp.data(), 0.95, iv); qv[1] = iv[0]; qv[3] = iv[1]; hdi(sample_size, temp.data(), 0.99, iv); qv[0] = iv[0]; qv[4] = iv[1];
        Rprintf("%3d%3d", ig + 1, jz + 1);
        for (int iq = 0; iq != 5; iq++) Rprintf("%12.4g", qv[iq]);
        Rprintf("\n");
        if (save_diagnose) { tests_out << std::setw(3) << ig + 1 << std::setw(3) << jz + 1;  for (int iq = 0; iq != 5; iq++) tests_out << std::setw(12) << qv[iq]; tests_out << std::endl; }
        jz++;
      }
    }
    Rprintf("mean rel. starting points per group [median, 96 and 99%% HDI]\n");
    if (save_diagnose) tests_out << "Start point per group" << std::endl;
    for (int ig = 0; ig != igroup; ig++) {
      int jz = 0;
      for (int iz = 0; iz != if2; iz++) if (dCOMP(2, iz)) {
        for (int j = 0; j != sample_size; j++) temp[j] = logit(avwtrans[2], dSAMPLE(j, jz + icomp[0] + icomp[1] + ig * icompg));
        gsl_sort(temp.data(), 1, sample_size);
        qv[2] = gsl_stats_median_from_sorted_data(temp.data(), 1, sample_size);
        double iv[2]; hdi(sample_size, temp.data(), 0.95, iv); qv[1] = iv[0]; qv[3] = iv[1]; hdi(sample_size, temp.data(), 0.99, iv); qv[0] = iv[0]; qv[4] = iv[1];
        Rprintf("%3d%3d", ig + 1, jz + 1);
        for (int iq = 0; iq != 5; iq++) Rprintf("%12.4g", qv[iq]);
        Rprintf("\n");
        if (save_diagnose) { tests_out << std::setw(3) << ig + 1 << std::setw(3) << jz + 1; for (int iq = 0; iq != 5; iq++) tests_out << std::setw(12) << qv[iq]; tests_out << std::endl; }
        jz++;
      }
    }
    
    Rprintf("SIGMA [median, 96 and 99%% HDI]\n");
    if (save_diagnose) tests_out << "SIG" << std::endl;
    int iz = isigoff + indi - 1;
    for (int ix = 0; ix != icompg; ix++) {
      for (int jz = ix; jz != icompg; jz++) {
        iz++;
        for (int j = 0; j != sample_size; j++) temp[j] = dSAMPLE(j, iz);
        gsl_sort(temp.data(), 1, sample_size);
        qv[2] = gsl_stats_median_from_sorted_data(temp.data(), 1, sample_size);
        double iv[2]; hdi(sample_size, temp.data(), 0.95, iv); qv[1] = iv[0]; qv[3] = iv[1]; hdi(sample_size, temp.data(), 0.99, iv); qv[0] = iv[0]; qv[4] = iv[1];
        Rprintf("%3d%3d", ix + 1, jz + 1);
        for (int iq = 0; iq != 5; iq++) Rprintf("%12.4g", qv[iq]);
        Rprintf("\n");
        if (save_diagnose) { tests_out << std::setw(3) << ix + 1 << std::setw(3) << jz + 1; for (int iq = 0; iq != 5; iq++) tests_out << std::setw(12) << qv[iq]; tests_out << std::endl; }
      }
    }
      
    Rprintf("Mean motor/encoding times per group [median, 96 and 99%% HDI]\n");
    if (save_diagnose) tests_out << "Mean motor/encoding times per group" << std::endl;
    iz = irmuoff;
    
    for (int ir = 0; ir != igroup * respno; ir++) {
      for (int j = 0; j != sample_size; j++) temp[j] = dSAMPLE(j, (iz + ir));
      gsl_sort(temp.data(), 1, sample_size);
      qv[2] = gsl_stats_median_from_sorted_data(temp.data(), 1, sample_size);
      double iv[2]; hdi(sample_size, temp.data(), 0.95, iv); qv[1] = iv[0]; qv[3] = iv[1]; hdi(sample_size, temp.data(), 0.99, iv); qv[0] = iv[0]; qv[4] = iv[1];
      for (int iq = 0; iq != 5; iq++) Rprintf("%12.4g", qv[iq]);
      Rprintf("\n");
      if (save_diagnose) { for (int iq = 0; iq != 5; iq++) tests_out << std::setw(12) << qv[iq]; tests_out << std::endl; }
    }
    
    Rprintf("GAMMA [median, 96 and 99%% HDI]\n");
    if (save_diagnose) tests_out << "Sig motor/encoding times" << std::endl;
    iz =  (icompg * (icompg + 1)) / 2 + isigoff + indi- 1;
    for (int ip = 0; ip != respno; ip++) {
      for (int jp = ip; jp != respno; jp++) {
        iz++;
        for (int j = 0; j != sample_size; j++) temp[j] = dSAMPLE(j, iz);
        gsl_sort(temp.data(), 1, sample_size);
        qv[2] = gsl_stats_median_from_sorted_data(temp.data(), 1, sample_size);
        double iv[2]; hdi(sample_size, temp.data(), 0.95, iv); qv[1] = iv[0]; qv[3] = iv[1]; hdi(sample_size, temp.data(), 0.99, iv); qv[0] = iv[0]; qv[4] = iv[1];
        Rprintf("%3d%3d", ip + 1, jp + 1);
        for (int iq = 0; iq != 5; iq++) Rprintf("%12.4g", qv[iq]);
        Rprintf("\n");
        if (save_diagnose) { tests_out << std::setw(3) << ip + 1 << std::setw(3) << jp + 1; for (int iq = 0; iq != 5; iq++) tests_out << std::setw(12) << qv[iq]; tests_out << std::endl; }
      }
    }
    
    Rprintf("Omega^2 [median, 96 and 99%% HDI]\n");
    if (save_diagnose) tests_out << "Residual variance" << std::endl;
    iz = n_all_parameters - 1;
    for (int j = 0; j != sample_size; j++) temp[j] = dSAMPLE(j, iz);
    gsl_sort(temp.data(), 1, sample_size);
    qv[2] = gsl_stats_median_from_sorted_data(temp.data(), 1, sample_size);
    double iv[2]; hdi(sample_size, temp.data(), 0.95, iv); qv[1] = iv[0]; qv[3] = iv[1]; hdi(sample_size, temp.data(), 0.99, iv); qv[0] = iv[0]; qv[4] = iv[1];
    for (int iq = 0; iq != 5; iq++) Rprintf("%12.4g", qv[iq]);
    Rprintf("\n");
    if (save_diagnose) {
      for (int iq = 0; iq != 5; iq++) {tests_out << std::setw(12) << qv[iq];} tests_out << std::endl;
    }
    
    
    // Daten zum Vergleich:
    double s = 0.0;
    std::vector<double> u(indi, 0.0);
    std::vector<int> nj(indi, 0);
    
    for (int i = 0; i != datenzahl; i++) { u[daten[i].person] += daten[i].rt / 1000.0; nj[daten[i].person]++; }
    for (int t = 0; t != indi; t++) { u[t] /= nj[t]; }
    
    for (int i = 0; i != datenzahl; i++) {s += gsl_pow_2(daten[i].rt / 1000.0 - u[daten[i].person]);} s = s / (datenzahl - 1);
    double grand = 0.0; for (int t = 0; t != indi; t++) grand += (u[t] * nj[t]) / datenzahl;
    double salph = 0.0; for (int t = 0; t != indi; t++) salph += gsl_pow_2(u[t] - grand); salph = salph / (indi - 1);
    Rprintf("RT: mean, variance, residual variance, due to Individual differences\n");
    Rprintf("%12g%12g%12g\n", grand, s, salph);
    if (save_diagnose) {
      tests_out << "RT: mean, variance, residual variance, Due to Individual differences " << std::endl;
      tests_out << std::setw(12) << grand << std::setw(12) << s << std::setw(12) << salph << std::endl;
    }
    
    R_CheckUserInterrupt();
  }
  
  //compute likelihoods for dic and log-likelihood
  void make_p_ind_cat(std::vector<double> rts, int t, int j, double* tavw, double mu, double sig, std::vector<double>& ps) {
    int itree = cat2tree[j];
    
    std::vector < std::vector<double>> p; p.clear();
    
    for (int k = 0; k != branch[j]; k++) {
      int pfadlength = dNDRIN(j, k);
      std::vector<double> a(pfadlength);
      std::vector<double> v(pfadlength);
      std::vector<double> w(pfadlength);
      std::vector<int> low_or_up(pfadlength);
      for (int ir = 0; ir != pfadlength; ir++) {
        int r = dDRIN(j, k, ir);
        low_or_up[ir] = dAR(j, k, r);
        a[ir] = dTAVW(t, 0, dTREE_AND_NODE2PAR(itree, r, 0));
        v[ir] = dTAVW(t, 1, dTREE_AND_NODE2PAR(itree, r, 1));
        w[ir] = dTAVW(t, 2, dTREE_AND_NODE2PAR(itree, r, 2));
      }
      
      std::vector<double> pbranch; pbranch.clear();
      convolution2(rts, pfadlength, low_or_up.data(), a.data(), v.data(), w.data(), mu, sig, pbranch);
      p.push_back(pbranch);
    }
    ps.clear();
    for (int x = 0; x != static_cast<int>(rts.size()); x++) {
      double temp = 0.0;
      for (int k = 0; k != branch[j]; k++) temp += p[k][x];
      ps.push_back(temp);
    }
  }
  
  //compute DIC and output likelihood for WAIC, etc. (optional)
  void dic(std::vector <trial> daten, double* sample) {
    // std::ofstream log_lik(LOGLIK);
    // if (log_lik_flag) {
    //   log_lik << std::setprecision(12);
    // }

    double dbar = 0.0, pv = 0.0; //, pd = 0.0
    
    std::vector<double> tavw(ifreemax * 3 * indi, 0.0);
    std::vector<double> lambdas((respno + 1) * indi, 0.0);
    
    std::vector<double> tavw_old(ifreemax * 3 * indi, 0.0);
    std::vector<double> lambdas_old((respno + 1) * indi, 0.0);
    
    
    std::vector<double> temp; temp.clear();
    
    std::vector<std::vector<double>> icdaten(indi*kerncat, temp);
    std::vector<std::vector<double>> icstore(indi*kerncat, temp);
    
    
    for (int x = 0; x != datenzahl; x++) {
      trial one = daten[x]; int t = one.person; int c = one.category;
      icdaten[t * kerncat + c].push_back(one.rt / 1000.0);
    }
    for (int t = 0; t != indi; t++) for (int c = 0; c != kerncat; c++) {
      std::sort(icdaten[t * kerncat + c].begin(), icdaten[t * kerncat + c].end());
    }
    
    double progress = 0.0;
    int ML_bar = 50;
    if (PROG_BAR_FLAG) {
      Rprintf("\nCalculating DIC.\nThis requires numerical integration and takes some time.\n");
      Rprintf("[");
      for (int i = 0; i < ML_bar; i++) Rprintf(" ");
      Rprintf("] 0%%");
    }
    
    double persample_old = 0.0;
    int myind = 0;
    for (int is = 0; is != sample_size; is++) {
      
      progress = 1.0*(is+1)/sample_size;
      
      double persample = 0.0;
      belege_ts(sample, is, tavw.data());
      belege_lambdas_mus(sample, is, lambdas.data());
      bool same = true;
      for (int i = 0; (i != 3 * ifreemax * indi) && (same); i++) same = (same) && (tavw[i] == tavw_old[i]);
      if (same) for (int i = 0; (i != (respno + 1) * indi) && (same); i++) same = (same) && (lambdas[i] == lambdas_old[i]);
      if (!same) {
        std::vector<double> icpersample(kerncat * indi, 0.0);
        
        /* prepare threads */
        int NThreads = DIC_CPUs;
        int maxThreads = std::thread::hardware_concurrency();
        if (maxThreads == 0) Rprintf("Could not find out number of threads. Taking 2 threads.\n");
        int suppThreads = maxThreads == 0 ? 2 : maxThreads;
        int AmntOfThreads = suppThreads > NThreads ? NThreads : suppThreads;
        int NperThread = indi / AmntOfThreads;
        AmntOfThreads = AmntOfThreads > indi ? indi : AmntOfThreads;
        NperThread = AmntOfThreads == indi ? 1 : NperThread;
        std::vector<std::thread> threads(AmntOfThreads-1);
        
        /* starting threads while ... */
        for (int ithread = 0; ithread < AmntOfThreads-1; ithread++) {
          threads[ithread] = std::thread([&, ithread]() {
            for (int t = ithread*NperThread; t < (ithread+1)*NperThread; t++) {
              for (int j = 0; j < kerncat; j++) {
                if (icdaten[t * kerncat + j].size() > 0) {
                  icstore[t * kerncat + j].clear();
                  std::vector<double> ps; ps.clear();
                  int r = cat2resp[j];
                  double mu = lambdas[t * respno + r]; double sig = lambdas[indi * respno + t];
                  double xsi = log(gsl_cdf_tdist_P(mu / sig, degf) * sig);
                  make_p_ind_cat(icdaten[t * kerncat + j], t, j, tavw.data(), mu, sig, ps);
                  for (int x = 0; x != static_cast<int>(icdaten[t * kerncat + j].size()); x++) {
                    double p = ps[x];
                    if ((p <= 0) || !(p == p))
                      Rprintf("DIC loglik Problem\n");
                    if (log_lik_flag) {
                      icstore[t * kerncat + j].push_back(log(p) - xsi);
                    }
                    icpersample[t * kerncat + j] += -2 * log(p);
                  }
                  icpersample[t * kerncat + j] += 2 * icdaten[t * kerncat + j].size() * xsi;
                }
              }
            }
          });
        }
        
        /* the main thread also runs */
        for (int t = (AmntOfThreads-1)*NperThread; t < indi; t++) {
          for (int j = 0; j < kerncat; j++) {
            if (icdaten[t * kerncat + j].size() > 0) {
              icstore[t * kerncat + j].clear();
              std::vector<double> ps; ps.clear();
              int r = cat2resp[j];
              double mu = lambdas[t * respno + r]; double sig = lambdas[indi * respno + t];
              double xsi = log(gsl_cdf_tdist_P(mu / sig, degf) * sig);
              make_p_ind_cat(icdaten[t * kerncat + j], t, j, tavw.data(), mu, sig, ps);
              for (int x = 0; x != static_cast<int>(icdaten[t * kerncat + j].size()); x++) {
                double p = ps[x];
                if ((p <= 0) || !(p == p))
                  Rprintf("DIC loglik Problem\n");
                if (log_lik_flag) {
                  icstore[t * kerncat + j].push_back(log(p) - xsi);
                }
                icpersample[t * kerncat + j] += -2 * log(p);
              }
              icpersample[t * kerncat + j] += 2 * icdaten[t * kerncat + j].size() * xsi;
            }
          }
        }
        
        for (int ithread = 0; ithread < AmntOfThreads-1; ithread++) {
          threads[ithread].join();
        }
        
        
        for (int tj = 0; tj != indi * kerncat; tj++)  persample += icpersample[tj];
      }
      if (log_lik_flag) {
        for (int t = 0; t < indi; t++) {
          for (int j = 0; j < kerncat; j++) {
            int tmp_K = static_cast<int>(icstore[t * kerncat + j].size());
            // Rprintf("tmp_K = %d\n", tmp_K);
            for (int k = 0; k != tmp_K; k++) {
              // log_lik << std::setw(20) << icstore[t * kerncat + j][k];
              
              loglik_vec[myind] = icstore[t * kerncat + j][k];
              myind++;
              // Rprintf("%10d\n" , index[t*kerncat + j][k] );
            }
          }
        }
      }

      if (same) persample = persample_old; else persample_old = persample;
      if (!same) {
        for (int i = 0; i != ifreemax * 3 * indi; i++) tavw_old[i] = tavw[i];
        for (int i = 0; i != (respno + 1) * indi; i++) lambdas_old[i] = lambdas[i];
      }
      //Rprintf("%20g%10d\n", persample, same);
      // log_lik << std::endl;
      dbar += persample;
      pv += gsl_pow_2(persample);
      //Rprintf("is %d\n", is);
      
      //#pragma optimize("", off)
      if (PROG_BAR_FLAG) {
        int pos = static_cast<int>(ML_bar * progress);
        Rprintf("\r[");
        for (int i = 0; i < ML_bar; i++) {
          if (i < pos) {
            Rprintf("=");
          } else if (i == pos) {
            Rprintf(">");
          } else Rprintf(" ");
        }
        Rprintf("] %d%%", static_cast<int>(progress * 100.0));
      }
      //#pragma optimize("", on)
      
      R_CheckUserInterrupt();
          
    }
    Rprintf("\n");
    double xn = sample_size * 1.0; dbar /= xn; pv /= xn;
    
    pv = pv - gsl_pow_2(dbar); pv = xn / (xn - 1.0) * pv; pv = 0.5 * pv;
    //std::cout << std::setprecision(8);
    if (save_diagnose) tests_out << std::setprecision(8);
    
    Rprintf("DIC, pv: \n");
    Rprintf("%15g%15g\n", pv + dbar, pv);
    if (save_diagnose) {
      tests_out << "DIC2,  pv: " << std::endl;
      tests_out << std::setw(15) << pv + dbar << std::setw(15) << pv << std::endl;
    }
    // log_lik.close();
    R_CheckUserInterrupt();
  }
  
  //simple posterior contrast test
  void test(double *t1, double *t2, std::string what) {
    double p = 0.0; double out[2] = { 0.0,0.0 };
    double r, count;
    for (int i = 0; i != sample_size; i++) {
      r = 1.0 / (i + 1);
      count = (t1[i] < t2[i]) ? 1.0 : 0.0;
      out[0] += (t1[i] - out[0])*r; out[1] += (t2[i] - out[1])*r; p += (count - p)*r;
    }
    Rprintf("\n%s\n", what.c_str());
    Rprintf("%12.4g%12.4g%12.4g\n", out[0], out[1], p);
    if (save_diagnose) {
      tests_out << std::endl;
      tests_out << what << std::endl;
      tests_out << std::setprecision(4);
      tests_out << std::setw(12) << out[0] << std::setw(12) << out[1] << std::setw(12) << p << std::endl;
    }
    double iv[2];
    gsl_vector_view vt1 = gsl_vector_view_array(t1, sample_size);
    gsl_vector_view vt2 = gsl_vector_view_array(t2, sample_size);
    gsl_vector_sub(&vt1.vector, &vt2.vector);
    gsl_sort(t1, 1, sample_size);
    hdi(sample_size, t1, .95, iv);
    Rprintf("95%% HDI\n");
    if (save_diagnose) tests_out << "95% HDI" << std::endl;
    for (int iq = 0; iq != 2; iq++) {Rprintf("%12.4g", iv[iq]);}
    Rprintf("\n");
    if (save_diagnose) {
      for (int iq = 0; iq != 2; iq++) {tests_out << std::setw(12) << iv[iq];} tests_out << std::endl;
    }
  }
  
  
  //computes path and category probabilities for one person
  void make_pij_for_individual(double* x, double* pij, double* pj) {
    // berechnet  pj, und pij fuer Individuum t, entspricht altem estimate
    
    double d0ha;
    
    for (int j = 0; j != kerncat; j++) {
      pj[j] = 0.0;
      for (int k = 0; k != branch[j]; k++) {
        dPIJ(j, k) = 1.0;
        for (int xr = 0; xr != dNDRIN(j, k); xr++) {
          int r = dDRIN(j, k, xr);
          int ix = dAR(j, k, r); int itree = cat2tree[j];
          //				int ia = dTREE_AND_NODE2PAR(itree, r, 0), iv = dTREE_AND_NODE2PAR(itree, r, 1), iw = dTREE_AND_NODE2PAR(itree, r, 2);
          int im = dTREE_AND_NODE2MAP(itree, r);
          d0ha = (ix > 0) ? x[im] : 1 - x[im];
          dPIJ(j, k) *= d0ha;
        }
        pj[j] += dPIJ(j, k);
      }
    }
    for (int j = 0; j != kerncat; j++)
      if (pj[j] != 0.0) for (int k = 0; k != branch[j]; k++) dPIJ(j, k) = dPIJ(j, k) / pj[j];
      else for (int k = 0; k != branch[j]; k++) dPIJ(j, k) = 1.0 / (1.0 * branch[j]);
  }
  
  //samples path from path probabilities
  int make_path_for_one_trial(int branchno, double* pij, gsl_rng* rst) {
    // int exit_status = 0;
    int help = 0; double temp;
    if (branchno > 1) {
      double u = oneuni(rst); temp = pij[help];
      while (u > temp) {
        help++;
        if (help > branchno - 1) { Rprintf("Achtung non-multinomial"); /*exit_status = -1;*/ }
        temp += pij[help];
      }
    }
    return(help);
    
  }
  
  //posterior tests category frequencies and mean category response times
  void aggregate(int n_all_parameters, int kerntree, int* idaten, const std::vector<trial>& daten, int* nks, int* jks, int* tree2cat, double* sample, gsl_rng* rst) {
    int exit_status = 1;
    int keig = kerncat * igroup; //kein = kerncat * indi, 
    
    std::vector<double> tavw(ifreemax * 3 * indi);
    
    std::vector<double> t1(sample_size, 0.0);
    std::vector<double> t2(sample_size, 0.0);
    std::vector<double> tt1(sample_size, 0.0);
    std::vector<double> tt2(sample_size, 0.0);
    
    std::vector<double> expe(kerncat);
    std::vector<int> rep(kerncat);
    
    std::vector<int> sobs(keig, 0);
    std::vector<double> sexp(keig);
    std::vector<int> srep(keig);
    
    //	double* tobs = 0; if (!(tobs = (double*)calloc(kerncat, sizeof(double)))) { Rprintf("Allocation failure\n"); exit_status = -1; }
    std::vector<double> texp(kerncat);
    std::vector<double> trep(kerncat);
    
    std::vector<double> stobs(keig, 0.0);
    std::vector<double> stexp(keig);
    std::vector<double> strep(keig);
    
    std::vector<double> pij(zweig * kerncat);
    std::vector<double> onepij(zweig);
    std::vector<double> x(no_patterns);
    
    std::vector<double> tdaten(indi * kerncat, 0.0);
    
    std::vector<int> nobs(keig, 0);
    std::vector<int> ntree(igroup * kerntree, 0);
    std::vector<int> nrep(keig);
    
    std::vector<double> d(kerncat);
    std::vector<double> x1(keig, 0.0);
    std::vector<double> xt1(keig, 0.0);
    std::vector<double> x2(keig * sample_size, 0.0);
    std::vector<double> xt2(keig * sample_size, 0.0);
    std::vector<unsigned int> drep(kerncat);
    
  #define dX1(IG,J) x1[IG*kerncat+J]
  #define dXT1(IG,J) xt1[IG*kerncat+J]
  #define dX2(IS,IG,J) x2[IS*igroup*kerncat+IG*kerncat+J]
  #define dXT2(IS,IG,J) xt2[IS*igroup*kerncat+IG*kerncat+J]
  #define dSOBS(IG,J) sobs[IG*kerncat+J]
  #define dSEXP(IG,J) sexp[IG*kerncat+J]
  #define dSREP(IG,J) srep[IG*kerncat+J]
  #define dSTOBS(IG,J) stobs[IG*kerncat+J]
  #define dSTEXP(IG,J) stexp[IG*kerncat+J]
  #define dSTREP(IG,J) strep[IG*kerncat+J]
  #define dNOBS(IG,J) nobs[IG*kerncat+J]
  #define dNTREE(IG,IT) ntree[IG*kerntree + IT]
  #define dNREP(IG,J) nrep[IG*kerncat+J]
  #define dTDATEN(T,J) tdaten[T*kerncat+J]
    
    
    for (int it = 0; it != datenzahl; it++) { trial one = daten[it]; dTDATEN(one.person, one.category) += one.rt / 1000.0; }
    
    for (int t = 0; t != indi; t++) {
      int ig = t2group[t];
      for (int j = 0; j != kerncat; j++) {
        dSOBS(ig, j) += dIDATEN(t, j);
        if (dIDATEN(t, j) > 0) {
          dTDATEN(t, j) /= dIDATEN(t, j);
          dSTOBS(ig, j) += dTDATEN(t, j); dNOBS(ig, j)++;
        }
      }
      for (int it = 0; it != kerntree; it++) dNTREE(ig, it) += dNKS(t, it);
    }
    
    for (int ig = 0; ig != igroup; ig++) for (int j = 0; j != kerncat; j++) dX1(ig, j) = dSOBS(ig, j) * 1.0 / dNTREE(ig, cat2tree[j]);
    for (int igc = 0; igc != keig; igc++) xt1[igc] = stobs[igc] / nobs[igc];
    
    ars_archiv ars_store;
    
    for (int is = 0; is != sample_size; is++) {
      for (int j = 0; j != keig; j++) { sexp[j] = 0.0; srep[j] = 0; stexp[j] = strep[j] = 0.0; nrep[j] = 0; }
      // compute exp pro Person und rep pro Person; aggregate, compute chi-square
      belege_ts(sample, is, tavw.data());
      bool change = (is==0)?true:dSAMPLE(is,0)!=dSAMPLE(is-1,0);
      if (change) {
        ars_store.hstore.clear(); ars_store.lowerstore.clear(); ars_store.startstore.clear(); ars_store.upperstore.clear(); ars_store.scalestore.clear(); ars_store.normstore.clear(); ars_store.sstore.clear();
        for (int t = 0; t != indi; t++) initialize_ars(t, tavw.data(), ars_store);
      }
      //		char xy; std::cin >> xy;
      for (int t = 0; t != indi; t++) {
        
        for (int im = 0; im != no_patterns; im++) {
          int ia = t * 3 * ifreemax + dCOMB(im, 0), iv = t * 3 * ifreemax + ifreemax + dCOMB(im, 1), iw = t * 3 * ifreemax + 2 * ifreemax + dCOMB(im, 2);
          x[im] = exp(logprob_upperbound(1, tavw[ia], tavw[iv], tavw[iw]));
        }
        make_pij_for_individual(x.data(), pij.data(), expe.data());
        
        for (int j = 0; j != kerncat; j++) {
          int r = cat2resp[j];
          double mu = dSAMPLE(is, irmuoff + t2group[t] * respno + r) + dSAMPLE(is, ilamoff + t * respno + r),
            sig = exp(dSAMPLE(is, isigoff + t));
          double z = -mu / sig;
          z = gsl_ran_tdist_pdf(0.0, degf) * degf * 1.0 / (degf - 1.0) / gsl_cdf_tdist_Q(z, degf) / pow(1.0 + gsl_pow_2(z) / degf, (degf - 1.0) / 2.0);
          texp[j] = mu + sig * z;
          
          for (int k = 0; k != branch[j]; k++) {
            double temp = 0.0;
            for (int xr = 0; xr != dNDRIN(j, k); xr++) {
              int r = dDRIN(j, k, xr);
              int itree = cat2tree[j];
              int ia = t * 3 * ifreemax + dTREE_AND_NODE2PAR(itree, r, 0), iv = t * 3 * ifreemax + ifreemax + dTREE_AND_NODE2PAR(itree, r, 1), iw = t * 3 * ifreemax + 2 * ifreemax + dTREE_AND_NODE2PAR(itree, r, 2);
              int pm = (dAR(j, k, r) > 0) ? 1 : 0;
              temp +=exp_mean(pm, tavw[ia], tavw[iv], tavw[iw]);
            }
            texp[j] += dPIJ(j, k) * temp;
          }
          
          if (gsl_isnan(texp[j]))	Rprintf("na dann");
          
          // compute exp und variance
        }
        
        for (int j = 0; j != kerncat; j++) rep[j] = 0;
        for (int it = 0; it != kerntree; it++) {
          int jksit = jks[it];
          for (int j = 0; j != jksit; j++) {
            d[j] = expe[dTREE2CAT(it, j)];
            drep[j] = 0;
          }
          gsl_ran_multinomial(rst, jks[it], dNKS(t, it), d.data(), drep.data());
          for (int j = 0; j != jksit; j++) rep[dTREE2CAT(it, j)] = drep[j];
        }
        int ig = t2group[t];
        for (int j = 0; j != kerncat; j++)  expe[j] = dNKS(t, cat2tree[j]) * expe[j];
        for (int j = 0; j != kerncat; j++) { dSEXP(ig, j) += expe[j];  dSREP(ig, j) += rep[j]; dX2(is, ig, j) += (rep[j] * 1.0); }
        
        //go on with times
        for (int j = 0; j != kerncat; j++) {
          int r = cat2resp[j];
          double mu = dSAMPLE(is, irmuoff + t2group[t] * respno + r) + dSAMPLE(is, ilamoff + t * respno + r),
            sig = exp(dSAMPLE(is, isigoff + t));
          trep[j] = 0.0;
          int repj = rep[j];
          for (int ir = 0; ir != repj; ir++) {
            double temp = 0.0;
            int brj = branch[j];
            for (int k = 0; k != brj; k++) onepij[k] = dPIJ(j, k);
            int ipath = make_path_for_one_trial(branch[j], onepij.data(), rst);
            int ndjip = dNDRIN(j, ipath);
            for (int xr = 0; xr != ndjip; xr++) {
              int r = dDRIN(j, ipath, xr);
              int itree = cat2tree[j];
              int ia = dTREE_AND_NODE2PAR(itree, r, 0), iv = dTREE_AND_NODE2PAR(itree, r, 1), iw = dTREE_AND_NODE2PAR(itree, r, 2);
              int m = dMAP(ia, iv, iw);
              ia = t * 3 * ifreemax + ia; iv = t * 3 * ifreemax + ifreemax + iv; iw = t * 3 * ifreemax + 2 * ifreemax + iw;
              int pm = (dAR(j, ipath, r) > 0) ? 1 : 0;
              temp += make_rwiener(t, m, pm, ars_store, GSL_POSINF, tavw[ia], tavw[iv], tavw[iw], rst);
            }
            double xtemp = 0;
            do {
              xtemp = mu + sig * gsl_ran_tdist(rst, degf);
            } while (xtemp < 0);
            temp += xtemp;
            trep[j] += temp;
          }
          if (rep[j] > 0) trep[j] /= rep[j]; else trep[j] = 0.0;
        }
        
        for (int j = 0; j != kerncat; j++) {
          dSTEXP(ig, j) += texp[j];
          if (trep[j] > 0) { dSTREP(ig, j) += trep[j]; dNREP(ig, j)++; }
        }
        
        
      }
      for (int ig = 0; ig != igroup; ig++) for (int j = 0; j != kerncat; j++)
      {
        t1[is] += gsl_pow_2(dSOBS(ig, j) - dSEXP(ig, j)) / dSEXP(ig, j);
        t2[is] += gsl_pow_2(dSREP(ig, j) - dSEXP(ig, j)) / dSEXP(ig, j);
        dX2(is, ig, j) /= dNTREE(ig, cat2tree[j]);
        
        dSTEXP(ig, j) /= ng[ig];
        if (dNREP(ig, j) > 0) dSTREP(ig, j) /= dNREP(ig, j);
        dXT2(is, ig, j) = dSTREP(ig, j);
      }
      
      for (int j = 0; j != keig; j++) {
        tt1[is] += gsl_pow_2(xt1[j] - stexp[j]) / stexp[j]; tt2[is] += gsl_pow_2(strep[j] - stexp[j]) / stexp[j];
      }
    }
    
    test(t1.data(), t2.data(), "Posterior predictive checks: frequencies");
    
    double qv[5];
    std::ofstream meansout(MEANSOUT);
    
    //std::cout << std::setprecision(4);
    for (int ig = 0; ig != igroup; ig++) {
      for (int j = 0; j != kerncat; j++) {
        Rprintf("%3d%12.4g", j, dX1(ig,j));
        int old_new = j % 2;
        int tartype = (j / 2) % 3;
        int bias = (j / 6);
        int correct = 0; correct = ((tartype == 0) && (old_new == 1)) || ((tartype > 0) && (old_new == 0));
        meansout << std::setw(3) << old_new << std::setw(3) << tartype << std::setw(3) << bias << std::setw(12) << dX1(ig, j);
        for (int is = 0; is != sample_size; is++) t2[is] = dX2(is, ig, j);
        gsl_sort(t2.data(), 1, sample_size);
        qv[2] = gsl_stats_median_from_sorted_data(t2.data(), 1, sample_size);
        double iv[2]; hdi(sample_size, t2.data(), 0.95, iv); qv[1] = iv[0]; qv[3] = iv[1]; hdi(sample_size, t2.data(), 0.99, iv); qv[0] = iv[0]; qv[4] = iv[1];
        Rprintf("%12.4g%12.4g%12.4g\n", qv[1], qv[2], qv[3]);
        meansout << std::setw(12) << qv[1] << std::setw(12) << qv[2] << std::setw(12) << qv[3] << std::setw(3) << correct << std::endl;
      }
    }
      
    // Das Ganze fuer die Zeiten
      
    test(tt1.data(), tt2.data(), "Posterior predictive checks: latencies");
    
    // int nq = 1;
    
    //std::cout << std::setprecision(4);
    for (int ig = 0; ig != igroup; ig++) {
      for (int j = 0; j != kerncat; j++) {
        Rprintf("%3d%12.4g", j, dX1(ig, j));
        int old_new = j % 2;
        int tartype = (j / 2) % 3;
        int bias = (j / 6);
        int correct = 0; correct = ((tartype == 0) && (old_new == 1)) || ((tartype > 0) && (old_new == 0));
        meansout << std::setw(3) << old_new << std::setw(3) << tartype << std::setw(3) << bias << std::setw(12) << dXT1(ig, j);
        for (int is = 0; is != sample_size; is++) t2[is] = dXT2(is, ig, j);
        gsl_sort(t2.data(), 1, sample_size);
        qv[2] = gsl_stats_median_from_sorted_data(t2.data(), 1, sample_size);
        double iv[2]; hdi(sample_size, t2.data(), 0.95, iv); qv[1] = iv[0]; qv[3] = iv[1]; hdi(sample_size, t2.data(), 0.99, iv); qv[0] = iv[0]; qv[4] = iv[1];
        Rprintf("%12.4g%12.4g%12.4g\n", qv[1], qv[2], qv[3]);
        meansout << std::setw(12) << qv[1] << std::setw(12) << qv[2] << std::setw(12) << qv[3] << std::setw(3) << correct << std::endl;
      }
    }
    meansout.close();
    
  }
  
  
  //contrasts between two groups in population means
  void groupwise(double* sample) {
    int exit_status = 0;
    
    
    std::vector<double> t1(sample_size);
    std::vector<double> t2(sample_size);
    
    for (int ip = 0; ip != ifreeg; ip++)
    {
      for (int is = 0; is != sample_size; is++) {
        //				par.push_back(gsl_cdf_ugaussian_P(dSAMPLE(is,ip))-gsl_cdf_ugaussian_P(dSAMPLE(is,ip+ifree)));
        t2[is] = logit(avwtrans[is_type(ip)], dSAMPLE(is, ip));
        t1[is] = logit(avwtrans[is_type(ip)], dSAMPLE(is, ip + ifreeg));
      }
      test(t1.data(), t2.data(), "group-tests mu");
    }
    
    for (int ir = 0; ir != respno; ir++)
    {
      for (int is = 0; is != sample_size; is++) {
        //				par.push_back(gsl_cdf_ugaussian_P(dSAMPLE(is,ip))-gsl_cdf_ugaussian_P(dSAMPLE(is,ip+ifree)));
        t2[is] = 1000 * dSAMPLE(is, irmuoff + ir);
        t1[is] = 1000 * dSAMPLE(is, irmuoff + ir + respno);
      }
      test(t1.data(), t2.data(), "group-tests mu");
    }
  }
  
  
  
  //main
  void diagnosis(const std::vector<trial> & daten, int *idaten, int kerntree, gsl_rng *rst) {
    std::vector<int> nks(indi * kerntree);
    std::vector<int> jks(kerntree);
    std::vector<int> tree2cat(kerntree * kerncat);
    
    
n_all_parameters = icompg*igroup + indi * icompg + (icompg * (icompg + 1)) / 2 + respno*igroup + (respno + 1) * indi + (respno * (respno + 1)) / 2 + 1;
    //                    ma,mv,mw  a,v,w            sig                           rmu     lambdas+sig_t          gam            omega
    
    std::vector<double> sample;
    lies_sample(n_all_parameters, sample);
    
    if (save_diagnose) tests_out.open(TESTSOUT);
    quantiles(daten, n_all_parameters, sample.data());
    
    // make nks
    for (int t = 0; t != indi; t++) for (int it = 0; it != kerntree; it++) dNKS(t, it) = 0;
    for (int t = 0; t != indi; t++) for (int j = 0; j != kerncat; j++) dNKS(t, cat2tree[j]) += dIDATEN(t, j);
    for (int it = 0; it != kerntree; it++) jks[it] = 0;
    for (int j = 0; j != kerncat; j++) { dTREE2CAT(cat2tree[j], jks[cat2tree[j]]) = j; jks[cat2tree[j]]++; }
    aggregate(n_all_parameters, kerntree, idaten, daten, nks.data(), jks.data(), tree2cat.data(), sample.data(), rst);
    
    if (DIC) dic(daten, sample.data());
    
    if (save_diagnose) tests_out.close();
    
  }

}
