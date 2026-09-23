#include "rts.h"
#include "gauss.h"

#include <cassert>
#include <atomic>

std::atomic<bool> cancel_flag{false};

void check_user_interrupt(void*) {
  R_CheckUserInterrupt();
}


//computes log(exp(xa) + exp(xb))
double logsum(double xa, double xb) {
  double temp;
  if (xa <= GSL_NEGINF) return xb;
  if (xb <= GSL_NEGINF) return xa;
  if (xa > xb) temp = xa + gsl_log1p(exp(xb - xa));
  else temp = xb + gsl_log1p(exp(xa - xb));
  return temp;
}

//generates one uniform variate
double oneuni(gsl_rng *rst) {
  return gsl_rng_uniform_pos(rst);
}

//generates one standard normal variate
double onenorm(gsl_rng *rst) {
  return gsl_ran_ugaussian(rst);
}

//truncnorm computes normal variate with sd = 1, mean = b, and value >= 0
double truncnorm(double b, gsl_rng* rst) {
  double temp;
  if (b >= 0.0) {
    do  temp = onenorm(rst);  while (temp < -b);
    temp += b;
  }
  else temp = gsl_ran_ugaussian_tail(rst, -b) + b;
  return temp;
}

int fopen_s(FILE **f, const char *name, const char *mode) {
  int ret = 0;
  assert(f);
  *f = fopen(name, mode);
  
  if (!*f)
    ret = errno;
  return ret;
}


namespace drtmpt {
  
  // compute log(1-exp(z))
  double log1pem1(double z) {
    //	if (z > 0) std::cout << "log1pem1";
    if (fabs(z) < 1.0e-2) return log(-gsl_expm1(z));
    else return gsl_log1p(-exp(z));
  }
  
  // computes log(exp(xa) - exp(xb)); xa should be larger than xb
  double logdiff(double xa, double xb) {
    // double result;
    if (xb >= xa) {
      //		std::cout << "logdiff";
      return(GSL_NEGINF);
    }
    if (xb <= GSL_NEGINF) return(xa);
    return xa + log1pem1(xb - xa);
  }
  
  //store parameters between parallelized sampling blocks of size ireps
  void push(int ithread, int n_value_store, int n_all_parameters, gsl_vector* hampar, double* tavw, double* tlams, double* ai, double* loglambdas, double* bi,  double* alltaus, double* rest,
            int trialno, int* paths, int* nips,  double liknorm[6], double activeeps, double epsm, double Hobjective, double* valuestore, double* parmon, double* parmonstore) {
    int offset = ithread * n_value_store;
    
    gsl_vector_view t1 = gsl_vector_view_array(valuestore, NOTHREADS * n_value_store);
    gsl_vector_view t2 = gsl_vector_subvector(&t1.vector, offset, (phase <= 2) ? nhamil : n_all_parameters);
    gsl_vector_memcpy(&t2.vector, hampar);
    offset += n_all_parameters;
    for (int type = 0; type != 3; type++) {
      int ift = ifree[type];
      for (int i = 0; i != ift; i++) if (dCOMP(type, i)) {
        for (int t = 0; t != indi; t++) {
          int iavw = t * 3 * ifreemax + type * ifreemax + i;
          valuestore[offset++] = tavw[iavw];
        }
      }
    }
    
    for (int i = 0; i!= icompg; i++) valuestore[offset++] = ai[i];
    
    for (int i = 0; i != indi; i++) valuestore[offset++] = loglambdas[i];
    
    for (int i = 0; i != respno; i++) valuestore[offset++] = bi[i];
    
    for (int i = 0; i!= respno * indi; i++) valuestore[offset++] = tlams[i];
    
    for (int i = 0; i != trialno; i++) valuestore[offset++] = paths[i];
    
    int i2n = indi * 2 * no_patterns;
    for (int i = 0; i != i2n; i++) valuestore[offset++] = nips[i];
    
    for (int i = 0; i != 6; i++) valuestore[offset++] = liknorm[i];
    
    gsl_vector_view t4 = gsl_vector_subvector(&t1.vector, offset, ntau);
    gsl_vector_view t5 = gsl_vector_view_array(alltaus, ntau);
    gsl_vector_memcpy(&t4.vector, &t5.vector);
    offset += ntau;
    gsl_vector_view t7 = gsl_vector_subvector(&t1.vector, offset, datenzahl);
    gsl_vector_view t8 = gsl_vector_view_array(rest, datenzahl);
    gsl_vector_memcpy(&t7.vector, &t8.vector);
    offset += datenzahl;
    
    valuestore[offset++] = activeeps;
    valuestore[offset++] = epsm;
    valuestore[offset++] = Hobjective;
    
    
    int znall = 2 * n_all_parameters;
    offset = ithread * znall;
    for (int i = 0; i != znall; i++) parmonstore[offset + i] = parmon[i];
    offset = 0;
  }
  
  //retrieve parameters between parallelized sampling blocks of size ireps
  void pop(int ithread, int n_value_store, int n_all_parameters, gsl_vector* hampar, double* tavw, double* tlams, double* ai, double* loglambdas, double* bi,  double* alltaus, double* rest,
           int trialno, int* paths, int* nips, double liknorm[6], double &activeeps, double &epsm, double &Hobjective, double* valuestore, double* parmon, double* parmonstore) {
    
    int offset = ithread * n_value_store;
    
    gsl_vector_view t1 = gsl_vector_view_array(valuestore, NOTHREADS * n_value_store);
    gsl_vector_view t2 = gsl_vector_subvector(&t1.vector, offset, (phase <= 2) ? nhamil : n_all_parameters);
    gsl_vector_memcpy(hampar, &t2.vector);
    offset += n_all_parameters;
    for (int type = 0; type != 3; type++) {
      int ift = ifree[type];
      for (int i = 0; i != ift; i++) if (dCOMP(type, i)) {
        for (int t = 0; t != indi; t++) {
          int iavw = t * 3 * ifreemax + type * ifreemax + i;
          tavw[iavw] = valuestore[offset++];
        }
      }
      else for (int t = 0; t != indi; t++) dTAVW(t, type, i) = dCONSTS(type, i);
    }
    
    for (int i = 0; i != icompg; i++) ai[i] = valuestore[offset++];
    
    for (int i = 0; i != indi; i++) loglambdas[i] = valuestore[offset++];
    
    for (int i = 0; i != respno; i++) bi[i] = valuestore[offset++];
    
    for (int i = 0; i != indi * respno; i++) tlams[i] = valuestore[offset++];
    
    for (int i = 0; i != trialno; i++) paths[i] = static_cast<int>(valuestore[offset++]);
    
    int i2n = indi * 2 * no_patterns;
    for (int i = 0; i != i2n; i++) nips[i] = static_cast<int>(valuestore[offset++]);
    
    for (int i = 0; i!= 6; i++) liknorm[i] = valuestore[offset++];
    
    gsl_vector_view t4 = gsl_vector_subvector(&t1.vector, offset, ntau);
    gsl_vector_view t5 = gsl_vector_view_array(alltaus, ntau);
    gsl_vector_memcpy(&t5.vector, &t4.vector);
    offset += ntau;
    gsl_vector_view t7 = gsl_vector_subvector(&t1.vector, offset, datenzahl);
    gsl_vector_view t8 = gsl_vector_view_array(rest, datenzahl);
    gsl_vector_memcpy(&t8.vector, &t7.vector);
    offset += datenzahl;
    
    activeeps = valuestore[offset++];
    epsm = valuestore[offset++];
    Hobjective = valuestore[offset++];
    
    
    int znall = 2 * n_all_parameters;
    offset = ithread * znall;
    for (int i = 0; i != znall; i++) parmon[i] = parmonstore[offset + i];
    offset = 0;
  }
  
  
  //show interim results after sampling blocks of size ireps for one thread
  void on_screen3(int n_all_parameters, double* xwbr, double* parmon, double* consts, double rmax, int imax, int irun, int update_block, int previous_samples) {
    
    int jz;
    Rprintf("\nThresholds\n");
    Rprintf("estim:");
    for (int ig = 0; ig != igroup; ig++) {
      for (int jp = 0; jp != kernpar; jp++) {
        jz = dKERN2FREE(0, jp);
        if (dCOMP(0, jz)) Rprintf((jp==0 ? (ig==0 ? "%11g" : "%17g") : "%15g"), logit(avwtrans[0], dPARMON(1, ig * icompg + dFREE2COMP(0, jz)))); else Rprintf((jp==0 ? (ig==0 ? "%11g" : "%17g") : "%15g"), dCONSTS(0, jz));
      }
      Rprintf("\n");
    }
    Rprintf("Rhat:");
    for (int ig = 0; ig != igroup; ig++) {
      for (int jp = 0; jp != kernpar; jp++) {
        jz = dKERN2FREE(0, jp);
        if (dCOMP(0, jz)) Rprintf((jp==0 ? (ig==0 ? "%12g" : "%18g") : "%15g"), dXWBR(3, ig * icompg + dFREE2COMP(0, jz))); else Rprintf((jp==0 ? (ig==0 ? "%12g" : "%18g") : "%15g"), 0.0);
      }
      Rprintf("\n");
    }
    Rprintf("--------\n");
    
    Rprintf("Drift\n");
    Rprintf("estim:");
    for (int ig = 0; ig != igroup; ig++) {
      for (int jp = 0; jp != kernpar; jp++) {
        jz = dKERN2FREE(1, jp);
        if (dCOMP(1, jz)) Rprintf((jp==0 ? (ig==0 ? "%11g" : "%17g") : "%15g"), logit(avwtrans[1], dPARMON(1, ig * icompg + icomp[0]  + dFREE2COMP(1, jz)))); else Rprintf((jp==0 ? (ig==0 ? "%11g" : "%17g") : "%15g"), dCONSTS(1, jz));
      }
      Rprintf("\n");
    }
    Rprintf("Rhat:");
    for (int ig = 0; ig != igroup; ig++) {
      for (int jp = 0; jp != kernpar; jp++) {
        jz = dKERN2FREE(1, jp);
        if (dCOMP(1, jz)) Rprintf((jp==0 ? (ig==0 ? "%12g" : "%18g") : "%15g"), dXWBR(3, ig * icompg   + icomp[0] + dFREE2COMP(1, jz))); else Rprintf((jp==0 ? (ig==0 ? "%12g" : "%18g") : "%15g"), 0.0);
      }
      Rprintf("\n");
    }
    Rprintf("--------\n");
    
    Rprintf("Bias\n");
    Rprintf("estim:");
    for (int ig = 0; ig != igroup; ig++) {
      for (int jp = 0; jp != kernpar; jp++) {
        jz = dKERN2FREE(2, jp);
        if (dCOMP(2, jz)) Rprintf((jp==0 ? (ig==0 ? "%11g" : "%17g") : "%15g"), logit(avwtrans[2], dPARMON(1, ig * icompg + icomp[0] + icomp[1] + dFREE2COMP(2, jz)))); else Rprintf((jp==0 ? (ig==0 ? "%11g" : "%17g") : "%15g"), dCONSTS(2, jz));
      }
      Rprintf("\n");
    }
    Rprintf("Rhat:");
    for (int ig = 0; ig != igroup; ig++) {
      for (int jp = 0; jp != kernpar; jp++) {
        jz = dKERN2FREE(2, jp);
        if (dCOMP(2, jz)) Rprintf((jp==0 ? (ig==0 ? "%12g" : "%18g") : "%15g"), dXWBR(3, ig * icompg + icomp[0] + icomp[1]  + dFREE2COMP(2, jz)));  else Rprintf((jp==0 ? (ig==0 ? "%12g" : "%18g") : "%15g"), 0.0);
      }
      Rprintf("\n");
    }
    Rprintf("--------\n");
    
    int js = jz = irmuoff;
    Rprintf("Motor-Time Means\n");
    Rprintf("estim:");
    for (int ig = 0; ig != igroup; ig++) {
      for (int ir = 0; ir != respno; ir++) {
        Rprintf((ir==0 ? (ig==0 ? "%11g" : "%17g") : "%15g"), dPARMON(1, jz)); jz++;
      }
      Rprintf("\n");
    }
    jz = js;
    Rprintf("Rhat:");
    for (int ig = 0; ig != igroup; ig++) {
      for (int ir = 0; ir != respno; ir++) {
        Rprintf((ir==0 ? (ig==0 ? "%12g" : "%18g") : "%15g"), dXWBR(3, jz)); jz++;
      }
      Rprintf("\n");
    }
    Rprintf("--------\n");
    
    Rprintf("Omega-Square\n");
    Rprintf("estim:");
    Rprintf("%11g", exp(dPARMON(1, n_all_parameters - 1)));
    Rprintf("\n");
    Rprintf("Rhat:");
    Rprintf("%12g", dXWBR(3, n_all_parameters - 1));
    Rprintf("\n");
    Rprintf("------------------------\n");
    
    if (goon) {
      //continuation mode (update_drtmpt): sampling starts immediately in phase 4,
      //so progress is based on blocks saved during this update call
      //(update_block = ioff + 1; ioff is incremented after each saved block),
      //not on the convergence counter used by fit_drtmpt
      const int update_iterations = update_block * ireps;
      const int actual_samples = previous_samples + update_iterations / THIN;
      double pct_temp = 100.0 * update_iterations / (1.0 * THIN * ADDITION / NOTHREADS);
      if (pct_temp > 100.0) pct_temp = 100.0;
      Rprintf("max(Rhats): %12g\n", rmax);
      Rprintf("     Phase: %10d/4\n", phase);
      Rprintf("Iterations: %12d [sampling: %g%%]\n", update_iterations, pct_temp);
      Rprintf("Actual samples:%9d\n", actual_samples);
    } else {
      if (rmax < RMAX && phase == 4) RMAX_reached += 1;
      else RMAX_reached = 0;
      double pct_temp = (RMAX_reached>1) ? (100.0*ireps*(RMAX_reached-1)/(1.0*(THIN*SAMPLE_SIZE/NOTHREADS))) : 0.0;

      Rprintf("max(Rhats): %12g\n", rmax);
      //std::cout << std::setw(15) << "rmax " << std::setw(15) << rmax << std::endl;
      Rprintf("     Phase: %10d/4\n", phase);
      //std::cout << std::setw(15) << "Phase " << std::setw(15) << phase << std::endl;
      if (phase == 4) Rprintf("Iterations: %12d [sampling: %g%%]\n", (irun + 1)*ireps, pct_temp);
      else            Rprintf("Iterations: %12d\n", (irun + 1) * ireps);
      if (RMAX_reached == 0 && phase == 4) Rprintf("Sampling starts when max(Rhats)<%g\n", RMAX);
      if (RMAX_reached == 1 && phase == 4) Rprintf("Sampling starts now.\n");
    }
    //std::cout << std::setw(15) << "Iterationen " << std::setw(15) << (irun + 1) * ireps << std::endl;
      
    Rprintf("__");
    if (kernpar > respno) {
      for (int i = 0; i < kernpar; i++) Rprintf("_______________");
    } else {
      for (int i = 0; i < respno; i++) Rprintf("_______________");
      // Rprintf("_______________");
    }
    Rprintf("\n");
    
    // int m = (irun + 1) * ireps;
    
  } // end on_screen3
  
  //go from mavw and avw to vector hampar
  void make_hampar_avw(double* mavw, double* avw, gsl_vector* hampar) {
    int jj = 0;
    for (int ig = 0; ig != igroup; ig++) {
      for (int type = 0; type != 3; type++) {
        int ift = ifree[type];
        for (int ip = 0; ip != ift; ip++)
          if (dCOMP(type, ip)) {
            gsl_vector_set(hampar, jj, dMAVW(ig, type, ip)); jj++;
          }
      }
    }
    for (int t = 0; t != indi; t++) {
      for (int type = 0; type != 3; type++) {
        int ift = ifree[type];
        for (int ip = 0; ip != ift; ip++)
          if (dCOMP(type, ip)) {
            gsl_vector_set(hampar, jj, dAVW(t, type, ip)); jj++;
          }
      }
    }
  }
  
  //go from vector hampar to mavw and avw
  void inv_make_hampar_avw(double* mavw, double* avw, gsl_vector* hampar) {
    int jj = 0;
    for (int ig = 0; ig != igroup; ig++) {
      for (int type = 0; type != 3; type++) {
        int ift = ifree[type];
        for (int ip = 0; ip != ift; ip++)
          if (dCOMP(type, ip)) {
            dMAVW(ig, type, ip) = gsl_vector_get(hampar, jj); jj++;
          }
      }
    }
    for (int t = 0; t != indi; t++) {
      for (int type = 0; type != 3; type++) {
        int ift = ifree[type];
        for (int ip = 0; ip != ift; ip++)
          if (dCOMP(type, ip)) {
            dAVW(t, type, ip) = gsl_vector_get(hampar, jj); jj++;
          }
      }
    }
  }
  
  //go from rmu und labmda to hampar
  void make_hampar_rmu_lambda(double* rmu, double* lambda, gsl_vector* hampar) {
    int igrre = igroup * respno, indreindi = indi * respno + indi;
    for (int ig = 0; ig != igrre; ig++) gsl_vector_set(hampar, irmuoff + ig, rmu[ig]);
    for (int ir = 0; ir != indreindi; ir++) gsl_vector_set(hampar, ilamoff + ir, lambda[ir]);
  }
  
  //go from hampar to rmu and lambda
  void inv_make_hampar_rmu_lambda(double* rmu, double* lambda, gsl_vector* hampar) {
    int igrre = igroup * respno, indreindi = indi * respno + indi;
    for (int ig = 0; ig != igrre; ig++) rmu[ig] = gsl_vector_get(hampar, irmuoff + ig);
    for (int ir = 0; ir != indreindi; ir++) lambda[ir] = gsl_vector_get(hampar, ilamoff + ir);
  }
  
  
  
  //initialize without maximum-likelihood estimation per individual
  void initialize_new0(std::vector <trial> daten, gsl_vector* hampar, double* tavw, double* tlams, double* ai, double* loglambdas, double* bi, int* paths, gsl_rng* rst) {
    
    for (int x = 0; x != datenzahl; x++) {
      paths[x] = gsl_rng_uniform_int(rst, branch[daten[x].category]);
    }

    std::vector<double> mavw(ifreemax * 3 * igroup);
    std::vector<double> avw(ifreemax * 3 * indi, 0.0);
    std::vector<double> rmu(respno * igroup);
    std::vector<double> lambdas(indi * (respno + 1));

    for (int type = 0; type != 3; type++) {
      int ift = ifree[type];
      for (int i = 0; i != ift; i++) {
        for (int ig = 0; ig != igroup; ig++) dMAVW(ig, type, i) = 0.0;
        for (int t = 0; t != indi; t++) {
          dAVW(t, type, i) = (dCOMP(type, i)) ?  onenorm(rst) : 0.0;
          dMAVW(t2group[t], type, i) += dAVW(t, type, i);
        }
        for (int ig = 0; ig != igroup; ig++) dMAVW(ig, type, i) = (dCOMP(type, i)) ? dMAVW(ig, type, i) / ng[ig] + 1.0 / sqrt(ng[ig]) * onenorm(rst) : invlogit(avwtrans[type], dCONSTS(type, i));
        for (int t = 0; t != indi; t++) {
          if (dCOMP(type, i))  dAVW(t, type, i) = (dAVW(t, type, i) - dMAVW(t2group[t], type, i));
          dTAVW(t, type, i) = logit(avwtrans[type], dMAVW(t2group[t], type, i) + dAVW(t, type, i));
        }
      }
    }

    std::vector<double> xb(datenzahl);
    for (int i = 0; i != respno; i++) {
      for (int ig = 0; ig != igroup; ig++) rmu[ig * respno + i] = 0.0;
      for (int t = 0; t != indi; t++) {
        int jj = 0;
        for (int x = 0; x != datenzahl; x++) if ((daten[x].person == t) && (cat2resp[daten[x].category] == i)) xb[jj++] = daten[x].rt / 1000.0;
        double temp = gsl_stats_mean(xb.data(), 1, jj);
        lambdas[t * respno + i] = 0.8 * temp + 0.01 * onenorm(rst);
        
        rmu[t2group[t] * respno + i] += lambdas[t * respno + i];
      }
      for (int ig = 0; ig != igroup; ig++) rmu[ig * respno + i] = rmu[ig * respno + i] / ng[ig] + 0.01 * onenorm(rst);
      for (int t = 0; t != indi; t++) {
        lambdas[t * respno + i] -= rmu[t2group[t] * respno + i];
      }
    }
    for (int t = 0; t != indi; t++) {
      int jj = 0;
      for (int x = 0; x != datenzahl; x++) if (daten[x].person == t) xb[jj++] = daten[x].rt / 1000.0;
      lambdas[indi * respno + t] = 0.9 * gsl_stats_sd(xb.data(), 1, jj) + 0.001 * oneuni(rst);
      loglambdas[t] = log(lambdas[indi * respno + t]);
    }
    for (int i = 0; i != icompg; i++) ai[i] = 1.0;
    for (int r = 0; r != respno; r++)  bi[r] = 1.0;
    
    for (int t = 0; t != indi; t++) for (int r = 0; r != respno; r++) tlams[t * respno + r] = rmu[t2group[t] * respno + r] + lambdas[t * respno + r];
    
    make_hampar_avw(mavw.data(), avw.data(), hampar);
    make_hampar_rmu_lambda(rmu.data(), lambdas.data(), hampar);
  }



  //initialize after maximum-likelihood estimation per individual
  void initialize_new1(std::vector <trial> daten, gsl_vector* hampar, double* tavw, double* tlams, double* ai, double* loglambdas, double* bi, int* paths, gsl_rng* rst) {
    
    std::vector<double> xx(no_patterns);
    std::vector<double> pj(kerncat);
    std::vector<double> pij(zweig * kerncat);
    std::vector<double> ppath(zweig);

    std::vector<double> mavw(ifreemax * 3 * igroup);
    std::vector<double> avw(ifreemax * 3 * indi, 0.0);
    std::vector<double> rmu(respno * igroup);
    std::vector<double> lambdas(indi * (respno + 1));
    
    inv_make_hampar_avw(mavw.data(), avw.data(), hampar);
    inv_make_hampar_rmu_lambda(rmu.data(), lambdas.data(), hampar);
    
    for (int x = 0; x != static_cast<int>(daten.size()); x++) {
      paths[x] = gsl_rng_uniform_int(rst, branch[daten[x].category]);
    }
    
    for (int type = 0; type != 3; type++) {
      int ift = ifree[type];
      for (int i = 0; i != ift; i++) {
        for (int ig = 0; ig != igroup; ig++) dMAVW(ig, type, i) = 0.0;
        for (int t = 0; t != indi; t++) {
          dAVW(t, type, i) = (dCOMP(type, i)) ? dAVW(t, type, i) : 0.0;
          if (dAVW(t, type, i) < -3.0) dAVW(t, type, i) = -3.0;
          if (dAVW(t, type, i) > 3.0) dAVW(t, type, i) = 3.0;
          dMAVW(t2group[t], type, i) += dAVW(t, type, i);
        }
        for (int ig = 0; ig != igroup; ig++) dMAVW(ig, type, i) = (dCOMP(type, i)) ? (dMAVW(ig, type, i) + onenorm(rst)) / (ng[ig] + 1) : invlogit(avwtrans[type], dCONSTS(type, i));
        for (int t = 0; t != indi; t++) {
          if (dCOMP(type, i)) {
            dAVW(t, type, i) = (dAVW(t, type, i) - dMAVW(t2group[t], type, i)) + 0.1 * onenorm(rst);
            dTAVW(t, type, i) = logit(avwtrans[type], dMAVW(t2group[t], type, i) + dAVW(t, type, i));
          }
          else dTAVW(t, type, i) = dCONSTS(type, i);
        }
      }
    }
    
    for (int t = 0; t != indi; t++) {
      for (int im = 0; im != no_patterns; im++) {
        int ia = t * 3 * ifreemax + dCOMB(im, 0), iv = t * 3 * ifreemax + ifreemax + dCOMB(im, 1), iw = t * 3 * ifreemax + 2 * ifreemax + dCOMB(im, 2);
        xx[im] = exp(logprob_upperbound(1, tavw[ia], tavw[iv], tavw[iw]));
      }
      make_pij_for_individual(xx.data(), pij.data(), pj.data());
      for (int x = 0; x != datenzahl; x++) if (daten[x].person == t) {
        int j = daten[x].category;
        int bj = branch[j];
        for (int k = 0; k != bj; k++) ppath[k] = dPIJ(j, k);
        if (branch[j] > 1) {
          std::vector<unsigned int> nn(zweig);
          gsl_ran_multinomial(rst, branch[j], 1, ppath.data(), nn.data());
          for (int k = 0; k != bj; k++) if (nn[k] > 0) { paths[x] = k; break;}
        }
        else paths[x] = 0;
      }
    }
    
    
    std::vector<double> temp_rest(2 * indi);
    for (int i = 0; i != indi; i++) { temp_rest[i] = lambdas[i * respno]; }
    for (int i = 0; i != indi; i++) { temp_rest[indi + i] = lambdas[indi * respno + i]; }
    for (int i = 0; i != indi; i++) {
      if (temp_rest[i] < -1.0) temp_rest[i] = -1.0;
      if (temp_rest[i] > 0.7) temp_rest[i] = 0.7;
      if (temp_rest[indi + i] > 0.3) temp_rest[indi + i] = 0.3;
    }
    
    for (int ig = 0; ig != igroup; ig++) {
      int ng = 0; double rsmu = 0.0;
      for (int t = 0; t != indi; t++) if (t2group[t] == ig) { rsmu += temp_rest[t]; ng++; }
      for (int ir = 0; ir != respno; ir++) rmu[ig * respno + ir] = (rsmu + 0.7 * oneuni(rst)) / (ng + 1);
    }
    for (int t = 0; t != indi; t++) for (int r = 0; r != respno; r++) lambdas[t * respno + r] = temp_rest[t] - rmu[t2group[t] * respno + r] + 0.01 * onenorm(rst);
    for (int t = 0; t != indi; t++) {
      lambdas[indi * respno + t] = (indi - 1) * 1.0 / indi * temp_rest[t + indi] + 1.0 / indi * 0.1 * oneuni(rst);
      loglambdas[t] = log(lambdas[indi * respno + t]);
    }

    for (int i = 0; i != icompg; i++) ai[i] = 1.0;
    for (int r = 0; r != respno; r++)  bi[r] = 1.0;
    
    for (int t = 0; t != indi; t++) for (int r = 0; r != respno; r++) tlams[t * respno + r] = rmu[t2group[t] * respno + r] + lambdas[t * respno + r];
    
    make_hampar_avw(mavw.data(), avw.data(), hampar);
    make_hampar_rmu_lambda(rmu.data(), lambdas.data(), hampar);
    
  }
  
  void initialize(int flag, const std::vector<trial> & daten, double xeps, double* parmonstore, int n_value_store, double* valuestore, std::vector<gsl_rng*>& rsts) {
    //flag = 0 initialize random; flag = 1 initialize with max lik
    gsl_rng* xst;   xst = gsl_rng_alloc(T_rng);
    double liknorm[6] = { 0 * 6 };
    
    gsl_vector* hampar = gsl_vector_alloc(nhamil);
    std::vector<double> tavw(ifreemax * 3 * indi);
    std::vector<double> tlams(indi * respno);
    std::vector<double> loglambdas(indi);
    std::vector<int> paths(datenzahl);
    std::vector<int> nips(no_patterns * 2 * indi);
    std::vector<double> ai(icompg);
    std::vector<double> bi(respno);
    double epsm = 0.0; double activeeps = xeps; double Hobjective = 0.0;
    std::vector<double> avw_temp;
    std::vector<double> lambdas_temp;
    std::vector<double> parmon(2 * n_all_parameters, 0.0);
    std::vector<double> alltaus(ntau, 0.0);
    std::vector<double> rest(datenzahl);
    
    if (flag == 1) {
      avw_temp.resize(ifreemax * 3 * indi);
      lambdas_temp.resize(indi * (respno + 1));
      if (generate_or_diagnose) tby_individuals(daten, avw_temp.data(), lambdas_temp.data(), rsts);
    }
    
    for (int ithread = 0; ithread != NOTHREADS; ithread++) {
      if (flag == 1) {
        int jj = igroup * icompg;
        for (int t = 0; t != indi; t++) {
          for (int type = 0; type != 3; type++) {
            int ift = ifree[type];
            for (int ip = 0; ip != ift; ip++)
              if (dCOMP(type, ip)) {
                gsl_vector_set(hampar, jj, avw_temp[t * 3 * ifreemax + type * ifreemax + ip]); jj++;
              }
          }
        }
        jj = (indi + igroup) * icompg + igroup * respno;
        for (int it = 0; it != indi; it++) for (int ir = 0; ir != respno; ir++)
          gsl_vector_set(hampar, jj + it * respno + ir, lambdas_temp[it]);
        jj += indi * respno;
        for (int it = 0; it != indi; it++) gsl_vector_set(hampar, jj + it, lambdas_temp[indi + it]);
      }
      
      gsl_rng_memcpy(xst, rsts[ithread]);

      if (flag == 0) initialize_new0(daten, hampar, tavw.data(), tlams.data(), ai.data(), loglambdas.data(), bi.data(), paths.data(), xst);
      else initialize_new1(daten, hampar, tavw.data(), tlams.data(), ai.data(), loglambdas.data(), bi.data(), paths.data(), xst);
      
      for (int x = 0; x != datenzahl; x++) {
        trial one = daten[x]; double rt = one.rt / 1000.0; int c = one.category; int t = one.person; int itree = one.tree;
        int path = paths[x]; int pfadlength = dNDRIN(c, path);
        std::vector<double> as(pfadlength);
        std::vector<double> vs(pfadlength);
        std::vector<double> ws(pfadlength);
        std::vector<double> taus(pfadlength);
        std::vector<int> pms(pfadlength);
        for (int ir = 0; ir != pfadlength; ir++) {
          int n = dDRIN(c, path, ir);
          int ia = dTREE_AND_NODE2PAR(itree, n, 0);
          int iv = dTREE_AND_NODE2PAR(itree, n, 1);
          int iw = dTREE_AND_NODE2PAR(itree, n, 2);
          as[ir] = dTAVW(t, 0, ia);
          vs[ir] = dTAVW(t, 1, iv);
          ws[ir] = dTAVW(t, 2, iw);
          pms[ir] = (dAR(c, path, n) + 1) / 2;
        }
        int icount = 0;
        NEW: rest[x] = rt;
        icount++;
        for (int ir = 0; ir != pfadlength; ir++) {
          taus[ir] = rwiener_diag(pms[ir], rt, as[ir], vs[ir], ws[ir], xst);
          rest[x] -= taus[ir];
          if (rest[x] <= 0) {
            if (icount < 100000) goto NEW;
          }
        }
        if (rest[x] > 0)
          for (int ir = 0; ir != pfadlength; ir++) {
            int n = dDRIN(c, path, ir);
            alltaus[dTAU_BY_NODE(x, n, pms[ir])] = taus[ir] * dAR(c, path, n);
          }
          else {
            rest[x] = rt * 2.0 / 3.0;
            for (int ir = 0; ir != pfadlength; ir++) {
              int n = dDRIN(c, path, ir);
              alltaus[dTAU_BY_NODE(x, n, pms[ir])] = (rt - rest[x]) / pfadlength * dAR(c, path, n);
            }
            
          }
      }
      
      make_nips(daten, paths.data(), nips.data());
      push(ithread, n_value_store, n_all_parameters, hampar, tavw.data(), tlams.data(), ai.data(), loglambdas.data(), bi.data(), alltaus.data(), rest.data(), datenzahl, paths.data(), nips.data(), liknorm, activeeps, epsm, Hobjective, valuestore, parmon.data(), parmonstore);
      gsl_rng_memcpy(rsts[ithread], xst);
    }
    
    gsl_rng_free(xst);
    gsl_vector_free(hampar);
    
  }
  
  
  //compute R-hat for convergence check
  void r_statistic(int ido, int n_all_parameters, int istream, int iter, double *parmon, double *xwbr, double &rmax, int &imax) {
    int drnall = 3 * n_all_parameters;
    if (ido == 1) for (int i = 0; i != drnall; i++) xwbr[i] = 0;
    double r = 1.0 / (istream + 1);
    
    gsl_vector* temp1 = gsl_vector_alloc(n_all_parameters);
    gsl_vector_view tpm = gsl_vector_view_array(parmon, 2 * n_all_parameters);
    gsl_vector_view tpm1 = gsl_vector_subvector(&tpm.vector, 0, n_all_parameters);
    gsl_vector_view tpm2 = gsl_vector_subvector(&tpm.vector, n_all_parameters, n_all_parameters);
    gsl_vector_memcpy(temp1, &tpm1.vector);
    
    gsl_vector_view txwbr = gsl_vector_view_array(xwbr, 3 * n_all_parameters);
    gsl_vector_view txwbr1 = gsl_vector_subvector(&txwbr.vector, 0, n_all_parameters);
    gsl_vector_view txwbr2 = gsl_vector_subvector(&txwbr.vector, n_all_parameters, n_all_parameters);
    gsl_vector_view txwbr3 = gsl_vector_subvector(&txwbr.vector, 2 * n_all_parameters, n_all_parameters);
    
    gsl_blas_daxpy(-1.0, &txwbr3.vector, temp1);
    gsl_vector* temp2 = gsl_vector_alloc(n_all_parameters);
    gsl_vector_memcpy(temp2, temp1);
    gsl_vector_mul(temp1, temp1);
    gsl_blas_daxpy(1.0 - r, temp1, &txwbr2.vector);
    gsl_blas_daxpy(r, temp2, &txwbr3.vector);
    gsl_vector_memcpy(temp2, &tpm2.vector);
    gsl_blas_daxpy(-1.0, &txwbr1.vector, temp2);
    gsl_blas_daxpy(r, temp2, &txwbr1.vector);
    
    if (ido == 3) {
      gsl_vector_memcpy(temp1, &txwbr1.vector);
      gsl_vector_memcpy(temp2, &txwbr2.vector);
      gsl_vector_scale(temp2, 1.0 / istream);
      gsl_blas_daxpy(1.0 / iter, temp1, temp2);
      gsl_vector_div(temp2, temp1);
      gsl_vector_scale(temp2, iter - 1.0);
      for (int i = 0; i != n_all_parameters; i++) gsl_vector_set(temp2, i, sqrt(gsl_vector_get(temp2, i)));
      gsl_vector_memcpy(&txwbr3.vector, temp2);
      imax = gsl_vector_max_index(&txwbr3.vector);
      rmax = dXWBR(3, imax);
    }
    gsl_vector_free(temp1); gsl_vector_free(temp2);
  }
  
  //type of index given on scale from 0 .. ifreeg-1
  int is_type(int ip) {
    int temp = -1;
    if (ip < 0) return(temp);
    if (ip < ifree[0]) return(0);
    if (ip < ifree[0] + ifree[1]) return(1);
    if (ip < ifree[0] + ifree[1] + ifree[2]) return(2);
    return(temp);
  }
  
  //number of index given on scale from 0..ifreeg-1, and type, to scale (type,ifree[type])
  int ind_type(int type, int ip) {
    int temp = -1;
    switch (type) {
    case 0: temp = ip; break;
    case 1: temp = ip - ifree[0]; break;
    case 2: temp = ip - ifree[0] - ifree[1]; break;
    }
    return temp;
  }
  
  //logit transform
  double logit(transform par, double  x) {
    double z = par.scale * x + par.loc;
    if (z < -700) return par.a;
    return par.a + par.range/ (1.0 + exp(-z));
  }
  
  //derivative of logit transform
  double dlogit(transform par, double  x) {
    double z = par.scale * x + par.loc;
    if (z < -700) return 0.0;
    if (z > 700) return 0.0;
    double temp = exp(-z);
    return par.range * par.scale * temp / gsl_pow_2(1.0 + temp);
  }
  
  //inverse logit transform
  double invlogit(transform par, double u) {
    u = (u - par.a) / par.range;
    u = log(u / (1 - u));
    u -= par.loc;
    u /= par.scale;
    return u;
  }
  
  // density of scaled t distribution without constants
  double log_tdist_pdf(double mu, double sig, int ddf, double x) {
    return -(ddf + 1.0) / 2.0 * gsl_log1p(gsl_pow_2((x-mu) / sig) / ddf);
  }
  
  //prepare logit transform so that lower bound is a, upper bound b, value of zero gives loc, value of 1 gives loc + scale
  transform prep_transform(double a, double b, double loc, double scale) {
    transform atrans;
    atrans.a = a;
    atrans.b = b;
    atrans.range = atrans.b - atrans.a;
    atrans.loc = (loc - atrans.a) / (atrans.b - atrans.a);
    atrans.loc = log(atrans.loc / (1 - atrans.loc));
    atrans.scale = (loc + scale - atrans.a) / (atrans.b - atrans.a);
    atrans.scale = log(atrans.scale / (1 - atrans.scale));
    atrans.scale = (atrans.scale - atrans.loc);
    return atrans;
  }
  
  //store required parameters and settings for sample increase after the algorithm has finished.
  void push_continue(int n_value_store, int irun, double* valuestore, double* parmonstore, std::vector<gsl_rng*>& rsts) {
    
    FILE* random;
    fopen_s(&random, RANDOM, "wb");
    for (std::size_t i = 0; i < rsts.size(); ++i) {
      gsl_rng_fwrite(random, rsts[i]);
    }
    // gsl_rng_fwrite(random, rst1);
    // gsl_rng_fwrite(random, rst2);
    // gsl_rng_fwrite(random, rst3);
    // gsl_rng_fwrite(random, rst4);
    fclose(random);
    
    std::ofstream contin; contin.open(CONTINUE);
    
    contin << std::setw(5) << irun << std::endl;
    for (int i = 0; i != NOTHREADS * n_value_store; i++) contin << std::setw(20) << valuestore[i];
    contin << std::endl;
    for (int i = 0; i != NOTHREADS * 2 * n_all_parameters; i++) contin << std::setw(20) << parmonstore[i];
    contin << std::endl;
    // supersig does not need to be stored - delete?
    //		for (int i = 0; i != NOTHREADS * n_all_parameters * n_all_parameters; i++) contin << setw(20) << supersig[i];
    //		contin << std::endl;
    for (int i = 0; i != n_all_parameters; i++) for (int j = 0; j != n_all_parameters; j++) contin << std::setw(20) << gsl_matrix_get(sigisqrt, i, j);
    contin << std::endl;
    for (int i = 0; i != n_all_parameters; i++) for (int j = 0; j != n_all_parameters; j++) contin << std::setw(20) << gsl_matrix_get(supsig, i, j);
    contin << std::endl;
    contin.close();
  }
  
  //retrieve required parameters and settings for sample increase after the algorithm has finished.
  void pop_continue(int n_value_store, int& irun, double* valuestore, double* parmonstore, std::vector<gsl_rng*>& rsts) {
    
    //	supersig = (double*)calloc(NOTHREADS * n_all_parameters * n_all_parameters, sizeof(double));
    FILE* random;
    fopen_s(&random, RANDOM, "rb");
    for (std::size_t i = 0; i < rsts.size(); ++i) {
      gsl_rng_fread(random, rsts[i]);
    }
    // gsl_rng_fread(random, rst1);
    // gsl_rng_fread(random, rst2);
    // gsl_rng_fread(random, rst3);
    // gsl_rng_fread(random, rst4);
    fclose(random);
    
    std::ifstream contin; contin.open(CONTINUE);
    
    contin >> irun;
    for (int i = 0; i != NOTHREADS * n_value_store; i++) contin >> valuestore[i];
    for (int i = 0; i != NOTHREADS * 2 * n_all_parameters; i++) contin >> parmonstore[i];
    //	for (int i = 0; i != NOTHREADS * n_all_parameters * n_all_parameters; i++) contin >> supersig[i];
    for (int i = 0; i != n_all_parameters; i++) for (int j = 0; j != n_all_parameters; j++) {
      double temp; contin >> temp;
      gsl_matrix_set(sigisqrt, i, j, temp);
    }
    for (int i = 0; i != n_all_parameters; i++) for (int j = 0; j != n_all_parameters; j++) {
      double temp; contin >> temp;
      gsl_matrix_set(supsig, i, j, temp);
    }
    contin.close();
  }
  
  //integrand for convolution transformed as in STAN
  int nstep2(unsigned dim, const double* x, void* p,
             unsigned fdim, double* retval)
  {
    struct my_params* params = (struct my_params*)p;
    // int pfadlength = (params->pfadlength);
    double* a = (params->a);
    double* v = (params->v);
    double* w = (params->w);
    int* low_or_up = (params->low_or_up);
    double mu = (params->mu);
    double sig = (params->sig);
    double t = (params->rt_rest);
    
    std::vector<double> tau(dim + 1);
    
    //Transformation aus Stan unit simplex
    
    double prod = 1.0, diff = t;
    unsigned int i;
    for (i = 0; i < dim + 1; i++) {
      tau[i] = diff;
      if (i < dim) {
        prod *= diff;  tau[i] *= x[i];
        prod *= exp(dwiener_d(tau[i] * low_or_up[i], a[i], v[i], w[i], accuracy));
      }
      diff -= tau[i];
    }
    
    prod *= gsl_ran_tdist_pdf((tau[dim] - mu) / sig, degf);
    retval[0] = prod;
    return 0;
  }
  
  //convolution
  void convolution2(const std::vector<double> & rts, int pfadlength, int* low_or_up, double* a, double* v, double* w, double mu, double sig, std::vector<double>& pbranch) {
    //	pbranch.clear();
    
    double val, err;
    double reltol = 1.0e-4, abstol = 0.0;
    
    std::vector<double> xmin(pfadlength);
    std::vector<double> xmax(pfadlength);
    
    for (int i = 0; i != pfadlength; ++i) {
      xmin[i] = 0.0;
      xmax[i] = 1.0;
    }
    
    int rtss = rts.size();
    for (int x = 0; x != rtss; x++) {
      double rt = rts[x];
      my_params params = { pfadlength, a, v, w, low_or_up, mu, sig, rt };
      hcubature(nstep2, &params, pfadlength, xmin.data(), xmax.data(), 0, abstol, reltol, &val, &err);
      pbranch.push_back(val);
    }
  }

}
