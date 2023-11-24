#ifndef RTS_H
#define RTS_H
// authors: Christoph Klauer and Raphael Hartmann

//#define _CRT_SECURE_CPP_OVERLOAD_STANDARD_NAMES 1
//#include <string.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <cfloat> // for DBL_MAX
#include <thread>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_multimin.h>

#include <R.h>
#include <Rinternals.h>


struct trial {
  int person, tree, category, item,group,rt;
};

extern const char *MODEL;
const double sqr2pi = 2.506628274631000502415765e0;
extern int kernpar;
extern int zweig;
extern int nodemax;
extern int kerncat;
extern int datenzahl;
extern int* ng;
extern int indi;
extern int *t2group;


namespace ertmpt {
  
  int mainx(int *k2f, int *f2k); //int argc, char *argv[]
  
  // #define RAUS "raus"
  //#define LOGLIK "log_likelihood"
  // #define log_lik_flag false
  // #define for_bridge_flag false
  
  //#define diagnose
  #define generate_and_save
  
  #ifdef generate_and_save
  const bool generate_or_diagnose = true;
  const bool save_diagnose=true;
  #endif
  
  
  #ifdef diagnose
  const bool generate_or_diagnose = false;
  const bool save_diagnose=true;
  #endif
  
  #define DEBUG false
  
  //#define PRIOR 1.00
  
  //#define adf 2
  // #define IREP 2200
  // #define BURNIN 200 // (must be smaller than IREP)
  // #define THIN 11 // (must be a divisor of IREP)
  // #define NOTHREADS 4
  // #define SAMPLE_SIZE 20000 //(a multiple of IREP*NOTHREADS/THIN)
  // #define RMAX 1.05
  
  #define T_rng gsl_rng_ranlxd1
  
  #define LNNORM_MAX_X 38.0
  #define LNNORM_MIN_X -1.00e9
  

  struct pfadinfo {
  int a;
  std::vector<int> r;
  std::vector<int> pfad_par;
  std::vector<int> pm;
  };
  
  
  struct point {
  	double x;
  	double h;
  	double dh;
  };
  
  
  extern int log_lik_flag;
  extern int for_bridge_flag;
  
  extern const char *DATA;
  extern const char *RAUS;
  extern const char *diagn_tests;
  extern const char *LOGLIK;
  
  extern int IREP;
  extern int BURNIN;
  extern int THIN;
  extern int NOTHREADS;
  extern int SAMPLE_SIZE;
  extern double RMAX;
  
  extern int nKERN;
  extern int nRESP;
  
  extern int *CatToResp;
  
  extern double *ConstProb;
  extern int *CompMinus;
  extern int *CompPlus;
  
  extern double *complete_sample;
  extern double *complete_bridge;
  
  //Globale Variablen
  
  extern int *cat2tree;
  extern int *branch;
  extern int *ar;
  extern int *nodes_per_tree;
  extern int *tree_and_node2par;
  extern bool *comp;
  extern int ilamfree, ifree,ipred;
  
  extern int *ndrin,*drin;
  extern int n_all_parameters;
  extern int *nppr;
  extern int n_bridge_parameters;
  
  extern int RMAX_reached;
  extern bool BURNIN_flag;
  
  extern int igroup;
  extern int ireps;
  extern int *cat2resp;
  extern int respno;
  extern int alphaoff;
  extern int sigalphaoff;
  extern int restparsno;
  extern int *free2kern;
  extern int *kern2free;
  extern double *consts;
  
  extern double pr_df_sigma_sqr;
  extern double pr_shape_omega_sqr;
  extern double pr_rate_omega_sqr;
  extern double pr_mean_mu_gamma;
  extern double pr_var_mu_gamma;
  extern double PRIOR;
  extern double pr_shape_exp_mu_beta;
  extern double pr_rate_exp_mu_beta;
  extern double pr_sf_scale_matrix_SIG;
  extern double pr_sf_scale_matrix_TAU;
  extern int pr_df_add_inv_wish;
  
  extern int *pfad_index;
  extern std::vector<pfadinfo> path_info;
  
  void make_idaten(std::vector<trial> daten,int *idaten);
  
  //void by_individuals(std::vector<trial> daten,int kerntree,double *beta, double *g2,double *likeli,gsl_rng *rst);
  
  void make_pij_for_individual(double *x, double *pij, double *pj);
  
  void make_mu(double *mu, double *lams, double *beta,  int *nnodes, double *z, gsl_rng *rst);
  
  void make_lams(double *mu, double *lams, double *beta,  int *nnodes, double *z, gsl_rng *rst);
  
  
  void bayesreg(int n,double *mean,double *sigma, double *out, gsl_rng *rst);
  
  void invwis(int cases, int nvar, double *xx, double *ssig, double *sigi, double eps, gsl_rng *rst);
  
  double truncexp(double lambda, double upper, gsl_rng *rst);
  
  double equation(int t, int ip, double *mu, double *lams, double *beta);
  
  double ars(double step, double &scale,double totallow, double n, double xp, double *beta, double *sigi, double *lambdas,double *lams,int tt, int iz, double start, gsl_rng *rst,
  	void gamma_prior(double scale, double norm, double n, double alpha, double p, double *beta, double *sigi, double *lambdas,double *lams,int t,  int iz, bool deriv, point &h));
  
  
  void r_statistic(int ido, int n_all_parameters, int istream, int iter, double *parmon, double *xwbr, double &rmax);
  
  double rexp(double x);
  
  void diagnosis(std::vector<trial> daten,int *idaten, int kerntree, gsl_rng *rst);
  
  //void make_index(bool *lambda_comp, int *index);
  
  int make_path_for_one_trial(int branchno, double *pij, double p, gsl_rng *rst);
  
  double oneexp(double lambda, gsl_rng *rst);
  
  void lies(std::vector<trial> &daten);
  
  void model_design(int kerntree,int *ar, int *branch, int *nodes_per_par, int *nodes_per_tree, int *tree_and_node2par);
  
  //void ddf(std::vector <trial> &daten, int kerntree, double *g2);
  
  double elogdiff(double xa, double xb);
  
  void make_parameters_for_all(double *mu, double *lams, double *beta, double *x_for_all);
  
  void make_pij_for_one_trial(trial one,double *x_for_all,  double *pij, double &pj);
  
  void make_zs_one_trial(trial one,int itrial, int ipath, double *mu, double *lams, double *beta,int *nz_position, double *z,gsl_rng *rst) ;
  
  double malpha(int t, int r,double* restpars,double* slams);
  
  void make_rtau(double *restpar, double *taui, double *slams, gsl_rng *rst);
  
  void make_rmu(std::vector<trial> daten,double *factor, double *rest,double *restpar, double* slams, gsl_rng *rst);
  
  void make_ralpha(std::vector<trial> daten,double *factor, double *rest, double *restpar,double *slams, double *taui, gsl_rng *rst);
  
  void make_rsigalpha(std::vector<trial> daten, double* factor, double *rest, double *restpar, double *slams, bool xflag, gsl_rng *rst);
  
  void make_rsig(std::vector<trial> daten, double *rest, double *restpar, gsl_rng *rst);
  
  void make_slams(std::vector<trial> daten,double* factor, double *rest, double *restpars, double *slams, gsl_rng *rst);
  
  double double_truncnorm(double lower, double upper, gsl_rng *rst);
  
  double logexgaussian(double lam, double mu, double sd, double t);
  
  void gibbs_times_new(std::vector<trial> daten,int *nnodes,int nz, int *nz_position,double *beta,int ntau, int *ntau_position,
  	gsl_rng *rst1, gsl_rng *rst2, gsl_rng *rst3, gsl_rng *rst4, gsl_rng *rst5, gsl_rng *rst6, gsl_rng *rst7, gsl_rng *rst8, gsl_rng *rst9, gsl_rng *rst10, gsl_rng *rst11, gsl_rng *rst12, gsl_rng *rst13, gsl_rng *rst14, gsl_rng *rst15, gsl_rng *rst16,
  	double *lambdas, double* restpars);
  
  void make_tij_for_one_trial_new(trial one,double *rhos, double *lambdas,double *restpars,double *pij);
  
  void tby_individuals(std::vector<trial> daten, int kerntree,double *beta,double *lambdas, double *restpars, gsl_rng *rst);
  
  double logdiff(double xa, double xb);
  
  void extract_pfadinfo(int *pfad_index, std::vector<pfadinfo> &path_info);
  
  void loggammagaussian(int n, double lam, double mu, double sd, double t, double& hplus, double& hminus);
  
  double  logf_tij(int a, std::vector<int> r, double *lams, double *loglams, double mu, double sd, double t);
  
  int gsl_linalg_tri_lower_invert(gsl_matrix * T);

}

void set_ns(std::vector<trial> daten, int &indi, int &kerntree, int &kerncat,int &igroup);
double logsum(double xa, double xb);
double lnnorm(double x);
double oneuni(gsl_rng *rst);
double onenorm(gsl_rng *rst);
double truncnorm(double b, gsl_rng *rst);

#endif
