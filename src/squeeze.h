/*
 *  Version: 2.0
 *  SQUEEZE 2 - Image reconstruction software for optical interferometry,
 *  based on Monte-Carlo Markov Chain algorithms.
 *
 *  Copyright (c) 2006-2013 Fabien Baron, John Monnier, Michael Ireland
 *
 *  New parallel tempering engine based on SQUEEZE 1 written by Pr. Fabien Baron (Georgia State University).
 *  Old simulated annealing engine based on MACIM written by Dr. Michael Ireland (Macquarie University) and
 *  Pr. John Monnier (University of Michigan).
 *
 *  SQUEEZE is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  SQUEEZE is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with Squeeze.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

// Minimization engines
#define ENGINE_SIMULATED_ANNEALING 1
#define ENGINE_PARALLEL_TEMPERING 2

// Regularizers
#define NREGULS 10
#define REG_MODELPARAM 0
#define REG_CENTERING 1
#define REG_PRIORIMAGE 2
#define REG_ENTROPY 3
#define REG_DARKENERGY 4
#define REG_TV 5
#define REG_SPOT 6
#define REG_LAP 7
#define REG_L0 8
#define REG_TRANSPECL2 9

// Mathematical constants
#define MAS_RAD          206264806.2
#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

#define STEPS_PER_OUTPUT 1
#define FRAC_COPYCAT     .1
#define FRAC_ANYWHERE    .05

#define DEFAULT_NITER    250
#define DEFAULT_DEPTH    500

#define BURN_IN_FRAC     5 /* At least 1/10 of the iterations must be used for burn-in (for mean image output) */
#define DEFAULT_CENT_MULT  1 /* Chi^2 changes by 1 when a flux element moves width/2 */

/* Old units: 25 Corresponded to 0.2 pix RMS at T=1.0 */

/* Initial Scaled Fit statistic (scaled by temperature), and
 the time in steps to change the temperature. The algorithm sets T
 so that chi^2/T = TARGET_SCALED_CHI2, with T >= 1.0
 NB Increase TARGET_SCALED_CHI2 or decrease number of elements if there is
 a problem with convergence.*/
#define TARGET_SCALED_CHI2   4.0
#define INITIAL_SCALED_CHI2_MULT  2
#define INITIAL_TMIN 1.0 /* Never start with a temperature less than 1 */
#define TEMP_CHANGE_TIME    1000.0
#define DEFAULT_TMIN 1.0
/* The damping time for parameter movement probability calculations */
#define PARAM_DAMPING_TIME 10.0
/* The e-folding time for adjusting the stepsize to give the desired
  mprob */
#define STEPSIZE_ADJUST_TIME 20.0
/* The target movement probability for chi^2 scaling */
#define TARGET_MPROB 0.3
/* The following define is for a model as part of the fit */
#define USE_MODEL
#define MAX_PARAMS       20
/* Number of times parameters are changed per element */
#define PARAMS_PER_ELT   2
/* Types of steps */
#define STEP_COPYCAT     3
#define STEP_ANYWHERE    2
#define STEP_MEDIUM      1
#define STEP_SMALL       0
#define STEP_MAX         3

/* Damping time for calculating probability of flux movement
 Currently this has to be relatively small so that natural fluctuations
 help with making different step types...*/
#define DAMPING_TIME     25.0
#define MPROB_LOW  0.2
#define MPROB_HIGH 0.45

/* For the MEM (or DEN) image, stop looking for mean chi^2=1.0 at this point.*/
#define MEM_DELTA_CHI2 0.02

/* The step type can not be changed until chi^2 has reduced to this
 multiple of chi^2 for a random, flat image. This also sets the
 starting temperature. */
#define FLAT_CHI2_MULT 3/* 0.25 */

/* Defines for dark energy */
#define DEN_ADD      0
#define DEN_SUBTRACT 1
#define DEN_INIT     65535

/* For regularization parameters (alpha, beta etc) this is the maximum possible value */
#define MAX_REG_PARAM 20


/* Max number of characters in a filename */
#define MAX_STRINGS   1024

#include <complex.h>
#include <stdbool.h>

/* Function prototypes for fred.c. Note that you have to include complex.h, and
 the complex number i is I.
 creal(vis_sig):  Error parallel to vis
 cimag(vis_sig):  Error perpendicular to vis
  ... same for bs_sig.*/
void * main_loop(void *index);
void printerror(int status);
void intHandler(int signum);
void printhelp(void);

bool read_commandline(int* argc, char** argv, bool* benchmark, bool* use_v2, bool* use_t3amp, bool* use_t3phi, bool* use_visamp, bool* use_visphi, bool* use_diffvis, bool* use_threadfits, bool* use_bandwidthsmearing, int* minimization_engine, bool* dumpchain,double* mas_pixel, unsigned short* axis_len, long* depth, long* niter, long* nelements, double* f_anywhere, double* f_copycat, int *nthreads, double* tempschedc, double* fov, double* chi2_temp, double* chi2_target, double* tmin, double* prob_auto, double* uvtol, char* output_filename, char* init_filename, char* prior_filename, double* v2s, double* v2a, double* t3amps, double* t3ampa, double* t3phia, double* t3phis, double* visamps, double* visampa, double* visphis, double* visphia, double* fluxs, double* cvfwhm, double* reg_param, double* init_param, double* wavmin, double* wavmax);

void print_diagnostics(int iThread, long current_iter, long nvis, long nv2, long nt3, long nt3phi, long nt3amp, long nvisamp, long nvisphi, double chi2v2, double chi2t3amp,double chi2t3phi,double chi2visphi,double chi2visamp, double lPosterior, double lPrior, double lLikelihood, const double* reg_param, const double* reg_value, const double* cent_xoffset, const double* cent_yoffset, long nelements, int nwavr, long niter, const double* temperature, double prob_movement, const double* params, const double* stepsize);

void compute_lLikelihood(double* likelihood, const double complex * __restrict mod_vis, double* __restrict res, double* __restrict mod_obs, double *chi2v2, double *chi2t3amp, double *chi2visamp, double *chi2t3phi, double *chi2visphi);
void compute_lPrior(double* lPrior, const long chan, const double* reg_param, const double* reg_value);

void vis_to_obs(const double complex *mod_vis, double *mod_obs);
void obs_to_res(const double *mod_obs, double *res);
double residuals_to_chi2(const double *res, double *chi2v2, double *chi2t3amp, double *chi2visamp, double *chi2t3phi, double *chi2visphi) ;

double get_flat_chi2(bool benchmark);
double fill_min_elts(long *min_elts, long depth, long threadnum);
static inline double dewrap(double diff) __attribute__((always_inline));
static inline double modsq(double complex input)  __attribute__((always_inline));

double fill_iframeburned(long *iframeburned, long depth, long threadnum, long nelements, long niter,  double *saved_lPosterior, double *saved_lLikelihood, double *saved_reg_value);
int find_reg_param(double *regparam, long *iframeburned, long depth, long niter, long ndf, long nelements);
int writeasfits(char *filename, double *image,
                long depth, long min_elt, double chi2, double temperature, long nelems, double* regpar , double* regval,
                long niter, unsigned short axis_len, double ndf, double tmin, double chi2_temp, double chi2_target, double mas_pixel, int nthreads, double *saved_params, double logZ, double logZ_err,
                char* init_filename, char* prior_filename);

void mcmc_annealing_image(char *file, double *image, long *iframeburned, long depth, long nelements,
                          unsigned short axis_len, double complex * __restrict xtransform, double complex * __restrict ytransform,
                          double *mn_chi2, unsigned short *saved_x, unsigned short *saved_y, double *saved_params, long niter, int nchanr);

void mcmc_tempering_image(char *file, double *image, long lowtempthread, long depth, long nelements,
                          unsigned short axis_len, double complex * __restrict xtransform, double complex * __restrict ytransform,
                          double *mn_chi2, unsigned short *saved_x, unsigned short *saved_y, double *saved_params, long niter, int nchanr);

void compute_logZ(const double* temperature , const unsigned short* iStoragetoThread, const double* lLikelihood_expectation, const double *lLikelihood_deviation, int nthreads, double* logZ, double* logZ_err);

void mcmc_fullchain(char* file, long nthreads, long niter, int nchanr, long nelements, unsigned short axis_len, unsigned short *saved_x, unsigned short *saved_y, double *saved_params, double *saved_lLikelihood, double *saved_lPrior, double *saved_lPosterior, double *temperature, unsigned short* iThreadtoStorage);

void compute_regularizers(double *reg_param, double *reg_value, double* image, double* prior_image, unsigned short* initial_x, unsigned short* initial_y, unsigned short* element_x,  unsigned short* element_y, int nwavr, unsigned short axis_len, long nelements, double* cent_xoffset, double* cent_yoffset, double fov, double cent_mult);

void compute_model_visibilities_fromelements(double complex* mod_vis, double complex* im_vis, double complex* param_vis, double* params, double* fluxratio_image, const unsigned short* element_x, const unsigned short* element_y, const double complex* xtransform, const double complex* ytransform, double* lPriorModel, long nparams, long nelements);

void compute_model_visibilities_fromimage(double complex* mod_vis, double complex* im_vis, double complex* param_vis, const double* params, double* fluxratio_image, const double *image, const double complex* xtransform, const double complex* ytransform, double* lPriorModel, long nparams, long nelements, unsigned short axis_len);

void initialize_image(int iThread, double* image, unsigned short* element_x, unsigned short* element_y, unsigned short* initial_x, unsigned short* initial_y, unsigned short axis_len, int nwavr, long nelements, char* init_filename);

/* Function prototype for extract_oifits.c*/
int extract_oifits(char* filename, bool use_v2, bool use_t3amp, bool use_t3phi, bool use_visamp, bool use_visphi,
                   double v2a, double v2s, double t3ampa, double t3amps, double t3phia, double t3phis,
                   double visampa, double visamps, double visphia, double visphis, double fluxs, double cwhm, double uvtol, double* wavmin, double *wavmax, double *timemin, double *timemax);
int write_best_oifits(char* filestring, double complex * mod_vis);

/* Function prototype for modelcode.c */
int model_vis(const double *params, double complex *modvis, double *logl, double *flux_frac);

/* regularizations */

double entropy(unsigned long s);
double den_full(const double* x, const double* pr, const double eps, const int nx, const int ny);
double den_change(const double *image, const unsigned short i, const unsigned short j, const unsigned short direction, const unsigned short axis_len);
double UDreg(const double* x, const double* pr, const double eps, const int nx, const int ny);
double TV(const double* x, const double* pr, const double eps, const int nx, const int ny);
double LAP(const double* x, const double* pr, const double eps, const int nx, const int ny);
double L0(const double* x, const double* pr, const double eps, const int nx, const int ny);
double entropy_full(const double* x, const double* pr, const double eps, const int nx, const int ny);
double L2(const double* x, const double* pr, const double eps, const int nx, const int ny);
double transpec(int nchanr, long axis_len, double *image);
double cent_change(int channel, double* cent_xoffset, double *cent_yoffset, long new_x, long new_y, long old_x, long old_y, unsigned short axis_len, double fov, double cent_mult);
/* helper functions */

double sinc(double x);
double mean(long *x, long n);
double stddev(long *x, long n);

inline void swapi(unsigned short* a, unsigned short* b);
inline void swapd(double* a, double* b);
double xatan2(double y, double x);
double xatan2_u1(double y, double x);
