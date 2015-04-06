/*
 *  SQUEEZE - Image reconstruction software for optical interferometry,
 *  based on Monte-Carlo Markov Chain algorithms.
 *
 *  Copyright (c) 2006-2014 Fabien Baron, John Monnier, Michael Ireland
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
#ifdef _OPENMP
#include <omp.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>
#include <signal.h>
#include <time.h>
#include <string.h>
#include "../lib/rngstreams/src/RngStream.h"
#include "squeeze.h"
#include "exchange.h" // includes cfitio

#define TEXT_COLOR_RED     "\x1b[31m"
#define TEXT_COLOR_GREEN   "\x1b[32m"
#define TEXT_COLOR_YELLOW  "\x1b[33m"
#define TEXT_COLOR_BLUE    "\x1b[34m"
#define TEXT_COLOR_MAGENTA "\x1b[35m"
#define TEXT_COLOR_CYAN    "\x1b[36m"
#define TEXT_COLOR_BLACK   "\x1b[0m"

static bool ctrlcpressed = FALSE; // CHECKS CTRL+C PRESSED DURING EXECUTION
int oi_hush_errors = 0; // flag for read_fits.c

/* Global OIFITS variables */
long nuv;
char *oifits_file;
long nvis, nv2, nt3, nt3phi, nt3amp, nt3amp_orphans, nt3phi_orphans, nvisamp, nvisphi, nvisamp_orphans, nvisphi_orphans;
long *visin, *v2in, *t3in1, *t3in2, *t3in3;
double *u, *v;
int ntimer;
bool use_diffvis = FALSE; // default for VIS tables = complex vis, not differential vis
long *dvisnwav;
long **dvisindx;

double* __restrict uv_time;
double* __restrict uv_lambda;
double* __restrict uv_dlambda;
int *uvwav2chan = NULL;
int *uvtime2chan = NULL;
double* __restrict data;
double* __restrict data_err;

/* Parametric Model Fitting*/
long nparams;
double init_params[MAX_PARAMS], init_stepsize[MAX_PARAMS];

/* SQUEEZE MAIN LOOP */
int main(int argc, char** argv) {

	// TBD: documentation for all variables coming soon.....
	signal(SIGINT, intHandler);

	int minimization_engine = ENGINE_SIMULATED_ANNEALING;
	int nchains = 0; // default number of parallel MCM chains
	int nthreads = 0; // default number of threads for computation
	int nmaxthreads = 0;
	int nthreadsperchain = 1;
	double cent_mult = 0.0, fov = 1; // centroid regularization
	double reg_param[NREGULS];
	long niter = DEFAULT_NITER;
	double mas_pixel;
	unsigned short axis_len;
	long nelements = 0;
	/* Important derived variables */
	double ndf, flat_chi2;
	double tmin = DEFAULT_TMIN, chi2_temp = TARGET_SCALED_CHI2, chi2_target = 0.0;
	double *prior_image = NULL;

	unsigned short *initial_x, *initial_y;

	double prob_auto = -1.0;

	double f_copycat = FRAC_COPYCAT, f_anywhere = FRAC_ANYWHERE;

	bool dumpchain = FALSE;
	bool benchmark = FALSE;
	//  const char *reg_names[NREGULS] = {"PARAM", "C", "PRI", "ENT", "DEN", "TV", "UD", "LAP", "L0", "TS"};
	long i, j, k, w, depth = DEFAULT_DEPTH, offset = 0;

	double* initial_image = NULL;
	char dummy_char[MAX_STRINGS], param_string[MAX_STRINGS];
	/* Variables needed for file i/o and creation of transform */
	double multx = 0, multy = 0, dtemp, ftot;
	double final_chi2, final_chi2v2, final_chi2t3amp, final_chi2t3phi, final_chi2visamp, final_chi2visphi;

	/* Stuff for fits file output */
	fitsfile *fptr; /* pointer to the FITS file, defined in fitsio.h */
	int status;

	/* initialize FITS image parameters */
	char output_filename[MAX_STRINGS] = "output";
	char prior_filename[MAX_STRINGS] = "";
	char init_filename[MAX_STRINGS] = "";

	double nullval = 0;
	int dummy_int;

	int nfound;
	double in_mas_pixel;
	long nwavi = 1, nwavp = 1;
	long* in_naxes;

	int nwavr;
	bool use_v2 = TRUE, use_t3amp = TRUE, use_t3phi = TRUE, use_visamp = TRUE, use_visphi = TRUE;
	bool use_tempfitswriting = FALSE, use_bandwidthsmearing = TRUE;
	double v2a = 0., t3ampa = 0., t3phia = 0., visampa = 0., visphia = 0.;
	double v2s = 1., t3amps = 1., t3phis = 1., visamps = 1., visphis = 1.;
	double cvfwhm = 0., uvtol = 1e3, fluxs = 1.;
	double tempschedc = 3.0;

	double *timemin = NULL;
	double *timemax = NULL;
	double *wavmin = NULL ;
	double *wavmax = NULL;

	printf(TEXT_COLOR_RED"SQUEEZE - Version %1.1f\n"TEXT_COLOR_BLACK, SQUEEZE_VERSION);

	/* Initialize regularization parameters */
	for (i = 0; i < NREGULS; i++)
		reg_param[i] = 0;
	/* Unless parameters are input, default to no model parameters (image only) */
	nparams = 0;

	/* Defults is to have the model parameters equal to 0 and fixed. */
	for (i = 0; i < MAX_PARAMS; i++)
		init_params[i] = 0.0;
	for (i = 0; i < MAX_PARAMS; i++)
		init_stepsize[i] = 0.0;

	mas_pixel = 0.0;
	axis_len = 0;

	// READ IN COMMAND LINE ARGUMENTS

	if (read_commandline(&argc, argv, &benchmark, &use_v2, &use_t3amp, &use_t3phi, &use_visamp, &use_visphi, &use_diffvis, &use_tempfitswriting,
			&use_bandwidthsmearing, &minimization_engine, &dumpchain, &mas_pixel, &axis_len, &depth, &niter, &nelements, &f_anywhere, &f_copycat, &nchains,
			&nthreads, &tempschedc, &fov, &chi2_temp, &chi2_target, &tmin, &prob_auto, &uvtol, &output_filename[0], &init_filename[0], &prior_filename[0], &v2s,
			     &v2a, &t3amps, &t3ampa, &t3phia, &t3phis, &visamps, &visampa, &visphis, &visphia, &fluxs, &cvfwhm, reg_param, init_params, wavmin, wavmax, &nwavr) == FALSE)
		return 0;

	printf("DEBUG wavmin pointer: %p \n", &wavmin[0]);

	// Check nchains and nthreads are consistent, and if not overwrite them
	// First check if nchains and nthreads have been set
	if (nchains == 0)
		nchains = 1;
#ifdef _OPENMP
	nmaxthreads = omp_get_max_threads()/2;
#else
	nmaxthreads = 1;
#endif

	if (nthreads == 0)
		nthreads = nmaxthreads; //use max available threads if not set

	if (nthreads > nmaxthreads) {
		nthreads = nmaxthreads;
		printf("Command line -- Requested number of threads: %d is greater than physically available: %d\n", nthreads, nmaxthreads);
	}

	// compute the resulting number of threads within each chain
	nthreadsperchain = nthreads / nchains;
	if (nthreadsperchain == 0)
		nthreadsperchain = 1;

	printf("Threading    -- Number of Chains: %d\n", nchains);
	printf("Threading    -- Number of OS threads: %d out of %d possible maximum\n", nthreads, nmaxthreads);
	printf("Threading    -- Number of OS threads per chain that may be used: %d\n", nthreadsperchain);

	// If we want to use parallel tempering but no threads are defined,
	if ((minimization_engine == ENGINE_PARALLEL_TEMPERING) && ((nchains == 1) || (nthreads == 1)))
	  {
	    printf("Threading    -- Parallel tempering requested but the requested number of chains is %d\n", nchains);
	    printf("Threading    -- Please restart SQUEEZE and set this number accordingly.\n");
	    return 0;
	  }

	// If no wavelength info was given, then we are in monochromatic mode
	if((wavmin == NULL) || (wavmax == NULL))
	  {
	    printf("Command line -- Monochromatic reconstruction\n");
	    nwavr = 1;
	    wavmin = malloc(sizeof(double));
	    wavmax = malloc(sizeof(double));
	  }

	fflush (stdout);

	/* Read in oifits file */

	if (import_single_epoch_oifits(argv[1], use_v2, use_t3amp, use_t3phi, use_visamp, use_visphi, v2a, v2s, t3ampa, t3amps, t3phia, t3phis, visampa, visamps, visphia,
				       visphis, fluxs, cvfwhm, uvtol, nwavr, wavmin, wavmax, timemin, timemax))
	  {
		printf("Error opening %s. \n", argv[1]);
		return 0;
	  }



	oifits_file = argv[1];

	if (nuv == 0) {
		printf("No usable data in OIFITS file.\nExiting...\n");
		ctrlcpressed = TRUE;
	}
	fflush(stdout);

	/* Initialise default values...(use multx,multy as temp variables)*/
	/* Find longest (multx) and shortest (multy) baselines */
	for (i = 0; i < nuv; i++) {
		dtemp = u[i] * u[i] + v[i] * v[i];
		if (i == 0) {
			multx = dtemp;
			multy = dtemp;
		} else {
			if (dtemp > multx)
				multx = dtemp;
			if (dtemp < multy)
				multy = dtemp;
		}
	}

	multx = sqrt(multx);
	multy = sqrt(multy);
	if (mas_pixel == 0.0)
		mas_pixel = MAS_RAD / multx / 6.; /* Default 6 pix per max baseline fringe */

	if (axis_len == 0) {
		axis_len = ceil(multx / multy * 6.); // we add a factor two for binaries which may not be centered
		if (axis_len > 1024)
			axis_len = 1024; /* This is a crazy image size...*/
	}

	/* Number of elements depends on degrees of freedom and size of image */
	/* Completely empirical formula */
	if (nelements == 0) {
		nelements = 2. * ceil(axis_len * pow(nv2 + nvisamp + nvisphi + nt3amp + nt3phi, 0.333));
		if (nelements < 500)
			nelements = 500;
	}

	printf("Reconst setup -- Baseline range:\t%ld - %ld wavelengths\n", (long) round(multy), (long) round(multx));
	printf("Reconst setup -- Pixel scale:   \t%lf mas/pixel\n", mas_pixel);
	printf("Reconst setup -- Image width:   \t%hu pixels\n", axis_len);

	/* Initial/Starting image (command line -i)*/
	/* This can be either a command such as random that generates a random starting point */
	/* Or a FITS image. Unlike the prior image, the starting image is then "converted" into elements */
	/* Note that a FITS starting image may include starting parameters for models */

	if (init_filename[0] != 0) {

		if (!strcmp(init_filename, "random")) {
			printf("\nInitial image -- Random image common to all parallel chains and wavelengths\n");
			nwavi = 1; // TBD : implement polychromatic init
			initial_x = malloc(nelements * sizeof(unsigned short));
			initial_y = malloc(nelements * sizeof(unsigned short));

			RngStream initrng = RngStream_CreateStream("init");
			for (i = 0; i < nelements; i++) {
				initial_x[i] = RngStream_RandInt(initrng, 0, 2147483647) % axis_len;
				initial_y[i] = RngStream_RandInt(initrng, 0, 2147483647) % axis_len;
			}
			RngStream_DeleteStream(&initrng);
		} else if (!strcmp(init_filename, "randomthr")) {
			printf("\nInitial image -- Random images unique to each parallel chain\n");
			nwavi = 1; // TBD : implement polychromatic init
			initial_x = malloc(nchains * nelements * sizeof(unsigned short));
			initial_y = malloc(nchains * nelements * sizeof(unsigned short));

			RngStream initrng = RngStream_CreateStream("init");
			for (j = 0; j < nchains; j++) {
				for (i = 0; i < nelements; i++) {
					initial_x[j * nelements + i] = RngStream_RandInt(initrng, 0, 2147483647) % axis_len;
					initial_y[j * nelements + i] = RngStream_RandInt(initrng, 0, 2147483647) % axis_len;
				}
			}
			RngStream_DeleteStream(&initrng);
		} else {
			status = 0;
			printf("\nInitial image -- Trying to open: '%s'\n", init_filename);
			if (fits_open_file(&fptr, init_filename, READONLY, &status))
				printerror(status);

			fits_read_key_lng(fptr, "NAXIS", &k, dummy_char, &status);

			if ((nwavr > 1) && (k == 3))
				in_naxes = malloc(3 * sizeof(long));
			else
				in_naxes = malloc(2 * sizeof(long));

			if (fits_read_keys_lng(fptr, "NAXIS", 1, 3, in_naxes, &nfound, &status))
				printerror(status);
			//printf("k:%ld nfound: %d in_naxes[0]: %ld in_naxes[1]: %ld in_naxes[2]: %ld \n",k,nfound,in_naxes[0], in_naxes[1], in_naxes[2]);

			if ((nwavr > 1) && (k == 3) && (nfound == 3)) {
				if (in_naxes[2] == nwavr) // polychromatic reconstruction using polychromatic init
					nwavi = in_naxes[2];
				else {
					printf("\nInitial image -- Error when reading NAXIS3\n");
					nwavi = 1;
				}
			} else
				nwavi = 1;

			printf("Initial image -- reading %ld x %ld pixels x %ld channel(s)\n", in_naxes[0], in_naxes[1], k > 2 ? in_naxes[2] : 1);

			if (in_naxes[0] != in_naxes[1]) {
				printf("Initial image -- input fits file must have a square array.");
				getchar();
			}

			/* In axis_len has changed, we have to offset the old versus new images.*/
			offset = (axis_len - in_naxes[0]) / 2;
			if (axis_len != in_naxes[0]) {
				printf("Initial image -- offset by : %ld pixels in x and y\n", offset);
			}

			/* read the SCALE keyword to get image scale */ //Disabled at the moment
			if (fits_read_key_dbl(fptr, "SCALE", &in_mas_pixel, dummy_char, &status)) {
				printf("Initial image -- no SCALE keyword in initial FITS image\n");
				printf("Initial image -- pixel scale assumed to be same as current image\n");
				in_mas_pixel = mas_pixel;
				status = 0;
			} else {
				if (fabs(mas_pixel - in_mas_pixel) / fabs(mas_pixel) > 1e-3) {
					printf("Initial image -- WARNING init image scale: %lf current %lf\n", in_mas_pixel, mas_pixel);
					getchar();
				}
			}

			/* read in model parameters, if there are any... */
			if (nparams > 0) {
				status = 0;
				for (i = 0; i < nparams; i++) {
					sprintf(param_string, "MNPARAM%1ld", i + 1);
					if (!fits_read_key_dbl(fptr, param_string, &init_params[i], dummy_char, &status)) {
						status = 0;
						printf("Initial image -- Could not find MNPARAM%1ld in header\n", i);
					}
				}
				status = 0;
			}

			initial_image = malloc(in_naxes[0] * in_naxes[0] * nwavi * sizeof(double));

			if (fits_read_img(fptr, TDOUBLE, 1, in_naxes[0] * in_naxes[0] * nwavi, &nullval, initial_image, &dummy_int, &status))
				printerror(status);

			// Check flux normalization
			for (w = 0; w < nwavi; w++) {
				ftot = 0.0;
				for (j = 0; j < in_naxes[0]; j++)
					for (i = 0; i < in_naxes[0]; i++)
						ftot += initial_image[i + j * in_naxes[0] + w * in_naxes[0] * in_naxes[0]];
				printf("Initial image -- checking total flux : %lf in channel %ld\n", ftot, w);
				// TODO: renormalize image
				if (fabs((ftot - 1.0)) > 1e-4) {
					printf("Initial image -- renormalizing channel %ld\n", w);
					for (j = 0; j < in_naxes[0]; j++)
						for (i = 0; i < in_naxes[0]; i++)
							initial_image[i + j * in_naxes[0] + w * in_naxes[0] * in_naxes[0]] /= ftot;
				}
			}

			initial_x = malloc(nwavi * nelements * sizeof(unsigned short));
			initial_y = malloc(nwavi * nelements * sizeof(unsigned short));

			for (w = 0; w < nwavi; w++) {
				/* 1) Dummy int counts through the elements
				 *  2) i,j count through the pixels
				 *  3) dtemp counts up flux in the input image.
				 *  Once we get to 1 element, we add this in */
				dummy_int = 0;
				if (offset < 0)
					i = -offset;
				else
					i = 0;
				j = i;
				dtemp = 0.0;
				/* i,j are in the units of the image (i.e. with in_naxes[0] the width)
				 *  initial_x and initial_y have to be set in units of axis_len. dtemp accumulates fractional flux. */

				while (dummy_int < nelements) {
					if (dtemp >= 1) {
						dtemp -= 1.0; /* We've used up one unit of flux */
						initial_x[w * nelements + dummy_int] = i + offset;
						initial_y[w * nelements + dummy_int] = j + offset;
						dummy_int++;
					} else {
						i++;
						/* If in_naxes[0] > axis_len, then we want to start mid-way through the image
						 *      In axis_len > in_naxes[0], then i,j must always be greater than zero */
						if (i == in_naxes[0] || ((offset < 0) && (i == axis_len - offset))) {
							if (offset < 0)
								i = -offset;
							else
								i = 0;
							j++;
							if (j == in_naxes[0] || (offset < 0 && j == axis_len - offset)) {
								if (offset < 0)
									j = -offset;
								else
									j = 0;
							}
						}
						dtemp += initial_image[w * in_naxes[0] * in_naxes[0] + i + j * in_naxes[0]] * nelements; /* Add the flux from the next pixel...*/
					}
				}
			}

			free(in_naxes);
			free(initial_image);
		}

	} else /* If we didn't define a starting image, the default is a dirac*/
	{
		nwavi = 1;
		initial_x = malloc(nelements * sizeof(unsigned short));
		initial_y = malloc(nelements * sizeof(unsigned short));
		for (i = 0; i < nelements; i++) {
			initial_x[i] = axis_len / 2;
			initial_y[i] = axis_len / 2;
		}
	}

	/* Defaults is to have the model parameters equal to 0 and fixed. */
	if (nparams == 0) {
		for (i = 0; i < MAX_PARAMS; i++)
			init_params[i] = 0.0;
		for (i = 0; i < MAX_PARAMS; i++)
			init_stepsize[i] = 0.0;
	}

	/* Initialize the prior image */
	if (prior_filename[0] != 0) {
		status = 0;
		printf("Prior image -- trying to open: '%s'\n", prior_filename);
		if (fits_open_file(&fptr, prior_filename, READONLY, &status))
			printerror(status);

		fits_read_key_lng(fptr, "NAXIS", &k, dummy_char, &status);

		if ((nwavr > 1) && (k == 3))
			in_naxes = malloc(3 * sizeof(long));
		else
			in_naxes = malloc(2 * sizeof(long));

		if (fits_read_keys_lng(fptr, "NAXIS", 1, 2, in_naxes, &nfound, &status))
			printerror(status);

		if ((nwavr > 1) && (k == 3)) {
			if (in_naxes[2] == nwavr) // polychromatic reconstruction using polychromatic init
				nwavp = in_naxes[2];
		} else
			nwavp = 1;

		if ((in_naxes[0] != in_naxes[1]) || (in_naxes[0] != axis_len)) {
			printf("Prior image -- Dimensions do not match other settings.");
			reg_param[REG_PRIORIMAGE] = 0;
		} else {
			prior_image = malloc(nwavp * axis_len * axis_len * sizeof(double));

			if (fits_read_img(fptr, TDOUBLE, 1, nwavp * in_naxes[0] * in_naxes[0], &nullval, prior_image, &dummy_int, &status))
				printerror(status);
			/* What we really mean by a prior is */
			for (w = 0; w < nwavp; w++)
				for (i = 0; i < axis_len; i++)
					for (j = 0; j < axis_len; j++)
						if (prior_image[i + axis_len * j + axis_len * axis_len * w] > 0)
							prior_image[i + axis_len * j + axis_len * axis_len * w] = -log(prior_image[i + axis_len * j + axis_len * axis_len * w]);
						else
							prior_image[i + axis_len * j + axis_len * axis_len * w] = 1e12;
			if (reg_param[REG_PRIORIMAGE] == 0.)
				reg_param[REG_PRIORIMAGE] = 1.;
		}
		free(in_naxes);
	}

	/* Set the derived parameters - number of degrees of freedom and the chi^2 for a random image... */
	/* Note that adding nparams to ndf is only valid if the parameters are free AND
	 *    there is some a-priori information for each parameter in reg_value[REG_MODELPARAM]... */
	ndf = (double) (nvisamp + nvisphi + nt3amp + nt3phi + nv2 + nparams);

	/* cent_mult should only be non-zero if there is no model to fit */
	if (nparams == 0) {
		cent_mult = DEFAULT_CENT_MULT;
		reg_param[REG_CENTERING] = 1.0;
		ndf += 2.;
	}

	printf("Reconst setup -- Degrees of freedom:\t%ld\n", (long) round(ndf));

	flat_chi2 = get_flat_chi2(benchmark);
	printf("Reconst setup -- Chi2r random image:\t%lf\n", flat_chi2 / ndf);

	/* Print out important parameters */
	printf("Reconst setup -- Number of elements:\t%ld\n", nelements);
	printf("Reconst setup -- Number of iterations:\t%ld\n", niter);
	if (depth > niter)
		depth = niter;
	printf("Reconst setup -- Depth of final image:\t%ld\n", depth);

	for (i = 0; i < nparams; i++)
		printf("Reconst setup -- Parametric model: parameter %2ld: %le, with stepsize: %lf \n", i, init_params[i], init_stepsize[i]);
	fflush(stdout);

	/* Now make big matrix - we'll just make this a big chunk of
	 *    memory. */
	double complex
	* __restrict xtransform = malloc(axis_len * nuv * sizeof(double complex));
	double complex
	* __restrict ytransform = malloc(axis_len * nuv * sizeof(double complex));
	for (k = 0; k < nuv; k++) {
		// dtemp = exp ( -M_PI * M_PI / 4.0 / log ( 2 ) * ( u[k] * u[k] + v[k] * v[k] ) / MAS_RAD / MAS_RAD  * ( cvfwhm * cvfwhm ) );

		if (use_bandwidthsmearing == TRUE) {
			for (j = 0; j < axis_len; j++) {
				xtransform[j * nuv + k] = sinc(uv_dlambda[k] / uv_lambda[k] * (j - axis_len / 2) * u[k] * mas_pixel / MAS_RAD)
						* cexp(I * (j - axis_len / 2) * 2.0 * M_PI * u[k] * mas_pixel / MAS_RAD);
				ytransform[j * nuv + k] = sinc(uv_dlambda[k] / uv_lambda[k] * (j - axis_len / 2) * v[k] * mas_pixel / MAS_RAD)
						* cexp(I * (j - axis_len / 2) * 2.0 * M_PI * -v[k] * mas_pixel / MAS_RAD);
			}
		} else {
			for (j = 0; j < axis_len; j++) {
				xtransform[j * nuv + k] = cexp(I * (j - axis_len / 2) * 2.0 * M_PI * u[k] * mas_pixel / MAS_RAD);
				ytransform[j * nuv + k] = cexp(I * (j - axis_len / 2) * 2.0 * M_PI * -v[k] * mas_pixel / MAS_RAD);
			}
		}
	}

	// Shared OpenMP memory

	unsigned short *burn_in_times = calloc(nchains, sizeof(unsigned short));
	double *lLikelihood_expectation = calloc(nchains, sizeof(double));
	double *lLikelihood_deviation = calloc(nchains, sizeof(double));
	double *saved_lLikelihood = calloc(nchains * niter, sizeof(double));
	double *saved_lPosterior = calloc(nchains * niter, sizeof(double));
	double *saved_lPrior = calloc(nchains * niter, sizeof(double));
	double *saved_reg_value = calloc(nchains * niter * NREGULS, sizeof(double));
	double *saved_params = calloc(nchains * niter * nparams, sizeof(double));
	unsigned short *saved_x = calloc(nchains * niter * nwavr * nelements, sizeof(unsigned short));
	unsigned short *saved_y = calloc(nchains * niter * nwavr * nelements, sizeof(unsigned short));
	double *temperature = calloc(nchains, sizeof(double));
	unsigned short *iChaintoStorage = calloc(nchains, sizeof(unsigned short)); // determines where each parallel MCMC data will be stored
	// iChaintoStorage[i] also gives the index in the sorted list of temperatures
	// of the i-st lowest temperature
	unsigned short *iStoragetoChain = calloc(nchains, sizeof(unsigned short)); // determines where each storage/temperature is ; converse of iChaintoStorage

	unsigned short *iMovedChain = calloc(nchains, sizeof(unsigned short));

	for (i = 0; i < nchains; i++)
		burn_in_times[i] = niter; // for ENGINE_SIMULATED_ANNEALING, unless T gets to tmin, burn-in is never achieved

#ifdef _OPENMP
	omp_set_num_threads(nchains);
	//omp_set_dynamic(0);
	//omp_set_nested(1); // we allow nested parallelism, typical application is if you have a low number of chains and lots of CPU cores
	// Start nchains MCMC
	//
#pragma omp parallel private(i,j,k,w) \
		shared(temperature, iChaintoStorage, iStoragetoChain, iMovedChain, burn_in_times, saved_x, saved_y, saved_lLikelihood, \
				saved_lPosterior, saved_lPrior, saved_params, saved_reg_value, minimization_engine, use_tempfitswriting,init_filename, \
				ctrlcpressed, f_anywhere, f_copycat, prob_auto, tmin, chi2_target, mas_pixel, niter, chi2_temp, flat_chi2, \
				axis_len, lLikelihood_expectation, lLikelihood_deviation , nwavr, nelements, nchains, tempschedc, \
				uvwav2chan, uvtime2chan, nuv,nv2,nt3amp,nt3phi,nvisamp,nvisphi,init_params,init_stepsize,initial_x,initial_y, \
				reg_param, prior_image,cent_mult,fov,nparams,xtransform,ytransform, ndf)
#endif
	{
		/* The current system state */

#ifdef _OPENMP
		int iChain = omp_get_thread_num();
		//	omp_set_nested(1);
#else
		int iChain = 0;
#endif
		//	omp_set_num_threads(1);
		//	printf("Subchains available %d\n", omp_get_max_threads());
		//getchar();
		if (minimization_engine == ENGINE_SIMULATED_ANNEALING)
			temperature[iChain] = 1.0;
		else {
			temperature[iChain] = pow(((double) nchains) / (1.0 * (nchains - iChain)), tempschedc); // T[0] = 1, T[nchains-1] = nchains^tempschedc
			printf("Reconst setup -- Chain %d has temperature %f\n", iChain, temperature[iChain]);
		}
#pragma omp barrier

		iChaintoStorage[iChain] = iChain; // initally, the storage unit for chain N is saved_xxx[N * ...],
		iStoragetoChain[iChain] = iChain; // initally,  storage[N] has temperature[N]
		iMovedChain[iChain] = 0; // no threads have been moved

		long chan = 0, rlong, xstep = 0, ystep = 0, steptype = STEP_MEDIUM;
		unsigned short new_x = 0, new_y = 0, old_x = 0, old_y = 0;
		long old_pos = 0, new_pos = 0;
		long current_elt;

		double lPosterior = 0, new_lPosterior, lPrior = 0, new_lPrior, transition_test;
		double lLikelihood = 0, new_lLikelihood;
		//double change_UDreg, change_TV;

		double prob_movement = MPROB_LOW;
		double chi2v2, chi2t3amp, chi2t3phi, chi2visamp, chi2visphi;
		double complex
		*dummy_cpointer = NULL;

		double* image = malloc(nwavr * axis_len * axis_len * sizeof(double));
		unsigned short* element_x = malloc(nwavr * nelements * sizeof(unsigned short));
		unsigned short* element_y = malloc(nwavr * nelements * sizeof(unsigned short));
		double complex
		*im_vis = malloc(nuv * sizeof(double complex));
		double complex
		*new_im_vis = malloc(nuv * sizeof(double complex));
		double complex
		*new_mod_vis = malloc(nuv * sizeof(double complex));
		double complex
		*mod_vis = malloc(nuv * sizeof(double complex));
		double complex
		*param_vis = malloc(nuv * sizeof(double complex));
		double complex
		*new_param_vis = malloc(nuv * sizeof(double complex));
		double *params = malloc(MAX_PARAMS * sizeof(double));
		double *new_params = malloc(MAX_PARAMS * sizeof(double));
		double *stepsize = malloc(MAX_PARAMS * sizeof(double));
		double *prob_pmovement = malloc(MAX_PARAMS * sizeof(double));
		double *res = malloc((nv2 + nt3amp + nt3phi + nvisamp + nvisphi) * sizeof(double)); // current residuals
		double *mod_obs = malloc((nv2 + nt3amp + nt3phi + nvisamp + nvisphi) * sizeof(double)); // current observables
		double *centroid_image_x = malloc(nwavr * sizeof(double));
		double *centroid_image_y = malloc(nwavr * sizeof(double));
		double *reg_value = malloc(nwavr * NREGULS * sizeof(double));
		double *new_reg_value = malloc(nwavr * NREGULS * sizeof(double));
		double *fluxratio_image = malloc(nuv * sizeof(double));
		double *new_fluxratio_image = malloc(nuv * sizeof(double));
		unsigned short chain1, chain2;
		double logZ = 0; // chose to have logZ to be a private variable
		double logZ_err = 0;

		// Randomization

		char rngname[80];
		sprintf(rngname, "rng%02d", iChain);
		RngStream rng = RngStream_CreateStream(rngname); // will be multithreaded

		/* Initialize the filename for temporarily saving*/
		char temp_filename[80];
		sprintf(temp_filename, "chain%02d", iChain);

		for (w = 0; w < nwavr; w++) {
			centroid_image_x[w] = 0.0;
			centroid_image_y[w] = 0.0;
		}

		for (w = 0; w < nwavr; w++)
			for (i = 0; i < NREGULS; i++) {
				reg_value[w * NREGULS + i] = 0;
				new_reg_value[w * NREGULS + i] = 0;
			}

		/* Initialise the model parameters, if any */
		for (i = 0; i < MAX_PARAMS; i++) {
			params[i] = init_params[i];
			new_params[i] = init_params[i];
			stepsize[i] = init_stepsize[i];
			prob_pmovement[i] = 0.5;
		}

		// Set up the model vs image chromatic flux ratios (one for each uv point !)
		for (i = 0; i < nuv; i++) {
			fluxratio_image[i] = 1.0;
			new_fluxratio_image[i] = 1.0;
		}

		//
		// INITIALIZE IMAGE
		//
		initialize_image(iChain, image, element_x, element_y, initial_x, initial_y, axis_len, nwavr, nelements, &init_filename[0]);

		//
		// COMPUTE INITIAL REGULARIZER VALUE
		//
		compute_regularizers(reg_param, reg_value, image, prior_image, (double) nelements, initial_x, initial_y, nwavr, axis_len, nelements, centroid_image_x,
				centroid_image_y, fov, cent_mult);

		//
		// COMPUTE INITIAL VISIBILITIES
		//
		compute_model_visibilities_fromelements(mod_vis, im_vis, param_vis, params, fluxratio_image, element_x, element_y, xtransform, ytransform,
				&reg_value[REG_MODELPARAM], nparams, nelements);

		//Compute initial values for prior, likelihood, and posterior
		compute_lLikelihood(&lLikelihood, mod_vis, res, mod_obs, &chi2v2, &chi2t3amp, &chi2visamp, &chi2t3phi, &chi2visphi);
		compute_lPrior(&lPrior, chan, reg_param, reg_value);
		lPosterior = lLikelihood + lPrior;

		if (minimization_engine == ENGINE_SIMULATED_ANNEALING) {
			/* Temperature is based on lLikelihood */
			temperature[iChain] = 2.0 * lLikelihood / ndf / INITIAL_SCALED_CHI2_MULT / chi2_temp;

			if (temperature[iChain] > FLAT_CHI2_MULT * flat_chi2 / ndf)
				temperature[iChain] = FLAT_CHI2_MULT * flat_chi2 / ndf;

			if (temperature[iChain] < INITIAL_TMIN) {
				temperature[iChain] = INITIAL_TMIN;
				burn_in_times[iChain] = niter / BURN_IN_FRAC;
			}
		}

		/*--------------------------
		 *    This is the main loop...
		 *    ----------------------------*/
		for (i = 0; i < niter * nwavr * nelements; i++) {
			// Reset new regularizer values
			for (w = 0; w < nwavr; w++) {
				for (j = 0; j < NREGULS; j++) {
					new_reg_value[w * NREGULS + j] = reg_value[w * NREGULS + j];
				}
			}

			if ((i % (nwavr * nelements)) == 0) {
				chain1 = iChaintoStorage[iChain];
				// Save (x,y) element positions and probabilities
#pragma omp critical(savestate) // not sure if this is really needed
				{
					// because we entered this section, this means that
					// i is a multiple of nwavr * nelements,
					// and the actual iteration number is i / (nwavr * nelements)
					for (w = 0; w < nwavr; w++) {
						for (j = 0; j < nelements; j++) {

							// Note: in parallel tempering mode, iChaintoStorage[iChain] tracks which chain correspond to saving slot
							//       in simulated annealing mode, iChaintoStorage[iChain] is just = iChain
							// i = (0...niter - 1) * nwavr * nelements
							saved_x[chain1 * nwavr * nelements * niter + w * nelements + i + j] = element_x[w * nelements + j];
							saved_y[chain1 * nwavr * nelements * niter + w * nelements + i + j] = element_y[w * nelements + j];
						}

					}

					saved_lLikelihood[chain1 * niter + i / (nwavr * nelements)] = lLikelihood;
					saved_lPrior[chain1 * niter + i / (nwavr * nelements)] = lPrior;
					saved_lPosterior[chain1 * niter + i / (nwavr * nelements)] = lPosterior;

					for (j = 0; j < nparams; j++)
						saved_params[chain1 * nparams * niter + i / (nwavr * nelements) * nparams + j] = params[j];

					for (j = 0; j < NREGULS; j++)
						saved_reg_value[chain1 * NREGULS * niter + i / (nwavr * nelements) * NREGULS + j] = reg_value[j];
				}

				if ((minimization_engine == ENGINE_PARALLEL_TEMPERING) && (i / (nwavr * nelements) > 0)) // note: prevent swapping states before the second iteration
						{

					// Marginal likelihood estimation -- Step 1 = compute likelihood expectation as MCMC averages
					lLikelihood_expectation[iChain] = 0;
					lLikelihood_deviation[iChain] = 0;
#pragma omp barrier   // synchronize chains needed here to prevent swapping while computing the averages
					if (i / (nwavr * nelements) > ceil(0.3 * niter)) // note: we need to check for burn-in info here instead
							{
						for (j = ceil(0.3 * niter); j < i / (nwavr * nelements); j++) {
							// we will average the likelihood for the current chain temperature
							lLikelihood_expectation[iChain] += saved_lLikelihood[iChaintoStorage[iChain] * niter + j];
						}
						lLikelihood_expectation[iChain] /= (i / (nwavr * nelements) - ceil(0.3 * niter));

						for (j = ceil(0.3 * niter); j < i / (nwavr * nelements); j++) {
							// we will average the likelihood for the current chain temperature
							lLikelihood_deviation[iChain] += (saved_lLikelihood[iChaintoStorage[iChain] * niter + j] - lLikelihood_expectation[iChain])
									* (saved_lLikelihood[iChaintoStorage[iChain] * niter + j] - lLikelihood_expectation[iChain]);
						}
						if (i / (nwavr * nelements) > ceil(0.3 * niter) + 1)
							lLikelihood_deviation[iChain] /= (i / (nwavr * nelements) - ceil(0.3 * niter) - 1);

					}

#pragma omp barrier   // synchronize chains
					//    printf("Chain: %d temperature: %lf iter: %ld lLikelihood_avg: %f \n",
					//    iChain, temperature[iChain], i / (nwavr * nelements), lLikelihood_expectation[iChain]);
#pragma omp single
					{
						// Marginal likelihood estimation -- Step 2 = compute log Z
						compute_logZ(temperature, iStoragetoChain, lLikelihood_expectation, lLikelihood_deviation, nchains, &logZ, &logZ_err);
						printf("log Z computed by chain: %d logZ: %f +/- %f \n", iChain, logZ, logZ_err);

					}

					//
					// Inter-chain exchanges
					//

					//  printf("Chain %d is ready to switch -- now idle at iteration i= %ld\n", iChain, i);
					fflush(stdout);
					iMovedChain[chain1] = 0; // set all chains as ready to switch
#pragma omp barrier
#pragma omp critical(chainswap)
					{
						// Inter-chain  Metropolis-Hastings temperature swap
						// 1) only one chain should try to swap temperature at the same time, hence the omp critical statement above
						// 2) we only swap adjacent temperatures. There is no warranty that iChain and iChain+1 have adjacent temperatures,
						//    so we actually do not try to switch these. Instead we keep track of the temperatures.
						//    iStoragetoChain[N] give us the index of the chain that currently has the Nth temperature (in the sorted list of increasing temperatures)
						//    so iStoragetoChain[0] gives use the lowest temperature, i.e. T=1 for typical parallel tempering
						//    and iStoragetoChain[1] would give us the next lowest temperature
						//    the following code is set so we only try to swap iStoragetoChain[iChain] with iStoragetoChain[iChain+1]
						// 3) we also track which chains have been swapped through iMovedChain, to prevent sequential swapping (0->1 then 1->2)

						if (iChain < (nchains - 1)) {
							chain1 = iStoragetoChain[iChain];
							chain2 = iStoragetoChain[iChain + 1];

							if ((iMovedChain[chain1] == 0) && (iMovedChain[chain2] == 0)) {
								if (RngStream_RandU01(rng) < 0.3) // base probability of swapping is 30%
										{

									transition_test = (1. / temperature[chain1] - 1. / temperature[chain2])
											* (saved_lPosterior[iChaintoStorage[chain1] * niter + i / (nwavr * nelements)]
													- saved_lPosterior[iChaintoStorage[chain2] * niter + i / (nwavr * nelements)]);

									if (log(RngStream_RandU01(rng)) < transition_test) {
										printf("Swap called from chain %d -- will be swapping chains %d : %f and %d : %f at iteration %ld\n", iChain, chain1,
												temperature[chain1], chain2, temperature[chain2], i / (nwavr * nelements));

										//                   for(j=0;j<nchains;j++)
										// printf("BEFORE Chain %ld \t iChaintoStorage %d \t iStoragetoChain %d \t temp[j] %lf \t temp[iStorage[j]] %lf \t temp[iChain[j] %lf\n", j, iChaintoStorage[j], iStoragetoChain[j], temperature[j], temperature[iStoragetoChain[j]],temperature[iChaintoStorage[j]] );

										//Swap temperatures and other relevant arrays
										swapd(&temperature[chain1], &temperature[chain2]);
										swapi(&iStoragetoChain[iChaintoStorage[chain1]], &iStoragetoChain[iChaintoStorage[chain2]]);
										swapi(&iChaintoStorage[chain1], &iChaintoStorage[chain2]);

										// printf("Swap done at chain %d -- has swapped chains %d : %f and %d : %f at iteration %ld\n",
										//         iChain, chain2, temperature[chain2], chain1,  temperature[chain1], i / (nwavr * nelements));

										// for(j=0;j<nchains;j++)
										//  printf("AFTER Chain %ld \t iChaintoStorage %d \t iStoragetoChain %d \t temp[j] %lf \t temp[iStorage[j]] %lf \t temp[iChain[j] %lf\n", j, iChaintoStorage[j], iStoragetoChain[j], temperature[j], temperature[iStoragetoChain[j]],temperature[iChaintoStorage[j]] );

										iMovedChain[chain1] = 1;
										iMovedChain[chain2] = 1;
									}
								}
							}
						}
					}
#pragma omp barrier   //re-synchronise Chains before continuing
					//fflush(stdout);
					//printf("%d Chain %d movement reset to 0 iter %ld\n", iChain, chain1, i / (nwavr * nelements) );
				}

			}

			//
			// DISPLAY DIAGNOSTICS
			//
			if ((i % (STEPS_PER_OUTPUT * nwavr * nelements)) == 0) {
				if (use_tempfitswriting == TRUE)
				  writeasfits(temp_filename, image, nwavr, 1, (iChain) * niter + i / (nelements * nwavr), 2.0 * lLikelihood / ndf, chi2v2, chi2t3amp, chi2t3phi, chi2visamp, chi2visphi,
					                temperature[iChain], nelements,	&reg_param[0], &reg_value[0], niter, axis_len, ndf, tmin,
					                chi2_temp, chi2_target, mas_pixel, nchains, 0, 0, "", "",
							&saved_params[iChaintoStorage[iChain] * nparams * niter], NULL);

				// PRINT DIAGNOSTICS
				print_diagnostics(iChain, (i / (nwavr * nelements) + 1), nvis, nv2, nt3, nt3phi, nt3amp, nvisamp, nvisphi, chi2v2, chi2t3amp, chi2t3phi,
						chi2visphi, chi2visamp, lPosterior, lPrior, lLikelihood, reg_param, reg_value, centroid_image_x, centroid_image_y, nelements, nwavr,
						niter, temperature, prob_movement, params, stepsize);

				if (prob_auto > 0)
					tmin = tmin * (1.0 - .5 * (prob_movement - prob_auto)); // BUG ? should this be here ?

			}

			//
			// Generate the main random number
			//
			rlong = RngStream_RandInt(rng, 0, 2147483647);

			//
			// Select wavelength channel based on iteration number or random number ???
			//
			if (nwavr > 1) {
				// alternate
				chan = rlong % nwavr;
				rlong /= nwavr;
				//printf("%ld", chan);
			} else
				chan = 0;

			current_elt = rlong % (nelements + nparams * PARAMS_PER_ELT);
			rlong /= (nelements + nparams * PARAMS_PER_ELT);

			if ((nparams == 0) || (current_elt < nelements)) /* Attempt image movement rather than parametric model movement */
			{

				// select step
				switch (steptype) {
				case STEP_SMALL: /* One step in either x or y direction */
					xstep = rlong % 4;
					rlong /= 4;
					if (xstep > 1)      // 2 or 3, we select a y move
							{
						ystep = (xstep - 2) * 2 - 1;
						xstep = 0;
					} else     // 0 or 1, we select a x movement
					{
						ystep = 0;
						xstep = xstep * 2 - 1;
					}
					break;
				case STEP_MEDIUM: /* Up to 2 steps in x and y direction */
					xstep = (rlong % 5) - 2;
					rlong /= 5;
					ystep = (rlong % 5) - 2;
					rlong /= 5;
					break;
				case STEP_ANYWHERE: /* Move flux anywhere in image */
					//xstep = (rlong % axis_len) - axis_len / 2;
					// rlong = RngStream_RandInt(rng, 0, 2147483647);
					// ystep = (rlong % axis_len) - axis_len / 2;
					// rlong /= axis_len;
					xstep = (RngStream_RandInt(rng, 0, 2147483647) % axis_len) - axis_len / 2;
					// rlong = RngStream_RandInt(rng, 0, 2147483647);
					ystep = (RngStream_RandInt(rng, 0, 2147483647) % axis_len) - axis_len / 2;
					// rlong /= axis_len;
					break;

				case STEP_COPYCAT: /* Copy the position of another element */
					xstep = element_x[chan * nelements + rlong % nelements] - element_x[chan * nelements + current_elt];
					ystep = element_y[chan * nelements + rlong % nelements] - element_y[chan * nelements + current_elt];
					rlong = RngStream_RandInt(rng, 0, 2147483647);
					break;
				}

				old_x = element_x[chan * nelements + current_elt];
				old_y = element_y[chan * nelements + current_elt];
				old_pos = chan * axis_len * axis_len + old_y * axis_len + old_x;

				new_x = (element_x[chan * nelements + current_elt] + xstep + axis_len) % axis_len;
				new_y = (element_y[chan * nelements + current_elt] + ystep + axis_len) % axis_len;
				new_pos = chan * axis_len * axis_len + new_y * axis_len + new_x;

				//
				// Regularization
				//

				if (reg_param[REG_ENTROPY] > 0.0)
					new_reg_value[chan * NREGULS + REG_ENTROPY] = reg_value[chan * NREGULS + REG_ENTROPY] - entropy(image[old_pos])
							+ entropy(image[old_pos] - 1);
				if (reg_param[REG_PRIORIMAGE] > 0.0)
					new_reg_value[chan * NREGULS + REG_PRIORIMAGE] = reg_value[chan * NREGULS + REG_PRIORIMAGE] - prior_image[old_pos];
				if (reg_param[REG_DARKENERGY] > 0.0)
					new_reg_value[chan * NREGULS + REG_DARKENERGY] = reg_value[chan * NREGULS + REG_DARKENERGY]
							- den_change(image, old_x, old_y, DEN_SUBTRACT, axis_len);
				if (reg_param[REG_SPOT] > 0.0)
					reg_value[chan * NREGULS + REG_SPOT] = UDreg(&image[chan * axis_len * axis_len], NULL, 0.0, axis_len, axis_len, (const double) nelements);
				if (reg_param[REG_TV] > 0.0)
					reg_value[chan * NREGULS + REG_TV] = TV(&image[chan * axis_len * axis_len], NULL, 0.0, axis_len, axis_len, (const double) nelements);
				if (reg_param[REG_LAP] > 0.0)
					reg_value[chan * NREGULS + REG_LAP] = LAP(&image[chan * axis_len * axis_len], NULL, 0.0, axis_len, axis_len, (const double) nelements);
				if (reg_param[REG_L0] > 0.0)
					reg_value[chan * NREGULS + REG_L0] = L0(&image[chan * axis_len * axis_len], NULL, 0.0, axis_len, axis_len, (const double) nelements);
				if (reg_param[REG_TRANSPECL2] > 0.0)
					reg_value[REG_TRANSPECL2] = transpec(nwavr, axis_len, image, (const double) nelements);

				image[old_pos]--;

				if (reg_param[REG_ENTROPY] > 0.0)
					new_reg_value[chan * NREGULS + REG_ENTROPY] += entropy(image[new_pos] + 1) - entropy(image[new_pos]);
				if (reg_param[REG_PRIORIMAGE] > 0.0)
					new_reg_value[chan * NREGULS + REG_PRIORIMAGE] += prior_image[new_pos];
				if (reg_param[REG_DARKENERGY] > 0.0)
					new_reg_value[chan * NREGULS + REG_DARKENERGY] += den_change(image, new_x, new_y, DEN_ADD, axis_len);

				image[new_pos]++;

				if (reg_param[REG_SPOT] > 0.0)
					new_reg_value[chan * NREGULS + REG_SPOT] = UDreg(&image[chan * axis_len * axis_len], NULL, 0.0, axis_len, axis_len,
							(const double) nelements);  // before
				if (reg_param[REG_TV] > 0.0)
					new_reg_value[chan * NREGULS + REG_TV] = TV(&image[chan * axis_len * axis_len], NULL, 0.0, axis_len, axis_len, (const double) nelements); // before
				if (reg_param[REG_LAP] > 0.0)
					new_reg_value[chan * NREGULS + REG_LAP] = LAP(&image[chan * axis_len * axis_len], NULL, 0.0, axis_len, axis_len, (const double) nelements); // before

				if (reg_param[REG_L0] > 0.0)
					new_reg_value[chan * NREGULS + REG_L0] = L0(&image[chan * axis_len * axis_len], NULL, 0.0, axis_len, axis_len, (const double) nelements); // before

				if (reg_param[REG_TRANSPECL2] > 0.0)
					new_reg_value[REG_TRANSPECL2] = transpec(nwavr, axis_len, image, (const double) nelements);

				// Go back to current state
				image[old_pos]++;
				image[new_pos]--;

				if (reg_param[REG_CENTERING] > 0.0)
					new_reg_value[chan * NREGULS + REG_CENTERING] = reg_value[chan * NREGULS + REG_CENTERING]
							+ fov * cent_change(chan, centroid_image_x, centroid_image_y, new_x, new_y, old_x, old_y, axis_len, fov, cent_mult);
				if (reg_param[REG_MODELPARAM] > 0.0)
					new_reg_value[REG_MODELPARAM] = reg_value[REG_MODELPARAM];

				/* Modify the visibilities */
				//#pragma omp parallel for simd
				for (j = 0; j < nuv; j++) {
					if (uvwav2chan[j] == chan) {
						new_im_vis[j] = im_vis[j]
								+ (xtransform[new_x * nuv + j] * ytransform[new_y * nuv + j] - xtransform[old_x * nuv + j] * ytransform[old_y * nuv + j])
										* fluxratio_image[j] / (double) nelements;
					} else
						new_im_vis[j] = im_vis[j];
					new_mod_vis[j] = new_im_vis[j] + param_vis[j];
				}

				prob_movement *= (1.0 - 1.0 / DAMPING_TIME);

			} else /* Attempt parametric model movement -- model visibilities are recomputed from scratch */
			{

				xstep = rlong % 2; /* The type of step to be taken */
				rlong /= 2;
				j = (current_elt - nelements) / PARAMS_PER_ELT;
				/* !!! This next line shouldn't be required. Given that it is, there must be a bug !!! */
				for (k = 0; k < nparams; k++)
					new_params[k] = params[k];
				new_params[j] = params[j] + stepsize[j] * (xstep * 2.0 - 1.0);
				model_vis(new_params, new_param_vis, &new_reg_value[REG_MODELPARAM], new_fluxratio_image);
				prob_pmovement[j] *= 1.0 - 1.0 / PARAM_DAMPING_TIME;
				stepsize[j] *= 1.0 + (prob_pmovement[j] - TARGET_MPROB) / STEPSIZE_ADJUST_TIME;
				for (j = 0; j < nuv; j++) {
					new_im_vis[j] = im_vis[j] * new_fluxratio_image[j] / fluxratio_image[j];
					new_mod_vis[j] = new_im_vis[j] + new_param_vis[j];
				}
			}

			//
			// Evaluate posterior probability
			//

			compute_lLikelihood(&new_lLikelihood, new_mod_vis, res, mod_obs, &chi2v2, &chi2t3amp, &chi2visamp, &chi2t3phi, &chi2visphi);
			compute_lPrior(&lPrior, chan, reg_param, reg_value);
			compute_lPrior(&new_lPrior, chan, reg_param, new_reg_value);
			new_lPosterior = new_lLikelihood + new_lPrior;

			// BUG: think how to rescale priors with fluxratio_image ?
			transition_test = (new_lLikelihood - lLikelihood) / temperature[iChain] + new_lPrior - lPrior;

			if ((double) (rlong % 1024) + 0.5 < 1024.0 * exp(-transition_test)) {
				// We accept the new state

				lLikelihood = new_lLikelihood;
				lPosterior = new_lPosterior;
				lPrior = new_lPrior;
				reg_value[chan * NREGULS + REG_CENTERING] = new_reg_value[chan * NREGULS + REG_CENTERING];

				/* Swap mod_vis and new_mod_vis*/
				dummy_cpointer = mod_vis;
				mod_vis = new_mod_vis;
				new_mod_vis = dummy_cpointer;

				/* Swap im_vis and new_im_vis */
				dummy_cpointer = im_vis;
				im_vis = new_im_vis;
				new_im_vis = dummy_cpointer;

				/* For flux movement, update the image and entropy */
				if (current_elt < nelements) {
					image[old_pos]--;
					element_x[chan * nelements + current_elt] = new_x;
					element_y[chan * nelements + current_elt] = new_y;
					image[new_pos]++;

					reg_value[chan * NREGULS + REG_ENTROPY] = new_reg_value[chan * NREGULS + REG_ENTROPY];
					reg_value[chan * NREGULS + REG_DARKENERGY] = new_reg_value[chan * NREGULS + REG_DARKENERGY];
					reg_value[chan * NREGULS + REG_PRIORIMAGE] = new_reg_value[chan * NREGULS + REG_PRIORIMAGE];
					reg_value[chan * NREGULS + REG_TV] = new_reg_value[chan * NREGULS + REG_TV];
					reg_value[chan * NREGULS + REG_LAP] = new_reg_value[chan * NREGULS + REG_LAP];
					reg_value[chan * NREGULS + REG_L0] = new_reg_value[chan * NREGULS + REG_L0];
					reg_value[chan * NREGULS + REG_SPOT] = new_reg_value[chan * NREGULS + REG_SPOT];
					reg_value[chan * NREGULS + REG_CENTERING] = new_reg_value[chan * NREGULS + REG_CENTERING];
					reg_value[REG_TRANSPECL2] = new_reg_value[REG_TRANSPECL2];
					prob_movement += 1.0 / DAMPING_TIME;
				}
				/* For parameter movement, update parameters */
				else {
					j = (current_elt - nelements) / PARAMS_PER_ELT;
					prob_pmovement[j] += 1.0 / PARAM_DAMPING_TIME;
					/* Swap param_vis and new_param_vis */
					dummy_cpointer = param_vis;
					param_vis = new_param_vis;
					new_param_vis = dummy_cpointer;
					/* Swap params and new_params */
					params[j] = new_params[j];
					/* Update fluxratio_image and reg_value[REG_MODELPARAM] */
					for (j = 0; j < nuv; j++)
						fluxratio_image[j] = new_fluxratio_image[j];

					reg_value[REG_MODELPARAM] = new_reg_value[REG_MODELPARAM];
				}

				if (minimization_engine == ENGINE_SIMULATED_ANNEALING) {
					/* If chi^2 has changed, change the temperature... */
					temperature[iChain] *= 1.0
							+ 1.0 / temperature[iChain] / TEMP_CHANGE_TIME * (2. * lLikelihood / ndf - chi2_temp * temperature[iChain])
									* (2. * lLikelihood / ndf - chi2_target) * ndf / (2. * lLikelihood);

					if ((temperature[iChain] < tmin) || (2. * lLikelihood / ndf < chi2_target)) // convergence
						if ((i / (nwavr * nelements) + niter / BURN_IN_FRAC) < burn_in_times[iChain])
							burn_in_times[iChain] = i / (nwavr * nelements) + niter / BURN_IN_FRAC;

					if (temperature[iChain] < tmin)
						temperature[iChain] = tmin;
					if (temperature[iChain] > FLAT_CHI2_MULT * flat_chi2 / ndf)
						temperature[iChain] = FLAT_CHI2_MULT * flat_chi2 / ndf;
				}

			} else {
				// The proposed new state was rejected
				if (current_elt < nelements) // flux movement was attempted and rejected
						{
					//if((xstep == 0) && (ystep == 0))
					//     printf("Rejected movement of 0 steps, delta_chi2 = %lf: \n", 2.*(new_lLikelihood - lLikelihood) );

					/* Reverse the centering change */
					if (reg_param[REG_CENTERING] > 0)
						cent_change(chan, centroid_image_x, centroid_image_y, old_x, old_y, new_x, new_y, axis_len, fov, cent_mult);
				} else     // parameter movement was attempted and rejected
				{
					// error checking
					if (stepsize[(current_elt - nelements) / PARAMS_PER_ELT] == 0)
						printf("Failed parameter movement of 0 stepsize! %ld %6.3lf Iter: %ld P: %6.4lf %6.3lf \n", (current_elt - nelements) / PARAMS_PER_ELT,
								2. * (new_lLikelihood - lLikelihood), i, params[0], new_params[0]);
				}

			}

			/* Now that prob_movement has changed, we may want to change the step type */
			if ((2. * lLikelihood) < FLAT_CHI2_MULT * flat_chi2) {
				if (steptype == STEP_COPYCAT)
					steptype = STEP_SMALL;

				if (steptype == STEP_ANYWHERE)
					steptype = STEP_SMALL;

				else if (prob_movement > MPROB_HIGH) {
					steptype++;
					if (steptype > STEP_MAX)
						steptype = STEP_MAX;
				} else if (prob_movement < MPROB_LOW) {
					steptype--;
					if (steptype < 0)
						steptype = 0;
				}
				if (rlong % (int) (1. / f_copycat) == 3)
					steptype = STEP_COPYCAT;
				//if ((steptype == STEP_SMALL) &&
				//     if (rlong % (int)(1./f_anywhere) == 2) steptype=STEP_ANYWHERE;

			} else {
				steptype = STEP_MEDIUM; /* If lLikelihood is too high, fix the step type.*/
			}

			if (rlong % (int) (1. / f_copycat) == 3)
				steptype = STEP_COPYCAT;
			//if ((steptype == STEP_SMALL) &&
			if (rlong % (int) (1. / f_anywhere) == 2)
				steptype = STEP_ANYWHERE;

			/* If we're on the smallest step size, we may want to change to steptype COPYCAT */

			if (ctrlcpressed == TRUE)
				break;
		} // end iterations

		//printf("End of chain %i\n", iChain);

		/* Write the fits file */
		if (ctrlcpressed == FALSE) {
			if (use_tempfitswriting == TRUE)
			  writeasfits(temp_filename, image, nwavr, 1, iChain * niter + i / (nelements * nwavr), 2.0 * lLikelihood / ndf,  chi2v2, chi2t3amp, chi2t3phi, chi2visamp, chi2visphi,
					    temperature[iChain], nelements, &reg_param[0], &reg_value[0], niter, axis_len, ndf, tmin, chi2_temp, chi2_target, mas_pixel, nchains, 0, 0, "", "",
						&saved_params[iChain * niter * nparams], NULL);
		}

		RngStream_DeleteStream(&rng);

		free(centroid_image_x);
		free(centroid_image_y);
		free(reg_value);
		free(new_reg_value);
		free(element_x);
		free(element_y);
		free (im_vis);
		free (new_im_vis);
		free (new_mod_vis);
		free (mod_vis);
		free (param_vis);
		free (new_param_vis);
		free(params);
		free(new_params);
		free(stepsize);
		free(prob_pmovement);
		free(res);
		free(mod_obs);
		free(image);

		free(fluxratio_image);
		free(new_fluxratio_image);

#pragma omp barrier    // Synchronize chains
	}

	//
	// End of parallel chains
	//

	//
	if (ctrlcpressed == FALSE) {

		// Determine the number of usable frames for statistics and image averaging
		// depth is the requested maximum number of usable frames
		long nburned = 0; // actual number of burnt frames (for parallel simulated annealing)
		// iframeburned will contain the index of the frames that have actually been burned,
		// these indexes are for saved_elements and thus go accross all chains
		long *iframeburned = malloc(nchains * niter * sizeof(long));
		double lowest_lLikelihood = 1e99;
		long lowest_lLikelihood_indx = 0;

		//
		// Also allocate memory for the final image and parameters
		//
		double* final_image = calloc(nwavr * axis_len * axis_len, sizeof(double));
		double* final_params = calloc(nparams, sizeof(double));
		double* final_params_std = calloc(nparams, sizeof(double));
		double *final_reg_value = calloc(nwavr * NREGULS, sizeof(double));
		double *centroid_image_x = calloc(nwavr, sizeof(double));
		double *centroid_image_y = calloc(nwavr, sizeof(double));

		if (final_image == NULL) {
			printf("Output -- Not enough memory to allocate final image\n");
		}

		if (minimization_engine != ENGINE_PARALLEL_TEMPERING) // CASE = (parallel) simulated annealing
				{
			for (i = 0; i < nchains; i++) {
				for (j = burn_in_times[i]; j < niter; j++) {
					iframeburned[nburned] = i * niter + j; // we note the location of the frame
					if (saved_lLikelihood[iframeburned[nburned]] < lowest_lLikelihood) {
						lowest_lLikelihood = saved_lLikelihood[iframeburned[nburned]];
						lowest_lLikelihood_indx = iframeburned[nburned];
					}
					nburned++;
				}
				printf("Output -- Burn-in time / burned frames for chain %ld: %d / %ld\n", i, burn_in_times[i], niter - burn_in_times[i]);
			}

			// If no chains got beyond their burn-in, then nburned =0 at this point
			// then we go back through and just use the depth number instead.
			// BUG - to check
			if (nburned == 0) {
				printf("Output -- DID NOT REACH BURN IN ! Writing final image based on depth parameter anyway\n");
				printf("Output -- DID NOT REACH BURN IN ! Be careful when interpreting the final image.\n");
				for (i = 0; i < nchains; i++)
					for (j = niter - depth; j < niter; j++) {
						burn_in_times[i] = niter - depth; // note: we previously ensured depth <= niter
						iframeburned[nburned] = i * nchains + j;
						nburned++;
					}
			}

			mcmc_annealing_results(output_filename, final_image, iframeburned, nburned, nelements, axis_len, xtransform, ytransform,
					       &final_chi2, &final_chi2v2, &final_chi2t3amp, &final_chi2t3phi, &final_chi2visamp, &final_chi2visphi, saved_x, saved_y,
					       saved_params, niter, nwavr, final_params, final_params_std, reg_param, final_reg_value, prior_image, initial_x, initial_y, centroid_image_x,
					       centroid_image_y, fov, cent_mult);

			// Recompute all the regularizers, likelihood, prior & posterior probability
			writeasfits(output_filename, final_image, nwavr, k, niter - burn_in_times[0] - 1,
                        final_chi2/ndf, final_chi2v2, final_chi2t3amp, final_chi2t3phi, final_chi2visamp, final_chi2visphi,
                        -1, nelements, reg_param, final_reg_value, niter, axis_len, ndf, tmin, chi2_temp, chi2_target,
                        mas_pixel, nchains, 0, 0, init_filename, prior_filename, final_params, final_params_std);


		} else // (minimization_engine == ENGINE_PARALLEL_TEMPERING)
		{
			// for(i = 0; i < nchains; i++)
			//      {
			//         printf("Chain: %ld\t temperature: %lf\t ChaintoStorage: %d Storage(%ld)toChain: %d\n", i, temperature[i], iChaintoStorage[i], i, iStoragetoChain[i]);
			//     }

			mcmc_tempering_results(output_filename, final_image, iStoragetoChain[0], nburned, nelements, axis_len, xtransform, ytransform, &final_chi2, saved_x, saved_y, saved_params, niter, nwavr);

			double logZ_final = 0, logZ_final_err = 0.;

			compute_logZ(temperature, iStoragetoChain, lLikelihood_expectation, lLikelihood_deviation, nchains, &logZ_final, &logZ_final_err);
			printf("Output -- Final logZ: %f +/- %f\n", logZ_final, logZ_final_err);

			writeasfits(output_filename, final_image, nwavr, k, niter - burn_in_times[0] - 1,
                    final_chi2/ndf, final_chi2v2, final_chi2t3amp, final_chi2t3phi, final_chi2visamp, final_chi2visphi,
                    -1, nelements, &reg_param[0], NULL, niter, axis_len, ndf,
					tmin, chi2_temp, chi2_target, mas_pixel, nchains, logZ_final, logZ_final_err, init_filename, prior_filename, final_params,
					final_params_std);
		}

		printf("Output -- Total number of realizations in final image: %ld\n", nburned);
		printf("Output -- Global chi2r of final image:\t %lf \n", final_chi2/ndf);
        if(nv2 > 0)
            printf("Output -- chi2r V2   of final image:\t %lf \n", final_chi2v2);
        if(nt3amp > 0)
            printf("Output -- chi2r T3A  of final image:\t %lf \n", final_chi2t3amp);
        if(nt3phi>0)
            printf("Output -- chi2r T3P  of final image:\t %lf \n", final_chi2t3phi);
        if(nvisamp>0)
            printf("Output -- chi2r VISA of final image:\t %lf \n", final_chi2visamp);
		if(nvisphi>0)
            printf("Output -- chi2r VISP of final image:\t %lf \n", final_chi2visphi);


		if (lowest_lLikelihood < 1e99)
			printf("Output -- Best single-frame chi2: %f obtained at iteration: %ld chain: %ld.\n", 2.0 * lowest_lLikelihood / ndf,
					lowest_lLikelihood_indx % niter, lowest_lLikelihood_indx / niter);

		if (dumpchain == TRUE)
			mcmc_fullchain(output_filename, nchains, niter, nwavr, nelements, axis_len, saved_x, saved_y, saved_params, saved_lLikelihood, saved_lPrior,
					saved_lPosterior, temperature, iChaintoStorage);

		free(final_params);
		free(final_params_std);
		free(final_reg_value);
		free(centroid_image_x);
		free(centroid_image_y);
		free(final_image);
		free(iframeburned);

	}

	fflush(stdout);

	free(uvwav2chan);
	free(uvtime2chan);
	free(wavmax);
	free(wavmin);

	free(saved_x);
	free(saved_y);
	free(saved_params);
	free(lLikelihood_expectation);
	free(lLikelihood_deviation);
	free(saved_lLikelihood);
	free(saved_lPrior);
	free(saved_lPosterior);
	free(saved_reg_value);
	free(burn_in_times);
	free(temperature);
	free(iChaintoStorage);
	free(iStoragetoChain);
	free(iMovedChain);
	free (xtransform);
	free (ytransform);
	free(visin);
	free(v2in);
	free(t3in1);
	free(t3in2);
	free(t3in3);
	free(u);
	free(v);
	free(uv_lambda);
	free(uv_dlambda);

	free(dvisnwav);
	if (use_diffvis == TRUE)
		for (i = 0; i < nvis; i++)
			free(dvisindx[i]);
	free(dvisindx);

	free(data);
	free(data_err);

	free(initial_x);
	free(initial_y);
	free(prior_image);
	return 0;
}

/*****************************************************/
/* Print out cfitsio error messages and exit program */
/*****************************************************/
void printerror(int status)
{
	if (status)
        {
		fits_report_error(stderr, status); /* print error report */
		exit(status); /* terminate the program, returning error status */
	    }
	return;
}

/*****************************************************/
/* Print out help text and exits program.            */
/*****************************************************/
void printhelp(void) {
	printf("SQUEEZE: an image reconstruction code for optical interferometry\n\n");
	printf("Usage: squeeze data.oifits -s scale -w width -o image.fits\n\n");
	printf("Options:\n");

	printf("\n***** MAIN IMAGE SETTINGS ***** \n");
	printf("  -s scale       : Size of a pixel in milli-arcseconds.\n");
	printf("  -w width       : Width in pixels.\n");
	printf("  -n iter        : Number of iterations per chain.\n");
	printf("  -e elements    : Number of elements per realisation.\n");
	printf("  -d depth       : Number of realizations to average to together for final mean image (default %i).\n", DEFAULT_DEPTH);
	printf("  -chains  N     : Number of simultaneous Markov Chains SQUEEZE will run.\n");
	printf("  -threads N     : Number of simultaneous threads SQUEEZE is allowed to use, has to be at least equal to nchains\n");
	printf("  -tempschedc c  : Temperature schedule power c for parallel tempering (default = 3).\n");
	printf("  -nobws         : Do not compute bandwidth smearing factors when computing visibilities.\n");

	printf("\n***** OUTPUT SETTINGS ***** \n");
	printf("  -o filename     : Squeeze outputs as a FITS image file.\n");
	printf("  -fullchain      : Output the full MCMC chain into the file output.fullchain .\n");
	printf("  -temporaryfits  : Enable continuous writing of chainxx.fits to monitor execution.\n");

	printf("\n***** SIMULTANEOUS MODEL FITTING SETTINGS ***** \n");
	printf("  -P p0 p1...    : Initial parameter input.\n");
	printf("  -S s0 s1...    : Initial parameter step sizes (NB: must come after -P option).\n");
	printf("     Current model implemented = polychromatic uniform disc \n");
	printf("     For this model the parameters are in the following order:\n");
	printf("         lambda_0 : the reference wavelength \n");
	printf("         fs_0     : the stellar flux fraction at the reference wavelength \n");
	printf("         diam     : diameter of the uniform disc\n");
	printf("         d_ind    : the flux power law index for the environment \n");

	printf("\n***** OIFITS IMPORT SETTINGS ***** \n");

	printf("  -wavchan         : Define polychromatic channels for reconstruction.\n");
	printf("                  Usage  : -wavchan 0 wavmin_0 wavmax_0 1 wavmin_1 wavmax_1 ...\n");
	printf("                  Then the wavelength channel i will have points with wavmin_i <= lambda < wavmax_i \n");
	printf("                  Example: -wavchan 0 1.2e-6 1.35e-6 1 1.35e-6 1.43e-6 2 1.6e-6 1.8e-6\n\n");

	printf("  -diffvis      : Switch forcing VIS tables to be treated as differential visibilities, not complex visibilities.\n");
	printf("  -novis        : Switch disabling all complex visibility data.\n");
	printf("  -novisamp     : Switch disabling visibility amplitudes.\n");
	printf("  -novisphi     : Switch disabling visibility phases.\n");
	printf("  -not3         : Switch disabling all bispectrum data.\n");
	printf("  -not3amp      : Switch disabling bispectrum amplitudes.\n");
	printf("  -not3phi      : Switch disabling bispectrum phases (closure phases).\n");
	printf("  -nov2         : Switch disabling all powerspectrum data.\n\n");
	printf("  -visamps mult : Scaling factor for visibility amplitude errors.\n");
	printf("  -visampa add  : Additive error to visibility amplitude errors.\n");
	printf("  -visphis mult : Scaling factor for visibility phase/differential phase errors.\n");
	printf("  -visphia add  : Additive error to visibility phase/differential phase errors.\n");
	printf("  -v2s mult     : Scaling factor for visibility-squared errors.\n");
	printf("  -v2a add      : Additive error to visibility-squared errors.\n");
	printf("  -t3amps mult  : Scaling factor for bispectrum amplitude errors.\n");
	printf("  -t3ampa add   : Additive error to bispectrum amplitude errors.\n");
	printf("  -t3phis mult  : Scaling factor for bispectrum amplitude errors.\n");
	printf("  -t3phia add   : Additive error to bispectrum amplitude errors.\n");
	printf("  -fs mult      : Flux scaling factor (zero-baseline visibility intercept).\n");
	printf("  -cv fwhm      : Convolve data and fits by gaussian FWHM mas (default 0.0).\n");
	printf("  -uvtol tol    : Consider all uv points to be the same within tolerance uvtol.\n");

	printf("\n***** REGULARIZATION & INIT SETTINGS ***** \n");
	printf("  -en param     : Entropy regularization multiplier.\n");
	printf("  -de param     : Dark energy regularization multiplier.\n");
	printf("  -ud param     : Uniform disc regularization multiplier.\n");
	printf("  -tv param     : Total variation regularization multiplier.\n");
	printf("  -la param     : Laplacian regularization multiplier.\n");
	printf("  -l0 param     : L0 sparsity norm multiplier.\n");
	printf("  -ts param     : Transpectral L2 regularization for polychromatic reconstructions.\n");
	printf("  -fv param     : Field of view regularizer.\n");
	printf("  -p file.fits  : Prior image. Log of this is the regularization.\n");
	printf("  -ps param     : Prior image regularization multiplier.\n");
	printf("  -i file.fits  : Initial image will be read from a FITS file.\n");
	printf("  -i random     : Initial image will be random, and common to all chain/wavelengths.\n");
	printf("  -i randomthr  : Initial image will be random, with a different image for each chain.\n");

	printf("\n***** CONVERGENCE SETTINGS ***** \n");
	printf("  -tm tmin      : Mininmum temperature (default 1.0 for MCMC).\n");
	printf("  -pa prob      : Adjust temperature to achieve given prob for simulated annealing (default = -1, unused).\n");
	printf("  -ct chi2_t    : Chi-squared / Temperature for sim annealing (default 4.0).\n");
	printf("  -fc target    : Attempt to fix chi-squared to this target (default 0.0).\n");

	printf("  -f_copy       : Fraction of steps that attempt a copycat move.\n");
	printf("  -f_any        : Fraction of step that attempt to move anywhere.\n");

	exit(0);
}

/**********************************************************************/
/* Calculate complex vis chi^2 taking into account the known_phases   */
/**********************************************************************/
void compute_lLikelihood(double* likelihood, const double complex * __restrict mod_vis, double* __restrict res, double* __restrict mod_obs, double *chi2v2, double *chi2t3amp, double *chi2visamp, double *chi2t3phi, double *chi2visphi)
{
	vis_to_obs(mod_vis, mod_obs);
	obs_to_res(mod_obs, res);
	*likelihood = 0.5 * residuals_to_chi2(res, chi2v2, chi2t3amp, chi2visamp, chi2t3phi, chi2visphi);
}

void vis_to_obs(const double complex * __restrict mod_vis, double* __restrict mod_obs)
{
	long i;
	double complex modt3;
	const long visampoffset = nv2 + nt3amp;
	const long t3phioffset = nv2 + nt3amp + nvisamp;
	const long visphioffset = nv2 + nt3amp + nvisamp + nt3phi;

	//#pragma omp for simd
	for(i = 0; i < nv2; i++)
	mod_obs[i] = modsq(mod_vis[v2in[i]]);

	//#pragma omp for simd
	for(i = 0; i < nt3; i++)
	{
		modt3 = mod_vis[ t3in1[i] ] * mod_vis[ t3in2[i] ] * conj(mod_vis[ t3in3[i] ]);

		if(nt3amp > 0) // as many instruments only have closure phases, useful check to gain cycles
		if(data_err[nv2 + i] > 0)
		mod_obs[nv2 + i] = cabs(modt3);

		if(nt3phi > 0)
		if(data_err[t3phioffset + i] > 0)
		//mod_obs[t3phioffset + i] = xatan2(cimag(modt3), creal(modt3));
		mod_obs[t3phioffset + i] = carg(modt3);
	}

	if(nvisamp > 0)
	  //#pragma omp for simd
	for(i = 0; i < nvis; i++)
	if(data_err[visampoffset + i] > 0)
	mod_obs[visampoffset + i] = cabs(mod_vis[ visin[i] ]);

	if(nvisphi > 0)
	  //#pragma omp for simd
	for(i = 0; i < nvis; i++)
	if(data_err[visphioffset + i] > 0)
	mod_obs[visphioffset + i] = carg(mod_vis[ visin[i] ]);

}

void obs_to_res(const double* __restrict mod_obs, double* __restrict res) {
	long i;
	// #pragma omp parallel for simd
	for (i = 0; i < nv2 + nt3amp + nvisamp; i++) {
		res[i] = (mod_obs[i] - data[i]) * data_err[i];
	}

	// #pragma omp parallel for simd
	for (i = nv2 + nt3amp + nvisamp; i < nv2 + nt3amp + nvisamp + nt3phi + nvisphi; i++)
		res[i] = dewrap(mod_obs[i] - data[i]) * data_err[i]; // TBD: improve wrapping
}

double residuals_to_chi2(const double *res, double *chi2v2, double *chi2t3amp, double *chi2visamp, double *chi2t3phi, double *chi2visphi) {
	long i;
	double temp1 = 0, temp2 = 0, temp3 = 0, temp4 = 0, temp5 = 0; // local accumulator (faster than using chi2)

	//#pragma omp simd reduction(+:temp1)
	for (i = 0; i < nv2; i++) {
		temp1 += res[i] * res[i];
	}
	*chi2v2 = temp1;

	//#pragma omp simd reduction(+:temp2)
	for (i = nv2; i < nv2 + nt3amp; i++) {
		temp2 += res[i] * res[i];
	}

	*chi2t3amp = temp2;

	//#pragma omp simd reduction(+:temp3)
	for (i = nv2 + nt3amp; i < nv2 + nt3amp + nvisamp; i++) {
		temp3 += res[i] * res[i];
	}

	*chi2visamp = temp3;

	//#pragma omp simd reduction(+:temp4)
	for (i = nv2 + nt3amp + nvisamp + nt3phi; i < nv2 + nt3amp + nvisamp + nt3phi + nvisphi; i++) {
		temp4 += res[i] * res[i];
	}

	*chi2visphi = temp4;

	//#pragma omp simd reduction(+:temp5)
	for (i = nv2 + nt3amp + nvisamp; i < nv2 + nt3amp + nvisamp + nt3phi; i++) {
		temp5 += res[i] * res[i];
	}

	*chi2t3phi = temp5;

	return temp1 + temp2 + temp3 + temp4 + temp5;
}

static inline double dewrap(double diff) //__attribute__((always_inline))
		{
	if (diff < -M_PI)
		diff += 2. * M_PI;
	if (diff > M_PI)
		diff -= 2. * M_PI;
	return diff;
}

static inline double modsq(double complex input) // __attribute__((always_inline))
{
	return creal(input) * creal(input) + cimag(input) * cimag(input);
}

/**********************************************************/
/* Calculate chi^2 for a flat image (for reference)       */
/**********************************************************/
double get_flat_chi2(bool benchmark) {
	long i, rlong, nbench;
	double dummy1, dummy2, dummy3, dummy4, dummy5, startTime, endTime;
	double complex
	*mod_vis = malloc(nuv * sizeof(double complex));
	double *res = malloc((nv2 + nt3amp + nt3phi + nvisamp + nvisphi) * sizeof(double)); // current residuals
	double *mod_obs = malloc((nv2 + nt3amp + nt3phi + nvisamp + nvisphi) * sizeof(double)); // current observables
	RngStream rngflat = RngStream_CreateStream("flatchi2");
	for (i = 0; i < nuv; i++) {
		rlong = RngStream_RandInt(rngflat, 0, 2147483647);
		mod_vis[i] = cexp(-10.0 + I * (rlong % 1000) * 0.006283);
	}

	if (benchmark == TRUE) {
		printf("Reconst setup -- Benchmarking the chi2...\n");
		nbench = 500000;
	}

	else
		nbench = 1;
	double flat_chi2 = 0;
	startTime = (double) clock() / CLOCKS_PER_SEC;
	for (i = 0; i < nbench; i++)
		compute_lLikelihood(&flat_chi2, mod_vis, res, mod_obs, &dummy1, &dummy2, &dummy3, &dummy4, &dummy5);

	endTime = (double) clock() / CLOCKS_PER_SEC;
	if (benchmark == TRUE)
		printf("Reconst setup -- %ld iterations in %f seconds = %f it/sec for this datafile\n", nbench, endTime - startTime,
				(double) nbench / (endTime - startTime));

	free(res);
	free(mod_obs);
	free (mod_vis);
	RngStream_DeleteStream(&rngflat);
	if (benchmark == TRUE)
		getchar();

	return 2. * flat_chi2;
}

#include "extract_oifits.c"
#include "./models/modelcode.c"
#include "regularizations.c"

/***********************************/
/* Write fits image cube           */
/***********************************/
int writeasfits(char *file, double *image, int nwavr, long depth, long min_elt, double chi2, double chi2v2, double chi2t3amp, double chi2t3phi, double chi2visamp, double chi2visphi,
		double temperature, long nelems, double* regpar, double* regval, long niter,
		unsigned short axis_len, double ndf, double tmin, double chi2_temp, double chi2_target, double mas_pixel, int nchains, double logZ, double logZ_err,
		char *init_filename, char *prior_filename, double* params, double* params_std)
 {
	int status = 0;
	int i, j, w;
	fitsfile *fptr; /* pointer to the FITS file, defined in fitsio.h */
	long fpixel;
	int bitpix = DOUBLE_IMG; /* 32-bit unsigned long pixel values       */
	char param_string[MAX_STRINGS];
	char filename[MAX_STRINGS];
	sprintf(filename, "!%s.fits", file);
	long naxis, naxes[4];
	double moment1 = 0.0, moment2 = 0.0;

	/* Prepare variables for fits file writing */
	naxes[0] = axis_len;
	naxes[1] = axis_len;
	naxes[2] = nwavr;
	naxis = 3;

	/* Save the image to a fits file...*/

	/* Delete old file if it already exists */

	//if(remove(filename) != 0)
	//printf("Error deleting temporary fits file\n");
	if (fits_create_file(&fptr, filename, &status)) /* create new FITS file */
		printerror(status);

	/* write the required keywords for the primary array image.     */
	/* Since bitpix = USHORT_IMG, this will cause cfitsio to create */
	/* a FITS image with BITPIX = 16 (signed short integers) with   */
	/* BSCALE = 1.0 and BZERO = 32768.  This is the convention that */
	/* FITS uses to store unsigned integers.  Note that the BSCALE  */
	/* and BZERO keywords will be automatically written by cfitsio  */
	/* in this case.                                                */

	if (fits_create_img(fptr, bitpix, naxis, naxes, &status))
		printerror(status);

	fpixel = 1; /* first pixel to write      */

	/* write the array of unsigned integers to the FITS file */
	if (fits_write_img(fptr, TDOUBLE, fpixel, naxes[0] * naxes[1] * naxes[2], image, &status))
		printerror(status);

	/* write another optional keyword to the header */
	/* Note that the ADDRESS of the value is passed in the routine */
	if (fits_write_key_dbl(fptr, "SCALE", mas_pixel, 3, "Milli-arcsecs per pixel", &status))
		printerror(status);

	/* write WCS keywords */
	if (fits_write_key_dbl(fptr, "CDELT1", -mas_pixel, 3, "Milli-arcsecs per pixel", &status))
		printerror(status);
	if (fits_write_key_dbl(fptr, "CDELT2", mas_pixel, 3, "Milli-arcsecs per pixel", &status))
		printerror(status);
	if (fits_write_key_dbl(fptr, "CRVAL1", 0.0, 3, "X-coordinate of ref pixel", &status))
		printerror(status);
	if (fits_write_key_dbl(fptr, "CRVAL2", 0.0, 3, "Y-coordinate of ref pixel", &status))
		printerror(status);
	if (fits_write_key_lng(fptr, "CRPIX1", naxes[0] / 2, "Ref pixel in X", &status))
		printerror(status);
	if (fits_write_key_lng(fptr, "CRPIX2", naxes[1] / 2, "Ref pixel in Y", &status))
		printerror(status);
	if (fits_write_key_str(fptr, "CTYPE1", "RA", "Name of X-coordinate", &status))
		printerror(status);
	if (fits_write_key_str(fptr, "CTYPE2", "DEC", "Name of Y-coordinate", &status))
		printerror(status);

	if (regpar != NULL) {

		for (i = 0; i < NREGULS; i++) {
			sprintf(param_string, "HYPER%1d", i);
			fits_write_key_dbl(fptr, param_string, regpar[i], 3, "Hyperparameter value", &status);
			printerror(status);
		}

	}

	if (regval != NULL) {
		for (w = 0; w < nwavr; w++) {
			for (i = 0; i < NREGULS; i++) {
				//printf("w: %d i: %d reg = %lf\n",w, i, regval[w * NREGULS + i]);
				sprintf(param_string, "REGUL%1dW%1d", i, w);
				fits_write_key_dbl(fptr, param_string, regval[w * NREGULS + i], 3, "Regularizer value", &status);
				printerror(status);
			}
		}
	}

	/* if ( fits_write_key_str ( fptr, "CUNIT1", "mas", "Unit of X-coordinate",
	 *    &status ) ) */
	/*   printerror ( status ); */
	/* if ( fits_write_key_str ( fptr, "CUNIT2", "mas", "Unit of Y-coordinate",
	 *    &status ) ) */
	/*   printerror ( status ); */
	//	getchar();
	/* Write model parameters */
	if ((min_elt >= 0) && (nparams > 0)) {

		if (params != NULL) {

			for (i = 0; i < nparams; i++) {
				sprintf(param_string, "MNPARAM%1d", i + 1);
				if (fits_write_key_dbl(fptr, param_string, params[i], 3, "Mean parameter for model", &status))
					printerror(status);
			}

			if (params_std != NULL) {
				for (i = 0; i < nparams; i++) {
					sprintf(param_string, "SDPARAM%1d", i + 1);
					if (fits_write_key_dbl(fptr, param_string, params_std[i], 3, "Stdev of parameter for model", &status))
						printerror(status);
				}
			}
		}
	}

	if (temperature >= 0) {
		if (fits_update_key(fptr, TDOUBLE, "TEMPER", &temperature, "Algorithm Temperature", &status))
			printerror(status);
	}

	if (fits_update_key(fptr, TDOUBLE, "CHI2", &chi2, "Chi-squared per degree of freedom", &status))
		printerror(status);
	if (fits_update_key(fptr, TDOUBLE, "CHI2V2", &chi2v2, "Chi-squared per degree of freedom", &status))
		printerror(status);
	if (fits_update_key(fptr, TDOUBLE, "CHI2T3A", &chi2t3amp, "Chi-squared per degree of freedom", &status))
		printerror(status);
	if (fits_update_key(fptr, TDOUBLE, "CHI2T3P", &chi2t3phi, "Chi-squared per degree of freedom", &status))
		printerror(status);
	if (fits_update_key(fptr, TDOUBLE, "CHI2VA", &chi2visamp, "Chi-squared per degree of freedom", &status))
		printerror(status);
	if (fits_update_key(fptr, TDOUBLE, "CHI2VP", &chi2visphi, "Chi-squared per degree of freedom", &status))
		printerror(status);
	if (fits_update_key(fptr, TDOUBLE, "TMIN", &tmin, "Setting for minimum temperature", &status))
		printerror(status);
	if (fits_update_key(fptr, TDOUBLE, "CHI2T", &chi2_temp, "Setting for chi-squared divided by temperature (convergence rate)", &status))
		printerror(status);
	if (fits_update_key(fptr, TDOUBLE, "CHI2TG", &chi2_target, "Setting for target chi-squared", &status))
		printerror(status);
	if (fits_update_key(fptr, TDOUBLE, "LOGZ", &logZ, "Marginal likelihood", &status))
		printerror(status);
	if (fits_update_key(fptr, TDOUBLE, "LOGZE", &logZ_err, "Error on Marginal likelihood", &status))
		printerror(status);
	if (fits_update_key(fptr, TSTRING, "INITIM", init_filename, "Initialization image.", &status))
		printerror(status);
	if (fits_update_key(fptr, TSTRING, "PRIORIM", prior_filename, "Prior image.", &status))
		printerror(status);
	if (fits_update_key(fptr, TSTRING, "OIFITS", oifits_file, "Input OIFITS file.", &status))
		printerror(status);
	if (fits_update_key(fptr, TLONG, "NITER", &niter, "Number of iterations per chain per element.", &status))
		printerror(status);
	if (fits_update_key(fptr, TINT, "NWAVR", &nwavr, "Number of spectral channels.", &status))
		printerror(status);
	if (fits_update_key(fptr, TLONG,    "DEPTH", &depth, "Number of realisations in image", &status))
		printerror(status);
	if (fits_update_key(fptr, TLONG,    "ELEMENTS", &nelems, "Number of elements per realisation", &status))
		printerror(status);
	if (fits_update_key(fptr, TINT,    "NCHAINS", &nchains, "Number of chains", &status))
		printerror(status);
	if (fits_update_key(fptr, TDOUBLE, "NDF", &ndf, "Number of degrees of freedom", &status))
		printerror(status);
	if (fits_update_key(fptr, TLONG,   "NV2", &nv2, "Number of V2", &status))
		printerror(status);
	if (fits_update_key(fptr, TLONG,   "NT3AMP", &nt3amp, "Number of T3AMP", &status))
		printerror(status);
	if (fits_update_key(fptr, TLONG,   "NT3PHI", &nt3phi, "Number of T3PHI", &status))
		printerror(status);
	if (fits_update_key(fptr, TLONG,   "NVISAMP", &nvisamp, "Number of VISAMP", &status))
		printerror(status);
	if (fits_update_key(fptr, TLONG,   "NVISPHI", &nvisphi, "Number of VISPHI", &status))
		printerror(status);

	if (fits_close_file(fptr, &status)) /* close the file */
		printerror(status);

	return status;
}

double mean(long *x, long n) {
	// returns average of vector x[0...n-1]
	double avg = 0;
	register long i;
	for (i = 0; i < n; i++)
		avg += x[i];
	if (avg != 0)
		avg /= (double) n;
	else
		avg = 0;
	return avg;
}

double stddev(long *x, long n) {
	// returns average of vector x[0...n-1]
	register long i;
	double avg = mean(x, n);
	double temp = 0;
	for (i = 0; i < n; i++)
		temp += (x[i] - avg) * (x[i] - avg);
	return sqrt(temp / (n - 1));
}

void compute_logZ(const double* temperature, const unsigned short* iStoragetoChain, const double* lLikelihood_expectation, const double *lLikelihood_deviation,
		int nchains, double* logZ, double* logZ_err) {
	long j;
	*logZ = 0;
	*logZ_err = (1. / temperature[iStoragetoChain[1]] - 1. / temperature[iStoragetoChain[0]])
			* (1. / temperature[iStoragetoChain[1]] - 1. / temperature[iStoragetoChain[0]]) * lLikelihood_deviation[iStoragetoChain[0]]
			+ (1. / temperature[iStoragetoChain[nchains - 1]] - 1. / temperature[iStoragetoChain[nchains - 2]])
					* (1. / temperature[iStoragetoChain[nchains - 1]] - 1. / temperature[iStoragetoChain[nchains - 2]])
					* lLikelihood_deviation[iStoragetoChain[nchains - 1]];

	//  for(j = 0; j < nchains; j++)
	// printf("temperature %lf lLike_avg[%ld] = %lf \n", temperature[iStoragetoChain[j]], j, lLikelihood_expectation[iStoragetoChain[j]]);

	for (j = 0; j < (nchains - 1); j++) {
		//      printf("logZ debug %ld %lf %lf %lf\n", j, (1. / temperature[iStoragetoChain[j]] - 1. / temperature[iStoragetoChain[j + 1]]),
		//                   (lLikelihood_expectation[iStoragetoChain[j]] + lLikelihood_expectation[iStoragetoChain[j + 1]]),
		//		   0.5 * (1. / temperature[iStoragetoChain[j]] - 1. / temperature[iStoragetoChain[j + 1]])
		//                   * (lLikelihood_expectation[iStoragetoChain[j]] + lLikelihood_expectation[iStoragetoChain[j + 1]]));

		*logZ += (1. / temperature[iStoragetoChain[j + 1]] - 1. / temperature[iStoragetoChain[j]])
				* (lLikelihood_expectation[iStoragetoChain[j]] + lLikelihood_expectation[iStoragetoChain[j + 1]]);
	}

	*logZ *= 0.5;

	for (j = 1; j < (nchains - 2); j++) {
		*logZ_err += (1. / temperature[iStoragetoChain[j + 1]] - 1. / temperature[iStoragetoChain[j - 1]])
				* (1. / temperature[iStoragetoChain[j + 1]] - 1. / temperature[iStoragetoChain[j - 1]]) * lLikelihood_deviation[iStoragetoChain[j]];
	}

	*logZ_err = sqrt(*logZ_err * 0.5);

}

//
void mcmc_fullchain(char *file, long nchains, long niter, int nwavr, long nelements, unsigned short axis_len, unsigned short* saved_x, unsigned short* saved_y,
		double* saved_params, double* saved_lLikelihood, double* saved_lPrior, double* saved_lPosterior, double* temperature, unsigned short* iChaintoStorage) {
	int i, n, w, t;

	// Dump all probabilities + temperatures
	char prob_filename[100];
	sprintf(prob_filename, "%s.fullprobs", file);
	FILE * pFile2 = fopen(prob_filename, "w");
	fprintf(pFile2, "%ld , %ld, 0. \n", nchains, niter);

	for (t = 0; t < (nchains / 3); t++) {
		fprintf(pFile2, "%lf , %lf, %lf\n", temperature[3 * t], temperature[3 * t + 1], temperature[3 * t + 2]);
	}

	for (t = 0; t < (nchains % 3); t++) {
		if (t == 0) {
			fprintf(pFile2, "%lf", temperature[3 * (nchains / 3) + t]);
		} else {
			fprintf(pFile2, " , %lf", temperature[3 * (nchains / 3) + t]);
		}
	}

	if ((nchains % 3) != 0) {
		for (t = 0; t < 3 - (nchains % 3); t++) {
			fprintf(pFile2, ", 0");
		}
		fprintf(pFile2, "\n");
	}

	for (t = 0; t < nchains; t++) {
		for (n = 0; n < niter; n++) {
			fprintf(pFile2, "%lf , %lf , %lf\n", saved_lLikelihood[iChaintoStorage[t] * niter + n], saved_lPrior[iChaintoStorage[t] * niter + n],
					saved_lPosterior[iChaintoStorage[t] * niter + n]);
		}
	}

	fclose(pFile2);
	printf("Output -- Probabilities output to output.fullprobs.\n");

	// Dump the content of saved_x, saved_y, saved_params into a file

	double *images = calloc(axis_len * axis_len * nwavr * niter * nchains, sizeof(double));
	if (images == NULL) {
		printf("Output -- Out of memory to write fullchain as a fits file.\n");
	} else {
		int status = 0;
		fitsfile *fptr;
		long fpixel = 1;
		int bitpix = DOUBLE_IMG;
		long naxis, naxes[5];
		char fullchain_filename[100];
		sprintf(fullchain_filename, "!%s_fullchain.fits.gz", file);

		for (t = 0; t < nchains; t++) {
			for (n = 0; n < niter; n++) {
				for (w = 0; w < nwavr; w++) {
					for (i = 0; i < nelements; i++) {
						images[t * niter * axis_len * axis_len * nwavr + n * axis_len * axis_len * nwavr + w * axis_len * axis_len
								+ saved_x[t * nwavr * nelements * niter + n * nwavr * nelements + w * nelements + i]
								+ axis_len * saved_y[t * nwavr * nelements * niter + n * nwavr * nelements + w * nelements + i]] += 1;
					}
				}
			}
		}

		naxes[0] = axis_len;
		naxes[1] = axis_len;
		naxes[2] = nwavr;
		naxes[3] = niter;
		naxes[4] = nchains;
		naxis = 5;

		if (fits_create_file(&fptr, fullchain_filename, &status))
			printerror(status);

		if (fits_create_img(fptr, bitpix, naxis, naxes, &status))
			printerror(status);
		if (fits_write_img(fptr, TDOUBLE, fpixel, naxes[0] * naxes[1] * naxes[2] * naxes[3] * naxes[4], images, &status))
			printerror(status);
		if (fits_update_key(fptr, TLONG, "NITER", &niter, "Number of iterations per chain per element.", &status))
			printerror(status);
		if (fits_update_key(fptr, TINT, "NCHAINS", &nchains, "Number of chains", &status))
			printerror(status);
		if (fits_update_key(fptr, TINT, "nwavR", &nwavr, "Number of spectral channels.", &status))
			printerror(status);
		if (fits_update_key(fptr, TLONG, "ELEMENTS", &nelements, "Number of elements per realization", &status))
			printerror(status);

		if (fits_close_file(fptr, &status))
			printerror(status);

		free(images);
		printf("Output -- Full MCMC chain output to %s.\n", fullchain_filename);

	}

}

/***************************************************************/
/* Fill an image from the saved_x and y, and a iframeburned vector */
/****************************************************************/
void mcmc_annealing_results(const char *file, double *final_image, const long *iframeburned, const long depth, const long nelements, const unsigned short axis_len,
			    const double complex * __restrict xtransform, const double complex * __restrict ytransform, double *final_chi2,
			    double *chi2v2, double *chi2t3amp, double *chi2t3phi, double *chi2visamp, double *chi2visphi,
			    const unsigned short *saved_x, const unsigned short *saved_y, const double *saved_params, const long niter,
			    const int nwavr, double* final_params, double* final_params_std,
			    const double* reg_param, double* final_reg_value, const double* prior_image, const unsigned short* initial_x, const unsigned short* initial_y,
			    double* centroid_image_x, double* centroid_image_y, const double fov, const double cent_mult)
// This averages the image obtained by MCMC over the iterations and chains
// the depth input should be the actual available depth, not the requested one
{
	int i, j, k, w;

	//
	// Compute annealed image
	//

	// Initialize image and parameters
	for(i = 0; i < nwavr * axis_len * axis_len; i++)
	final_image[i] = 0;

	// (iframeburned[k] / niter) = chain number
	// (iframeburned[k] % niter) = iteration number
	for(k = 0; k < depth; k++)
	{
		for(w = 0; w < nwavr; w++)
		{
			for(i = 0; i < nelements; i++)
			{

				final_image[w * axis_len * axis_len
				+ saved_y[(iframeburned[k] / niter) * nwavr * nelements * niter + (iframeburned[k] % niter) * nwavr * nelements + w * nelements + i] * axis_len
				+ saved_x[(iframeburned[k] / niter) * nwavr * nelements * niter + (iframeburned[k] % niter) * nwavr * nelements + w * nelements + i]] += 1.0 / (nelements * depth);

				//printf("%d ", saved_x[ (iframeburned[k] / niter) * nwavr * nelements * niter + (iframeburned[k] % niter) * nwavr * nelements + w * nelements + i]);
			}
		}
	}

	//
	// Compute annealed parametric model and errors
	//

	if(nparams > 0)
	{

		for(j = 0; j < nparams; j++)
		final_params[j] = 0;

		for(k = 0; k < depth; k++)
		{
			for(j = 0; j < nparams; j++)
			final_params[j] += saved_params[iframeburned[k] * nparams + j];
		}

		for(j = 0; j < nparams; j++)
		{
			final_params[j] /= (double) depth;
			final_params_std[j] = -1.;
		}
	}

	double complex *mod_vis = malloc(nuv * sizeof(double complex));
	double complex *im_vis = malloc(nuv * sizeof(double complex));
	double complex *param_vis = malloc(nuv * sizeof(double complex));
	double *final_fluxratio_image = malloc(nuv * sizeof(double));
	double lPriorModel = 0;
	compute_model_visibilities_fromimage(mod_vis, im_vis, param_vis, final_params,
		                             final_fluxratio_image, final_image, xtransform, ytransform, &lPriorModel, nparams, nelements, axis_len);

	double *res = malloc((nv2 + nt3amp + nt3phi + nvisamp + nvisphi) * sizeof(double)); // current residuals
	double *mod_obs = malloc((nv2 + nt3amp + nt3phi + nvisamp + nvisphi) * sizeof(double));// current observables
	for(i = 0; i < (nv2 + nt3amp + nt3phi + nvisamp + nvisphi); i++)
	{
		res[i] = 0;
		mod_obs[i] = 0;
	}
	// compute chi2s on final image
	compute_lLikelihood(final_chi2, mod_vis, res, mod_obs, chi2v2, chi2t3amp, chi2visamp, chi2t3phi, chi2visphi);
	*final_chi2 *= 2.;
    if(nv2>0)
        *chi2v2 /= (double) nv2;
    if(nt3amp>0)
        *chi2t3amp /= (double) nt3amp;
    if(nt3phi>0)
        *chi2t3phi /= (double) nt3phi;
    if(nvisamp>0)
        *chi2visamp /= (double) nvisamp;
    if(nvisphi>0)
        *chi2visphi /= (double) nvisphi;

	// recompute regularizers on final image
	compute_regularizers(reg_param, final_reg_value, final_image, prior_image, 1., initial_x, initial_y, nwavr, axis_len, nelements, centroid_image_x, centroid_image_y, fov, cent_mult);

	// Output obs and residuals in separate file
	char data_filename[MAX_STRINGS];
	sprintf(data_filename, "%s.data", file);
	FILE * pFile = fopen(data_filename, "w");
	fprintf(pFile, "%lf %lf %lf %lf\n", (double)nuv, (double)nv2 , (double)nt3amp , (double)nvisamp);
	fprintf(pFile, "%lf %lf %lf %lf\n", (double)nt3phi , (double)nt3, (double)nvisphi, (double)nwavr);
	for(i = 0; i < nuv; i++)
	fprintf(pFile, "%lf %lf %lf %lf\n", u[i], v[i], uv_lambda[i] * 1E6, uv_dlambda[i] * 1E6);
	for(i = 0; i < nt3; i++)
	fprintf(pFile, "%lf %lf %lf %lf\n", (double)(t3in1[i]), (double)(t3in2[i]), (double)(t3in3[i]), 0.0);
	for(i = 0; i < nv2 + nt3amp + nvisamp + nt3phi + nvisphi; i++)
	fprintf(pFile, "%lf %lf %lf %lf\n", mod_obs[i], data[i], data_err[i], res[i]);
	fclose(pFile);

	free(res);
	free(mod_obs);

	free(final_fluxratio_image);
	free(mod_vis);
	free(im_vis);
	free(param_vis);

}

void mcmc_tempering_results(char* file, double *image, long lowtempchain, long depth, long nelements, unsigned short axis_len,
		double complex * __restrict xtransform, double complex * __restrict ytransform, double *final_chi2,
		unsigned short *saved_x, unsigned short *saved_y, double *saved_params, long niter, int nwavr)
// note: depth is the actual available depth
{
	// TODO: rework this to use standard image2vis, etc. functions

	int i, j, k, w, ix, iy;
	double lPriorModel = 0;
	double avg_params[MAX_PARAMS];
	double complex *mod_vis = malloc(nuv * sizeof(double complex));
	double complex *im_vis = malloc(nuv * sizeof(double complex));
	double complex *param_vis = malloc(nuv * sizeof(double complex));
	double *avg_fluxratio_image = malloc(nuv * sizeof(double));

	for(j = 0; j < MAX_PARAMS; j++)
	avg_params[j] = 0;

	for(j = 0; j < nuv; j++)
	param_vis[j] = 0.0;

	// Compute mean image
	// Initialize image and parameters
	for(i = 0; i < nwavr * axis_len * axis_len; i++)
	image[i] = 0.;

	//    printf("MCMC tempering: %ld %ld %ld %d\n", niter, nwavr,  nelements, axis_len);

	depth = ceil(0.7 * niter);// 30% frames for burn-in
	// lowtempchain is
	for(k = niter - depth; k < niter; k++)
	{
		for(w = 0; w < nwavr; w++)
		{
			for(i = 0; i < nelements; i++)
			{
				image[w * axis_len * axis_len
				+ saved_y[lowtempchain * nwavr * nelements * niter + k * nwavr * nelements + w * nelements + i] * axis_len
				+ saved_x[lowtempchain * nwavr * nelements * niter + k * nwavr * nelements + w * nelements + i] ] += 1.0 / (nelements * depth);
			}
		}
	}

	// The following three blocks are the default image + param to vis routine
	// TODO: move all these to a function as they are common to all mcmc_image routines

	// Compute image visibilities
	for(j = 0; j < nuv; j++)
	{
		im_vis[j] = 0;
		for(ix = 0; ix < axis_len; ix++)
		  for(iy = 0; iy < axis_len; iy++)
		    im_vis[j] += image[uvwav2chan[j] * axis_len * axis_len + iy * axis_len + ix] * xtransform[ix * nuv + j] * ytransform[iy * nuv + j];
	}

	// Compute model visibilities
	if(nparams > 0)
	{
		for(k = 0; k < depth; k++)
		{
			for(j = 0; j < nparams; j++)
			avg_params[j] += saved_params[0 * nparams + j];
		}

		for(j = 0; j < MAX_PARAMS; j++)
		avg_params[j] /= (double) depth;
		// Note: we don't use the history of the flux ratio between image and model
		// as it is recomputed from the averaged parameterile.field1[0]s
		model_vis(&avg_params[0], param_vis, &lPriorModel, avg_fluxratio_image);
	}

	// Total visibilities
	for(j = 0; j < nuv; j++)
	{
		if(nparams > 0)
		mod_vis[j] = param_vis[j] + im_vis[j] * avg_fluxratio_image[j];
		else
		mod_vis[j] = im_vis[j];
	}

	double dummy1 = 0, dummy2 = 0, dummy3 = 0, dummy4 = 0, dummy5 = 0;
	double *res = malloc((nv2 + nt3amp + nt3phi + nvisamp + nvisphi) * sizeof(double)); // current residuals
	double *mod_obs = malloc((nv2 + nt3amp + nt3phi + nvisamp + nvisphi) * sizeof(double));// current observables
	for(i = 0; i < (nv2 + nt3amp + nt3phi + nvisamp + nvisphi); i++)
	{
		res[i] = 0;
		mod_obs[i] = 0;
	}
	compute_lLikelihood(final_chi2, mod_vis, res, mod_obs, &dummy1, &dummy2, &dummy3, &dummy4, &dummy5);

	// Output obs and residuals in separate file
	char data_filename[100];
	sprintf(data_filename, "%s.data", file);
	FILE * pFile = fopen(data_filename, "w");
	fprintf(pFile, "%lf %lf %lf %lf\n", (double)nuv, (double)nv2 , (double)nt3amp , (double)nvisamp);
	fprintf(pFile, "%lf %lf %lf %lf\n", (double)nt3phi , (double)nt3, (double)nvisphi, (double)nwavr);
	for(i = 0; i < nuv; i++)
	fprintf(pFile, "%lf %lf %lf %lf\n", u[i], v[i], uv_lambda[i] * 1E6, uv_dlambda[i] * 1E6);
	for(i = 0; i < nt3; i++)
	fprintf(pFile, "%lf %lf %lf %lf\n", (double)(t3in1[i]), (double)(t3in2[i]), (double)(t3in3[i]), 0.0);
	for(i = 0; i < nv2 + nt3amp + nvisamp + nt3phi + nvisphi; i++)
	fprintf(pFile, "%lf %lf %lf %lf\n", mod_obs[i], data[i], data_err[i], res[i]);
	fclose(pFile);

	free(res);
	free(mod_obs);

	free(avg_fluxratio_image);
	free(mod_vis);
	free(im_vis);
	free(param_vis);

}

inline void swapi(unsigned short* a, unsigned short* b) {
	unsigned short tmp = *a;
	*a = *b;
	*b = tmp;
}

inline void swapd(double* a, double* b) {
	double tmp = *a;
	*a = *b;
	*b = tmp;
}

double sinc(double x) {
	return sin(x + 1e-15) / (x + 1e-15);
}

void intHandler(int signum) {
	printf("\nExiting at next opportunity\n");
	if (ctrlcpressed == TRUE)
		exit(0);
	else
		ctrlcpressed = TRUE;

}

void compute_regularizers(const double *reg_param, double *reg_value, const double* image, const double* prior_image, const double fluxscaling,
		const unsigned short* initial_x, const unsigned short* initial_y, const int nwavr, const unsigned short axis_len, const long nelements,
		double* centroid_image_x, double* centroid_image_y, const double fov, const double cent_mult) {

	long w, i;

	//	for(i=0;i<NREGULS;i++)
	//  printf("CRi: %lf %lf\n", reg_param[i], reg_value[i]);

	for (w = 0; w < nwavr; w++) {
		if (reg_param[REG_DARKENERGY] > 0)
			reg_value[w * NREGULS + REG_DARKENERGY] = den_full(&image[w * axis_len * axis_len], NULL, 0.0, axis_len, axis_len, fluxscaling);

		if (reg_param[REG_ENTROPY] > 0.0)
			reg_value[w * NREGULS + REG_ENTROPY] = entropy_full(&image[w * axis_len * axis_len], NULL, 0.0, axis_len, axis_len, fluxscaling);

		if (reg_param[REG_SPOT] > 0.0)
			reg_value[w * NREGULS + REG_SPOT] = UDreg(&image[w * axis_len * axis_len], NULL, 0.0, axis_len, axis_len, fluxscaling);

		if (reg_param[REG_TV] > 0.0)
			reg_value[w * NREGULS + REG_TV] = TV(&image[w * axis_len * axis_len], NULL, 0.0, axis_len, axis_len, fluxscaling);

		if (reg_param[REG_LAP] > 0.0)
			reg_value[w * NREGULS + REG_LAP] = LAP(&image[w * axis_len * axis_len], NULL, 0.0, axis_len, axis_len, fluxscaling);

		if (reg_param[REG_L0] > 0.0)
			reg_value[w * NREGULS + REG_L0] = L0(&image[w * axis_len * axis_len], NULL, 0.0, axis_len, axis_len, fluxscaling);

		if (reg_param[REG_PRIORIMAGE] > 0)
			reg_value[w * NREGULS + REG_PRIORIMAGE] = reg_prior_image(&image[w * axis_len * axis_len], &prior_image[w * axis_len * axis_len], 0.0, axis_len,
					axis_len, fluxscaling);

		if (reg_param[REG_CENTERING] > 0) {
			reg_value[w * NREGULS + REG_CENTERING] = 0;
			for (i = 0; i < nelements; i++)
				reg_value[w * NREGULS + REG_CENTERING] += cent_change(w, centroid_image_x, centroid_image_y, initial_x[i], initial_y[i], axis_len / 2,
						axis_len / 2, axis_len, fov, cent_mult);
		}
		if (reg_param[REG_MODELPARAM] > 0)
			reg_value[w * NREGULS + REG_MODELPARAM] = 0.;

	}

	if (reg_param[REG_TRANSPECL2] > 0.0)
		reg_value[REG_TRANSPECL2] = transpec(nwavr, axis_len, image, fluxscaling);

}

void compute_model_visibilities_fromelements(double complex* mod_vis, double complex* im_vis, double complex* param_vis,
		double* params, double* fluxratio_image, const unsigned short* element_x, const unsigned short* element_y,
		const double complex* xtransform, const double complex* ytransform, double* lPriorModel, long nparams, long nelements)
{
	long i, j;
	if(nparams > 0)
	model_vis(params, param_vis, lPriorModel, fluxratio_image);
	else
	for(j = 0; j < nuv; j++)
	param_vis[j] = 0.0;

	// Compute visibilities from scratch
	for(j = 0; j < nuv; j++)
	{
		im_vis[j] = 0;
		for(i = 0; i < nelements; i++)
		  im_vis[j] += xtransform[element_x[ uvwav2chan[j] * nelements + i] * nuv + j] * ytransform[element_y[ uvwav2chan[j] * nelements + i] * nuv + j];
		im_vis[j] *= fluxratio_image[j] / (double) nelements; // Will add SED here
		mod_vis[j] = param_vis[j] + im_vis[j];
	}
}

void compute_model_visibilities_fromimage(double complex* mod_vis, double complex* im_vis, double complex* param_vis,
		const double* params, double* fluxratio_image, const double *image, const double complex* xtransform, const double complex* ytransform,
		double* lPriorModel, long nparams, long nelements, unsigned short axis_len)
{
	// Note: input image should be normalized

	long ix, iy, j;
	if(nparams > 0)
	model_vis(params, param_vis, lPriorModel, fluxratio_image);
	else
	for(j = 0; j < nuv; j++)
	param_vis[j] = 0.0;

	// Compute image visibilities
	for(j = 0; j < nuv; j++)
	{
		im_vis[j] = 0;
		for(ix = 0; ix < axis_len; ix++)
		for(iy = 0; iy < axis_len; iy++)
		im_vis[j] += image[uvwav2chan[j] * axis_len * axis_len + iy * axis_len + ix] * xtransform[ix * nuv + j] * ytransform[iy * nuv + j];
	}

	for(j = 0; j < nuv; j++)
	{
		if(nparams > 0)
		{

			im_vis[j] *= fluxratio_image[j];
			mod_vis[j] = param_vis[j] + im_vis[j];
		}
		else mod_vis[j] = im_vis[j];
	}
}

void initialize_image(int iChain, double* image, unsigned short* element_x, unsigned short* element_y, unsigned short* initial_x, unsigned short* initial_y,
		unsigned short axis_len, int nwavr, long nelements, char* init_filename) {
	long i, k, w;
	for (i = 0; i < nwavr * axis_len * axis_len; i++)
		image[i] = 0;

	// Initial image channels
	if (!strcmp(init_filename, "randomthr"))
		k = iChain;
	else
		k = 0;
	for (w = 0; w < nwavr; w++) {
		for (i = 0; i < nelements; i++) {
			element_x[w * nelements + i] = initial_x[k * nelements + i];
			element_y[w * nelements + i] = initial_y[k * nelements + i];
			image[w * axis_len * axis_len + initial_y[i] * axis_len + initial_x[i]]++;
		}
	}
}

bool read_commandline(int* argc, char** argv, bool* benchmark, bool* use_v2, bool* use_t3amp, bool* use_t3phi, bool* use_visamp, bool* use_visphi,
		bool* use_diffvis, bool* use_tempfitswriting, bool* use_bandwidthsmearing, int* minimization_engine, bool* dumpchain, double* mas_pixel,
		unsigned short* axis_len, long* depth, long* niter, long* nelements, double* f_anywhere, double* f_copycat, int *nchains, int* nthreads,
		double* tempschedc, double* fov, double* chi2_temp, double* chi2_target, double* tmin, double* prob_auto, double* uvtol, char* output_filename,
		char* init_filename, char* prior_filename, double* v2s, double* v2a, double* t3amps, double* t3ampa, double* t3phia, double* t3phis, double* visamps,
		double* visampa, double* visphis, double* visphia, double* fluxs, double* cvfwhm, double* reg_param, double* init_param, double* wavmin,
	        double* wavmax, int* nwavr) 
{
	long i, j, k;

	
	/* Read in command line info... */
	if (*argc < 2) // need at least a filename to do something
	  {
	    printhelp();
	    return 0;
	  }

	// Determine the number of oifits files to use 
	// For this we determine the position of the first option
	int nfiles = 1;
	for(i=1;i<*argc;i++)
	  {
	    if(argv[i][0] == '-')
	    {
	      nfiles = i-1;
	      break;
	    }
	    if(i==(*argc-1)) 
	      nfiles = i;
	  }

	printf("Command line -- %d oifits to read\n", nfiles); 

	if ((strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "-help") == 0 || strcmp(argv[1], "--help") == 0) || (strcmp(argv[1], "--help") == 0)) {
		printhelp();
		return 0;
	}

	for (i = 2; i < *argc; i++) {
		/* First the options without arguments */
		if ((strcmp(argv[i], "-h") == 0) || (strcmp(argv[i], "-help") == 0) || (strcmp(argv[i], "--help") == 0)) {
			printhelp();
			return 0;
		} else if (strcmp(argv[i], "-benchmark") == 0) {
			*benchmark = TRUE; // run benchmark
		} else if (strcmp(argv[i], "-nov2") == 0) {
			*use_v2 = FALSE; // disable the v2
		} else if (strcmp(argv[i], "-not3amp") == 0) {
			*use_t3amp = FALSE; // disable the t3amp
		} else if (strcmp(argv[i], "-not3phi") == 0) {
			*use_t3phi = FALSE; // disable the t3phi
		} else if (strcmp(argv[i], "-not3") == 0) {
			*use_t3amp = FALSE; // disable the t3amp
			*use_t3phi = FALSE; // disable the t3phi
		} else if (strcmp(argv[i], "-novisamp") == 0) {
			*use_visamp = FALSE; // disable the visamp
		} else if (strcmp(argv[i], "-novisphi") == 0) {
			*use_visphi = FALSE; // disable the visphi
		} else if (strcmp(argv[i], "-novis") == 0) {
			*use_visamp = FALSE; // disable the t3amp
			*use_visphi = FALSE; // disable the t3phi
		} else if (strcmp(argv[i], "-diffvis") == 0) {
			*use_diffvis = TRUE; // interpret VIS tables as differential visibilities, not complex visibilities
		} else if (strcmp(argv[i], "-temporaryfits") == 0) {
			*use_tempfitswriting = TRUE; // enable writing chainxx.fits to disc
		} else if (strcmp(argv[i], "-nobws") == 0) {
			*use_bandwidthsmearing = FALSE; // disable bandwidth smearing
		} else if (strcmp(argv[i], "-tempering") == 0) {
			*minimization_engine = ENGINE_PARALLEL_TEMPERING;
			printf("Command line -- Using parallel tempering\n");
		} else if (strcmp(argv[i], "-fullchain") == 0) {
			*dumpchain = TRUE; // write the full MCMC chain in output.fullchain
		} else if (*argc > i + 1) /* Now the options with arguments*/
		{
			if (strcmp(argv[i], "-s") == 0)
				sscanf(argv[i + 1], "%lf", mas_pixel);
			else if (strcmp(argv[i], "-w") == 0)
				sscanf(argv[i + 1], "%hu", axis_len);
			else if (strcmp(argv[i], "-e") == 0)
				sscanf(argv[i + 1], "%ld", nelements);
			else if (strcmp(argv[i], "-n") == 0)
				sscanf(argv[i + 1], "%ld", niter);
			else if (strcmp(argv[i], "-ps") == 0)
				sscanf(argv[i + 1], "%lf", &reg_param[REG_PRIORIMAGE]);
			else if (strcmp(argv[i], "-de") == 0)
				sscanf(argv[i + 1], "%lf", &reg_param[REG_DARKENERGY]);
			else if (strcmp(argv[i], "-en") == 0)
				sscanf(argv[i + 1], "%lf", &reg_param[REG_ENTROPY]);
			else if (strcmp(argv[i], "-ud") == 0)
				sscanf(argv[i + 1], "%lf", &reg_param[REG_SPOT]);
			else if (strcmp(argv[i], "-tv") == 0)
				sscanf(argv[i + 1], "%lf", &reg_param[REG_TV]);
			else if (strcmp(argv[i], "-la") == 0)
				sscanf(argv[i + 1], "%lf", &reg_param[REG_LAP]);
			else if (strcmp(argv[i], "-l0") == 0)
				sscanf(argv[i + 1], "%lf", &reg_param[REG_L0]);
			else if (strcmp(argv[i], "-ts") == 0)
				sscanf(argv[i + 1], "%lf", &reg_param[REG_TRANSPECL2]);
			else if (strcmp(argv[i], "-f_any") == 0)
				sscanf(argv[i + 1], "%lf", f_anywhere);
			else if (strcmp(argv[i], "-f_copy") == 0)
				sscanf(argv[i + 1], "%lf", f_copycat);
			else if (strcmp(argv[i], "-d") == 0)
				sscanf(argv[i + 1], "%ld", depth);
			else if (strcmp(argv[i], "-chains") == 0)
				sscanf(argv[i + 1], "%d", nchains);
			else if (strcmp(argv[i], "-threads") == 0)
				sscanf(argv[i + 1], "%d", nthreads);
			else if (strcmp(argv[i], "-tempschedc") == 0)
				sscanf(argv[i + 1], "%lf", tempschedc);
			else if (strcmp(argv[i], "-fv") == 0)
				sscanf(argv[i + 1], "%lf", fov);
			else if (strcmp(argv[i], "-ct") == 0)
				sscanf(argv[i + 1], "%lf", chi2_temp);
			else if (strcmp(argv[i], "-fc") == 0)
				sscanf(argv[i + 1], "%lf", chi2_target);
			else if (strcmp(argv[i], "-tm") == 0)
				sscanf(argv[i + 1], "%lf", tmin);
			else if (strcmp(argv[i], "-pa") == 0)
				sscanf(argv[i + 1], "%lf", prob_auto);
			else if (strcmp(argv[i], "-uvtol") == 0)
				sscanf(argv[i + 1], "%lf", uvtol);
			else if (strcmp(argv[i], "-o") == 0) {
				sscanf(argv[i + 1], "%s", output_filename);
				if (strstr(output_filename, ".fits") != NULL) {
					output_filename[strlen(output_filename) - 5] = 0;
				}
			} else if (strcmp(argv[i], "-i") == 0)
				sscanf(argv[i + 1], "%s", init_filename);
			else if (strcmp(argv[i], "-p") == 0)
				sscanf(argv[i + 1], "%s", prior_filename);
			else if (strcmp(argv[i], "-v2s") == 0) {
				sscanf(argv[i + 1], "%lf", v2s);
			} else if (strcmp(argv[i], "-v2a") == 0) {
				sscanf(argv[i + 1], "%lf", v2a);
			} else if (strcmp(argv[i], "-t3amps") == 0) {
				sscanf(argv[i + 1], "%lf", t3amps);
			} else if (strcmp(argv[i], "-t3ampa") == 0) {
				sscanf(argv[i + 1], "%lf", t3ampa);
			} else if (strcmp(argv[i], "-t3phia") == 0) {
				sscanf(argv[i + 1], "%lf", t3phia);
			} else if (strcmp(argv[i], "-t3phis") == 0) {
				sscanf(argv[i + 1], "%lf", t3phis);
			} else if (strcmp(argv[i], "-visamps") == 0) {
				sscanf(argv[i + 1], "%lf", visamps);
			} else if (strcmp(argv[i], "-visampa") == 0) {
				sscanf(argv[i + 1], "%lf", visampa);
			} else if (strcmp(argv[i], "-visphis") == 0) {
				sscanf(argv[i + 1], "%lf", visphis);
			} else if (strcmp(argv[i], "-visphia") == 0) {
				sscanf(argv[i + 1], "%lf", visphia);
			} else if (strcmp(argv[i], "-fs") == 0) {
				sscanf(argv[i + 1], "%lf", fluxs);
			} else if (strcmp(argv[i], "-cv") == 0) {
				sscanf(argv[i + 1], "%lf", cvfwhm);
			} else if (strcmp(argv[i], "-P") == 0) {
				reg_param[REG_MODELPARAM] = 1.0;
				for (j = i + 1; j < *argc; j++) {
					/* If the next parameter is "-" followed by a non-numeric character,
					 *            then we have reached the end of the parameter list. */
					if (argv[j][0] == 45)
						if (argv[j][1] > 64)
							break;
					nparams++;
				}
				for (j = 0; j < nparams; j++)
					sscanf(argv[i + 1 + j], "%lf", &(init_params[j]));
				i += nparams - 1;
			} else if (strcmp(argv[i], "-S") == 0) {
				/* NB the way this is set up, the -P option must come before a -S option */
				for (j = 0; j < nparams; j++)
					sscanf(argv[i + 1 + j], "%lf", &(init_stepsize[j]));
				i += nparams - 1;
			} else if (strcmp(argv[i], "-wavchan") == 0) 
			  {
				int wavchaninfo = 0;
				for (j = i + 1; j < *argc; j++) {
					/* If the next parameter is "-" followed by a non-numeric character,
					 *            then we have reached the end of the parameter list. */
					if (argv[j][0] == 45)
						if (argv[j][1] > 64)
							break;
					wavchaninfo++;
				}
				printf("Command line  -- set to polychromatic reconstruction\n");

				*nwavr = wavchaninfo / 3;
				printf("DEBUG nwavr = %d\n", *nwavr);
				wavmin = malloc( *nwavr * sizeof(double));
				wavmax = malloc( *nwavr * sizeof(double));
				printf("DEBUG wavmin pointer: %p \n", &wavmin[0]);
				for (j = 0; j < *nwavr; j++) {
					sscanf(argv[i + 1 + 3 * j], "%ld", &k);
					sscanf(argv[i + 1 + 3 * j + 1], "%lf", &wavmin[k]);
					sscanf(argv[i + 1 + 3 * j + 2], "%lf", &wavmax[k]);
					if (k != j) {
						printf("Command line  -- Misformed channel information \n");
						getchar();
						break;
					} else {
						printf("Command line  -- channel %ld = %lg m to\t %lg m\n", k, wavmin[k], wavmax[k]);
					}
				}
				i += wavchaninfo - 1;
			} else {
				printf("Command line  -- invalid option %s. Type squeeze -h for help.\n", argv[i]);
				return FALSE;
			}
			i++;
		} else {
			printf("Command line  -- invalid option. Type squeeze -h for help.\n");
			return FALSE;
		}

	}
	return TRUE;
}

void print_diagnostics(int iChain, long current_iter, long nvis, long nv2, long nt3, long nt3phi, long nt3amp, long nvisamp, long nvisphi, double chi2v2,
		double chi2t3amp, double chi2t3phi, double chi2visphi, double chi2visamp, double lPosterior, double lPrior, double lLikelihood, const double* reg_param,
		const double* reg_value, const double* centroid_image_x, const double* centroid_image_y, long nelements, int nwavr, long niter,
		const double* temperature, double prob_movement, const double* params, const double* stepsize) {
	// Note: current_iter = i / (nwavr * nelements) + 1

	long w, j;
	char diagnostics[250];
	int diagnostics_used = 0;

	/* Print output to screen (or wherever stdout is piped to) */

	// compute reduced chi2
	if (nv2 > 0)
		chi2v2 /= (double) nv2;
	if (nt3amp > 0)
		chi2t3amp /= (double) nt3amp;
	if (nt3phi > 0)
		chi2t3phi /= (double) nt3phi;
	if (nvisphi > 0)
		chi2visphi /= (double) nvisphi;
	if (nvisamp > 0)
		chi2visamp /= (double) nvisamp;

	for (w = 0; w < nwavr; w++) // diagnostics
			{

		if (nwavr > 1)
			diagnostics_used = snprintf(diagnostics, 250, "Chain: %d Chan: %ld lPost:%8.1f lPrior:%8.1f lLike:%9.1f ", iChain, w, lPosterior, lPrior,
					lLikelihood);
		else
			diagnostics_used = snprintf(diagnostics, 250, "Chain: %d lPost:%8.1f lPrior:%8.1f lLike:%9.1f ", iChain, lPosterior, lPrior, lLikelihood);

		if (nv2 > 0)
			diagnostics_used += snprintf(diagnostics + diagnostics_used, 250 - diagnostics_used,
			TEXT_COLOR_RED "V2:%5.2f " TEXT_COLOR_BLACK, chi2v2);
		if (nt3amp > 0)
			diagnostics_used += snprintf(diagnostics + diagnostics_used, 250 - diagnostics_used,
			TEXT_COLOR_BLUE "T3A:%5.2f " TEXT_COLOR_BLACK, chi2t3amp);
		if (nt3phi > 0)
			diagnostics_used += snprintf(diagnostics + diagnostics_used, 250 - diagnostics_used,
			TEXT_COLOR_GREEN "T3P:%5.2f " TEXT_COLOR_BLACK, chi2t3phi);
		if (nvisamp > 0)
			diagnostics_used += snprintf(diagnostics + diagnostics_used, 250 - diagnostics_used, "VA:%5.2f ", chi2visamp);
		if (nvisphi > 0)
			diagnostics_used += snprintf(diagnostics + diagnostics_used, 250 - diagnostics_used, "VP:%5.2f ", chi2visphi);
		if (reg_param[REG_PRIORIMAGE] > 0)
			diagnostics_used += snprintf(diagnostics + diagnostics_used, 250 - diagnostics_used, "PRI:%5.2f ",
					reg_param[REG_PRIORIMAGE] * reg_value[w * NREGULS + REG_PRIORIMAGE]);
		if (reg_param[REG_ENTROPY] > 0)
			diagnostics_used += snprintf(diagnostics + diagnostics_used, 250 - diagnostics_used, "ENT:%5.2f ",
					reg_param[REG_ENTROPY] * reg_value[w * NREGULS + REG_ENTROPY]);
		if (reg_param[REG_DARKENERGY] > 0)
			diagnostics_used += snprintf(diagnostics + diagnostics_used, 250 - diagnostics_used, "DEN:%5.2f ",
					reg_param[REG_DARKENERGY] * reg_value[w * NREGULS + REG_DARKENERGY]);
		if (reg_param[REG_SPOT] > 0)
			diagnostics_used += snprintf(diagnostics + diagnostics_used, 250 - diagnostics_used, "UD:%5.2f ",
					reg_param[REG_SPOT] * reg_value[w * NREGULS + REG_SPOT]);
		if (reg_param[REG_TV] > 0)
			diagnostics_used += snprintf(diagnostics + diagnostics_used, 250 - diagnostics_used, "TV:%5.2f ",
					reg_param[REG_TV] * reg_value[w * NREGULS + REG_TV]);
		if (reg_param[REG_LAP] > 0)
			diagnostics_used += snprintf(diagnostics + diagnostics_used, 250 - diagnostics_used, "LAP:%5.2f ",
					reg_param[REG_LAP] * reg_value[w * NREGULS + REG_LAP]);
		if (reg_param[REG_L0] > 0)
			diagnostics_used += snprintf(diagnostics + diagnostics_used, 250 - diagnostics_used, "L0:%5.2f ",
					reg_param[REG_L0] * reg_value[w * NREGULS + REG_L0]);
		if (reg_param[REG_TRANSPECL2] > 0)
			diagnostics_used += snprintf(diagnostics + diagnostics_used, 250 - diagnostics_used, "TS:%5.2f ",
					reg_param[REG_TRANSPECL2] * reg_value[REG_TRANSPECL2]);
		if (reg_param[REG_CENTERING] > 0)
			diagnostics_used += snprintf(diagnostics + diagnostics_used, 250 - diagnostics_used, "C:%5.1f XY:(%5.2f,%5.2f) ",
					reg_param[REG_CENTERING] * reg_value[w * NREGULS + REG_CENTERING], centroid_image_x[w] / nelements, centroid_image_y[w] / nelements);

		diagnostics_used += snprintf(diagnostics + diagnostics_used, 250 - diagnostics_used, "E: %5ld MPr: %4.2f T: %5.2f Iter: %4ld of %4ld", nelements,
				prob_movement, temperature[iChain], current_iter, niter);

		puts(diagnostics);
	}
	fflush (stdout);

	if (nparams > 0) {
		printf("Chain: %d Model Parameters: ", iChain);
		for (j = 0; j < nparams; j++)
			printf("P[%ld]: %7.5g +/- %7.5g  ", j, params[j], stepsize[j]);
		printf("\n");
	}
}

void compute_lPrior(double* lPrior, const long chan, const double* reg_param, const double* reg_value) {
	*lPrior = reg_param[REG_PRIORIMAGE] * reg_value[chan * NREGULS + REG_PRIORIMAGE]
	        + reg_param[REG_CENTERING]  * reg_value[chan * NREGULS + REG_CENTERING]
		+ reg_param[REG_ENTROPY]    * reg_value[chan * NREGULS + REG_ENTROPY]
                + reg_param[REG_MODELPARAM] * reg_value[chan * NREGULS + REG_MODELPARAM]
		+ reg_param[REG_DARKENERGY] * reg_value[chan * NREGULS + REG_DARKENERGY]
                + reg_param[REG_SPOT]       * reg_value[chan * NREGULS + REG_SPOT]
		+ reg_param[REG_TV]         * reg_value[chan * NREGULS + REG_TV]
                + reg_param[REG_LAP]        * reg_value[chan * NREGULS + REG_LAP]
		+ reg_param[REG_L0]         * reg_value[chan * NREGULS + REG_L0]
                + reg_param[REG_TRANSPECL2] * reg_value[REG_TRANSPECL2];

}

