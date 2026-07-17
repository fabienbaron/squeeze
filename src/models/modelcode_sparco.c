/*   modelcode_sparco.c:

SPARCO-style chromatic model (after Kluska et al. 2012), with an additional
over-resolved background component.

The image (the "environment") is reconstructed grey; all of the chromatism lives
here, in the parametric model. Each uv point gets

           f_star (l/l0)^-4 V_star + f_bg (l/l0)^bg_ind * 0 + f_env (l/l0)^env_ind V_env
   V_tot = ---------------------------------------------------------------------------
                  f_star (l/l0)^-4 + f_bg (l/l0)^bg_ind + f_env (l/l0)^env_ind

with f_env = 1 - f_star - f_bg. The background is over-resolved (V = 0), so it
contributes flux but no visibility: it only dilutes. model_vis() returns the star
term in modvis[] and hands the environment's flux fraction back in flux_frac[],
which is what scales the reconstructed image.

The star is a uniform disc of diameter params[1]; set it to 0 for a point source.

 params(0) = f_star   : stellar flux fraction at lambda_0
 params(1) = UD       : uniform-disc diameter of the star at lambda_0, in mas (0 = point source)
 params(2) = env_indx : flux power-law index of the environment (the image)
 params(3) = lambda_0 : reference wavelength (m), must be > 0
 params(4) = f_bg     : background flux fraction at lambda_0
 params(5) = bg_indx  : flux power-law index of the background

 Globals: nuv, u, v, uv_lambda */

#define NPARAMS_SPARCO 6
static const char *model_param_names_sparco[NPARAMS_SPARCO] = {"flux_star   ", "UD        ", "env_indx  ", "lambda_0 ", "flux_bg  ", "bg_indx "};

static int model_vis_sparco(const double *params, double complex *modvis, double *lPriorModel, double *flux_frac) {
  int status = 0;
  long i;
  double lPriorParams[NPARAMS_SPARCO];
  double f_primary_ref = params[0];
  double ud_primary = params[1];
  double env_ind = params[2];
  double lambda_ref = params[3];
  double f_bg_ref = params[4];
  double bg_ind = params[5];
  double tempd, f_primary, f_env, f_bg;
  double complex vis_primary;

  for (i = 0; i < nuv; i++) {
    // Compute visibilities for unity fluxes
    if (ud_primary > 0) {
      tempd = M_PI * ud_primary / MAS_RAD * sqrt(u[i] * u[i] + v[i] * v[i]) + 1e-15;
      vis_primary = 2.0 * j1(tempd) / tempd + 0 * I;
    } else
      vis_primary = 1.0 + 0 * I;
    // Compute fluxes
    f_primary = f_primary_ref * pow(uv_lambda[i] / lambda_ref, -4.0); // Stellar flux primary
    f_bg = f_bg_ref * pow(uv_lambda[i] / lambda_ref, bg_ind);
    f_env = (1.0 - f_primary_ref - f_bg) * pow(uv_lambda[i] / lambda_ref, env_ind); // environment flux
    flux_frac[i] = f_env / (f_primary + f_env + f_bg); // Flux ratio Disc/Total flux = Flux ratio Image/(Image + Model)
    // Visibilities for the model
    modvis[i] = f_primary * vis_primary / (f_primary + f_env + f_bg);
  }

  // Flux ratios (strictly positive and <=1 )
  if (params[0] >= 0)
    lPriorParams[0] = 0;
  else
    lPriorParams[0] = 1e99;

  if (params[4] >= 0)
    lPriorParams[4] = 0;
  else
    lPriorParams[4] = 1e99;

  // Diameters
  if (params[1] >= 0)
    lPriorParams[1] = 0;
  else
    lPriorParams[1] = 1e99;

  // Environment index -- could be many things
  if (params[2] > -5 && params[2] <= 5)
    lPriorParams[2] = 0;
  else
    lPriorParams[2] = 1e99;

  // Background index -- same
  if (params[5] > -5 && params[5] <= 5)
    lPriorParams[5] = 0;
  else
    lPriorParams[5] = 1e99;

  // Note: no priors on reference wavelength, most likely will be fixed anyway.
  // It is validated (> 0) at startup instead, since we divide by it.
  lPriorParams[3] = 0;

  *lPriorModel = lPriorParams[0] + lPriorParams[1] + lPriorParams[2] + lPriorParams[3] + lPriorParams[4] + lPriorParams[5];

  return (status != 0);
}
