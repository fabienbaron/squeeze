/*   modelcode_ud bwsmearing.c:

A polychromatic binary model with primary uniform disc star at the center of the image,
the secondary as an offset uniform disc with bandwidth smearing parameter.
The stars follow the lambda^-4 law for their fluxes, while the environment follows lambda^env_ind (user-defined)

with parameters :
      delta_ra: secondary delta RA in mas from image center (EAST => +)
      delta_dec: secondary delta DEC in mas from image center (NORTH => +)
      f_primary_ref : the stellar flux fraction at lambda_ref
      f_secondary_ref : the stellar flux fraction at lambda_ref
      diam_primary : the size of the primary, uniform disc at lambda_ref in mas
      diam_secondary : the size of the primary, uniform disc at lambda_ref in mas
      bws: bandwidth smearing parameter
      lambda_ref : reference wavelength (H band)
      env_ind : the flux power law index for the environment (== the image)

 params(0) = primary flux
 params(1) = UD size of primary (mas) [located at origin]
 params(2) = environment index
 params(3) = reference wavelength for above parameters
 params(4) = background flux
 params(5) = background power law
 Globals: nparams, nbaselines, u, v */

const char *model_param_names[6] = {"flux_star   ", "UD        ", "env_indx  ", "lambda_0 ", "flux_bg  ", "bg_indx "};

extern double j1(double);
int model_vis(const double *params, double complex *modvis, double *lPriorModel, double *flux_frac) {
  int status = 0;
  long i;
  double lPriorParams[6];
  double f_primary_ref =  params[0];
  double ud_primary = params[1];
  double env_ind = params[2];
  double lambda_ref = params[3];
  double f_bg_ref = params[4];
  double bg_ind = params[5];
  double tempd, f_primary, f_env, f_bg;
  double complex vis_primary;

  for(i = 0; i < nuv; i++)
        {
            // Compute visibilities for unity fluxes
            if (ud_primary > 0)
            {
              tempd = M_PI * ud_primary / MAS_RAD * sqrt(u[i] * u[i] + v[i] * v[i]) + 1e-15;
              vis_primary = 2.0 * j1(tempd) / tempd + 0*I;
            }
            else
                 vis_primary = 1.0 + 0*I;
            // Compute fluxes
            f_primary = f_primary_ref * pow(uv_lambda[i] / lambda_ref, -4.0); // Stellar flux primary
            f_bg = f_bg_ref * pow(uv_lambda[i] / lambda_ref, bg_ind);
            f_env = (1.0 - f_primary_ref - f_bg) * pow(uv_lambda[i] / lambda_ref, env_ind); // environment flux
            flux_frac[i] = f_env / (f_primary + f_env + f_bg ); // Flux ratio Disc/Total flux = Flux ratio Image/(Image + Model)
            // Visibilities for the model
            modvis[i] = f_primary * vis_primary / (f_primary + f_env + f_bg );
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

  // Note: no priors on reference wavelength, most likely will be fixed anyway
  lPriorParams[3] = 0;

  *lPriorModel = lPriorParams[0] + lPriorParams[1] + lPriorParams[2]
               + lPriorParams[3] + lPriorParams[4] + lPriorParams[5];

  return (status != 0);
}
