/*   modelcode.c: parametric model registry and dispatcher.

The model is chosen at RUNTIME with -modeltype <name>. Editing and recompiling
this file to switch model is no longer necessary.

To add a model:
  1. Write models/modelcode_<name>.c defining
       #define NPARAMS_<NAME> n
       static const char *model_param_names_<name>[NPARAMS_<NAME>] = {...};
       static int model_vis_<name>(const double *params, double complex *modvis,
                                   double *lPriorModel, double *flux_frac);
     model_vis_<name>() fills modvis[] with the model visibilities and flux_frac[]
     with the fraction of the total flux carried by the image ("environment").
  2. #include it below and add a row to squeeze_models[].
The parameter count is checked against -P at startup, and the names are used when
reporting fitted parameters.

Globals available to a model: nuv, u, v, uv_lambda, uv_dlambda, MAS_RAD.
*/

extern double j1(double); // Bessel J1, from libm; used by the uniform-disc models

#include "modelcode_sparco.c"
#include "modelcode_binary_bwsmearing.c"

const model_desc squeeze_models[] = {
    {"sparco", NPARAMS_SPARCO, model_param_names_sparco, model_vis_sparco},
    {"binary", NPARAMS_BINARY, model_param_names_binary, model_vis_binary},
};
const int n_squeeze_models = sizeof(squeeze_models) / sizeof(squeeze_models[0]);

int modeltype = 0; // default: sparco

int model_vis(const double *params, double complex *modvis, double *lPriorModel, double *flux_frac) {
  return squeeze_models[modeltype].model_vis_fn(params, modvis, lPriorModel, flux_frac);
}
