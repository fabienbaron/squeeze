/*   modelcode_ud bwsmearing.c: An offset UD with bw smearing parameter
  JDM (MACIM) + FB (SQUEEZE port)

 model_UD: Uniform Disk model centered at (0,0)
 params: visib at uv=0, diameter in milliarcseconds, unresolved flux
 inputs: u,v, are in the wavelengths.
 logl is the a priori -log likelihood of the parameter combination.
------------------------------
To change this: ANY model is possible here, as long as model complex
visibilities can be returned based on a combinarion of up to MAX_PARAMS
(default 20 in macim.h) parameters.

params(0) = delta RA in mas from phase center (EAST => +)
params(1) = delta DEC in mas from phase center (NORTH => +)
params(2) = Brightness ratio of the primary
params(3) = UD size of primary (mas)
params(4) = Fractional Bandwidth

Globals: nparams, nbaselines, u, v

*/
extern double j1(double);
int model_vis(const double *params, double complex *modvis, double *lPriorModel, double *flux_frac_0)
{
    int status = 0;
    long i;
    double tempd, vis_primary = 1.0, vis_bw = 1.0, delta_ra, delta_dec;
    double lPriorParams[5];


    delta_dec = params[1] / MAS_RAD;
    delta_ra  = params[0] / MAS_RAD;

    for(i = 0; i < nuv; i++)
        {

            if(params[3] > 0)
                {
                    tempd = M_PI * params[3] / MAS_RAD * sqrt(u[i] * u[i] + v[i] * v[i]) + 1e-15;
                    vis_primary = 2.0 * j1(tempd) / tempd;
                }
            else
                vis_primary = 1.0;

            if(params[4] > 0)
                {
                    tempd = M_PI * params[4] * (u[i] * delta_ra + v[i] * delta_dec) + 1e-15;
                    vis_bw = sin(tempd) / tempd;
                }
            else
                vis_bw = 1.0;

            modvis[i] = vis_bw * vis_primary * cexp(-2.0 * I * M_PI * (u[i] * delta_ra + v[i] * delta_dec));

            flux_frac_0[i] = 1. - params[2];
        }

    // Enforce positivity and bounds of parameters

    // No constrain on parameters RA and DEC

    lPriorParams[0] = 0;
    lPriorParams[1] = 0;

    // Flux ratio of the primary (strictly positive)
    if(params[2] > 0 && params[2] <= 1)
        lPriorParams[2] = 0;
    else
        lPriorParams[2] = 1e99;

    // Diameter (strictly positive)
    if(params[3] > 0)
        lPriorParams[3] = 0;
    else
        lPriorParams[3] = 1e99;

    // Fractional bandwidth
    if(params[4] > 0 && params[4] <= 1)
        lPriorParams[4] = 0;
    else
        lPriorParams[4] = 1e99;

    *lPriorModel = lPriorParams[0] + lPriorParams[1] + lPriorParams[2] + lPriorParams[3] + lPriorParams[4] ;

    return (status != 0);
}
