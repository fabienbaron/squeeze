/*   modelcode_ud bwsmearing.c: An offset UD with bw smearing parameter
// params(0) = delta RA in mas from phase center (EAST => +)
// params(1) = delta DEC in mas from phase center (NORTH => +)
// params(2) = secondary brightness
// params(3) = UD size of secondary (mas)
// params(4) = Fractional Bandwidth
// params(5) = UD size of primary (mas) [located at origin]
// params(6) = primarybrightness
/* Globals: nparams, nbaselines, u, v*/
int model_vis(const double *params, double complex *modvis, double *lPriorModel, double *flux_frac)
{
    int status = 0;
    long i;
    double tempd, vis_secondary = 1.0, vis_primary = 1.0, vis_bw = 1.0, delta_ra, delta_dec;
    double complex phase_factor;
    delta_dec = params[1] / (double)206264806.2;
    delta_ra  = params[0] / (double)206264806.2;

    for(i = 0; i < nuv; i++)
        {
            if(params[3] > 0)
                {
                    tempd = PI * params[3] / MAS_RAD * sqrt(u[i] * u[i] + v[i] * v[i]) + 1e-15;
                    vis_secondary = 2.0 * j1(tempd) / tempd;
                }
            else
              vis_secondary = 1.0;

            if(params[5] > 0)
                {
                    tempd = PI * params[5] / MAS_RAD * sqrt(u[i] * u[i] + v[i] * v[i]) + 1e-15;
                    vis_primary = 2.0 * j1(tempd) / tempd;
                }
            else
              vis_primary = 1.0;


            if(params[4] > 0)
                {
                    tempd = PI * params[4] * (u[i] * delta_ra + v[i] * delta_dec) + 1e-15;
                    vis_bw = sin(tempd) / tempd;
                }
            else
              vis_bw = 1.0;

            vis_secondary *= params[2] * vis_bw;
            vis_primary *= params[6];
            phase_factor = cexp(-2.0 * PI * (u[i] * delta_ra + v[i] * delta_dec);
            modvis[i] = vis_primary + vis_secondary * phase_factor ;
        }

    *flux_frac = 1.0 - params[2] - params[6];

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
