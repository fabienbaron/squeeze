/*
   Version: 2.0
   SQUEEZE 2 - Image reconstruction software for optical interferometry,
   based on Monte-Carlo Markov Chain algorithms.

   Copyright (c) 2006-2013 Fabien Baron, John Monnier, Michael Ireland, Jacques Kluska

   Based on MACIM written by Dr. Michael Ireland (Macquarie University) and
   Pr. John Monnier (University of Michigan).
   Also based on SQUEEZE 1.2 written by Pr. Fabien Baron (Georgia State University).

   SQUEEZE is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   SQUEEZE is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Squeeze.  If not, see <http://www.gnu.org/licenses/>.
------------------------------
 FB model code based on collaboration with Jacques Kluska on SPARCO.

 The chromatism is defined as follows (Kluska et al. 2012) :
        fs0 (lambda/ lambda_0)^-4 V_star + (1-fs0)*(lambda/lambda_0)^d_ind * V_env
V_tot = --------------------------------------------------------------------------
        fs0 (lambda/ lambda_0)^-4 + (1-fs0)*(lambda/lambda_0)^d_ind

with parameters :
       lambda_0 : reference wavelenght (H band)
       fs0 : the stellar flux fraction at lambda_0
       diam : the size of the uniform disc at lambda_0 in mas
       d_ind : the flux power law index for the environement (here the image)

Notes:
       V_env : the visibilities of the environment (here the image)
       V_star : visibilities of the star (here a UD)


------------------------------
Globals: nparams, nbaselines, u, v*/
extern double j1(double);
int model_vis(const double *params, double complex *modvis, double *lPriorModel, double *flux_frac_0)
{
    int status = 0;
    long i;
    double tempd, diam;
    double lPriorParams[4];
    double fs0, d_ind, fs, lambda_0, fd;

    lambda_0 = params[0];
    fs0 = params[1];
    diam = params[2];
    d_ind = params[3];


    // Compute flux ratio
    for(i = 0; i < nuv; i++)
        {
            fs = fs0 * pow((uv_lambda[i] / lambda_0), -4.); // Stellar flux
            fd = (1. - fs0) * pow((uv_lambda[i] / lambda_0), d_ind); // Disc flux
            flux_frac_0[i] = fd / (fs + fd); // Flux ratio Disc/Total flux = Flux ratio Image/(Image + Model)
        }

    // Compute model visibilities
    for(i = 0; i < nuv; i++)
        {
            tempd = M_PI * diam / MAS_RAD * sqrt(u[i] * u[i] + v[i] * v[i]) + 1e-15;
            //    tempd = PI*params[4]*(u[i]*delta_ra+v[i]*delta_dec)+1e-15;
            //    vis_bw=sin(tempd)/tempd;
            modvis[i] = 2. * (1. - flux_frac_0[i]) * j1(tempd) / tempd;
        }

    // Enforce positivity and bounds of parameters

    // Reference wavelength
    if(params[0] >= 0)
        lPriorParams[0] = 0;
    if(params[0] < 0)
        lPriorParams[0] = 1e99;

    // Flux ratio (positive, <= 1)
    if(params[1] >= 0 && params[1] <= 1)
        lPriorParams[1] = 0;
    if(params[1] < 0 || params[1] > 1)
        lPriorParams[1] = 1e99;

    // Diameter (positive)
    if(params[2] >= 0)
        lPriorParams[2] = 0;
    else
        lPriorParams[2] = 1e99;


    if(abs(params[3]) >= 10)
        lPriorParams[3] = 1e99;
    else
        lPriorParams[3] = 0;


    *lPriorModel = lPriorParams[0] + lPriorParams[1] + lPriorParams[2] +lPriorParams[3] ;

    return (status != 0);
}

