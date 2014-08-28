/*   modelcode_ud bwsmearing.c: An offset UD with bw smearing parameter
  JDM

  MACIM - software for the creation of images consistent with data sets
  from optical interferometry, using a Monte-Carlo Markov Chain algorithm.

  Copyright (c) 2006 California Institute of Technology.  Written by
  Dr. Michael Ireland with contributions from Prof. John Monnier (University
  of Michigan). For comments or questions about this software, please contact
  the author at mireland@gps.caltech.edu.

  This program is free software; you can redistribute it and/or modify it
  under the terms of the GNU General Public License as  published by the
  Free Software Foundation; either version 2 of the License, or (at your
  option) any later version.

  This program is provided "as is" and distributed in the hope that it
  will be useful, but WITHOUT ANY WARRANTY; without even the implied
  warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  In no
  event shall California Institute of Technology be liable to any party
  for direct, indirect, special, incidental or consequential damages,
  including lost profits, arising out of the use of this software and its
  documentation, even if the California Institute of Technology has been
  advised of the possibility of such damage.   The California Institute of
  Technology has no obligation to provide maintenance, support, updates,
  enhancements or modifications.  See the GNU General Public License for
  more details.

  You should have received a copy of the GNU General Public License along
  with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.

------------------------------
 JDM test code for imager program.
 Using JYoung code.

 model_UD: Uniform Disk model centered at (0,0)
 params: visib at uv=0, diameter in milliarcseconds, unresolved flux
 inputs: u,v, are in the wavelengths.
 logl is the a priori -log likelihood of the parameter combination.
------------------------------
To change this: ANY model is possible here, as long as model complex
visibilities can be returned based on a combinarion of up to MAX_PARAMS
(default 20 in macim.h) parameters.
*/

// params(0) = delta RA in mas from phase center (EAST => +)
// params(1) = delta DEC in mas from phase center (NORTH => +)
// params(2) = secondary brightness
// params(3) = UD size of secondary (mas)
// params(4) = Fractional Bandwidth
// params(5) = UD size of primary (mas) [located at origin]
// params(6) = primarybrightness
/* Globals: nparams, nbaselines, u, v*/

int model_vis(double *params, double complex *modvis, double *logl, double *flux_frac)
{
    int status = 0;
    long i;
    double tempd, vis_secondary = 1.0, vis_primary = 1.0, vis_bw = 1.0, delta_ra, delta_dec;
    double complex phase_factor;
    double l0, l1, l2, l3, l4, l5, l6, l_low;
    double p0, p1, p2, p3, p4, p5, p6; //priors
    double sig0, sig1, sig2, sig3, sig4, sig5, sig6; // width of priors

// Enter priors:
// Algol 2009Aug10

    p0 = 14.8;
//p0=15.6;
    sig0 = .1 ;

    p1 = -3.0;
//p1=-3.7;
    sig1 = .1;

    p2 = .1;
    sig2 = .03; //.01;

    p3 = .6;
    sig3 = .005;

    p4 = 1. / 45;
    sig4 = p4 * .01;

    p5 = .95;
    sig5 = .01; //.01;

    p6 = .8;
    sig6 = .003;


    delta_dec = params[1] / (double)206264806.2;
    delta_ra  = params[0] / (double)206264806.2;

    for(i = 0; i < nbaselines; i++)
        {
            if(params[3] > 0)
                {
                    tempd = PI * params[3] / MAS_RAD * sqrt(u[i] * u[i] + v[i] * v[i]) + 1e-15;
                    vis_secondary = 2 * j1(tempd) / tempd;
                }
            else vis_secondary = 1.0;

            if(params[5] > 0)
                {
                    tempd = PI * params[5] / MAS_RAD * sqrt(u[i] * u[i] + v[i] * v[i]) + 1e-15;
                    vis_primary = 2 * j1(tempd) / tempd;
//      printf("uv: vis1: %f %f\n",sqrt(u[i]*u[i]+v[i]*v[i]), vis_primary);
                }
            else vis_primary = 1.0;


            if(params[4] > 0)
                {
                    tempd = PI * params[4] * (u[i] * delta_ra + v[i] * delta_dec) + 1e-15;
                    vis_bw = sin(tempd) / tempd;
                    //   printf("uv: vis_bw: %f %f\n",sqrt(u[i]*u[i]+v[i]*v[i]), vis_bw);
                }
            else vis_bw = 1.0;

            vis_secondary *= params[2] * vis_bw;
            vis_primary *= params[6];
//printf("i,vis1,vis2 %i  %f %f \n",i,vis_primary,vis_secondary);
//        printf("delta dec and ra un urad: %lf %lf\n", delta_dec*1000000., delta_ra*1000000.);

            phase_factor = cos(-2.0 * PI * (u[i] * delta_ra + v[i] * delta_dec)) +
                           I * sin(-2.0 * PI * (u[i] * delta_ra + v[i] * delta_dec));
            //printf("Sin and Cos parts: %lf\n", u[i]*delta_ra+v[i]*delta_dec);
            // printf(" %lf %lf %lf \n",u[i],v[i],tempd);
            modvis[i] = vis_primary + vis_secondary * phase_factor ;
//    printf(" %ld %lf %lf: \n",i,tempd,modvis[i]);
        }
//  printf("Primary,Secondary vis, phase factor, params[2], params[3]: %lf %lf %lf %lf %lf\n",
//   vis_primary, vis_secondary, creal(phase_factor), params[2], params[3]);
    write_best_oifits("oldmod.bispectrum", modvis);
    /* Do custom likelihood */

// For vis0, allow equal probability between 0,1
    l0 = exp(-.5 * pow((params[0] - p0) / sig0 , 2));
    l1 = exp(-.5 * pow((params[1] - p1) / sig1 , 2));
    l2 = exp(-.5 * pow((params[2] - p2) / sig2 , 2));
    l3 = exp(-.5 * pow((params[3] - p3) / sig3 , 2));
    l4 = exp(-.5 * pow((params[4] - p4) / sig4 , 2));
    l5 = exp(-.5 * pow((params[5] - p5) / sig5 , 2));
    l6 = exp(-.5 * pow((params[6] - p6) / sig6 , 2));

    l_low = 1e-5;
    if(params[2] < 0 || params[2] > 1) l2 = l_low;
    if(params[3] < 0) l3 = l_low;
    if(params[4] < 0 || params[4] > 1) l4 = l_low;
    if(params[5] < 0) l5 = l_low;
    if(params[6] < 0 || params[6] > 1) l6 = l_low;



    *logl = -1 * (log(l0) + log(l1) + log(l2) + log(l3) + log(l4) + log(l5) + log(l6));
//printf("logl %f\n",*logl);
//printf("logx %f %f %f %f %f %f %f\n",l0,l1,l2,l3,l4,l5,l6);
//printf("l6,p6,sig6,params[6] %f %f %f %f \n",l6,p6,sig6,params[6]);

    *flux_frac = 1.0 - params[2] - params[6];


    return (status != 0);
}
