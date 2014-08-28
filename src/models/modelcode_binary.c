/*   modelcode_ud.c: A central uniform disk model input.

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

// params(0) = Separation (mas)
// params(1) = Position Angle (degs,
//          E of N, pointing from primary -> secondary)
// params(2) = Primary brightness
// params(3) = secondary brightness
// params(4) = UD size of primary (mas)
// params(5) = UD size of secondary (mas)
//
//
/* Globals: nparams, nbaselines, u, v*/

int model_vis(double *params, double complex *modvis, double *logl, double *flux_frac)
{
    int status = 0;
    long i;
    double tempd, vis_primary = 1.0, vis_secondary = 1.0, delta_ra, delta_dec;
    double complex phase_factor;
    double l0, l1;


    for(i = 0; i < nbaselines; i++)
        {
            if(params[4] > 0)
                {
                    tempd = PI * params[4] / MAS_RAD * sqrt(u[i] * u[i] + v[i] * v[i]) + 1e-15;
                    vis_primary = 2 * j1(tempd) / tempd;
                }
            else vis_primary = 1.0;
            if(params[5] > 0)
                {
                    tempd = PI * params[5] / MAS_RAD * sqrt(u[i] * u[i] + v[i] * v[i]) + 1e-15;
                    vis_secondary = 2 * j1(tempd) / tempd;
                }
            else vis_secondary = 1.0;
            vis_primary *= params[2];
            vis_secondary *= params[3];
            delta_dec = params[0] * sin((params[1] + 90.) * PI / 180.) / (double)206264806.2;
            delta_ra  = params[0] * cos((params[1] + 90.) * PI / 180.) / (double)(-206264806.2);
            //    printf("delta dec and ra un urad: %lf %lf\n", delta_dec*1000000., delta_ra*1000000.);
            phase_factor = cos(-2.0 * PI * (u[i] * delta_ra + v[i] * delta_dec)) +
                           I * sin(-2.0 * PI * (u[i] * delta_ra + v[i] * delta_dec));
            //printf("Sin and Cos parts: %lf\n", u[i]*delta_ra+v[i]*delta_dec);
            // printf(" %lf %lf %lf \n",u[i],v[i],tempd);
            modvis[i] = vis_primary + vis_secondary * phase_factor;
//    printf(" %ld %lf %lf: \n",i,tempd,modvis[i]);
        }
//  printf("Primary,Secondary vis, phase factor, params[2], params[3]: %lf %lf %lf %lf %lf\n",
//   vis_primary, vis_secondary, creal(phase_factor), params[2], params[3]);
    write_best_oifits("oldmod.bispectrum", modvis);
    /* Do custom likelihood */

// For vis0, allow equal probability between 0,1
    if(params[2] >= 0 && params[2] <= 1) l0 = 1.;
    if(params[2] < 0 || params[2] > 1) l0 = 1e-6;
    if(params[3] >= 0 && params[3] <= 1) l1 = 1.;
    if(params[3] < 0 || params[3] > 1) l1 = 1e-6;

    *logl = -1 * log(l0 * l1);
    *flux_frac = 1.0 - params[2] - params[3];


    return (status != 0);
}
