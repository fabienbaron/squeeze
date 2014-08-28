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
// params(2) = Brightness
// params(3) = UD size (mas)
//
/* Globals: nparams, nbaselines, u, v*/

int model_vis(double *params, double complex *modvis, double *logl, double *flux_frac)
{
    int status = 0;
    long i;
    double tempd, vis_primary = 1.0, delta_ra, delta_dec;
    double complex phase_factor;
    double l0, l1, l2, l3;
    double p0, p1, p2, p3; //priors
    double sig0, sig1, sig2, sig3; // width of priors

// Enter priors:
    p0 = 0.0;
    sig0 = 3.0 ;

    p1 = 0.0;
    sig1 = 3.;

    p2 = .16;
    sig2 = 0.4;

    p3   = 0.5;
    sig3 = 1.0;

    delta_dec = params[1] / MAS_RAD;
    delta_ra  = params[0] / MAS_RAD;

    for(i = 0; i < nbaselines; i++)
        {
            if(params[3] > 0)
                {
                    tempd = PI * params[3] / MAS_RAD * sqrt(u[i] * u[i] + v[i] * v[i]) + 1e-15;
                    vis_primary = 2 * j1(tempd) / tempd;
                }
            else vis_primary = 1.0;

            vis_primary *= params[2];
            phase_factor = cexp(-2.0 * PI * (u[i] * delta_ra + v[i] * delta_dec)) ;
            modvis[i] = vis_primary * phase_factor;
        }

// For vis0, allow equal probability between 0,1
    l0 = exp(-.5 * pow((params[0] - p0) / sig0 , 2));
    l1 = exp(-.5 * pow((params[1] - p1) / sig1 , 2));
    l2 = exp(-.5 * pow((params[2] - p2) / sig2 , 2));
    l3 = exp(-.5 * pow((params[3] - p3) / sig3 , 2));

    if(params[2] < -0.3 || params[2] > 0.)
        l2 = 1e-6;

    if(params[3] < 0. || params[3] > 1.0)
        l3 = 1e-6;

    if((sqrt(params[0]*params[0] + params[1] * params[1]) + params[3])  > 1.375)
        {
            l0 = 1e-7; l1 = 1e-7;
        }


    *logl = -1.0 * log(l0 * l1 * l2 * l3);
    *flux_frac = 1.0 - params[2];

    return (status != 0);

}
