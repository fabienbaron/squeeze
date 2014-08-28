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

/* Globals: nparams, nbaselines, u, v*/

int model_vis(double *params, double complex *modvis, double *logl, double *flux_frac)
{
    int status = 0;
    long i;
    double tempd;
    double l0, l1;


    for(i = 0; i < nbaselines; i++)
        {
            tempd = PI * params[1] / MAS_RAD * sqrt(u[i] * u[i] + v[i] * v[i]) + 1e-15;
            modvis[i] = 2 * params[0] * j1(tempd) / tempd;
        }

    /* Do custom likelihood */

// For vis0, allow equal probability between 0,1

    if(params[0] >= 0 && params[0] <= 1) l0 = 1.;
    if(params[0] < 0 || params[0] > 1) l0 = 1e-6;

// for mas, prefer sizes around 0.5, but with broad distribution
    l1 = exp(-0.5 * pow((params[1]) / 0.5, 2));
    if(params[1] < 0)  l1 = 1e-6;

    *logl = -1 * log(l0 * l1);

    *flux_frac = 1.0 - params[0] - params[2];


    return (status != 0);
}
