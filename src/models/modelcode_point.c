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

 model_UD: point source centered at (0,0)
 params: unresolved flux
 inputs: u,v, are in the wavelengths.
 logl is the a priori -log likelihood of the parameter combination.
------------------------------
To change this: ANY model is possible here, as long as model complex
visibilities can be returned based on a combinarion of up to MAX_PARAMS
(default 20 in macim.h) parameters.
*/

/* Globals: nparams, nbaselines, u, v*/

int model_vis(const double *params, double complex *modvis, double *logl, double *flux_frac)
{
    int status = 0;
    long i;
    double l0;

    for(i = 0; i < nbaselines; i++)
        {
            modvis[i] = 1.;
        }

// For flux allow equal probability between 0,1
    if(params[0] >= 0 && params[0] < 1.0) l0 = 1.;
    if(params[0] < 0 || params[0] >= 1.0) l0 = 1e-8;

    *logl = -1 * log(l0);
    *flux_frac = 1.0 - params[0];
    return (status != 0);
}
