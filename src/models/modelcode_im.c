/*   modelcode_im.c: An (incomplete) routine for adding an arbitrary fixed
     known image that takes up a certain fraction of the total flux.

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
 MJI test for a model where an arbitrary input .fits file is scales and
 added to the image.
------------------------------
To change this: ANY model is possible here, as long as model complex
visibilities can be returned based on a combinarion of up to MAX_PARAMS
(default 20 in macim.h) parameters.
*/

/* Globals: nparams, nbaselines, u, v*/

int first_run = 1;
double complex *unscaled_modvis;

int model_vis(double *params, double complex *modvis, double *logl, double *flux_frac)
{
    int status = 0;
    float *image;
    long i, j, k;
    double tempd, ftot;
    double l0, l1;
    /* Stuff for fits file output */
    fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */

    /* initialize FITS image parameters */
    char filename[80] = "output.fits";
    char mem_filename[80] = "MEM";
    char den_filename[80] = "DEN";
    long in_naxes[2];
    float nullval  = 0;
    int dummy_int;
    int nfound;

    if(first_run == 1)
        {
            unscaled_modvis = malloc(nbaselines * sizeof(double complex));
            if(fits_open_file(&fptr, "model_in.fits", READONLY, &status))
                printerror(status);
            /* read the NAXIS1 and NAXIS2 keyword to get image size */
            if(fits_read_keys_lng(fptr, "NAXIS", 1, 2, in_naxes, &nfound, &status))
                printerror(status);
            image = malloc((long)in_naxes[0] * (long)in_naxes[0] * sizeof(float));
            /* The image is now a float image, with flux that adds to 1.0 */
            if(fits_read_img(fptr, TFLOAT, 1, (long)in_naxes[0] * (long)in_naxes[0], &nullval,
                             image, &dummy_int, &status))
                printerror(status);
            ftot = 0.0;
            for(j = 0; j < in_naxes[0]; j++) for(i = 0; i < in_naxes[0]; i++) ftot += image[i + j * in_naxes[0]];
            printf("Total flux in model input image: %lf\n", ftot);
            for(j = 0; j < nbaselines; j++)
                {
                    unscaled_modvis[j] = 0;
                    for(k = 0; k < in_naxes[0]; k++) for(i = 0; i < in_naxes[0]; i++)
                            {
                                unscaled_modvis[j] += xtransform[i * nbaselines + j] *
                                                      ytransform[k * nbaselines + j] / ftot * image[i + k * in_naxes[0]];
                            }
                }
            printf("Filled image\n");
            free(image);
            image = NULL;
            first_run = 0;
        }
    for(j = 0; j < nbaselines; j++)
        modvis[j] = unscaled_modvis[j] * params[0];
    if(params[0] >= 0 && params[0] <= 1) l0 = 1.;
    if(params[0] < 0 || params[0] > 1) l0 = 1e-6;
    *logl = -1 * log(l0);
    *flux_frac = 1.0 - params[0];


    return (status != 0);
}
