#include <math.h>
#include <stdio.h>
/***************************************************************/
/* Change the centroid statistic                               */
/***************************************************************/
double cent_change(const int channel, double *centroid_image_x, double *centroid_image_y, const long new_x, const long new_y, const long old_x,
                   const long old_y, const unsigned short axis_len, const double fov, const double cent_mult) {
  double old_sumsqr, new_sumsqr;
  old_sumsqr = (centroid_image_x[channel] * centroid_image_x[channel] + centroid_image_y[channel] * centroid_image_y[channel]);
  centroid_image_x[channel] += (double)(new_x - old_x);
  centroid_image_y[channel] += (double)(new_y - old_y);
  const double centpos = 0.5 * ((double)axis_len); // rationale: choosing starting dirac has to fall on pixel for even number of pixels
  new_sumsqr = centroid_image_x[channel] * centroid_image_x[channel] + centroid_image_y[channel] * centroid_image_y[channel];
  return ((new_sumsqr - old_sumsqr) * cent_mult +
          (((double)new_x - centpos) * ((double)new_x - centpos) - ((double)old_x - centpos) * ((double)old_x - centpos) +
           ((double)new_y - centpos) * ((double)new_y - centpos) - ((double)old_y - centpos) * ((double)old_y - centpos)) *
              fov * fov) *
         4.0 / (double)(axis_len * axis_len);
}

double entropy(const double s) {
  if (s < 1e-7)
    return 0.0;
  else
    return lgamma(s);
}

double entropy_full(const double *x, const double *pr, const double eps, const int nx, const int ny, const double flux) {
  register int i;
  double reg = 0;
  for (i = 0; i < nx * ny; ++i) {
    if (x[i] > 0)
      reg += entropy(x[i]);
  }
  return reg;
}

double den_full(const double *x, const double *pr, const double eps, const int nx, const int ny, const double flux) {
  register int i, j;
  double reg = 2. * (nx + ny);
  for (i = 0; i < nx; ++i)
    for (j = 0; j < ny; ++j)
      reg += den_change(x, i, j, DEN_INIT, nx);
  return 0.5 * reg;
}

/*********************************************************/
/* A function that calculates the change in dark energy. */
/* DEN_ADD is 0, DEN_SUBTRACT is 1                       */
/*********************************************************/
double den_change(const double *image, const unsigned short i, const unsigned short j, const unsigned short direction, const unsigned short axis_len) {
  double delta_den = 0.0;
  double edge_amount = 1.0;

  const int pos = j * axis_len + i;

  if ((image[pos] == direction) || ((image[pos] == 0) && (direction == DEN_INIT))) {
    if (i == 0)
      delta_den += edge_amount;
    else if (image[pos - 1] == 0)
      delta_den += 1.0;

    if (j == 0)
      delta_den += edge_amount;
    else if (image[pos - axis_len] == 0)
      delta_den += 1.0;

    if (i == axis_len - 1)
      delta_den += edge_amount;
    else if (image[pos + 1] == 0)
      delta_den += 1.0;

    if (j == axis_len - 1)
      delta_den += edge_amount;
    else if (image[pos + axis_len] == 0)
      delta_den += 1.0;
  }
  return delta_den;
}

double transpec(const int nchan, const long imwidth, const double *image, const double flux) {
  long i, w;
  double temp1, temp2;

  temp1 = 0;
  for (i = 0; i < imwidth * imwidth; ++i) {
    temp2 = 0;
    for (w = 0; w < nchan; w++) {
      temp2 += image[w * imwidth * imwidth + i] * image[w * imwidth * imwidth + i];
    }

    temp1 += sqrt(temp2);
  }

  return temp1;
}

double tvsqspec(const int nchan, const long imwidth, const double *image, const double flux) {
  long i, w;
  double temp1, temp2;
  temp1 = 0;
  for (w = 1; w < nchan-1; w++)
  {
    for (i = 0; i < imwidth * imwidth; ++i)
    temp1 += (image[w * imwidth * imwidth + i] - image[(w-1) * imwidth * imwidth + i])*(image[w * imwidth * imwidth + i] - image[(w-1) * imwidth * imwidth + i]);
  }
  return temp1;
}

double tvsqspec_diff(long oldpos, long newpos, long chan, const int nchan, const long imsize, const double *image){
  double temp = 0.0;
  if (chan == nchan)
  {
    temp = -(image[chan * imsize + oldpos] - image[(chan-1) * imsize + oldpos])*(image[chan * imsize + oldpos] - image[(chan-1) * imsize + oldpos])+(image[chan * imsize + newpos] - image[(chan-1) * imsize + newpos])*(image[chan * imsize + newpos] - image[(chan-1) * imsize + newpos]);
    return temp;
  }

  if (chan == 0)
  {
    temp = -(image[imsize + oldpos] - image[oldpos])*(image[imsize + oldpos] - image[oldpos])+(image[imsize + newpos] - image[newpos])*(image[imsize + newpos] - image[newpos]);
    return temp;
  }

  // two sides
  temp = -(image[chan * imsize + oldpos] - image[(chan-1) * imsize + oldpos])*(image[chan * imsize + oldpos] - image[(chan-1) * imsize + oldpos])+(image[chan * imsize + newpos] - image[(chan-1) * imsize + newpos])*(image[chan * imsize + newpos] - image[(chan-1) * imsize + newpos])-(image[(chan+1) * imsize + oldpos] - image[chan * imsize + oldpos])*(image[(chan+1) * imsize + oldpos] - image[chan * imsize + oldpos])+(image[(chan+1) * imsize + newpos] - image[chan * imsize + newpos])*(image[(chan+1) * imsize + newpos] - image[chan * imsize + newpos]);
  return temp;
}




double transpec_diffpoint(long pos, long chan, double diff, const int nchan, const long imsize, const double *image) {
  long w;
  double temp = 0;
  for (w = 0; w < nchan; w++) {
    if (w == chan)
      temp += (image[w * imsize + pos] + diff) * (image[w * imsize + pos] + diff);
    else
      temp += (image[w * imsize + pos]) * (image[w * imsize + pos]);
  }
  return sqrt(temp);
}

double TV(const double *x, const double *pr, const double eps, const int nx, const int ny, const double flux) {
  // This regularizer is the TOTVAR (p=1.0) on the local gradient
  // the gradient is evaluated by backward difference, e is the threshold
  // by definition, L1g(x) = (sum( |grad im| ))
  // Note: all image borders are ignored

  // also ignore if pixel borders an area of prior with logp < -13.8
  // corresponding to a valuei n the prior image of about 1-e6

  register int i, j, off;
  double dx, dy, pixreg;
  //  double pr_threshold=-13.8;

  double L1g = -sqrt(eps) * (double)(nx * ny);

  // pr_threshold=-1e9; // turn off.
  //  printf("Max available here is : %d %d", omp_get_max_threads(), omp_get_num_threads());
  // printf("nx ny %i %i\n",nx,ny);
  // Compute the norm of the local image gradient on each point
  //#pragma omp for
  for (j = 1; j < ny; ++j) {
    off = nx * j;
    for (i = 1; i < nx; ++i) {
      if (i > 0) {
        dx = x[i + off] - x[i - 1 + off];
        // if ( (pr[i+off] <=pr_threshold) || (pr[i-1+off] <=pr_threshold)) dx=0.0;
      } else {
        dx = x[1 + off] - x[off];
        //  if ( (pr[1+off] <=pr_threshold) || (pr[off] <=pr_threshold)) dx=0.0;
      }

      if (j > 0) {
        dy = x[i + off] - x[i + off - nx];
        // if ( (pr[i+off] <=pr_threshold) || (pr[i+off-nx] <=pr_threshold)) dy=0.0;
      } else {
        dy = x[i + nx] - x[i];
        //  if ( (pr[i+nx] <=pr_threshold) || (pr[i] <=pr_threshold)) dy=0.0;
      }

      pixreg = sqrt(dx * dx + dy * dy + eps * eps);

      L1g += pixreg;
    }
  }
  return L1g / flux;
}

double UDreg(const double *x, const double *pr, const double eps, const int nx, const int ny, const double flux) {
  // This regularizer is the Lp norm (p=0.5) on the local gradient
  // It is somehow similar to the total variation, which is the L1 norm
  // the gradient is evaluated by backward difference, e is the threshold
  // by definition, L05g(x) = (sum( |grad im|^.5 ))^2
  // Note: all image borders are ignored

  // also ignore if pixel borders an area of prior with logp < -13.8
  // corresponding to a value in the prior image of about 1-e6

  register int i, j, off;
  double dx, dy, pixreg;
  //  double pr_threshold=-13.8;

  double L05g = -sqrt(eps) * (double)(nx * ny);

  // pr_threshold=-1e9; // turn off.

  // printf("nx ny %i %i\n",nx,ny);
  // Compute the norm of the local image gradient on each point
  for (j = 1; j < ny; ++j) {
    off = nx * j;
    for (i = 1; i < nx; ++i) {
      if (i > 0) {
        dx = x[i + off] - x[i - 1 + off];
      } else {
        dx = x[1 + off] - x[off];
      }

      if (j > 0) {
        dy = x[i + off] - x[i + off - nx];
      } else {
        dy = x[i + nx] - x[i];
      }
      pixreg = sqrt(sqrt(dx * dx + dy * dy + eps * eps));
      //  printf("i j dx dy x %i %i %5.2f %5.2f %5.2f \n",i,j,dx,dy,x[i*nx+j]);
      L05g += pixreg;
    }
  }
  //   return L05g;

  return L05g * L05g / flux;
}

double L0(const double *x, const double *pr, const double eps, const int nx, const int ny, const double flux) {
  register int i;
  double L0l = 0;
  for (i = 0; i < nx * ny; ++i) {
    if (x[i] > 0.0)
      L0l += 1.;
  }
  return L0l;
}

double entropy_xlogx(const double *x, const double *pr, const double eps, const int nx, const int ny, const double flux) {
  register int i;
  double ent = 0;
  for (i = 0; i < nx * ny; ++i) {
    if (x[i] > 1e-10)
      ent += x[i] * log(x[i]) - x[i];
  }
  return ent;
}

double L0_diff(const double *x, const long new_pos, const long old_pos) {
  // difference due to old_pos now having one less element: x[old_pos]==1 -> l0+= -1 else 0
  // difference due to new_pos now having one more element: x[new_pos]==0 -> l0+= 1 else 0
  return (x[new_pos] < 1) ? 1 : 0 - ((x[old_pos] > 0) && (x[old_pos] < 2)) ? 1 : 0;
}

double L1(const double *x, const double *pr, const double eps, const int nx, const int ny, const double flux) {
  register int i;
  double L1l = 0;
  for (i = 0; i < nx * ny; ++i) {
    L1l += fabs(x[i]);
  }
  return L1l / flux;
}

// If you're sure flux >= 0 this will be faster
double L1pos(const double *x, const double *pr, const double eps, const int nx, const int ny, const double flux) {
  register int i;
  double L1l = 0;
  for (i = 0; i < nx * ny; ++i) {
    L1l += x[i];
  }
  return L1l / flux;
}

double L2(const double *x, const double *pr, const double eps, const int nx, const int ny, const double flux) {
  register int i;
  double L2l = 0;
  for (i = 0; i < nx * ny; ++i) {
    L2l += x[i] * x[i];
  }
  return sqrt(L2l) / flux;
}

double L2sq(const double *x, const double *pr, const double eps, const int nx, const int ny, const double flux) {
  register int i;
  double L2l = 0;
  for (i = 0; i < nx * ny; ++i) {
    L2l += x[i] * x[i];
  }
  return L2l / flux;
}

double LAP(const double *x, const double *pr, const double eps, const int nx, const int ny, const double flux) {
  // L1 norm of the Laplacian
  // Contrary to totvar, allows slope (good for stellar contours)
  // Laplacian via 5 point stencil  L = x[ i - 1, j] + x[i+1, j] + x[i, j-1] + x[i, j+1] - 4 * x[i,j]

  register int i, j, off;
  double reg = 0;

  // for boundaries, we assume zeroes outside of the image
  for (j = 1; j < ny - 1; ++j) {
    off = nx * j;
    for (i = 1; i < nx - 1; ++i)
      reg += fabs(x[i - 1 + off] + x[i + 1 + off] + x[i + off - nx] + x[i + off + nx] - 4. * x[i + off]);
    // case i = 0
    reg += fabs(x[1 + off] + x[off - nx] + x[off + nx] - 4. * x[off]);
    // case i = nx -1
    reg += fabs(x[nx - 2 + off] + x[nx - 1 + off - nx] + x[nx - 1 + off + nx] - 4. * x[nx - 1 + off]);
  }

  for (i = 1; i < nx - 1; ++i) {
    // case j = 0
    off = 0;
    reg += fabs(x[i - 1] + x[i + 1] + x[i + off + nx] - 4. * x[i + off]);
    // case j = nx -1
    off = nx * (nx - 1);
    reg += fabs(x[i - 1 + off] + x[i + 1 + off] + x[i + off - nx] - 4. * x[i + off]);
  }

  return reg / flux;
}

double reg_prior_image(const double *x, const double *pr, const double eps, const int nx, const int ny, const double flux) {
  register int i;
  double rpi = 0;

  for (i = 0; i < nx * ny; ++i)
    if (x[i] > 0)
      rpi += pr[i];

  // For ref, initial methods was:
  //      for (i = 0; i < nelements; ++i)
  //    reg_value[w * NREGULS + REG_PRIORIMAGE] += prior_image[element_y[w * nelements + i] * axis_len + element_x[w * nelements + i]];

  return rpi;
}

//
// Wavelets
//

// CDF (5,3), used in JPEG2000 -- Note: this can be sped up a lot in the future
void fwt53(double *wav, const double *x, const int nx, const int ny) {
  const double a0 = -0.5;
  const double a1 = 0.25;
  const double s0 = sqrt(2.0);
  const double s1 = sqrt(2.0) * 0.5;
  register int i, j, off;

  double *tempx = malloc(nx * ny * sizeof(double));
  memcpy(tempx, x, nx * ny * sizeof(double));

  for (j = 0; j < ny; ++j) {
    off = nx * j;
    // Predict 1
    for (i = 1; i < nx - 1; i += 2) {
      tempx[off + i] += a0 * (tempx[off + i - 1] + tempx[off + i + 1]);
    }
    tempx[off + nx - 1] += 2 * a0 * tempx[off + nx - 2];

    // Update 1
    for (i = 2; i < nx; i += 2) {
      tempx[off + i] += a1 * (tempx[off + i - 1] + tempx[off + i + 1]);
    }
    tempx[off] += 2 * a1 * tempx[off + 1];
  }

  for (j = 0; j < ny; ++j) {
    off = nx * j;
    for (i = 0; i < nx; ++i) {
      if (i % 2 == 0)
        wav[(i / 2) * ny + j] = s0 * tempx[off + i];
      else
        wav[(nx / 2 + i / 2) * ny + j] = s1 * tempx[off + i];
    }
  }
  free(tempx);
}

// CDF (9,7), used in JPEG2000 -- Note: this can be sped up a lot in the future
void fwt97(double *wav, const double *x, const int nx, const int ny) {
  const double a0 = -1.586134342;
  const double a1 = -0.05298011854;
  const double a2 = 0.8829110762;
  const double a3 = 0.4435068522;
  const double s0 = 0.81289306611596146; // 1/1.230174104914
  const double s1 = 0.61508705245700002; // 0.5 * 1.230174104914
  register int i, j;
  int off;
  double *tempx = malloc(nx * ny * sizeof(double));
  memcpy(tempx, x, nx * ny * sizeof(double));

  for (j = 0; j < ny; ++j) {
    off = nx * j;
    // Predict 1
    for (i = 1; i < nx - 1; i += 2) {
      tempx[i + off] += a0 * (tempx[off + i - 1] + tempx[off + i + 1]);
    }
    tempx[off + nx - 1] += 2 * a0 * tempx[off + nx - 2];

    // Update 1
    for (i = 2; i < nx; i += 2) {
      tempx[off + i] += a1 * (tempx[off + i - 1] + tempx[off + i + 1]);
    }
    tempx[off] += 2 * a1 * tempx[off + 1];

    // Predict 2
    for (i = 1; i < nx - 1; i += 2) {
      tempx[off + i] += a2 * (tempx[off + i - 1] + tempx[off + i + 1]);
    }
    tempx[off + nx - 1] += 2 * a2 * tempx[off + nx - 2];

    // Update 2
    for (i = 2; i < nx; i += 2) {
      tempx[off + i] += a3 * (tempx[off + i - 1] + tempx[off + i + 1]);
    }
    tempx[off] += 2 * a3 * tempx[off + 1];
  }

  // Deinterlace, transpose and scale
  for (j = 0; j < ny; ++j) {
    off = nx * j;
    for (i = 0; i < nx; ++i) {
      if (i % 2 == 0)
        wav[(i / 2) * ny + j] = s0 * tempx[off + i];
      else
        wav[(nx / 2 + i / 2) * ny + j] = s1 * tempx[off + i];
    }
  }
  free(tempx);
}

void fwt97_2D(double *wav, const double *x, const int nx, const int ny, const int levels) {
  int i;
  for (i = 0; i < levels; ++i) {
    fwt97(wav, x, nx, ny);   // do on rows
    fwt97(wav, wav, nx, ny); // do on cols using the result
  }
}

void fwt53_2D(double *wav, const double *x, const int nx, const int ny, const int levels) {
  int i;
  for (i = 0; i < levels; ++i) {
    fwt53(wav, x, nx, ny);   // do on rows
    fwt53(wav, wav, nx, ny); // do on cols using the previous result
  }
}

double L0_CDF97(const double *x, const double *pr, const double eps, const int nx, const int ny, const double flux) {
  double *wav = malloc(nx * ny * sizeof(double));
  fwt97_2D(wav, x, nx, ny, 1);
  for (int i = 0; i < nx * ny; ++i)
    wav[i] = fabs(wav[i]);
  double reg = L0(wav, NULL, 0, nx, ny, 1.);
  free(wav);
  return reg;
}

double L0_CDF53(const double *x, const double *pr, const double eps, const int nx, const int ny, const double flux) {
  double *wav = malloc(nx * ny * sizeof(double));
  fwt53_2D(wav, x, nx, ny, 1);
  for (int i = 0; i < nx * ny; ++i)
    wav[i] = fabs(wav[i]);
  double reg = L0(wav, NULL, 0, nx, ny, 1.);
  free(wav);
  return reg;
}

double L1_CDF53(const double *x, const double *pr, const double eps, const int nx, const int ny, const double flux) {
  double *wav = malloc(nx * ny * sizeof(double));
  fwt53_2D(wav, x, nx, ny, 1);
  double reg = L1(wav, NULL, 0, nx, ny, 1.);
  free(wav);
  return reg;
}

double L1_CDF97(const double *x, const double *pr, const double eps, const int nx, const int ny, const double flux) {
  double *wav = malloc(nx * ny * sizeof(double));
  fwt97_2D(wav, x, nx, ny, 1);
  double reg = L1(wav, NULL, 0, nx, ny, 1.);
  free(wav);
  return reg;
}

void atrous_set(int idx) // set the wavelet coefficients for the a trous algorithm
{
  static float d1_lin3[3] = {0.25, 0.5, 0.25};
  static float d1_spl5[5] = {0.0625, 0.25, 0.375, 0.25, 0.0625};
  static float d2_lin3[3 * 3];
  static float d2_spl5[5 * 5];
  int i, j;

  switch (idx) {
  case 1: // linear interpolation filter
    atrous_1d_filter.ncof = 3;
    atrous_1d_filter.cc = d1_lin3;
    atrous_2d_filter.ncof = 3;
    atrous_2d_filter.cc = d2_lin3;
    break;
  case 2: // B_3-spline interpolation
    atrous_1d_filter.ncof = 5;
    atrous_1d_filter.cc = d1_spl5;
    atrous_2d_filter.ncof = 5;
    atrous_2d_filter.cc = d2_spl5;
    break;
  default:
    printf("unknown value in my_at_set\n");
    break;
  }

  // compute the two dimensional convolution mask
  for (i = 0; i < atrous_2d_filter.ncof; i++) {
    for (j = 0; j < atrous_2d_filter.ncof; j++) {
      atrous_2d_filter.cc[i * atrous_2d_filter.ncof + j] = atrous_1d_filter.cc[i] * atrous_1d_filter.cc[j];
      atrous_1d_filter.ioff = atrous_1d_filter.ncof / 2 + atrous_1d_filter.ncof % 2 - 1;
      atrous_2d_filter.ioff = atrous_2d_filter.ncof / 2 + atrous_2d_filter.ncof % 2 - 1;
    }
  }
}
int mpower(int basis, int exponent) {
  int result = 1;
  int expo;
  expo = exponent;
  while (expo-- > 0) {
    result *= basis;
  }
  return result;
}

void atrous_fwd(const double *x, double *wav, const int nx, const int ny, const int nscales) {
  int i, j, i2, j2, i3, j3, k, nf, nc, istep;
  nf = atrous_2d_filter.ncof;
  nc = atrous_2d_filter.ioff;
  int n = nx * ny;

  // forward transform
  float *storage = malloc(sizeof(float) * n);
  for (i = 0; i < n; i++) // copy initial image - set as first scale
    storage[i] = x[i];

  for (k = nscales - 1; k > 0; k--) { // determine step size for convolution depending on resolution scale

    istep = powf(2, nscales - 1 - k); // copy image to scale, starting with large k (small scales)

    for (i = 0; i < n; i++)
      wav[k * n + i] = x[i];

    // smooth image by convolution with filter
    for (i = 0; i < nx; i++) {
      for (j = 0; j < ny; j++) {
        storage[j * nx + i] = 0.;
        for (i2 = 0; i2 < nf; i2++) {
          i3 = i + (i2 - nc) * istep;
          // wrap around edges, using periodic boundary conditions
          i3 = (i3 >= 0 ? (i3 < nx ? i3 : i3 - nx) : i3 + nx);
          for (j2 = 0; j2 < nf; j2++) {
            j3 = j + (j2 - nc) * istep;
            // wrap around edges, using periodic boundary conditions
            j3 = (j3 >= 0 ? (j3 < ny ? j3 : j3 - ny) : j3 + ny);
            storage[j * nx + i] += atrous_2d_filter.cc[i2 * nf + j2] * wav[k * n + j3 * nx + i3];
          }
        }
      }
    }
    for (i = 0; i < n; i++) // construct detail coefficients by subtraction of smoothed image
      wav[k * n + i] -= storage[i];

    for (i = 0; i < n; i++)
      wav[i] = storage[i];
  }
  free(storage);
}

/* two dimensional a trous transform */
int a_trous(double **a, int nx, int ny, int nscales, int isign) {
  int i, j, i2, j2, i3, j3, k, nf, nc, istep;
  double *image, *image1;
  double **b;
  // char line[80];

  nf = atrous_2d_filter.ncof;
  nc = atrous_2d_filter.ioff;

  // printf("Starting a trous algorithm\n");
  // printf("n:%i, c:%i, nx:%i, ny:%i\n", nf, nc, nx, ny);

  if (isign > 0) {
    /* forward transform */
    image = malloc(sizeof(double) * nx * ny);
    for (i = 0; i < nx * ny; i++) /* copy image */
      image[i] = a[0][i];

    for (k = nscales - 1; k > 0; k--) { /* determine step size for convolution depending on resolution scale */

      istep = mpower(2, nscales - 1 - k); /* copy image to scale, starting with large k (small scales) */

      for (i = 0; i < nx * ny; i++)
        a[k][i] = image[i];

      /* smooth image by convolution with filter */
      for (i = 0; i < nx; i++) {
        for (j = 0; j < ny; j++) {
          image[j * nx + i] = 0.;
          for (i2 = 0; i2 < nf; i2++) {
            i3 = i + (i2 - nc) * istep;
            /* wrap around edges, using periodic boundary conditions */
            i3 = (i3 >= 0 ? (i3 < nx ? i3 : i3 - nx) : i3 + nx);
            for (j2 = 0; j2 < nf; j2++) {
              j3 = j + (j2 - nc) * istep;
              /* wrap around edges, using periodic boundary conditions */
              j3 = (j3 >= 0 ? (j3 < ny ? j3 : j3 - ny) : j3 + ny);
              image[j * nx + i] += atrous_2d_filter.cc[i2 * nf + j2] * a[k][j3 * nx + i3];
            }
          }
        }
      }
      for (i = 0; i < nx * ny; i++) /* construct detail coefficients by subtraction of smoothed image */
        a[k][i] -= image[i];

      for (i = 0; i < nx * ny; i++)
        a[0][i] = image[i];
    }
    free(image);
  } else {
    if (isign == -1) {
      /* inverse transform */
      /* the inverse transform of the a trous algorithm is a simple sum
       of the smoothed image plus the detail images at all scales */
      for (k = 1; k < nscales; k++) {
        for (i = 0; i < nx * ny; i++) {
          a[0][i] += a[k][i];
          a[k][i] = 0.;
        }
      }
    } else {
      /* Transposed transform */
      b = malloc(sizeof(double *) * nscales);
      image1 = malloc(sizeof(double) * nx * ny * nscales);
      for (k = 0; k < nscales; k++)
        b[k] = image1 + k * nx * ny;
      for (k = nscales - 1; k > 0; k--) {
        for (i = 0; i < nx * ny; i++)
          b[0][i] = a[k][i];
        a_trous(b, nx, ny, nscales + 1 - k, 1);
        for (i = 0; i < nx * ny; i++)
          a[k][i] = b[1][i];
      }
      for (i = 0; i < nx * ny; i++)
        b[0][i] = a[0][i];
      a_trous(b, nx, ny, nscales, 1);
      for (i = 0; i < nx * ny; i++)
        a[0][i] = b[0][i];
      for (k = 1; k < nscales; k++) {
        for (i = 0; i < nx * ny; i++) {
          a[0][i] += a[k][i];
          a[k][i] = 0.;
        }
      }
      free(b);
      free(image1);
    }
  }
  return 0;
}

double L0_ATROUS(const double *x, const double *pr, const double eps, const int nx, const int ny, const double flux) {
  const int nscales = 4;
  double *wav = malloc(nscales * nx * ny * sizeof(double));
  atrous_fwd(x, wav, nx, ny, nscales);
  // Take modulus for L1(pos) or L0
  for (int i = 0; i < nscales * nx * ny; ++i)
    wav[i] = fabs(wav[i]);
  double reg = L0(wav, NULL, 0, nscales * nx, ny, 1.);
  free(wav);
  return reg;
}

double L1_ATROUS(const double *x, const double *pr, const double eps, const int nx, const int ny, const double flux) {
  const int nscales = 4;
  double *wav = malloc(nscales * nx * ny * sizeof(double));
  atrous_fwd(x, wav, nx, ny, nscales);
  double reg = entropy_xlogx(wav, NULL, 0, nscales * nx, ny, 1.);
  free(wav);
  return reg;
}
