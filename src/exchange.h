/* $Id: exchange.h,v 1.5 2007-11-01 12:56:16 jsy1001 Exp $ */

/**
 * @file exchange.h
 * @ingroup oitable
 *
 * Data structure definitions and function prototypes for table-level
 * operations on OIFITS data.
 *
 * Copyright (C) 2007 John Young
 *
 *
 * This file is part of OIFITSlib.
 *
 * OIFITSlib is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * OIFITSlib is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with OIFITSlib.  If not, see
 * http://www.gnu.org/licenses/
 */

/**
 * @mainpage
 *
 * OIFITSlib is a C library for input/output, merging, filtering and
 * checking of optical/IR interferometry datasets in the OIFITS
 * exchange format.
 *
 * OIFITS is a standard for exchanging calibrated, time-averaged data
 * from astronomical optical interferometers, based on the FITS
 * Standard. OIFITS may be used to combine data from multiple
 * interferometer arrays for joint analysis and/or image
 * reconstruction.
 *
 * The OIFITS standard is described in Pauls et al. (2005) PASP 117,
 * 1255. A PDF reprint of this paper and other OIFITS resources are
 * available from the OIFITS website
 * http://www.mrao.cam.ac.uk/~jsy1001/exchange/
 *
 * OIFITSlib incorporates the following modules:
 * - @ref oitable
 * - @ref oifile
 * - @ref oimerge
 * - @ref oifilter
 * - @ref oicheck
 */

/**
 * @defgroup oitable  Table-level OIFITS I/O
 *
 * This module is derived from the "OIFITS example software in C". It
 * provides the same Application Programming Interface (API) as the
 * example software, with the addition of a set of free_oi_*()
 * functions.
 *
 * A higher-level API, containing functions to read and write an
 * entire OIFITS file, is provided by the @ref oifile module.
 *
 * @{
 */

#ifndef EXCHANGE_H
#define EXCHANGE_H

#include "../lib/cfitsio/fitsio.h"

typedef char BOOL;
typedef double DATA;
typedef int STATUS;


extern int oi_hush_errors; /**< If TRUE, don't report I/O errors to stderr */


/*
 * Data structures
 */
/* NB must allow for final null when dimensioning character arrays */

/** Array element. Corresponds to one row of an OI_ARRAY FITS table. */
typedef struct
{
  char tel_name[17];
  char sta_name[17];
  int sta_index;
  float diameter;
  double staxyz[3];
} element;

/** Data for OI_ARRAY FITS table */
typedef struct
{
  int revision;
  char arrname[FLEN_VALUE];
  char frame[FLEN_VALUE];
  double arrayx, arrayy, arrayz;
  int nelement;
  element *elem;
} oi_array;

/** Info on an observing target.
 *
 * Corresponds to one row of an OI_TARGET FITS table.
 */
typedef struct
{
  int target_id;
  char target[17];
  double raep0, decep0;
  float equinox;
  double ra_err, dec_err;
  double sysvel;
  char veltyp[9], veldef[9];
  double pmra, pmdec;
  double pmra_err, pmdec_err;
  float parallax, para_err;
  char spectyp[17];
} target;

/** Data for OI_TARGET FITS table */
typedef struct
{
  int revision;
  int ntarget;
  target *targ;
} oi_target;

/** Data for OI_WAVELENGTH FITS table */
typedef struct
{
  int revision;
  char insname[FLEN_VALUE];
  int nwave;
  float *eff_wave;
  float *eff_band;
} oi_wavelength;

/** Complex visibility record. Corresponds to one row of an OI_VIS
    FITS table. */
typedef struct
{
  int target_id;
  double time;
  double mjd;
  double int_time;
  DATA *visamp, *visamperr;
  DATA *visphi, *visphierr;
  double ucoord, vcoord;
  int sta_index[2];
  BOOL *flag;
} oi_vis_record;

/** Data for OI_VIS FITS table */
typedef struct
{
  int revision;
  char date_obs[FLEN_VALUE];
  char arrname[FLEN_VALUE]; /**< empty string "" means not specified */
  char insname[FLEN_VALUE];
  long numrec;
  int nwave;
  oi_vis_record *record;
} oi_vis;

/** Visibility squared record. Corresponds to one row of an OI_VIS2
    FITS table. */
typedef struct
{
  int target_id;
  double time;
  double mjd;
  double int_time;
  DATA *vis2data, *vis2err;
  double ucoord, vcoord;
  int sta_index[2];
  BOOL *flag;
} oi_vis2_record;

/** Data for OI_VIS2 FITS table */
typedef struct
{
  int revision;
  char date_obs[FLEN_VALUE];
  char arrname[FLEN_VALUE]; /**< empty string "" means not specified */
  char insname[FLEN_VALUE];
  long numrec;
  int nwave;
  oi_vis2_record *record;
} oi_vis2;

/** Triple product record. Corresponds to one row of an OI_T3 FITS table. */
typedef struct
{
  int target_id;
  double time;
  double mjd;
  double int_time;
  DATA *t3amp, *t3amperr;
  DATA *t3phi, *t3phierr;
  double u1coord, v1coord, u2coord, v2coord;
  int sta_index[3];
  BOOL *flag;
} oi_t3_record;

/** Data for OI_T3 FITS table */
typedef struct
{
  int revision;
  char date_obs[FLEN_VALUE];
  char arrname[FLEN_VALUE]; /**< empty string "" means not specified */
  char insname[FLEN_VALUE];
  long numrec;
  int nwave;
  oi_t3_record *record;
} oi_t3;


/*
 * Function prototypes
 */

/* Functions from write_fits.c */
STATUS write_oi_array(fitsfile *fptr, oi_array array, int extver,
                      STATUS *pStatus);
STATUS write_oi_target(fitsfile *fptr, oi_target targets, STATUS *pStatus);
STATUS write_oi_wavelength(fitsfile *fptr, oi_wavelength wave, int extver,
                           STATUS *pStatus);
STATUS write_oi_vis(fitsfile *fptr, oi_vis vis, int extver, STATUS *pStatus);
STATUS write_oi_vis2(fitsfile *fptr, oi_vis2 vis2, int extver,
                     STATUS *pStatus);
STATUS write_oi_t3(fitsfile *fptr, oi_t3 t3, int extver, STATUS *pStatus);
/* Functions from read_fits.c */
STATUS read_oi_target(fitsfile *fptr, oi_target *pTargets, STATUS *pStatus);
STATUS read_oi_array(fitsfile *fptr, char *arrname, oi_array *pArray,
                     STATUS *pStatus);
STATUS read_next_oi_array(fitsfile *fptr, oi_array *pArray, STATUS *pStatus);
STATUS read_oi_wavelength(fitsfile *fptr, char *insname, oi_wavelength *pWave,
                          STATUS *pStatus);
STATUS read_next_oi_wavelength(fitsfile *fptr, oi_wavelength *pWave,
                               STATUS *pStatus);
STATUS read_next_oi_vis(fitsfile *fptr, oi_vis *pVis, STATUS *pStatus);
STATUS read_next_oi_vis2(fitsfile *fptr, oi_vis2 *pVis2, STATUS *pStatus);
STATUS read_next_oi_t3(fitsfile *fptr, oi_t3 *pT3, STATUS *pStatus);
/* Functions from free_fits.c */
void free_oi_array(oi_array *pArray);
void free_oi_target(oi_target *pTargets);
void free_oi_wavelength(oi_wavelength *pWave);
void free_oi_vis(oi_vis *pVis);
void free_oi_vis2(oi_vis2 *pVis2);
void free_oi_t3(oi_t3 *pT3);

#endif /* #ifndef EXCHANGE_H */

/** @} */
