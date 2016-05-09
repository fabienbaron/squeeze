/**
 * @file
 * @ingroup oitable
 * Data structure definitions and function prototypes for table-level
 * operations on OIFITS data.
 *
 * Copyright (C) 2007, 2015 John Young
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

#include <fitsio.h>
#include <complex.h>

typedef char BOOL;
typedef double DATA;
typedef int STATUS;


extern int oi_hush_errors; /**< If TRUE, don't report I/O errors to stderr */


/*
 * Data structures
 */
/* NB must allow for final null when dimensioning character arrays */

/** FITS primary HDU header keywords. */
typedef struct {
  /* mandatory keywords */
  char origin[FLEN_VALUE];
  char date_obs[FLEN_VALUE];
  char content[FLEN_VALUE];
  char telescop[FLEN_VALUE];
  char instrume[FLEN_VALUE];
  char observer[FLEN_VALUE];
  char insmode[FLEN_VALUE];
  char object[FLEN_VALUE];
  /* optional keywords */
  char referenc[FLEN_VALUE];
  char author[FLEN_VALUE];
  char prog_id[FLEN_VALUE];
  char procsoft[FLEN_VALUE];
  char obstech[FLEN_VALUE];
} oi_header;

/** Array element. Corresponds to one row of an OI_ARRAY FITS table. */
typedef struct {
  char tel_name[17];
  char sta_name[17];
  int sta_index;
  float diameter;
  double staxyz[3];
  double fov;
  char fovtype[7];
} element;

/** Data for OI_ARRAY FITS table */
typedef struct {
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
typedef struct {
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
  char category[4];
} target;

/** Data for OI_TARGET FITS table */
typedef struct {
  int revision;
  int ntarget;
  target *targ;
  BOOL usecategory;  /**< is target::category being used? */
} oi_target;

/** Data for OI_WAVELENGTH FITS table */
typedef struct {
  int revision;
  char insname[FLEN_VALUE];
  int nwave;
  float *eff_wave;
  float *eff_band;
} oi_wavelength;

/** Data for OI_CORR FITS table (new in OIFITS2) */
typedef struct {
  int revision;
  char corrname[FLEN_VALUE];
  int ndata;
  int ncorr;  /**< number of non-zero correlations */
  int *iindx;
  int *jindx;
  double *corr;
} oi_corr;

/** Polarization record. Corresponds to one row of an OI_INSPOL FITS table. */
typedef struct {
  int target_id;
  char insname[FLEN_VALUE];
  double mjd_obs;
  double mjd_end;
  float complex *lxx;
  float complex *lyy;
  float complex *lxy;
  float complex *lyx;
  int sta_index;
} oi_inspol_record;

/** Data for OI_INSPOL FITS table (new in OIFITS2) */
typedef struct {
  int revision;
  char date_obs[FLEN_VALUE];
  int npol;
  char arrname[FLEN_VALUE];
  char orient[FLEN_VALUE];
  char model[FLEN_VALUE];
  long numrec;
  int nwave;
  oi_inspol_record *record;
} oi_inspol;

/** Complex visibility record. Corresponds to one row of an OI_VIS
    FITS table. */
typedef struct {
  int target_id;
  double time;
  double mjd;
  double int_time;
  DATA *visamp, *visamperr;
  int corrindx_visamp;
  DATA *visphi, *visphierr;
  int corrindx_visphi;
  BOOL *visrefmap;
  DATA *rvis, *rviserr;
  int corrindx_rvis;
  DATA *ivis, *iviserr;
  int corrindx_ivis;
  double ucoord, vcoord;
  int sta_index[2];
  BOOL *flag;
} oi_vis_record;

/** Data for OI_VIS FITS table */
typedef struct {
  int revision;
  char date_obs[FLEN_VALUE];
  char arrname[FLEN_VALUE];  /**< empty string "" means not specified */
  char insname[FLEN_VALUE];
  char corrname[FLEN_VALUE]; /**< empty string "" means not specified */
  char amptyp[FLEN_VALUE];   /**< empty string "" means not specified */
  char phityp[FLEN_VALUE];   /**< empty string "" means not specified */
  int amporder;              /**< -1 means not specified */
  int phiorder;              /**< -1 means not specified */
  BOOL usevisrefmap;         /**< is oi_vis_record::visrefmap being used? */
  BOOL usecomplex;           /**< are oi_vis_record::rvis etc. being used? */
  char complexunit[FLEN_VALUE];  /**< TUNITn for RVIS/RVISERR/IVIS/IVISERR */
  char ampunit[FLEN_VALUE];  /**< TUNITn for VISAMP/VISAMPERR */
  long numrec;
  int nwave;
  oi_vis_record *record;
} oi_vis;

/** Visibility squared record. Corresponds to one row of an OI_VIS2
    FITS table. */
typedef struct {
  int target_id;
  double time;
  double mjd;
  double int_time;
  DATA *vis2data, *vis2err;
  int corrindx_vis2data;
  double ucoord, vcoord;
  int sta_index[2];
  BOOL *flag;
} oi_vis2_record;

/** Data for OI_VIS2 FITS table */
typedef struct {
  int revision;
  char date_obs[FLEN_VALUE];
  char arrname[FLEN_VALUE]; /**< empty string "" means not specified */
  char insname[FLEN_VALUE];
  char corrname[FLEN_VALUE]; /**< empty string "" means not specified */
  long numrec;
  int nwave;
  oi_vis2_record *record;
} oi_vis2;

/** Triple product record. Corresponds to one row of an OI_T3 FITS table. */
typedef struct {
  int target_id;
  double time;
  double mjd;
  double int_time;
  DATA *t3amp, *t3amperr;
  int corrindx_t3amp;
  DATA *t3phi, *t3phierr;
  int corrindx_t3phi;
  double u1coord, v1coord, u2coord, v2coord;
  int sta_index[3];
  BOOL *flag;
} oi_t3_record;

/** Data for OI_T3 FITS table */
typedef struct {
  int revision;
  char date_obs[FLEN_VALUE];
  char arrname[FLEN_VALUE]; /**< empty string "" means not specified */
  char insname[FLEN_VALUE];
  char corrname[FLEN_VALUE]; /**< empty string "" means not specified */
  long numrec;
  int nwave;
  oi_t3_record *record;
} oi_t3;

/** Flux record. Corresponds to one row of an OI_FLUX FITS table. */
typedef struct {
  int target_id;
  /* no TIME in this table */
  double mjd;
  double int_time;
  DATA *fluxdata;
  DATA *fluxerr;
  int corrindx_fluxdata;
  int sta_index;  /**< -1 means not specified */
  BOOL *flag;
} oi_flux_record;

/** Data for OI_FLUX FITS table (new in OIFITS2; formerly OI_SPECTRUM) */
typedef struct {
  int revision;
  char date_obs[FLEN_VALUE];
  char arrname[FLEN_VALUE]; /**< empty string "" means not specified */
  char insname[FLEN_VALUE];
  char corrname[FLEN_VALUE]; /**< empty string "" means not specified */
  double fov;
  char fovtype[FLEN_VALUE]; /**< empty string "" means not specified */
  char calstat;  /**< first character of FITS keyword */
  char fluxunit[FLEN_VALUE];  /**< TUNITn for FLUXDATA/FLUXERR */
  long numrec;
  int nwave;
  oi_flux_record *record;
} oi_flux;


/*
 * Function prototypes
 */

/* Functions from write_fits.c */
STATUS write_oi_header(fitsfile *fptr, oi_header header, STATUS *pStatus);
STATUS write_oi_array(fitsfile *fptr, oi_array array, int extver,
                      STATUS *pStatus);
STATUS write_oi_target(fitsfile *fptr, oi_target targets, STATUS *pStatus);
STATUS write_oi_wavelength(fitsfile *fptr, oi_wavelength wave, int extver,
                           STATUS *pStatus);
STATUS write_oi_corr(fitsfile *fptr, oi_corr corr, int extver,
                     STATUS *pStatus);
STATUS write_oi_inspol(fitsfile *fptr, oi_inspol inspol, int extver,
                       STATUS *pStatus);
STATUS write_oi_vis(fitsfile *fptr, oi_vis vis, int extver, STATUS *pStatus);
STATUS write_oi_vis2(fitsfile *fptr, oi_vis2 vis2, int extver,
                     STATUS *pStatus);
STATUS write_oi_t3(fitsfile *fptr, oi_t3 t3, int extver, STATUS *pStatus);
STATUS write_oi_flux(fitsfile *fptr, oi_flux flux, int extver,
                     STATUS *pStatus);
/* Functions from read_fits.c */
STATUS read_oi_header(fitsfile *fptr, oi_header *pHeader, STATUS *pStatus);
STATUS read_oi_target(fitsfile *fptr, oi_target *pTargets, STATUS *pStatus);
STATUS read_oi_array(fitsfile *fptr, char *arrname, oi_array *pArray,
                     STATUS *pStatus);
STATUS read_next_oi_array(fitsfile *fptr, oi_array *pArray, STATUS *pStatus);
STATUS read_oi_wavelength(fitsfile *fptr, char *insname, oi_wavelength *pWave,
                          STATUS *pStatus);
STATUS read_next_oi_wavelength(fitsfile *fptr, oi_wavelength *pWave,
                               STATUS *pStatus);
STATUS read_oi_corr(fitsfile *fptr, char *corrname, oi_corr *pCorr,
                    STATUS *pStatus);
STATUS read_next_oi_corr(fitsfile *fptr, oi_corr *pCorr, STATUS *pStatus);
STATUS read_next_oi_inspol(fitsfile *fptr, oi_inspol *pInspol, STATUS *pStatus);
STATUS read_next_oi_vis(fitsfile *fptr, oi_vis *pVis, STATUS *pStatus);
STATUS read_next_oi_vis2(fitsfile *fptr, oi_vis2 *pVis2, STATUS *pStatus);
STATUS read_next_oi_t3(fitsfile *fptr, oi_t3 *pT3, STATUS *pStatus);
STATUS read_next_oi_flux(fitsfile *fptr, oi_flux *pFlux, STATUS *pStatus);
/* Functions from alloc_fits.c */
void alloc_oi_array(oi_array *pArray, int nelement);
void alloc_oi_target(oi_target *pTargets, int ntarget);
void alloc_oi_wavelength(oi_wavelength *pWave, int nwave);
void alloc_oi_corr(oi_corr *pCorr, int ncorr);
void alloc_oi_inspol(oi_inspol *pInspol, long numrec, int nwave);
void alloc_oi_vis(oi_vis *pVis, long numrec, int nwave);
void alloc_oi_vis2(oi_vis2 *pVis2, long numrec, int nwave);
void alloc_oi_t3(oi_t3 *pT3, long numrec, int nwave);
void alloc_oi_flux(oi_flux *pFlux, long numrec, int nwave);
/* Functions from free_fits.c */
void free_oi_array(oi_array *pArray);
void free_oi_target(oi_target *pTargets);
void free_oi_wavelength(oi_wavelength *pWave);
void free_oi_corr(oi_corr *pCorr);
void free_oi_inspol(oi_inspol *pInspol);
void free_oi_vis(oi_vis *pVis);
void free_oi_vis2(oi_vis2 *pVis2);
void free_oi_t3(oi_t3 *pT3);
void free_oi_flux(oi_flux *pFlux);

#endif /* #ifndef EXCHANGE_H */

/** @} */
