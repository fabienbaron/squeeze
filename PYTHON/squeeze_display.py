#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 28 11:29:42 2017

Procedure to monitor the output of Squeeze with Python

@author: montarges
"""

#pylint: disable-msg=E1101,C0103

import os
import argparse
import time
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

#Retrieving band via argument
parser = argparse.ArgumentParser()
parser.add_argument("--d", dest="directory", help="Directory to monitor", required=False)
parser.add_argument("--p", dest="power", help="Exponent for image display", required=False, \
                    type=float)
args = parser.parse_args()
dir_mon = args.directory
power = args.power

#If no directory given, default to parent one
if dir_mon is None:
    dir_mon = os.path.pardir
print "\nMonitoring chain00.fits in "+dir_mon+"\n"

#Rgulation methods
reg_names = np.array(['PARAM', 'CENT', 'IMPRIOR', 'ENT', 'DEN', 'TV', 'SPOT', 'LAP', 'L0', 'TRANSPEC'])
colsty = np.array(["b", "g", "r", "c", "m", "orange", "peru", "thistle", "lightslategray"])
linesty = np.array(["--", "-.", ":", ".", "v", "^", "s", "+"])

#Init
ltime = datetime.now()
mtime = ltime
chi2 = np.array([-1])

loop1 = True
nofig = True

while loop1:

    loop2 = True

    #Waiting for update of chain00.fits
    while loop2:

        #First, does the file exists ?
        try:
            mtime = os.path.getmtime(os.path.join(dir_mon, "chain00.fits"))
            mtime = datetime.fromtimestamp(mtime)
        except OSError:
            print "\nTrouble opening file (ignore if occasional message only)"

        #Update or not ?
        if mtime > ltime:
            break
        else:
            time.sleep(0.3)

    #Update of the monitored fits, let's go !
    ltime = mtime
    try:
        with (fits.open(os.path.join(dir_mon, "chain00.fits"))) as hdu:
            hdr = hdu[0].header
            img = hdu[0].data

        #Exponent scaling ?
        if power is not None:
            img = img**power

        #Multi wavelength ?
        if len(img.shape) == 2:
            nchanr = 1
        elif len(img.shape) > 2:
            nchanr = img.shape[0]
            if nchanr > 10:
                nchanr = 5

        #Creating figure if not already there + regulation parameters
        if nofig:

            fig, all_axes = plt.subplots(1, nchanr)
            fig.set_size_inches(4*nchanr, 4)
            nofig = False
            cbar = [None]*nchanr

            reg_params = np.zeros(len(reg_names))
            reg_vals = np.zeros([nchanr, len(reg_names)])
            reg_active = np.zeros(len(reg_names))

        #Header information
        ndf = hdr['NDF']
        newchi2 = hdr['CHI2']
        scale = hdr['scale']
        nx = hdr['NAXIS1']
        ny = hdr['NAXIS2']
        fovx = nx*scale
        fovy = ny*scale

        #Plotting current image
        time_str = mtime.strftime("%a. %d %b %Y, %H:%M:%S")
        fig.suptitle("Current/Mean image\n"+time_str)
        for n, ax in enumerate(all_axes):
            plot_im = ax.matshow(img[n, :, :], cmap='hot', \
            extent=[fovx/2., -fovx/2., -fovy/2., fovy/2.], origin='lower')
            ax.set_title("Channel "+str(n+1)+"\n")
            ax.set_xlabel("Rel. R.A. (mas)")
        all_axes[0].set_ylabel("Rel. Dec. (mas)")

        #Getting regularisation output
        for i in range(len(reg_names)):
            reg_params[i] = hdr['HYPER'+str(i)]
            if reg_params[i] > 0:
                for j in range(nchanr):
                    reg_vals[j, i] = hdr['REG'+str(i)+'W'+str(j)]

        #Detmining which regularizers are active
        reg_active = (reg_params>0)
        valid = False
        if reg_active.any():
            valid = True
            newregs = np.zeros([nchanr, reg_active.sum()])
            for i in range (nchanr):
                newregs[i, :] = reg_vals[i, reg_active]*reg_params[reg_active]/(2.*ndf)

        #Concatenating data + initializing figure at 1st iter
        if chi2[0] == -1:
            chi2 = np.array([newchi2])
            regs = newregs
            regs = np.reshape(newregs, [nchanr, reg_active.sum(), 1])
            fig_reg, ax_reg = plt.subplots(1, 1)
            ax_reg.set_yscale('log')
            ax_reg.set_xlabel('Trial')
            ax_reg.set_ylabel(r'$\chi^2$ and Reguls/(2*ndf)')
            ax_reg.grid()
            linesty = linesty[reg_active]
            colsty = colsty[reg_active]
            reg_names = reg_names[reg_active]
            first_plot = True
        else:
            fig_reg.suptitle(time_str)
            chi2 = np.append(chi2, newchi2)
            newregs = np.reshape(newregs, [nchanr, reg_active.sum(), 1])
            regs = np.append(regs, newregs, axis=2)

        #Plot
        ax_reg.plot(range(chi2.size), chi2, 'k', label=r'$\chi^2$')
        if valid:
            for i in range(reg_active.sum()):
                    for j in range(nchanr):
                        if j == 0:
                            ax_reg.plot(range(chi2.size), regs[j, i, :], linesty[i], \
                            color=colsty[i], label=reg_names[i])
                        else:
                            ax_reg.plot(range(chi2.size), regs[j, i, :], linesty[i], \
                            color=colsty[i])

        #Legend
        if first_plot:
            ax_reg.legend(loc=1)
            first_plot = False

        plt.pause(0.01)

    #In case the .fits was being written
    except IOError:
        print "\nTrouble opening file (ignore if occasional message only)"
