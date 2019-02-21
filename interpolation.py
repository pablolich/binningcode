'''This file will perform the first rebinning of the wavelengths, including interpolation'''
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import pdb
import pandas as pd

def interpolation(spec, lam, vel, bmin, bmax, datanew):
    
    datanew = datanew.drop(np.where(vel ==0 )[0]) #Renounce to some spectrum data becouse the velocities are 0.
    vel = vel[vel!=0]    
    velcenter = vel[10]


    #Create matrix of corrected lambdas. Only need the dopler velocities are needed to dedopplershift. Get these by subtracting 
    #velcenter to the vector of velocities.

    vdop = vel - velcenter
    vdop[vdop == -velcenter] = 0


    #This lambda is the lambda of the center, since it hasn't been corrected
    cor = velcenter/vel
    lamcor = np.transpose(np.repeat(lam, len(vdop)).reshape(len(lam), len(vdop)))
    for i in np.arange(len(vel)-1): #Correct each row by its velocity in the galaxy.
        lamcor[i,:] = lamcor[i,:]*cor[i]
    
    datanewcor =datanew.mul(cor**(-1), axis = 0)#Array of fluxes that have been summed up to match the 21 positions
                                                #where we have velocities. It has  been corrected inversly to assure
                                                #photon conservation. 

    #Now that we have deredshifted wavelength and corrected flux, we rebin it and rescale it. 
    from pysynphot import observation
    from pysynphot import spectrum 
    def rebin_spec(wave, specin, wavnew):
        spec = spectrum.ArraySourceSpectrum(wave=wave, flux=specin)
        f = np.ones(len(wave))
        filt = spectrum.ArraySpectralElement(wave, f, waveunits='angstrom')
        obs = observation.Observation(spec, filt, binset=wavnew, force='taper') 
        return obs.binflux

    fluxnew = np.zeros((len(vdop), len(lam)))
    plt.figure(1)
    plt.xlim((6660, 6780))
    for i in np.arange(0, len(vdop)):
        fluxnew[i, :] = rebin_spec(lamcor[i,:], datanewcor.iloc[i].values, lam)
        plt.plot(lam, fluxnew[i,:])

    #Now we write the deredshifted fluxes to a fits file so that we  can visualize it with ds9 and load it to make the rebinning
    #Fixing unprintable characters
    #you have to do it twice, otherwise the poor guy doesn't get it  
    
    #del hdr['PARAM0']
    #del hdr['PARAM61']
    #del hdr['PARAM62']
    #del hdr['PARAM63']
    
    #fits.writeto('slincen_rf0283bspecd.fits', fluxnew, hdr, overwrite = True)
    #fits.writeto('slincen_rf0283bspecdatanew.fits',datanew, hdr, overwrite = True)
    #hdulmod = fits.open('slincen_rf0283bspecd.fits')  

    return fluxnew
