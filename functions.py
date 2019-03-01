from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import ipdb
import pandas as pd
from interpolation import interpolation as inter
from astropy.table import Table
import sys
import warnings

def binningopt(input_path, flux_error, velocity, names, line, center, AGN, opt, output_path, spec = None, lam = None, bmin = None, bmax = None, tol = 20.0):

   #########################################################################################################

   '''This function determines the best binning according to the user specifications. It takes 8 to 13 arguments.
   
   0. path. String with the path where the files flux_error and velocities are kept.

   1. flux_error. String with the name of the file where fluxes and errors are kept. 

   4. velocities. Array with the velocities of each line at each position of the galaxy. Each column corresponds to
      one line. Each row to one position of the galaxy.

   8. names. Default: None.  Lsit of strings containing the names of the provided lines, ie. ['hb', 'ha', 'NII']. It only
      needs to be specified if the input data is in a csv format, since there are no line names by default there. 

   3. line. String specifying which line to bin on. It has to match one of the names in the header of the input fits
      file, or one of the elements in the array names if the input is a csv file.

   5. center. Integer with the row number corresponding to the center of the galaxy.

   6. AGN. Boolean specifiying wether the center of the galaxy is or isn't an AGN.

   7. opt. Float, int, or list. Depending on the data type, the code will bin differently. 
   
   9. spec. Default: None. String with the name of the file containing the spectrum of the galaxy. Supported filetypes are .csv and 
      .fits.

   10. lam. Default: None.  Array with wavelengths. This array can be spaced anyhow (it could be given in a 
      logarithic scale). 
                          
   11,12. bmin, bmax. Ranges of the preliminary binning. By default, the code prebins in columns of pixels by 3 THIS          WILL CHANGE. IT WILL BECOME ARGUMENT SET  TO 3 BY DEFAULT, JUST LIKE TOL. THE USER MAY BE ABLE TO  MODIFY           IT WITH BMIN AND BMAX. THEY CANNOT REMAIN IN OPTIONAL ARGUMENTS, BECAUSE THEY COULD BE SPECIFIED WITHOUT A          SPECTRUM, AND THE CODE WOULDNT WORK BECAUSE ONLY THAT HAS BEEN SPECIFIED, AND THE CODE IS CONCIEVED TO REC          IEVE A WHOLE BUNCH OF PT ARGS IF YOU WANT TO SPECIFY THE SPECTRUM (ACTUALLY ONLY 9 AND 10).
                          
   13- tol. Default: 20(%). Float specifying how much can the first and last row  differ from the S/N target. 

   #### BINNING OPTIONS ####

   Depending on what dtype is the 6th parameter, the code offers different options:
 
   OPTION 1. If the input is 0, the function bins up the best rows that don't contribute to decrease s/n. So 
   basically we bin up the best lines. It returns the s/n reached and the rows binned.
   
   OPTION 2. If you input a list of sets of rows to bin, it will return anotherlist with s/n reached in each new 
   position, and other two empty parameters. Note that the list of rows to bin should not have repeated lines: ie. 
   opt  =  [0,3,4,8,9,12] and **not opt = [0,3,3,8,8,12]**, because in that case we would be including the rows 3
   and 8 twice. In the former case, the code would bin [0,3], [4,8], [9,12]. Also note that when you write [0,3] 
   you are telling the program to sum from row 0 (1st row), to row 3 (4th row). 
   
   OPTION 3. If you in input a S_N target (either a float or int), it will return a list (rows that need to be 
   binned so that the S_N in everybinned rows is equal or higher than the target. It also returns a list with the 
   new s/n in each postition after the binning. Both lists are sotrted as the original rows, from the left of the
   galaxy to the right. '''

   #########################################################################################################b
   #Load files with fluxes, errors and velocities.
   err = ["" for x in range(len(names))] #Create an array with the names of lines and errors, alternatively.
   for i in np.arange(len(names)):
       err[i] = names[i]+'_err'
   
   namecol = np.asarray(zip(names, err)).flatten()
   data = pd.read_csv(input_path + '/' + flux_error, header =None, names = namecol)
   vel = np.genfromtxt(input_path + '/' + velocity) #Just an array.
   if not vel.any():
      raise ValueError("All the elements in %s are zero, can't compute anything!" % velocity )
   #Select the specified line from provided data, remove zeros and rename indices and redefine center:
   flux = data[line]
   if not flux.any():
      raise ValueError("All the elements in flux are zero, can't compute anything!")
   fc = flux[center]
   flux  = flux[flux != 0].reset_index(drop = True)
   error = data[line+'_err']
   ec = error[center]
   error = error[error !=0].reset_index(drop = True)

   #check if there are negative fluxes 
   boolean = flux < 0
   if  boolean.all():
      raise ValueError("Flux values are negative at all positions. No action taken.")
   if boolean.any():
      print('WARNING: There are some negative fluxes. They have not been taken into account')
      
   if fc != 0:
      center = np.where(flux == fc)[0]
      velcenter = (vel[vel!=0])[center]
      icor = vel/velcenter #Correcting flux and error from the fits.
      icor = icor[icor!=0]
      flux = flux * icor
      error = error * icor
   else:
      warnings.warn("The flux in the center row is 0. So is the velocity.\n Spectrum (if provided) can't be aligned ")

   #Set some initial values
   specnew = None
   binnedspec = None
   listb = None
   qnew = None
   dropf = False #A flag that tells us if rows have been dropped due to weak s/n
   dropl = False

   if spec is not None and fc != 0: #If the user specifies a spectrum, the alignment is computed here
      #First we redimensionalize the spectrum so it has the same dimension as the array of velocities, and
      #therefore it can be corrected.Doing this takes spec to  datanew, which later on, once corrections
      #have been performed, becomes specnew.
      spect = pd.DataFrame((fits.open('slincen_rf0283bspec.fits'))[0].data)	
      siz = spect.shape[1]
      datanew = pd.DataFrame(0, index = np.arange(len(bmin)), columns = np.arange(siz))
      q= np.asarray(zip(bmin,bmax)) #q is an array with 2 columns containing the row ranges determining the binning of the 
                                    #spectrum.
      q = np.reshape([ int(x) for x in q.reshape(21*2, 1) ], (21,2)) #Transform to integers to prevent warnings.
      for i in np.arange(0,len(bmin)):
         datanew.iloc[i] = (spect.iloc[q[i][0]:q[i][1]]).sum(axis =0)
         if i == 10:
            datanew.iloc[i] = spect.iloc[172]
      specnew = inter(spect, lam, vel, bmin, bmax, datanew) #Aligned spectrum.
      q = q[vel!=0]

   #Subtracting center from fits data and spectrum (if provided) if AGN is declared:
   if AGN and fc != 0:
      flux = flux[flux!=flux[center].values[0]].reset_index(drop=True)
      error = error[error !=error[center].values[0]].reset_index(drop=True)
      print('Center subtracted due to present AGN')
      if not flux.any():
         raise ValueError("All the elements in flux are zero, can't compute anything!")
      if spec is not None:
         specnew = np.delete(specnew, (center), axis = 0)
         print('Center of spectra also subtracted')

  

   #######   OPTION 1: Binning of the best rows   #######

    
   if opt == 0:
      #The first step is order from higer to lower s/n, and get a vector with
      #the positions of the ordered s/n in the original array
      ind = np.asarray(np.argsort(flux/error)[::-1])
      S_N = np.zeros(len(flux))
      binnedflux = np.zeros(len(flux))
      binnederror = np.zeros(len(error))
      s = n = 0
      t = 0
      listb = []
      for i in ind:
         s = s + flux[i]
         n = np.sqrt(n**2 + error[i]**2)
         binnedflux[t] = s
         binnederror[t] = n
         S_N[t]  = binnedflux[t]/binnederror[t]
         if t != 0 and S_N[t] <= S_N[t-1]:
            listb = ind[0:t]
            binnedflux = binnedflux[binnedflux != 0]; binnedflux = binnedflux[-2]
            binnederror = binnederror[binnederror != 0]; binnederror = binnederror[-2]
            S_N = (binnedflux/binnederror)
            break
         t += 1

      if t ==len(flux): #When all measurements increase S_N when added, the previous for will not break.
         listb = ind[0:t]
         if len(flux) == 1:
            S_N = float((binnedflux/binnederror))
         else:
             binnedflux = binnedflux[binnedflux != 0]; binnedflux = binnedflux[-1] #to avoid selecting the last element
                                                                                   #which contributes to worsen the S/N.
             binnederror = binnederror[binnederror != 0]; binnederror = binnederror[-1]
             S_N = (binnedflux/binnederror)

      #Bin the aligned spectrum according to the user specifications.
      if spec is not None and fc != 0:
          binnedspec = np.zeros((1,spect.shape[1])) #Inintial value binned spectra
          for i in listb:
            binnedspec = binnedspec + specnew[i, :]



   #######   OPTION 2: Binning of given sets of rows   #######


   if type(opt) == list: 
      S_N = np.zeros(len(opt))
      binnedflux = np.zeros(len(opt))
      binnederror = np.zeros(len(opt))
      if spec is not None and fc != 0:
         binnedspec = np.zeros((len(opt)/2, spect.shape[1]))
      for i in np.arange(0, len(opt)/2):
         binnedflux[i] = sum(flux[opt[2*i]:opt[2*i+1]])
         binnederror[i] = np.sqrt(sum(error[opt[2*i]:opt[2*i+1]]**2))
         S_N[i] = binnedflux[i]/binnederror[i]
         S_N = np.nan_to_num(S_N)
         if spec is not None and fc != 0:
            binnedspec[i, :] = sum(specnew[opt[2*i]:opt[2*i+1]]) #Bin the aligned spectrum according to opt.
      
      binnedflux = np.trim_zeros(binnedflux)
      binnederror = np.trim_zeros(binnederror); S_N = np.trim_zeros(S_N)
      


   #######   OPTION 3: Binning to a S/N target   #######



   if opt !=0 and  type(opt) != list:
      #Set initial values
      binnedflux = np.zeros(len(flux))
      binnederror = np.zeros(len(flux))
      s_n = None
      S_N = np.zeros(len(flux))
      binnedflux[center] = flux[center]
      binnederror[center] = error[center]
      S_N[center] = binnedflux[center]/binnederror[center]
      i  = 1
      listb = np.zeros(2*len(flux))
      nums = np.arange(0,len(flux))
      j = 0 
      z = 0 
      up  = center + i 
      down = center - i


      ### Binning the center ###

      while S_N[center] < opt:  #The center is binned symmetrically
         try:
            
            binnedflux[center] = binnedflux[center] + flux.values[center+i] + flux.values[center-i]
            binnederror[center]=np.sqrt(binnederror[center]**2 + error.values[center+i]**2 + error.values[center-i]**2)
            S_N[center] = np.nan_to_num(binnedflux[center]/binnederror[center])
            up = center + i+1 
            down = center - i-1
            i += 1
         except IndexError:
            raise IndexError('Can not reach the desired target, try a smaller S/N')
         

      #We prepare the start points for the binning of the sides
      up0 = up 
      down0 = down

      # Assign the binned center rows to the output
      listb[j], listb[j+1] = down+1, up-1
      j = j+2
      z = z+1


      ### Binning the left/upper side ###

      while down0 > -1:
         while S_N[center-z] < opt: 
            binnedflux[center-z]= binnedflux[center-z]+ flux.values[down]
            binnederror[center-z] = np.sqrt(binnederror[center-z]**2 + error.values[down]**2)
            S_N[center-z] = binnedflux[center-z]/binnederror[center-z]
            down = down - 1
            if down == -1:
               print 'I have binned everyting on the first group and the S/N reached is:', S_N[center-z]
               break

         listb[j], listb[j+1] = down0, down + 1
         j = j+2
         z = z+1
         down0 = down
         
      #If the first element doesn't meet the s/n target within the tolerance, it is dropped.
      if S_N[S_N!=0][0]<opt-tol/100*opt:
         S_N[np.where(S_N == S_N[S_N!=0][0])] = 0
         binnedflux[np.where(binnedflux == binnedflux[binnedflux!=0][0])] = 0
         binnederror[np.where(binnederror == binnederror[binnederror!=0][0])] = 0
         print 'The first row has been removed'
         dropf = True

      z = 1 # Reset z to go from the center to the right side
         
      
      ### Binning the right/lower side ###

      while up0 < len(flux):
         while S_N[center+z] < opt:
            binnedflux[center+z] = binnedflux[center+z] + flux.values[up]
            binnederror[center+z] = np.sqrt(binnederror[center+z]**2 + error.values[up]**2)
            S_N[center+z] = np.nan_to_num(binnedflux[center+z]/binnederror[center+z])
            up = up + 1
            if up  == len(flux):
               print 'I have binned everyting on the last  group and reached a S/N of:', S_N[center+z]
               break

         listb[j], listb[j+1] =  up0, up-1 #Note that if two consecutive elements in list 
                                           # are the same, it means that that row by itself has already  a S/N
                                           #higher than the target, so no binning is required. 
         j = j+2
         z = z +1
         up0 = up 

      #If the last element doesn't meet the s/n target within the tolerance, it is dropped.
      if S_N[np.nonzero(S_N)][-1]<opt - tol/100*opt:
         S_N[np.nonzero(S_N)[-1][-1]] = 0
         binnedflux[np.nonzero(binnedflux)[-1][-1]] = 0
         binnederror[np.nonzero(binnederror)[-1][-1]] = 0
         print 'The last row has been removed'
         dropl = True
      
      #Making readable and usable outputs before adjusting the endpoints:
      listb = np.trim_zeros(listb).tolist()
      listb = sorted(listb)
      S_N = np.trim_zeros(S_N).tolist()
      if dropf:
         listb = listb[2:] #remove first row
      if dropl:
         listb = listb[:-2] #remove last row
      for i in np.arange(0,len(listb)): #converting elements in list to int.
         listb[i] = int(listb[i])

      #Bin the aligned spectrum according to the user specifications.
      if spec is not None and fc != 0:
         si = 0
         binnedspec = np.zeros((len(listb)/2,specnew.shape[1]))
         for i in np.arange(0, len(listb),2):
            binnedspec[si,:] = sum(specnew[listb[i]:listb[i + 1]+1])
            si +=1

      listb = np.reshape(listb, (len(listb)/2, 2)) #Making it more readable.

   if spec is not None and fc != 0:
      #Transform to DataFrame so that it is more readable.
      specnew = pd.DataFrame(specnew, columns = np.arange(siz))
      binnedspec = pd.DataFrame(binnedspec, columns = np.arange(siz))
      #fits.writeto('slincen_rf0283bspecdb.fits', binnedspec, hdr, overwrite = True)
      qnew = np.zeros((len(listb), 2))
      for i in np.arange(len(listb)):
         qnew[i, 0], qnew[i,1] = q[listb[i,0], 0], q[listb[i,1], 1]
   try:
       import os
       os.mkdir(output_path + '/results')
   except:
       print 'Existing results directory. Files are being stored at: (not yet)', output_path + '/results'
   return(listb # Array with the optimal binning row ranges.
          , S_N # Array S/N reached after the binning.
          , specnew # DataFrame containing the aligned spectrum.
          , binnedflux # Array with binned fitted fluxes.
          , binnederror  # Array with binned fitted errors.
          , binnedspec # Array with the binned spectrum
          , qnew # Array with the optimal binning in the original rows.
          )

