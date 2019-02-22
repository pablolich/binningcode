                              ####### WRAPPER CODE TO TEST functions.py AND interpolation.py ########
import numpy as np
from astropy.io import fits
from astropy.table import Table
import os #to delete files once they have been used.
import ipdb
import pandas as pd 
import sys 
from functions import binningopt
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

#Load the fits file as a table so that it is easier to read.
t = Table.read('merged_bztre.fits', format = 'fits')

#Create a csv file that my code accepts, with all the lines and errors.
selec = np.asarray([np.asarray([i[0]=='F' and i[-1]!='R' for i in np.asarray(t.keys())])])[0] #selecting names of lines from table.
lines = np.asarray(t.keys())[selec] 


selerr = np.asarray([np.asarray([i[0]=='F' and i[-6:]!='SYSERR' for i in np.asarray(t.keys())])])[0]
linerr = np.asarray(t.keys())[selerr] #To extract the data (fluxes and errors), create selerr
data = t[linerr.tolist()][:] #Table with all lines and errors (cols) (not systematic errors) of all galaxies (rows).
AGN_flags = np.random.randint(2, size = len(t))
l1 = []
SN1 = []
l2 = []
SN2 = []
SN3 = []
# Apply the code to all galaxies

#NO SPECTRUM

totalbins = np.zeros(len(lines))
for i in np.arange(len(t)): # Runs through every galaxy
    csv_file = np.zeros((len(linerr),len(data.columns[0][0])))    
    print 'galaxy #:', i
    for j in np.arange(len(data.columns)): #Create ith csv file with all lines and errors for the ith galaxy.
        csv_file[j,:] = data.columns[j][i]

    csv_file_df = pd.DataFrame(csv_file.transpose(), columns = linerr)
    vel = pd.DataFrame(t['L8_PARAM'][i,:,0]) #Velocity of the ith galaxy at each position ([:]).
    vel.to_csv('vel.csv', header = False, index=False) #save data to be loaded by the function.
    csv_file_df.to_csv('gal.csv', header=False, index=False)
    
    #OPTION 1
    try:
        listb, S_N,  specnew, binnedflux, binnederror, binnedspec, qnew = binningopt(flux_error = 'gal.csv', 
                                                                                     velocity = 'vel.csv',
                                                                                     names = lines, 
                                                                                     line = lines[3],
                                                                                     center = 10,
                                                                                     AGN = AGN_flags[i],
                                                                                     opt = 0)
        
        l1.append(listb)
        SN1.append(S_N)
        S_N_lines = [] 
        for k in np.arange(len(lines)): #For every galaxy loop through all the lines.
            l = t[linerr[2*k]][i] #list with the fluxes 
            l = l[l!=0] #subtract 0
            if len(l) == 0:
                S_N_lines.append('Error') 
                pass
            else:
                signal = sum(l[l1[i]])
                err = t[linerr[2*k+1]][i]
                err = err[err!=0]
                noise = np.sqrt(sum(err[l1[i]]**2))
                S_N_lines.append(signal/noise)
        totalbins = np.vstack((totalbins,S_N_lines))
    except (ValueError, IndexError) as e:
        print(e)
        l1.append('Error')
        SN1.append('Error')
        totalbins = np.vstack((totalbins, np.zeros(len(lines))))

nerrs = len([i for i in l1 if i == 'Error'])
totalbins = totalbins[1:,:]
totalbins = pd.DataFrame(np.transpose(totalbins), index = lines, columns = t['NAME'])

#Bin all the lines according to l1.

'''

    #OPTION 2
    try:
        listb, S_N,  specnew, binnedflux, binnederror, binnedspec, qnew = binningopt(flux_error = 'gal.csv', 
                                                                                     velocity = 'vel.csv',
                                                                                     names = lines, 
                                                                                     line = lines[2],
                                                                                     center = 10,
                                                                                     AGN = AGN_flags[i],
                                                                                     opt = [8,11])
        SN3.append(S_N)
    except (ValueError, IndexError) as e:
        print(e)
        SN3.append('Error')    

    #OPTION 3
    try:
        listb, S_N,  specnew, binnedflux, binnederror, binnedspec, qnew = binningopt(flux_error = 'gal.csv', 
                                                                                     velocity = 'vel.csv',
                                                                                     names = lines, 
                                                                                     line = lines[2],
                                                                                     center = 10,
                                                                                     AGN = AGN_flags[i],
                                                                                     opt = 50)
        l2.append(listb)
        SN2.append(S_N)
    except (ValueError, IndexError) as e:
        print(e)
        l2.append('Error')
        SN2.append('Error')

'''
#Apply the binning to the rest of the lines.

#Add spectrum, lam, bmin, bmax as parameters.
