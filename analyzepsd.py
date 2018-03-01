# -*- coding: utf-8 -*-
"""

@author: Mathieu Gibeault - Coda Consulting
"""

import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as mplcm
import matplotlib.colors as colors
import os
import ntpath
from glob import glob
import ConfigParser
import collections
from scipy.optimize import curve_fit


class AnalyzePSD:
    
    def __init__(self):
        pass
    

    def process(self,DataPath, outputPath, include = [], exclude = []):
        
        if not os.path.isdir(outputPath):
            print("Output path does not exist: " + outputPath)
            return
        files = getfilelist(DataPath,'T*.PSD', include, exclude)
        if not files:
            print('No PSD Data file found in: ' + DataPath)
        else:
            self.PSD_params = {}
            for filepath in files:
                print("Analysing file: " + filepath)
                try:
                    params = compute_PSD_params(filepath)
                    chipID = getchipID(filepath)
                    if chipID in self.PSD_params:
                        self.PSD_params[chipID] = self.PSD_params[chipID].append(pd.DataFrame([params],index=[ntpath.basename(filepath)]))
                    else:
                        self.PSD_params[chipID] = pd.DataFrame([params],index=[ntpath.basename(filepath)])
                except:
                    print("Warning: Analysis FAILD: " + filepath)
            self.write_params_to_file(outputPath)
    
    def write_params_to_file(self,outputPath):
        for chipID, data in self.PSD_params.items():
            path = os.path.join(outputPath,'PSD_Params_' + chipID + '.csv')
            data.to_csv(path)
        

        

###############################################################################
#
#  Module Functions
#
###############################################################################

def getfilelist(path, searchString, include, exclude):  
    files = [y for x in os.walk(path) for y in glob(os.path.join(x[0], searchString))]
    results = []
    if include:
        for chip in include:
            results.extend([s for s in files if chip in s])
    else:
        results = files
    if exclude:
        for chip in exclude:
            results = [ x for x in results if chip not in x]
    return results
    
def getchipID(path):
    path = path.replace('__','_')
    filenamestrings = ntpath.basename(path).split('_')
    return filenamestrings[2]   

def fitfunc(f, f0, alpha, fstar, offset):
    return np.log10((f0/f)**alpha + alpha*(f0/fstar)**(1+alpha)*(f/f0) + offset)


def corrected_L(f, S, f0, alpha, fstar, offset, df, B):
    integrand = S - alpha*(f0/fstar)**(1+alpha)*(f/f0) - offset
    return np.sqrt(np.sum(integrand)*df/B)

def old_L(f, S, df, B):
    return np.sqrt(np.sum(S)*df/B)

def compute_PSD_params(psdfile):

    B = 100000
    N = 10000
    psd = pd.read_csv(psdfile,sep='\t',names=['f','S','integral','norm'])
    f = psd['f'].values[1:N]
    df = f[1]-f[0]
    S = psd['norm'].values[1:N]
    popt, pcov = curve_fit(fitfunc, f, np.log10(S), p0=[1.0,1,1000.0, 0.0001], sigma=np.sqrt(np.arange(1,N)+np.sqrt(3)/3), maxfev=10000)
    #print popt
    #print np.sqrt(np.diag(pcov))


    f0 = popt[0]
    alpha_2 = popt[1]
    fStar = popt[2]
    Offset = popt[3]
    
    #L = old_L(f[:100], S[:100], df, B)
    L_2 = corrected_L(f[:100], S[:100], f0, alpha_2, fStar, Offset, df, B)
    
    return {'L_2': L_2, 'alpha_2': alpha_2, 'f0': f0, 'fStar': fStar, 'Offset' : Offset }
    
    
#    plt.loglog(f,S)
#    plt.loglog(f, 10**fitfunc(f, f0, alpha, fstar, offset))
#    plt.show()

#    print 'Adjusted L: {0:.6f}'.format(corrected_L(f[:100], S[:100], f0, alpha, fstar, offset, df, B))
#    print 'Old L: {0:.6f}'.format(old_L(f[:100], S[:100], df, B))
#    print 'New alpha: {0:.6f}'.format(popt[1])