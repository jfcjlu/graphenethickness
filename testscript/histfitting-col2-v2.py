# -*- coding: utf-8 -*-
"""
Created on Fri Aug  2 16:40:59 2019

@author: lenovo
"""

import pandas as pd
import numpy as np
import os
from itertools import chain
import matplotlib.pylab as plt
from scipy.optimize import curve_fit

def gauss1p(x,y0,a1,m1,s1):
    return y0+a1*np.exp(-((x-m1)/s1)**2)


FOLDER = '.\sampledata'
SAMPLE = 'analysissample'
FILElist = ['analysis-1', 'analysis-2', 'analysis-3','analysis-4']
ampfactor = 1.0   # transform the unit to nm



binstart=-1
binmid = 0.5
binend=3
binsize = 0.02
binrange = np.arange(binstart,binend+binsize, binsize)
plotfitx = np.arange(binstart,binend+binsize, 0.001)
binnumber = int((binend-binstart)/binsize)
binlownumber = int((binmid-binstart)/binsize)
bincenter = np.arange(binstart+binsize/2,binend, binsize)

summaryfolder = FOLDER + '\\' + 'indivlinesummary'

if not os.path.exists(summaryfolder):
    os.mkdir(summaryfolder)

HEIGHTSUM = summaryfolder + '\\' + SAMPLE + '.txt' 

if os.path.exists(HEIGHTSUM):
    os.remove(HEIGHTSUM) 

with open(HEIGHTSUM,'a+', encoding='utf-8') as sumf:
    sumf.write('linenumber'+'\t'+'bottomheight'+'\t'+'upheight'+'\t'+'height'+'\n')
    for FILE in FILElist:
        PATH = FOLDER + '\\' + FILE + '.txt' 
        outhist = summaryfolder+'\\' + FILE + '_hist.txt'
        outfig =  summaryfolder+'\\' + FILE  + '_fit.png'
        
        if os.path.exists(outhist):
            os.remove(outhist) 
        
        with open(PATH) as orgf:
            countline = len(orgf.readlines())
        
        train = pd.read_csv(PATH, sep = '\t', header = 0, index_col = False, engine='python') 
        datalist = train.values.tolist()
        
        heightlist=[]
        
        for i in range(0, countline-1):
            heightlist.append(datalist[i][1]*ampfactor)
            
#        plt.hist(heightlist, binrange)
#        plt.show()
            
        A=np.histogram(heightlist, bins=binrange)
        
        
        histx = bincenter 
        histy = A[0]
        
        
        with open(outhist,'a+', encoding='utf-8') as histf:
            for j in range(0, binnumber-1):
                histf.write(str(histx[j])+'\t'+str(histy[j])+'\n')

        firsthistx = histx[0:binlownumber-1]
        secondhistx = histx[binlownumber:binnumber-1]
        firsthisty = histy[0:binlownumber-1]
        secondhisty = histy[binlownumber:binnumber-1]
                
        popt1, pcov1 = curve_fit(gauss1p, firsthistx, firsthisty)
        popt2, pcov2 = curve_fit(gauss1p, secondhistx, secondhisty)
        
        bh = popt1[2]
        uph = popt2[2]
        
        plotfity1 = popt1[0]+popt1[1]*np.exp(-((plotfitx-popt1[2])/popt1[3])**2)
        plotfity2 = popt2[0]+popt2[1]*np.exp(-((plotfitx-popt2[2])/popt2[3])**2)
        
        plt.figure(1)
        plt.hist(heightlist, binrange,color='b')
        plt.plot(plotfitx, plotfity1,'r')
        plt.plot(plotfitx, plotfity2,'k')
#        plt.show()
        plt.savefig(outfig)
        plt.close(1)
        

        sumf.write(FILE+'\t'+ str(bh) +'\t'+ str(uph) +'\t'+str(uph-bh)+'\n')
        
            