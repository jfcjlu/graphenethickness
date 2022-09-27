# -*- coding: utf-8 -*-
"""
Created on Fri Aug  2 16:40:59 2019

@author: lenovo
"""

import pandas as pd
import numpy as np
import os
import time
from itertools import chain
import matplotlib.pylab as plt
from scipy.optimize import curve_fit

def gauss1p(x,y0,a1,m1,s1):
    return y0+a1*np.exp(-((x-m1)/s1)**2)


numofgr = 50
parentFOLDER = '.\sampledata'
FILE='analysis-whole' 

 
headerlines = 4
datasufix= '.txt'
ampfactor = 1.0e9   # transform the unit to nm

#p001=[9.60225704900395,100000,-0.02,0.2]
#p002=[27.998298554631248,17586,0.89,0.2]


binstart=-10
binmid = 0.1
heightcritical = 0.2
binend= 30 
binsize = 0.02
binrange = np.arange(binstart,binend+binsize, binsize)
plotfitx = np.arange(binstart,binend+binsize, binsize)
binnumber = int((binend-binstart)/binsize)
binlownumber = int((binmid-binstart)/binsize)
bincenter = np.arange(binstart+binsize/2,binend, binsize)


summaryfolder = parentFOLDER + '\\' + 'planesummary'
if not os.path.exists(summaryfolder):
    os.mkdir(summaryfolder)

Ssummaryfile = summaryfolder + '\\' 'Smallplane_summaryfile_indivplane.txt' 
localtime = time.asctime( time.localtime(time.time()) ) 

with open(Ssummaryfile,'a+', encoding='utf-8') as Sallsumf:
    Sallsumf.write('the-name-of-sample'+'\t'+'bottomheight(nm)'+'\t'+'upheight(nm)'+'\t'+'height(nm)'+\
                    '\t'+'bottomstddev'+'\t'+'upstddev'+'\t'+'upstddev-bottomstddevM'+'\t'+'bottom_W'+'\t'+'up_W'+'\n')
    Sallsumf.write(str(localtime) + '\n')
    currentfolder = summaryfolder + '\\' + FILE
    if not os.path.exists(currentfolder):
        os.mkdir(currentfolder)
    HEIGHTSUMfile = currentfolder + '\\' + FILE + '_smallplanesummary.txt' 
    if os.path.exists(HEIGHTSUMfile):
        os.remove(HEIGHTSUMfile)
    CONFIGfile = currentfolder + '\\' + FILE + '_configurationforallplane.txt' 
    if os.path.exists(CONFIGfile):
        os.remove(CONFIGfile)
    with open(HEIGHTSUMfile,'a+', encoding='utf-8') as smallsumf:
        with open(CONFIGfile,'a+', encoding='utf-8') as configsumf:
            localtime = time.asctime( time.localtime(time.time()) )
            smallsumf.write("processed time:" +'\t' +str(localtime) + '\n') 
            configsumf.write("processed time:" +'\t' +str(localtime) + '\n')
            PATH = parentFOLDER + '\\' + FILE + datasufix 
            outhist = currentfolder + '\\' + FILE + '_smallplanehist.txt'
            outfig = currentfolder + '\\' + FILE + '_smallplanefit.png'                
            
            train = pd.read_csv(PATH, sep = '\t', header = headerlines, index_col = False, engine='python')
            datalist = train.values.tolist() 
            clist= list(chain(*datalist))        
            datalistsmall = []        
            with open(PATH) as orgf:
                countline = len(orgf.readlines())        
            for i in range(0,countline-headerlines-1):
                heightlist =  datalist[i][:]
                heightlstlength = len(heightlist)
                heightlistnm = []
                for j in range(0, heightlstlength-1):
                    heightlistnm.append(heightlist[j]*ampfactor)
                countgr = 0
                for j in range(0, heightlstlength-1):
                    if heightlistnm[j] > heightcritical : countgr += 1
                if countgr > numofgr:
                    datalistsmall.append(heightlist) 
                    
            clist= list(chain(*datalistsmall))
            if len(clist) >0:
                ONECOLUMEOUT = currentfolder + '\\' + FILE + '_smallplaneoutto1colume.txt'
                if os.path.exists(ONECOLUMEOUT):
                    os.remove(ONECOLUMEOUT)        
                with open(ONECOLUMEOUT,'a+' ,encoding='UTF-8') as outf:
                    for line in clist:
                        outf.write(str(line) + '\n')       
                with open(ONECOLUMEOUT, 'r' ,encoding='UTF-8') as orgf:
                    countline = len(orgf.readlines())        
                train = pd.read_csv(ONECOLUMEOUT, sep = '\t', header = None, index_col = False, engine='python') 
                datalist = train.values.tolist()        
                heightlist=[]        
                for i in range(0, countline-1):
                    heightlist.append(datalist[i][0]*ampfactor)        
                heightlist.sort()
                numberofheight = len(heightlist) 
                theestimatestart = heightlist[int(0.05*numberofheight)] - 0.95
                theestimateend = heightlist[int(0.95*numberofheight)] + 1.95
                if theestimatestart > 0:
                    binstart = round(theestimatestart)
                else:
                    binstart = round(theestimatestart)-1
                if theestimateend > 0:
                    binend = round(theestimateend)
                else:
                    binend = round(theestimateend)-1
                binmid = (theestimatestart + theestimateend)/2
                binrange = np.arange(binstart,binend+binsize, binsize)
                plotfitx = np.arange(binstart,binend+binsize, binsize)
                binnumber = int((binend-binstart)/binsize)
                binlownumber = int((0.8*binmid-binstart)/binsize)
                bincenter = np.arange(binstart+binsize/2,binend, binsize)                                         
               
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
                        
                try:
                    popt1, pcov1 = curve_fit(gauss1p, firsthistx, firsthisty,p0=None)
                    popt2, pcov2 = curve_fit(gauss1p, secondhistx, secondhisty,p0=None)
                    
                    if popt1[1] > 0 and popt2[1] >0 \
                    and (popt1[2] > binstart and popt1[2] < binmid) \
                    and (popt2[2] > binmid and popt2[2] < binend):
                        bh = popt1[2]
                        uph = popt2[2]
                        bw = abs(popt1[3])/np.sqrt(2.0)
                        upw = abs(popt2[3])/np.sqrt(2.0)
                                       
                        plotfity1 = popt1[0]+popt1[1]*np.exp(-((plotfitx-popt1[2])/popt1[3])**2)
                        plotfity2 = popt2[0]+popt2[1]*np.exp(-((plotfitx-popt2[2])/popt2[3])**2)
                        
                        plt.figure(1)
                        plt.hist(heightlist, binrange,color='b')
                        plt.plot(plotfitx, plotfity1,'r')
                        plt.plot(plotfitx, plotfity2,'k')
                        plt.xlabel('height(nm)',fontsize=10)
                        plt.ylabel('counts',fontsize=10)
                        plt.xlim(-0.9 , 1.9)
                        plt.xticks(fontsize=10)
                        plt.yticks(fontsize=10)
                        plt.text(0.2*(bh+uph), 0.5*np.max(histy), 'Sample: ' + FILE + 'smallplane'+ '\n'\
                                 +'lowerheight='+ str(round(bh,4)) +' nm'+'\n'\
                                 +'upperheight='+ str(round(uph,4)) +' nm'+'\n'\
                                 +'height=' + str(round((uph-bh),4)) +' nm'+'\n'\
                                 +'lowersigma=' + str(round(bw,4)) +' nm'+'\n'\
                                 +'uppersigma=' + str(round(upw,4)) +' nm'+'\n'\
                                 ,fontsize=10)
                        
                        plt.show()
                        plt.savefig(outfig)       
    #                    plt.close(1)
                        Sallsumf.write(FILE + '\t' + str(round(bh,4)) +'\t'+ str(round(uph,4))+'\t'\
                                    +str(round((uph-bh),4))+'\t'+str(round(bw,4))+'\t'\
                                    +str(round(upw,4))+'\t'+str(round((upw-bw),4))+'\t'\
                                    +str(round(popt1[3],4))+'\t'+str(round(popt2[3],4))+'\n')
                        smallsumf.write('Sample: '+ '\t' + currentfolder + '\\' + FILE+'\t'+ 'smallplane' + '\n' \
                                   + 'bottomheight(nm):' + '\t' + str(round(bh,4)) +'\n'\
                                   + 'upperheight(nm):' + '\t' + str(round(uph,4)) +'\n'\
                                   + 'height:(nm)' + '\t' + str(round((uph-bh),4))+'\n'\
                                   + 'loweruncertaimty(nm):' + '\t' + str(round(bw,4))+ '\n'\
                                   + 'upperuncertaimty(nm):' + '\t'+ str(round(upw,4)) + '\n'\
                                   +'fitting equation:' + '\t' + 'y0+a1*np.exp(-((x-m1)/s1)**2)' + '\n' \
                                   +'fitting parameters:' + '\n' \
                                   + 'p01:' +'\n'\
                                   + str(popt1[0]) + '\t' + str(popt1[1]) + '\t' +str(popt1[2]) + '\t' + str(popt1[3]) + '\n'\
                                   + 'p02:'+'\n'\
                                   + str(popt2[0]) + '\t' + str(popt2[1]) + '\t' +str(popt2[2]) + '\t' + str(popt2[3]) + '\n')
                         
                                   
                        configsumf.write('FOLDER:' + '\n'\
                                   + "r'" + currentfolder + "'" +'\n'\
                                   + 'FILE:' + '\n'\
                                   + "'"+ FILE + "'" +'\n'\
                                   + 'the first data line is binstart binend binsize binmid' + '\n'\
                                   + 'the second data line is p01' + '\n'\
                                   + 'the third line is p02' + '\n'\
                                   + 'binstart' + '\t' + 'binend' + '\t' + 'binsize' + '\t' + 'binmid' + '\n'\
                                   + str(binstart)+'\t'+ str(binend)+'\t' + str(binsize)+'\t'+ str(binmid)+'\n'\
                                   + str(popt1[0]) + '\t' + str(popt1[1]) + '\t' +str(popt1[2]) + '\t' + str(popt1[3]) + '\n'\
                                   + str(popt2[0]) + '\t' + str(popt2[1]) + '\t' +str(popt2[2]) + '\t' + str(popt2[3]) + '\n')
                    else:
                        print(FILE + '\t' + 'fitting error'+ 'binmid=' + str(binmid))
                        plt.figure(1)
                        plt.hist(heightlist, binrange,color='b')
                        plt.xlabel('height(nm)',fontsize=10)
                        plt.ylabel('counts',fontsize=10)
                        plt.xlim(-0.9 , 1.9)
                        plt.xticks(fontsize=10)
                        plt.yticks(fontsize=10)
                        plt.savefig(outfig)       
                        plt.close(1)
                except RuntimeError:
                    print(FILE + '\t' + 'runtime error')
                    plt.figure(1)
                    plt.hist(heightlist, binrange,color='b')
                    plt.xlabel('height(nm)',fontsize=10)
                    plt.ylabel('counts',fontsize=10)
                    plt.xlim(-0.9 , 1.9)
                    plt.xticks(fontsize=10)
                    plt.yticks(fontsize=10)
                    plt.savefig(outfig)       
                    plt.close(1)
            else:
                print (FILE + '\t' + 'the amount of data is not enough for statistical analysis')



    
                    
                    

        
            