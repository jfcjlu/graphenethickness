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
ampfactor = 1.0e9   # transform the unit to nm

#p001=[0.6265687536355697,3840.434103326893,0.000180485352038041,0.069418342187439]
#p002=[21.224134056723685,2589.1156259708364,0.8616628037743324,0.08105487211651517]


currentfolder = parentFOLDER + '\\' + FILE + '_alllineanlysis'
currentsumfolder = parentFOLDER + '\\' + FILE + '_allline_sum'
summaryfolder = parentFOLDER + '\\' + 'allinessummary'
Psummaryfile = summaryfolder + '\\' 'P-alllines-summaryfile_indivplane.txt' 
outsummary = currentsumfolder + '\\' + FILE + '_alllines_outsummary_indivplane.txt'    
outsummaryexcept = currentsumfolder + '\\' + FILE +'_allines_outsummaryexcept_indivplane.txt'
alllineheightfig =  currentsumfolder + '\\' + FILE +'_allines_outsummary_indivplane.png'
    
if os.path.exists(outsummaryexcept):
    os.remove(outsummaryexcept) 
if os.path.exists(outsummary):
    os.remove(outsummary)

if not os.path.exists(summaryfolder):
    os.mkdir(summaryfolder)


binstart=-10
binmid = 0.3
heightgr = 0.9*binmid
binend= 30 
binsize = 0.02
binrange = np.arange(binstart,binend+binsize, binsize)
plotfitx = np.arange(binstart,binend+binsize, binsize)
binnumber = int((binend-binstart)/binsize)
binlownumber = int((binmid-binstart)/binsize)
bincenter = np.arange(binstart+binsize/2,binend, binsize)

localtime = time.asctime( time.localtime(time.time()) )

with open(Psummaryfile,'a+', encoding='utf-8') as Pallsumf:
       
#    Pallsumf.write('the name of sample'+ '\t' + 'average height (nm)' + '\n')
    Pallsumf.write(str(localtime) + '\n')
    exceptionlines = []
    PATH = parentFOLDER + '\\' + FILE + '.txt'

    if not os.path.exists(currentfolder):
        os.mkdir(currentfolder) 
    if not os.path.exists(currentsumfolder):
        os.mkdir(currentsumfolder)
    
                
    train = pd.read_csv(PATH, sep = '\t', header = 4, index_col = False, engine='python')
    datalist = train.values.tolist() 
    clist= list(chain(*datalist))          
    ONECOLUMEOUT = currentfolder + '\\' + FILE + '_wholeplaneoutto1colume.txt'
    if os.path.exists(ONECOLUMEOUT):
        os.remove(ONECOLUMEOUT)        
    with open(ONECOLUMEOUT,'a+' ,encoding='UTF-8') as outf:
        for line in clist:
            outf.write(str(line) + '\n')        
    with open(ONECOLUMEOUT, 'r' ,encoding='UTF-8') as orgf:
        countline = len(orgf.readlines())        

    
    train = pd.read_csv(PATH, sep = '\t', header = 4, index_col = False, engine='python')
    datalistA = train.values.tolist() 
    clistA= list(chain(*datalistA))
    print('length of clistA' + '\t' + str(len(clistA)))        
    datalistsmall = []        
    with open(PATH) as orgf:
        countline = len(orgf.readlines())        
    for i in range(0,countline-5):
        heightlisti =  datalistA[i][:]
        heightlstlength = len(heightlisti)
        heightlistnm = []
        for j in range(0, heightlstlength-1):
            heightlistnm.append(heightlisti[j]*1.0e9)
        countgr = 0
        for j in range(0, heightlstlength-1):
            if heightlistnm[j] > heightgr : countgr += 1
        if countgr > numofgr:
            datalistsmall.append(heightlisti)        
    clist= list(chain(*datalistsmall))
    print('length of clist' + '\t' + str(len(clist)))
    
    if len(clist) >0:         
        ONECOLUMEOUT = currentfolder + '\\' + FILE + '_wholeplaneoutto1colume.txt'
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
        theestimateend = heightlist[int(0.95*numberofheight)] + 0.95
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
        binlownumber = int((0.9*binmid-binstart)/binsize)
        bincenter = np.arange(binstart+binsize/2,binend, binsize)
        heightgr = 0.9*binmid

        
        
        
        A=np.histogram(heightlist, bins=binrange)      
        histx = bincenter 
        histy = A[0]    
        firsthistx = histx[0:binlownumber-1]
        secondhistx = histx[binlownumber:binnumber-1]
        firsthisty = histy[0:binlownumber-1]
        secondhisty = histy[binlownumber:binnumber-1]
        
        try:            
            popt1, pcov1 = curve_fit(gauss1p, firsthistx, firsthisty, p0=None)
            popt2, pcov2 = curve_fit(gauss1p, secondhistx, secondhisty, p0=None)
            
            p01 = popt1
            p02 = popt2
                
            if popt1[1] > 0 and popt2[1]>0 \
            and (popt1[2] > binstart and popt1[2] < binmid) \
            and (popt2[2] > binmid and popt2[2] < binend):
                bh = popt1[2]
                uph = popt2[2]
                bw = abs(popt1[3])/np.sqrt(2.0)
                upw = abs(popt2[3])/np.sqrt(2.0)
                               
                plotfity1 = popt1[0]+popt1[1]*np.exp(-((plotfitx-popt1[2])/popt1[3])**2)
                plotfity2 = popt2[0]+popt2[1]*np.exp(-((plotfitx-popt2[2])/popt2[3])**2)
                
                outfig = currentfolder + '\\' + FILE + '_wholeplanefit.png'
                plt.figure(1)
                plt.hist(heightlist, binrange,color='b')
                plt.plot(plotfitx, plotfity1,'r')
                plt.plot(plotfitx, plotfity2,'k')
                plt.xlabel('height(nm)',fontsize=10)
                plt.ylabel('counts',fontsize=10)
                plt.xlim(-0.8 , 1.8)
                plt.xticks(fontsize=10)
                plt.yticks(fontsize=10)
                plt.text(0.2*(bh+uph), 0.5*np.max(histy), 'Sample: ' + FILE + 'wholeplane'+ '\n'\
                         +'lowerheight='+ str(round(bh,4)) +' nm'+'\n'\
                         +'upperheight='+ str(round(uph,4)) +' nm'+'\n'\
                         +'height=' + str(round((uph-bh),4)) +' nm'+'\n'\
                         +'lowersigma=' + str(round(bw,4)) +' nm'+'\n'\
                         +'uppersigma=' + str(round(upw,4)) +' nm'+'\n'\
                         ,fontsize=10)
                
        #        plt.show()
                plt.savefig(outfig)       
                plt.close(1)
                
                
                 
                
                with open(PATH) as orgf:
                    countline = len(orgf.readlines())
        
                train = pd.read_csv(PATH, sep = '\t', header = 4, index_col = False, engine='python') 
                datalist = train.values.tolist()
                
                with open(outsummary,'a+', encoding='utf-8') as sumf:
                    localtime = time.asctime( time.localtime(time.time()) )
                    linenumlst = []
                    lowerhtlist = []
                    upperhtlist = []
                    histheightlist = []
                    with open(outsummaryexcept,'a+', encoding='utf-8') as sumexcpf: 
                        sumf.write('linenumber'+'\t'+'bottomheight(nm)'+'\t'+'upheight(nm)'+'\t'+'height(nm)'+\
                                   '\t'+'bottomstddev'+'\t'+'upstddev'+'\t'+'upstddev-bottomstddevM'+'\t'+'bottom_W'+'\t'+'up_W'+'\n')
                        sumexcpf.write('linenumber'+'\t'+'counts'+ 'error'+'\n')
                        
                        for i in range(0,countline-5):
                            heightlist =  datalist[i][:]
                            heightlstlength = len(heightlist)
                            heightlistnm = []
                            for j in range(0, heightlstlength-1):
                                heightlistnm.append(heightlist[j]*1.0e9)
                            countgr = 0
                            A=np.histogram(heightlistnm, bins=binrange)
                            histx = bincenter
                            histy = A[0]
                            uperpeakhistlist = []
                            binlowlimit4uper = int((p02[2]-3.6*abs(p02[3])-binstart)/binsize)
                            binuplimit4uper = int((p02[2]+3.6*abs(p02[3])-binstart)/binsize)
                            uperpeakhistlist = histy[binlowlimit4uper : binuplimit4uper]
                            countgr = np.sum(uperpeakhistlist)
                            
        #                    for j in range(0, uphtlstlength-1):
        #                        if uperpeaklist[j] > heightgr : countgr += 1
                            
                            if countgr > numofgr:
                                A=np.histogram(heightlistnm, bins=binrange)
                                histx = bincenter
                                histy = A[0]
                                outhist = currentfolder + '\\' + FILE +'_'+'line'+str(i)+'_hist.txt'
                                outfig = currentfolder + '\\' + FILE +'_'+'line'+str(i) +'_fit.png'
                                if os.path.exists(outhist):  os.remove(outhist)
                                with open(outhist,'a+', encoding='utf-8') as histf:
                                    for j in range(0, binnumber-1):histf.write(str(histx[j])+'\t'+str(histy[j])+'\n')
                                firsthistx = histx[0:binlownumber-1]
                                secondhistx = histx[binlownumber:binnumber-1]
                                firsthisty = histy[0:binlownumber-1]
                                secondhisty = histy[binlownumber:binnumber-1]
                                try:
                                    popt1, pcov1 = curve_fit(gauss1p, firsthistx, firsthisty, p0 = [0, np.max(firsthisty), p01[2], p01[3]])
                                    popt2, pcov2 = curve_fit(gauss1p, secondhistx, secondhisty,p0 = [0, np.max(secondhisty), p02[2], p02[3]])
                                    if popt1[1] > 0 and popt2[1] >0 and \
                                    (popt2[2] > p02[2]- 1.5*abs(p02[3]) and popt2[2] < p02[2]+ 1.5*abs(p02[3])) and\
                                    (popt1[2] > p01[2]- 1.5*abs(p01[3]) and popt1[2] < p02[2]+ 1.5*abs(p01[3])):
                                        bh = popt1[2]
                                        uph = popt2[2]
                                        bw = abs(popt1[3])/np.sqrt(2.0)
                                        upw = abs(popt2[3])/np.sqrt(2.0)
                                        plotfity1 = popt1[0]+popt1[1]*np.exp(-((plotfitx-popt1[2])/popt1[3])**2)
                                        plotfity2 = popt2[0]+popt2[1]*np.exp(-((plotfitx-popt2[2])/popt2[3])**2)
                                        plt.figure(1)
                                        plt.hist(heightlistnm, binrange,color='b')
                                        plt.plot(histx,histy,color='violet')
                                        plt.plot(plotfitx, plotfity1,'r')
                                        plt.plot(plotfitx, plotfity2,'k')
                                        plt.xlabel('height(nm)',fontsize=10)
                                        plt.ylabel('counts',fontsize=10)
                                        plt.xlim(-0.8 , 1.8)
                                        plt.xticks(fontsize=10)
                                        plt.yticks(fontsize=10)
                                        plt.title('line '+str(i)+' '+'height='+str(uph-bh)+'nm')
                                        plt.savefig(outfig)
                                        plt.close(1)
                                        linenumlst.append(i)
                                        lowerhtlist.append(bh)
                                        upperhtlist.append(uph)
                                        histheightlist.append(round((uph-bh),4))
                                        sumf.write(str(i)+'\t'+ str(round(bh,4)) +'\t'+ str(round(uph,4)) \
                                                   +'\t'+str(round((uph-bh),4))+'\t'+str(round(bw,4))+'\t'\
                                                   +str(round(upw,4))+'\t'+str(round((upw-bw),4))+'\t'\
                                                   +str(round(popt1[3],4))+'\t'+str(round(popt2[3],4))+'\n')
                                        print('line='+str(i)+'\t'+'stepheight='+str(uph-bh))
                                    else:
                                        try:
                                            binlownumberEXCPT = int((p02[2]-3.6*abs(p02[3])-binstart)/binsize)
                                            secondhistxn = histx[binlownumberEXCPT:binnumber-1]
                                            secondhistyn = histy[binlownumberEXCPT:binnumber-1]
                                            popt1, pcov1 = curve_fit(gauss1p, firsthistx, firsthisty, p0 = [0, np.max(firsthisty), p01[2], p01[3]])
                                            popt2, pcov2 = curve_fit(gauss1p, secondhistxn, secondhistyn,p0 = [0, np.max(secondhistyn), p02[2], p02[3]])
                                            if popt1[1] > 0 and popt2[1] >0 and popt2[2] > p02[2]-1.5*abs(p02[3]) and popt2[2] < p02[2]+ 1.5*abs(p02[3]):
                                                bh = popt1[2]
                                                uph = popt2[2]
                                                bw = abs(popt1[3])/np.sqrt(2.0)
                                                upw = abs(popt2[3])/np.sqrt(2.0)
                                                plotfity1 = popt1[0]+popt1[1]*np.exp(-((plotfitx-popt1[2])/popt1[3])**2)
                                                plotfity2 = popt2[0]+popt2[1]*np.exp(-((plotfitx-popt2[2])/popt2[3])**2)
                                                plt.figure(1)
                                                plt.hist(heightlistnm, binrange,color='b')
                                                plt.plot(histx,histy,color='violet')
                                                plt.plot(plotfitx, plotfity1,'r')
                                                plt.plot(plotfitx, plotfity2,'k')
                                                plt.xlabel('height(nm)',fontsize=10)
                                                plt.ylabel('counts',fontsize=10)
                                                plt.xlim(-0.8 , 1.8)
                                                plt.xticks(fontsize=10)
                                                plt.yticks(fontsize=10)
                                                plt.title('line '+str(i)+' '+'height='+str(uph-bh)+'nm')
                                                plt.savefig(outfig)
                                                plt.close(1)
                                                linenumlst.append(i)
                                                lowerhtlist.append(bh)
                                                upperhtlist.append(uph)
                                                histheightlist.append(round((uph-bh),4))
                                                sumf.write(str(i)+'\t'+ str(round(bh,4)) +'\t'+ str(round(uph,4)) \
                                                           +'\t'+str(round((uph-bh),4))+'\t'+str(round(bw,4))+'\t'\
                                                           +str(round(upw,4))+'\t'+str(round((upw-bw),4))+'\t'\
                                                           +str(round(popt1[3],4))+'\t'+str(round(popt2[3],4))+'\n')
                                                print('line='+str(i)+'\t'+'stepheight='+str(uph-bh))
                                            else:
                                                plt.figure(2)
                                                plt.hist(heightlistnm, binrange,color='b')
                                                plt.xlabel('height(nm)',fontsize=10)
                                                plt.ylabel('counts',fontsize=10)
                                                plt.xlim(-0.8 , 1.8)
                                                plt.xticks(fontsize=10)
                                                plt.yticks(fontsize=10)
                                                plt.title('line '+str(i))
                                                plt.savefig(outfig)
                                                plt.close(2)
                                                exceptionlines.append(str(i))
                                                sumexcpf.write(str(i)+'\t'+ str(countgr)+ '\t' + 'fitting error'+'\n')
                                                print('line='+str(i)+'\t'+'countgr='+str(countgr)+ '\t' + 'fitting error')        
                                        except RuntimeError:
                                            A=np.histogram(heightlistnm, bins=binrange)
                                            histx = bincenter
                                            histy = A[0]
                                            outhist = currentfolder + '\\' + FILE +'_'+'line'+str(i)+'_hist.txt'
                                            outfig = currentfolder + '\\' + FILE +'_'+'line'+str(i) +'_fit.png'
                                            if os.path.exists(outhist):  os.remove(outhist)
                                            with open(outhist,'a+', encoding='utf-8') as histf:
                                                for j in range(0, binnumber-1):histf.write(str(histx[j])+'\t'+str(histy[j])+'\n')
                                            plt.figure(2)
                                            plt.hist(heightlistnm, binrange,color='b')
                                            plt.xlabel('height(nm)',fontsize=10)
                                            plt.ylabel('counts',fontsize=10)
                                            plt.xlim(-0.8 , 1.8)
                                            plt.xticks(fontsize=10)
                                            plt.yticks(fontsize=10)
                                            plt.title('line '+str(i))
                                            plt.savefig(outfig)
                                            plt.close(2)
                                            exceptionlines.append(str(i))
                                            sumexcpf.write(str(i)+'\t'+ str(countgr)+ '\t' + 'run time error'+'\n')
                                            print('line='+str(i)+'\t'+'countgr='+str(countgr)+ '\t' + 'run time error')
        #                            else:
        #                                plt.figure(2)
        #                                plt.hist(heightlistnm, binrange,color='b')
        #                                plt.xlabel('height(nm)',fontsize=10)
        #                                plt.ylabel('counts',fontsize=10)
        #                                plt.xlim(-0.8 , 1.8)
        #                                plt.xticks(fontsize=10)
        #                                plt.yticks(fontsize=10)
        #                                plt.title('line '+str(i))
        #                                plt.savefig(outfig)
        #                                plt.close(2)
        #                                sumexcpf.write(str(i)+'\t'+ str(countgr)+ '\t' + 'fitting error'+'\n')
        #                                print('line='+str(i)+'\t'+'countgr='+str(countgr)+ '\t' + 'fitting error')             
                                        
                                except RuntimeError:
                                    A=np.histogram(heightlistnm, bins=binrange)
                                    histx = bincenter
                                    histy = A[0]
                                    outhist = currentfolder + '\\' + FILE +'_'+'line'+str(i)+'_hist.txt'
                                    outfig = currentfolder + '\\' + FILE +'_'+'line'+str(i) +'_fit.png'
                                    if os.path.exists(outhist):  os.remove(outhist)
                                    with open(outhist,'a+', encoding='utf-8') as histf:
                                        for j in range(0, binnumber-1):histf.write(str(histx[j])+'\t'+str(histy[j])+'\n')
                                    plt.figure(2)
                                    plt.hist(heightlistnm, binrange,color='b')
                                    plt.xlabel('height(nm)',fontsize=10)
                                    plt.ylabel('counts',fontsize=10)
                                    plt.xlim(-0.8 , 1.8)
                                    plt.xticks(fontsize=10)
                                    plt.yticks(fontsize=10)
                                    plt.title('line '+str(i))
                                    plt.savefig(outfig)
                                    plt.close(2)
                                    exceptionlines.append(str(i))
                                    sumexcpf.write(str(i)+'\t'+ str(countgr)+ '\t' + 'run time error'+'\n')
                                    print('line='+str(i)+'\t'+'countgr='+str(countgr)+ '\t' + 'run time error')                
                            else:
                                A=np.histogram(heightlistnm, bins=binrange)
                                histx = bincenter
                                histy = A[0]
                                outhist = currentfolder + '\\' + FILE +'_'+'line'+str(i)+'_hist.txt'
                                outfig = currentfolder + '\\' + FILE +'_'+'line'+str(i) +'_fit.png'
                                if os.path.exists(outhist):  os.remove(outhist)
                                with open(outhist,'a+', encoding='utf-8') as histf:
                                    for j in range(0, binnumber-1):histf.write(str(histx[j])+'\t'+str(histy[j])+'\n')
                                plt.figure(3)
                                plt.hist(heightlistnm, binrange,color='b')
                                plt.xlabel('height(nm)',fontsize=10)
                                plt.ylabel('counts',fontsize=10)
                                plt.xlim(-0.8 , 1.8)
                                plt.xticks(fontsize=10)
                                plt.yticks(fontsize=10)
                                plt.title('line '+str(i))
                                plt.savefig(outfig) 
                                plt.close(3)
                                print('line='+str(i)+'\t'+'countgr='+str(countgr))
                plt.figure()
                plt.subplot(131)
                plt.plot(linenumlst, lowerhtlist, 'bs-')
                plt.xlabel('line number',fontsize=20)
                plt.ylabel('bottom height (nm)',fontsize=20)
                plt.xticks(fontsize=20)
                plt.yticks(fontsize=20)
                plt.subplot(132)
                plt.plot(linenumlst, upperhtlist, 'r^-')
                plt.xlabel('line number',fontsize=20)
                plt.ylabel('up height (nm)',fontsize=20)
                plt.xticks(fontsize=20)
                plt.yticks(fontsize=20)
                plt.subplot(133)
                plt.plot(linenumlst, histheightlist, 'ko-')
                plt.xlabel('line number',fontsize=20)
                plt.ylabel('step height (nm)',fontsize=20)
                plt.xticks(fontsize=20)
                plt.yticks(fontsize=20)
                fig=plt.gcf()
                fig.subplots_adjust(left=0.15, bottom=0.15, right=0.9, top=0.9, wspace=0.35, hspace=0.35)
                fig.set_size_inches(20, 10.5)
                plt.show
                plt.savefig(alllineheightfig)
                plt.close()
                Pallsumf.write(FILE + '\t' + str(round(np.mean(histheightlist),4)) + '\t' + '\t'  +'exception lines: ' + str(exceptionlines) +'\n')
#                ALLavgheight.append(round(np.mean(histheightlist),4))
            else:
                print(FILE + '\t' + 'plane Runtime error')
                
        except RuntimeError:
            print(FILE + '\t' + 'plane Runtime error')
    else:
        print(FILE + '\t' + 'not enough data')

            


                

        
            