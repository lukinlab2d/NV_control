import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import time
import os
from scipy.optimize import curve_fit

NO_MS_EQUALS_1 = 0
Q_FINAL = 1
THREE_PI_HALF_FINAL = 2
REF_MINUS_SIG  = 3

def sinusoid(t, A, Tpi, phi,C):
    return A*np.cos(np.pi/Tpi*t + phi) + C

def sinusoidDecay(t, A, Tpi, phi,C, T2):
    return A*np.cos(np.pi/Tpi*t + phi)*np.exp(-t/T2) + C

def fitSinusoid(xdata, ydata, guess=None):
    lowerBounds = (0,0,-np.pi, -np.inf)
    upperBounds = (np.inf, np.inf, np.pi, np.inf)
    popt, pcov = curve_fit(sinusoid, xdata, ydata, p0=guess, bounds=(lowerBounds, upperBounds))
    perr = np.sqrt(np.diag(pcov))
    xfit = np.linspace(xdata[0], xdata[-1], 1001)
    yfit = sinusoid(xfit, *popt)
    return xfit, yfit, popt, perr

def fitSinusoidDecay(xdata, ydata, guess=None):
    lowerBounds = (0,0,-np.pi, -np.inf,0)
    upperBounds = (np.inf, np.inf, np.pi, np.inf, np.inf)
    popt, pcov = curve_fit(sinusoidDecay, xdata, ydata, p0=guess, bounds=(lowerBounds, upperBounds))
    perr = np.sqrt(np.diag(pcov))
    xfit = np.linspace(xdata[0], xdata[-1], 1001)
    yfit = sinusoidDecay(xfit, *popt)
    return xfit, yfit, popt, perr


def readData(datafile, type=None, typeNorm=0, guess=None):
    readfile = np.loadtxt(datafile)
    
    x_s = [xAxisAndResult[0] for xAxisAndResult in readfile]
    ref = [xAxisAndResult[1] for xAxisAndResult in readfile]
    sig = [xAxisAndResult[2] for xAxisAndResult in readfile]

    x_s = np.array(x_s)
    if type == 'XY8': x_s = 8*x_s

    fig,ax = plt.subplots()
    ax.plot(x_s, sig, 'o-', label='sig', color='C0')
    ax.plot(x_s, ref, 'o-', label='ref', color = 'C1')
    ax.legend(loc='best')
    ax.set_xlabel(r"$\tau$ (ns)")

    sig = np.array(sig); ref = np.array(ref); sigOverRef = sig/ref; 

    fig2,ax = plt.subplots()
    if typeNorm == NO_MS_EQUALS_1:
        ax.plot(x_s, sigOverRef, 'o-', label='sig/ref', color='C0')
    elif typeNorm == Q_FINAL or typeNorm == THREE_PI_HALF_FINAL:
        ax.plot(x_s, (sig-ref)/(sig+ref), 'o-', label='(sig-ref)/(sig+ref)', color='C0')
    elif typeNorm == REF_MINUS_SIG:
        ax.plot(x_s, -(sig-ref)/(sig+ref), 'o-', label='(sig-ref)/(sig+ref)', color='C0')
    ax.legend(loc='best')
    ax.set_xlabel(r"$\tau$ (ns)")

    if type == 'Rabi':
        xfit, yfit, popt, perr = fitSinusoid(x_s, sigOverRef, guess=guess)
        print(popt)
        ax.plot(xfit, yfit, color='C1')
        ax.set_title('$\pi$-pulse = %.2f $\pm$ %.2f ns' % (popt[1], perr[1]))
    
    elif type == 'RabiDecay':
        xfit, yfit, popt, perr = fitSinusoidDecay(x_s, sigOverRef, guess=guess)
        print(popt)
        ax.plot(xfit, yfit, color='C1')
        ax.set_title('$\pi$-pulse = %.2f $\pm$ %.2f ns' % (popt[1], perr[1]))

    
    plt.show()

def readDataSoftAvg(dateFolder, indices, type=None, typeNorm=0):
    sig = []; ref = []
    for idx in indices:
        for root, dirs, files in os.walk(dateFolder):
            for dataFolder in dirs:
                idxString = str(int(idx)) + '_'
                if idxString in str(dataFolder):
                    dataFolderPath = dateFolder + dataFolder + '/'
                    for root2, dirs2, files2 in os.walk(dataFolderPath):
                        for name in files2:
                            if name.endswith((".dat")):
                                datafile = dataFolderPath + name
                                print(datafile)

                                readfile = np.loadtxt(datafile)
                                # print(readfile)
                                x_s = [xAxisAndResult[0] for xAxisAndResult in readfile]
                                ref += [xAxisAndResult[1] for xAxisAndResult in readfile]
                                sig += [xAxisAndResult[2] for xAxisAndResult in readfile]

                                sig = np.array(sig); ref = np.array(ref); sigOverRef = sig/ref; x_s = np.array(x_s)
   
    fig,ax = plt.subplots()
    ax.plot(x_s, sig, 'o-', label='sig', color='C0')
    ax.plot(x_s, ref, 'o-', label='ref', color = 'C1')
    ax.legend(loc='best')
    ax.set_xlabel(r"$\tau$ (ns)")

    
    fig2,ax = plt.subplots()
    if typeNorm == NO_MS_EQUALS_1:
        ax.plot(x_s, sigOverRef, 'o-', label='sig/ref', color='C0')
    elif typeNorm == Q_FINAL or typeNorm == THREE_PI_HALF_FINAL:
        ax.plot(x_s, (sig-ref)/(sig+ref), 'o-', label='(sig-ref)/(sig+ref)', color='C0')
    ax.legend(loc='best')
    ax.set_xlabel(r"$\tau$ (ns)")

    if type == 'Rabi':
        guess=(0.2, 50, 0, 0.9)
        xfit, yfit, popt, perr = fitSinusoid(x_s, sigOverRef, guess=guess)
        print(popt)
        ax.plot(xfit, yfit, color='C1')
        # ax.plot(xfit, sinusoid(xfit, *guess), color='C2')
        ax.set_title('$\pi$-pulse = %.2f $\pm$ %.2f ns' % (popt[1], perr[1]))
    
    plt.show()

def readDataConcatenate(dateFolder, indices, type=None, typeNorm=0):
    x_s_all = []; sig_all = []; ref_all = []
    for idx in indices:
        for root, dirs, files in os.walk(dateFolder):
            for dataFolder in dirs:
                idxString = str(int(idx)) + '_'
                if idxString in str(dataFolder):
                    dataFolderPath = dateFolder + dataFolder + '/'
                    for root2, dirs2, files2 in os.walk(dataFolderPath):
                        for name in files2:
                            if name.endswith((".dat")):
                                datafile = dataFolderPath + name
                                print(datafile)

                                readfile = np.loadtxt(datafile)
                                # print(readfile)
                                x_s = [xAxisAndResult[0] for xAxisAndResult in readfile]
                                ref = [xAxisAndResult[1] for xAxisAndResult in readfile]
                                sig = [xAxisAndResult[2] for xAxisAndResult in readfile]

                                sig = np.array(sig); ref = np.array(ref); sigOverRef = sig/ref; x_s = np.array(x_s)
                                x_s_all = np.concatenate((x_s_all, x_s))
                                sig_all = np.concatenate((sig_all, sig))
                                ref_all = np.concatenate((ref_all, ref))
   
    fig,ax = plt.subplots()
    ax.plot(x_s_all, sig_all, 'o-', label='sig', color='C0')
    ax.plot(x_s_all, ref_all, 'o-', label='ref', color = 'C1')
    ax.legend(loc='best')
    ax.set_xlabel(r"$\tau$ (ns)")

    
    fig2,ax = plt.subplots()
    if typeNorm == NO_MS_EQUALS_1:
        ax.plot(x_s_all, sig_all/ref_all, 'o-', label='sig/ref', color='C0')
    elif typeNorm == Q_FINAL or typeNorm == THREE_PI_HALF_FINAL:
        ax.plot(x_s_all, (sig_all-ref_all)/(sig_all+ref_all), 'o-', label='(sig-ref)/(sig+ref)', color='C0')
    ax.legend(loc='best')
    ax.set_xlabel(r"$\tau$ (ns)")

    plt.show()

def readDataSigMinusRef(datafile):
    readfile = np.loadtxt(datafile)
    # print(readfile)
    x_s = [xAxisAndResult[0] for xAxisAndResult in readfile]
    ref = [xAxisAndResult[1] for xAxisAndResult in readfile]
    sig = [xAxisAndResult[2] for xAxisAndResult in readfile]
    sigOverRef = [xAxisAndResult[3] for xAxisAndResult in readfile]

    # fig,ax = plt.subplots()
    # ax.plot(x_s, sig, label='sig')
    # ax.plot(x_s, ref, label='ref')
    # ax.legend(loc='best')
    # ax.set_xlabel(r"$\tau$ (ns)")
    # ax.set_title(datafile[43:46])

    sig = np.array(sig); ref = np.array(ref)
    fig,ax = plt.subplots()
    ax.plot(x_s, sig-ref, label='sig-ref')
    # ax.plot(x_s, ref/sig-np.average(ref/sig), label='sig-ref')
    ax.legend(loc='best')
    ax.set_xlabel(r"$\tau$ (ns)")
    return fig

def readDataNoPlot(datafile):
    readfile = np.loadtxt(datafile)
    # print(readfile)
    x_s = [xAxisAndResult[0] for xAxisAndResult in readfile]
    ref = [xAxisAndResult[1] for xAxisAndResult in readfile]
    sig = [xAxisAndResult[2] for xAxisAndResult in readfile]
    
    return x_s, sig, ref

if __name__ == '__main__':
    # datafile = 'C:/Users/lukin2dmaterials/data/2022-11-11/#021_ODMR_CW_17-40-40/ODMRObject_sig_set.dat'

    # datafile = 'C:/Users/lukin2dmaterials/data/2022-12-13/#028_Rabi_13-41-39/RabiObject_sig_set.dat'
    # guess=(0.2, 40, 0, 0.9,60)
    # fig2 = readData(datafile, type='RabiDecay', guess=guess)

    # datafile = 'C:/Users/lukin2dmaterials/data/2022-12-01/#014_T2R_18-02-17/T2RObject_sig_set.dat'
    # fig2 = readData(datafile, type='None', typeNorm=1)

    # datafile = 'C:/Users/lukin2dmaterials/data/2022-12-13/#010_T2E_03-56-23/T2EObject_sig_set.dat'
    # fig2 = readData(datafile, type='None', typeNorm=1)

    # datafile = 'C:/Users/lukin2dmaterials/data/2022-12-02/#008_T1_01-45-23/T1Object_sig_set.dat'
    # fig2 = readData(datafile, type='None', typeNorm=0)

    # datafile = 'C:/Users/lukin2dmaterials/data/2022-12-14/#025_XY8_14-16-22/XY8Object_sig_set.dat'
    # fig2 = readData(datafile, type='XY8', typeNorm=1)

    # datafile = 'C:/Users/lukin2dmaterials/data/2022-12-02/#015_CalibrateReadoutLength_12-42-8/CalibrateReadoutLengthObject_sig_set.dat'
    # fig2 = readData(datafile, type='None', typeNorm=0)

    # dateFolder = 'C:/Users/lukin2dmaterials/data/2022-12-02/'
    # indices = np.linspace(15,17,3)
    # readDataSoftAvg(dateFolder, indices, type=None, typeNorm=1)

    # dateFolder = 'C:/Users/lukin2dmaterials/data/2022-12-03/'
    # indices = (2,9)
    # readDataConcatenate(dateFolder, indices, type=None, typeNorm=0)

    # plt.show()

    mainFolder = 'C:/Users/lukin2dmaterials/data/2022-12-19/'
    for dataFolder in os.listdir(mainFolder):
        # print(dataFolder)
        if 'XY8' in dataFolder:
            idx = int(dataFolder[1:4])
            if idx >= 0:
                datafile = mainFolder + dataFolder +'/XY8Object_sig_set.dat'
                fig = readData(datafile, type='XY8', typeNorm=1)

    # mainFolder = 'C:/Users/lukin2dmaterials/data/2022-12-19/'
    # for dataFolder in os.listdir(mainFolder):
    #     # print(dataFolder)
    #     if 'T2E' in dataFolder:
    #         idx = int(dataFolder[1:4])
    #         if idx >= 0:
    #             datafile = mainFolder + dataFolder +'/T2EObject_sig_set.dat'
    #             fig = readData(datafile, type='T2E', typeNorm=1)
    plt.show()

    # #-----------------------------------------------------------------------
    # # Loop

    # dataArr = []; dataArr2 = []; dataArr3 = []
    # for dataFolder in os.listdir('C:/Users/lukin2dmaterials/data/2022-11-29/'):
    #     idx = int(dataFolder[1:4])
    #     if idx in np.linspace(3,9,7):
    #         datafile = 'C:/Users/lukin2dmaterials/data/2022-11-29/' + dataFolder +'/T2RObject_sig_set.dat'
    #         x_s, sig, ref = readDataNoPlot(datafile)
    #         sig = np.array(sig); ref = np.array(ref)
    #         contrast = sig/ref
    #         dataArr.append(contrast)
    #         print(len(contrast))
    #     if idx in np.linspace(10,11,2):
    #         datafile = 'C:/Users/lukin2dmaterials/data/2022-11-29/' + dataFolder +'/T2RObject_sig_set.dat'
    #         x_s2, sig, ref = readDataNoPlot(datafile)
    #         sig = np.array(sig); ref = np.array(ref)
    #         contrast = sig/ref
    #         dataArr2.append(contrast)
    #         print(len(contrast))
    # for dataFolder in os.listdir('C:/Users/lukin2dmaterials/data/2022-11-30/'):
    #     idx = int(dataFolder[1:4])
    #     if idx in np.linspace(1,4,4):
    #         datafile = 'C:/Users/lukin2dmaterials/data/2022-11-30/' + dataFolder +'/T2RObject_sig_set.dat'
    #         x_s3, sig, ref = readDataNoPlot(datafile)
    #         sig = np.array(sig); ref = np.array(ref)
    #         contrast = sig/ref
    #         dataArr3.append(contrast)
    #         print(len(contrast))
    # dataArr = np.array(dataArr)
    # dataArr2 = np.array(dataArr2)
    # dataArr3 = np.array(dataArr3)

    # fig, ax = plt.subplots()
    # y = np.linspace(2867, 2870, 7)
    # X, Y = np.meshgrid(x_s, y)
    # plot = ax.pcolormesh(X,Y, dataArr, cmap='inferno')
    # fig.colorbar(plot, orientation='vertical')

    # y2 = np.linspace(2870.5, 2871, 2)
    # X2, Y2 = np.meshgrid(x_s2, y2)
    # plot2 = ax.pcolormesh(X2,Y2, dataArr2, cmap='inferno')

    # y3 = np.linspace(2871.5, 2873, 4)
    # X3, Y3 = np.meshgrid(x_s3, y3)
    # plot3 = ax.pcolormesh(X3,Y3, dataArr3, cmap='inferno')
    # plt.show()
    