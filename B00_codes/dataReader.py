import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import time
import os

# read dat file
# datafile = 'C:/Users/lukin2dmaterials/data/2022-10-26/#014_ODMR_CW_02-57-50/ODMRObject_sig_set.dat'
# datafile = 'C:/Users/lukin2dmaterials/data/2022-10-26/#020_ODMR_CW_10-32-32/ODMRObject_sig_set.dat'
# datafile = 'C:/Users/lukin2dmaterials/data/2022-10-26/#024_Rabi_12-32-23/RabiObject_sig_set.dat'
def readData(datafile):
    readfile = np.loadtxt(datafile)
    # print(readfile)
    x_s = [xAxisAndResult[0] for xAxisAndResult in readfile]
    ref = [xAxisAndResult[1] for xAxisAndResult in readfile]
    sig = [xAxisAndResult[2] for xAxisAndResult in readfile]

    fig,ax = plt.subplots()
    ax.plot(x_s, sig, 'o', label='sig', color='C0')
    ax.plot(x_s, ref, 'o', label='ref', color = 'C1')
    ax.legend(loc='best')
    ax.set_xlabel(r"$\tau$ (ns)")

    sig = np.array(sig); ref = np.array(ref); sigOverRef = sig/ref
    fig2,ax = plt.subplots()
    ax.plot(x_s, sigOverRef, 'o', label='sig/ref', color='C0')
    # ax.plot(x_s, ref/sig-np.average(ref/sig), label='sig-ref')
    ax.legend(loc='best')
    ax.set_xlabel(r"$\tau$ (ns)")
    # print(x_s[np.argmin(sigOverRef)])
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
    # datafile = 'C:/Users/lukin2dmaterials/data/2022-11-11/#030_Rabi_19-22-18/RabiObject_sig_set.dat'
    datafile = 'C:/Users/lukin2dmaterials/data/2022-11-20/#001_T2R_02-53-02/T2RObject_sig_set.dat'
    fig2 = readData(datafile)
    plt.show()

    # mainFolder = 'C:/Users/lukin2dmaterials/data/2022-11-06/'
    # for dataFolder in os.listdir(mainFolder):
    #     idx = int(dataFolder[1:4])
    #     if idx in np.linspace(10,11,2):
    #         # print(idx)
    #         datafile = mainFolder + dataFolder +'/RabiObject_sig_set.dat'
    #         # fig = readData(datafile)
    #         fig = readDataSigMinusRef(datafile)
    # plt.show()

    # dataArr = []
    # for dataFolder in os.listdir('C:/Users/lukin2dmaterials/data/2022-10-28/'):
    #     idx = int(dataFolder[1:4])
    #     if idx in np.linspace(37,67,31):
    #         datafile = 'C:/Users/lukin2dmaterials/data/2022-10-28/' + dataFolder +'/RabiObject_sig_set.dat'
    #         x_s, sig, ref = readDataNoPlot(datafile)
    #         sig = np.array(sig); ref = np.array(ref)
    #         contrast = sig/ref - np.mean(sig/ref)
    #         dataArr.append(contrast)
    # dataArr = np.array(dataArr)
    # # print(data)
    # plt.imshow(dataArr, cmap = "inferno")
    # plt.show()
    