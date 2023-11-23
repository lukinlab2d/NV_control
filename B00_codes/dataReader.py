import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import time
import os
from scipy.optimize import curve_fit
import matplotlib.colors as colors
import math
import scipy
from scipy import integrate
from scipy.special import jn
from scipy.special import iv
import pandas as pd

NO_MS_EQUALS_1 = 0
Q_FINAL = 1
THREE_PI_HALF_FINAL = 2
REF_MINUS_SIG  = 3

def sinusoid(t, A, Tpi, phi,C):
    return A*np.cos(np.pi/Tpi*t + phi) + C
def cosThree(t, A, f1, p1, B,f2,p2, C, f3,p3,D):
    return A*np.cos(2*np.pi*f1*t + p1) + B*np.cos(2*np.pi*f2*t + p2) + C*np.cos(2*np.pi*f3*t + p3) + D
def cosFour(t, A, f1, p1, B,f2,p2, C, f3,p3, D,f4,p4, E):
    return A*np.cos(2*np.pi*f1*t + p1) + B*np.cos(2*np.pi*f2*t + p2) + C*np.cos(2*np.pi*f3*t + p3) + D*np.cos(2*np.pi*f4*t + p4) + E
def linear(x, a,b):
    return a*x+b
def sinusoidDecay(t, A, Tpi, phi,C, T2):
    return A*np.cos(np.pi/Tpi*t + phi)*np.exp(-t/T2) + C
def saturation(P,I0,Ps):
    return I0*(P/Ps)/(1 + P/Ps)
def saturationQuad(P,a,Ps):
    return a*(P**2)/(1 + P/Ps)
def decay(t, A ,C, T2):
    return A*np.exp(-t/T2) + C
def strDecay(t, A, T2, n, C):
    return A*np.exp(-(t/T2)**n) + C
def strDecaySinusoid(t, A, T2, n, Tosc, phi, B, C):
    return A*np.exp(-(t/T2)**n)*np.cos(2*np.pi/Tosc*t + phi) + B*np.exp(-(t/T2)**n) + C
def lor(f, A, f0, g, C):
    return A*g/((f-f0)**2 + (g/2)**2) + C
def lorThree(f, A0, f0, g0, A1, f1, g1, A2, f2, g2, C):
    return A0*g0/((f-f0)**2 + (g0/2)**2) + A1*g1/((f-f1)**2 + (g1/2)**2) + A2*g2/((f-f2)**2 + (g2/2)**2) + C
def lorFour(f, A0, f0, g0, A1, f1, g1, A2, f2, g2, A3, f3, g3, C):
    return A0*g0/((f-f0)**2 + (g0/2)**2) + A1*g1/((f-f1)**2 + (g1/2)**2) + A2*g2/((f-f2)**2 + (g2/2)**2) + A3*g3/((f-f3)**2 + (g3/2)**2) + C
def twoPois(x,A,mu,B,nu):
    if x <= 170:
        return A*np.exp(-mu)*mu**x/scipy.special.factorial(x) + B*np.exp(-nu)*nu**x/scipy.special.factorial(x)
    else: #Stirling approximation
        aterm = A*np.exp(-mu)/np.sqrt(2*np.pi*x) * (mu*np.e/x)**x
        bterm = B*np.exp(-nu)/np.sqrt(2*np.pi*x) * (nu*np.e/x)**x
        return aterm+bterm
def pois(x,A,mu):
    if x <= 170:
        return A*np.exp(-mu)*(mu**(x/2))*(mu**(x/2))/scipy.special.factorial(x)
    else: #Stirling approximation
        return A*np.exp(-mu)/np.sqrt(2*np.pi*x) * (mu*np.e/x)**x
def itgrand_m_even(tau, n, gm, g0, nm, n0, tr):
    A = np.sqrt(g0*gm*tau/(tr-tau))
    B = np.exp(-gm*tau - g0*(tr-tau))
    navg = nm*tau + n0*(tr-tau)
    C = pois(n,1,navg)
    D = iv(1, 2*np.sqrt(g0*gm*tau*(tr-tau)))
    return A*B*C*D
def itgrand_m_odd(tau, n, gm, g0, nm, n0, tr):
    A = gm
    B = np.exp(-gm*tau - g0*(tr-tau))
    navg = nm*tau + n0*(tr-tau)
    C = pois(n,1,navg)
    D = iv(0, 2*np.sqrt(g0*gm*tau*(tr-tau)))
    return A*B*C*D
def itgrand_m_0(n, gm, nm, tr):
    A = np.exp(-gm*tr)
    B = pois(n,1,nm*tr)
    return A*B
def pm(n, gm, g0, nm, n0, tr):
    pm0 = itgrand_m_0(n, gm, nm, tr)
    pm_even, err_pm_even = integrate.quad(itgrand_m_even,0,tr,args=(n, gm, g0, nm, n0, tr))
    pm_odd,  err_pm_odd  = integrate.quad(itgrand_m_odd, 0,tr,args=(n, gm, g0, nm, n0, tr))
    return pm0 + pm_even + pm_odd
def p0(n, gm, g0, nm, n0, tr):
    return pm(n, g0, gm, n0, nm, tr)
def p(n, gm, g0, nm, n0, tr):
    A0 = gm/(g0+gm); Am = g0/(g0+gm)
    return A0*p0(n, gm, g0, nm, n0, tr) + Am*pm(n, gm, g0, nm, n0, tr)
def ps(ns, gm, g0, nm, n0, tr):
    p_arr = []
    for n in ns:
        result = p(n, gm, g0, nm, n0, tr)
        p_arr.append(result)
    return np.array(p_arr)
def p0s(ns, gm, g0, nm, n0, tr):
    p_arr = []
    for n in ns:
        result = p0(n, gm, g0, nm, n0, tr)
        p_arr.append(result)
    return np.array(p_arr)
def pms(ns, gm, g0, nm, n0, tr):
    p_arr = []
    for n in ns:
        result = pm(n, gm, g0, nm, n0, tr)
        p_arr.append(result)
    return np.array(p_arr)

def fitTwoPois(xdata, ydata, guess=None):
    lowerBounds = (0,0,0,0)
    upperBounds = (np.inf, np.inf, np.pi, np.inf)
    # fig, ax = plt.subplots()
    # ax.plot(xdata, twoPois(xdata, *guess))
    # ax.set_xlim((0,xdata[-1]))
    # ax.set_yscale('log')
    # print(ydata)
    popt, pcov = curve_fit(twoPois, xdata, ydata, p0=guess, bounds=(lowerBounds, upperBounds))
    perr = np.sqrt(np.diag(pcov))
    xfit = xdata
    yfit = twoPois(xfit, *popt)
    return xfit, yfit, popt, perr

def fitSinusoid(xdata, ydata, guess=None):
    lowerBounds = (0,0,-np.pi, -np.inf)
    upperBounds = (np.inf, np.inf, np.pi, np.inf)
    popt, pcov = curve_fit(sinusoid, xdata, ydata, p0=guess, bounds=(lowerBounds, upperBounds))
    perr = np.sqrt(np.diag(pcov))
    xfit = np.linspace(xdata[0], xdata[-1], 1001)
    yfit = sinusoid(xfit, *popt)
    return xfit, yfit, popt, perr

def fitCosThree(xdata, ydata, guess=None, lowerBounds=None, upperBounds=None):
    if lowerBounds is None: lowerBounds = (0, 0,-np.pi, 0, 0,-np.pi, 0, 0,-np.pi, -np.inf)
    if upperBounds is None: upperBounds = (1, np.inf, np.pi, 1, np.inf, np.pi, 1, np.inf, np.pi, np.inf)
    popt, pcov = curve_fit(cosThree, xdata, ydata, p0=guess, bounds=(lowerBounds, upperBounds))
    perr = np.sqrt(np.diag(pcov))
    xfit = np.linspace(xdata[0], xdata[-1], 1001)
    yfit = cosThree(xfit, *popt)
    return xfit, yfit, popt, perr

def fitCosFour(xdata, ydata, guess=None, lowerBounds=None, upperBounds=None):
    if lowerBounds is None: lowerBounds = (0, 0,-np.pi, 0, 0,-np.pi, 0, 0,-np.pi, 0, 0,-np.pi,-np.inf)
    if upperBounds is None: upperBounds = (1, np.inf, np.pi, 1, np.inf, np.pi, 1, np.inf, np.pi, 1, np.inf, np.pi, np.inf)
    popt, pcov = curve_fit(cosFour, xdata, ydata, p0=guess, bounds=(lowerBounds, upperBounds))
    perr = np.sqrt(np.diag(pcov))
    xfit = np.linspace(xdata[0], xdata[-1], 1001)
    yfit = cosFour(xfit, *popt)
    return xfit, yfit, popt, perr

def fitSinusoidDecay(xdata, ydata, guess=None):
    lowerBounds = (0,0,-np.pi, -np.inf,0)
    upperBounds = (np.inf, np.inf, np.pi, np.inf, np.inf)
    popt, pcov = curve_fit(sinusoidDecay, xdata, ydata, p0=guess, bounds=(lowerBounds, upperBounds))
    perr = np.sqrt(np.diag(pcov))
    xfit = np.linspace(xdata[0], xdata[-1], 1001)
    yfit = sinusoidDecay(xfit, *popt)
    return xfit, yfit, popt, perr

def fitDecay(xdata, ydata, guess=None):
    lowerBounds = (-np.inf, -np.inf,0)
    upperBounds = (np.inf, np.inf, 1e8)
    popt, pcov = curve_fit(decay, xdata, ydata, p0=guess, bounds=(lowerBounds, upperBounds))
    perr = np.sqrt(np.diag(pcov))
    xfit = np.linspace(xdata[0], xdata[-1], 1001)
    yfit = decay(xfit, *popt)
    return xfit, yfit, popt, perr

def fitSaturation(xdata, ydata, guess=None):
    lowerBounds = (0,0)
    upperBounds = (np.inf, np.inf)
    popt, pcov = curve_fit(saturation, xdata, ydata, p0=guess, bounds=(lowerBounds, upperBounds))
    perr = np.sqrt(np.diag(pcov))
    xfit = np.linspace(xdata[0], xdata[-1], 1001)
    yfit = saturation(xfit, *popt)
    return xfit, yfit, popt, perr

def fitSaturationQuad(xdata, ydata, guess=None):
    lowerBounds = (0,0)
    upperBounds = (np.inf, np.inf)
    popt, pcov = curve_fit(saturationQuad, xdata, ydata, p0=guess, bounds=(lowerBounds, upperBounds))
    perr = np.sqrt(np.diag(pcov))
    xfit = np.linspace(xdata[0], xdata[-1], 1001)
    yfit = saturationQuad(xfit, *popt)
    return xfit, yfit, popt, perr

def fitLinear(xdata, ydata, guess=None):
    popt, pcov = curve_fit(linear, xdata, ydata, p0=guess)
    perr = np.sqrt(np.diag(pcov))
    xfit = np.linspace(xdata[0], xdata[-1], 1001)
    yfit = linear(xfit, *popt)
    return xfit, yfit, popt, perr

def fitStrDecay(xdata, ydata, guess=None, upperBounds=None, lowerBounds=None):
    if lowerBounds is None: lowerBounds = (-np.pi, 0, 0, -np.inf)
    if upperBounds is None: upperBounds = (np.inf, np.inf, 5, np.inf)
    popt, pcov = curve_fit(strDecay, xdata, ydata, p0=guess, bounds=(lowerBounds, upperBounds))
    perr = np.sqrt(np.diag(pcov))

    xfit = np.linspace(xdata[0], xdata[-1], 1001)
    yfit = strDecay(xfit, *popt)
    return xfit, yfit, popt, perr

def fitStrDecaySinusoid(xdata, ydata, guess=None, upperBounds=None, lowerBounds=None):
    if lowerBounds is None: lowerBounds = (-np.pi, 0, 0,      0,      -np.inf, -np.inf, -np.inf)
    if upperBounds is None: upperBounds = (np.inf, np.inf, 5, np.inf, np.inf,  np.inf, np.inf)
    popt, pcov = curve_fit(strDecaySinusoid, xdata, ydata, p0=guess, bounds=(lowerBounds, upperBounds))
    perr = np.sqrt(np.diag(pcov))

    xfit = np.linspace(xdata[0], xdata[-1], 1001)
    yfit = strDecaySinusoid(xfit, *popt)
    return xfit, yfit, popt, perr

def fitLor(xdata, ydata, guess=None, upperBounds=None, lowerBounds=None):
    if lowerBounds is None: lowerBounds = (-np.inf,0,0,0)
    if upperBounds is None: upperBounds = (0, np.inf, np.inf, np.inf)
    popt, pcov = curve_fit(lor, xdata, ydata, p0=guess, bounds=(lowerBounds, upperBounds))
    perr = np.sqrt(np.diag(pcov))

    xfit = np.linspace(xdata[0], xdata[-1], 1001)
    yfit = lor(xfit, *popt)
    return xfit, yfit, popt, perr

def fitLorThree(xdata, ydata, guess=None, upperBounds=None, lowerBounds=None):
    if lowerBounds is None: lowerBounds = (-np.inf,0,0,-np.inf,0,0,-np.inf,0,0,0)
    if upperBounds is None: upperBounds = (0, np.inf, np.inf, 0, np.inf, np.inf,0, np.inf, np.inf,np.inf)
    popt, pcov = curve_fit(lorThree, xdata, ydata, p0=guess, bounds=(lowerBounds, upperBounds))
    perr = np.sqrt(np.diag(pcov))

    xfit = np.linspace(xdata[0], xdata[-1], 1001)
    yfit = lorThree(xfit, *popt)
    return xfit, yfit, popt, perr

def fitLorFour(xdata, ydata, guess=None, upperBounds=None, lowerBounds=None):
    if lowerBounds is None: lowerBounds = (-np.inf,0,0,  -np.inf,0,0,  -np.inf,0,0,  0)
    if upperBounds is None: upperBounds = (0, np.inf, np.inf,   0, np.inf, np.inf,  0, np.inf, np.inf,  0, np.inf, np.inf,  np.inf)
    popt, pcov = curve_fit(lorFour, xdata, ydata, p0=guess, bounds=(lowerBounds, upperBounds))
    perr = np.sqrt(np.diag(pcov))

    xfit = np.linspace(xdata[0], xdata[-1], 1001)
    yfit = lorFour(xfit, *popt)
    return xfit, yfit, popt, perr

def fitBlinkTwoPois(xdata, ydata, guess=None):
    lowerBounds = (0,0,0,0,guess[-1]*0.99999)
    upperBounds = (np.inf, np.inf, np.inf, np.inf, guess[-1]*1.00001)
    popt, pcov = curve_fit(ps, xdata, ydata, p0=guess, bounds=(lowerBounds, upperBounds))
    perr = np.sqrt(np.diag(pcov))
    xfit = xdata
    yfit = ps(xfit, *popt)
    return xfit, yfit, popt, perr





def readDataFullData(datafile, num_of_bins=10, binwidth=0, plot_hist_every=5, 
                     ifDataSavedAsCountRate=False, ifLogColor=False, ifSubtractRef=False, ifPlotRef=False,
                     ifRealTimeMonitor=False, ifPlot=True, num_of_iter=0):
    readfile = np.loadtxt(datafile)
    # print(readfile)
    tau = [xAxisAndResult[0] for xAxisAndResult in readfile]
    ref = [xAxisAndResult[2] for xAxisAndResult in readfile]
    sig = [xAxisAndResult[4] for xAxisAndResult in readfile]

    tau = np.array(tau); sig = np.array(sig); ref = np.array(ref)
    if ifRealTimeMonitor: 
        num_of_iter_same_tau = num_of_iter
    else: 
        num_of_iter_same_tau = np.count_nonzero(tau == np.min(tau))
    num_of_different_tau = int(len(tau)/num_of_iter_same_tau)

    sig = np.reshape(sig,(num_of_different_tau, num_of_iter_same_tau))
    ref = np.reshape(ref,(num_of_different_tau, num_of_iter_same_tau))

    if math.isnan(np.max(sig)): sig[-1] = np.zeros(num_of_iter_same_tau)
    if math.isnan(np.max(ref)): ref[-1] = np.zeros(num_of_iter_same_tau)

    if ifSubtractRef: sig = sig-ref

    xplot = np.linspace(1,num_of_iter_same_tau,num_of_iter_same_tau)
    yplot = tau[::num_of_iter_same_tau]
    Xplot, Yplot = np.meshgrid(xplot,yplot)

    if ifDataSavedAsCountRate:
        for i in range(num_of_different_tau):
            sig[i] = sig[i]*1e3*yplot[i]/1e9
            ref[i] = sig[i]*1e3*yplot[i]/1e9

    if ifPlot:
        fig,ax = plt.subplots()
        plot = ax.pcolormesh(Xplot, Yplot/1e3, sig, cmap='inferno')
        ax.set_xlabel(r"Iterations")
        ax.set_ylabel(r"Readout time ($\mu$s)")
        ax.set_title("Signal")
        fig.colorbar(plot, orientation='vertical')

    if ifPlotRef:
        fig,ax = plt.subplots()
        plot = ax.pcolormesh(Xplot, Yplot/1e3, ref, cmap='inferno')
        ax.set_xlabel(r"Iterations")
        ax.set_ylabel(r"Readout time ($\mu$s)")
        ax.set_title("Reference")
        fig.colorbar(plot, orientation='vertical')

    if binwidth == 0: num_of_bins = num_of_bins
    else: num_of_bins = int((np.max(sig)-np.min(sig))/binwidth)
    # print("num of bins = " + str(num_of_bins))
    bins1DArray = np.linspace(np.min(sig), np.max(sig), num_of_bins+1)
    hist2DArray_sig = np.zeros((num_of_different_tau,num_of_bins))
    hist2DArray_ref = np.zeros((num_of_different_tau,num_of_bins))
    

    for i in range(num_of_different_tau):
        # Create a histogram
        hist, bins = np.histogram(sig[i], bins=bins1DArray)
        hist2DArray_sig[i] = hist
        if np.max(hist) == num_of_iter_same_tau: hist2DArray_sig[i] = np.zeros(len(hist)) # get rid of unfinished line scan

        hist, bins = np.histogram(ref[i], bins=bins1DArray)
        hist2DArray_ref[i] = hist
        # if np.max(hist) == num_of_iter_same_tau: hist2DArray_ref[i] = np.zeros(len(hist)) # get rid of unfinished line scan

        if np.mod(i, plot_hist_every) == 0:
            if ifPlot:
                fig,ax = plt.subplots(figsize=(2.5,2.5))
                ax.hist(sig[i], bins=bins1DArray, edgecolor='black', density=True)
                ax.set_xlabel(r"Count")
                ax.set_ylabel(r"Probability")
                if ifSubtractRef: text = " ms. Sig-ref. "
                else: text = " ms. Signal. "
                ax.set_title(r"$\tau$ = " + str(np.round(yplot[i]/1e6,2)) + text + str(np.round(num_of_iter_same_tau/1e5,1)) + "e5 iters", 
                            fontsize=10)
            if ifPlotRef:
                fig,ax = plt.subplots(figsize=(2.5,2.5))
                ax.hist(ref[i], bins=bins1DArray, edgecolor='black', density=True)
                ax.set_xlabel(r"Count")
                ax.set_ylabel(r"Probability")
                ax.set_title(r"$\tau$ = " + str(np.round(yplot[i]/1e6,2)) + " ms. Ref. " + str(np.round(num_of_iter_same_tau/1e5,1)) + "e5 iters", 
                            fontsize=10)

    # Plot the histogram
    binsMidPoint = (bins1DArray[0:-1] + bins1DArray[1:])/2
    binPlot, Yplot = np.meshgrid(binsMidPoint,yplot)
    norm = colors.LogNorm(vmin=hist2DArray_sig.min()+0.01, vmax=hist2DArray_sig.max()) # Create a logarithmic color scale

    if ifPlot:
        fig,ax = plt.subplots()
        if ifLogColor: plot = ax.pcolormesh(binPlot, Yplot/1e3, hist2DArray_sig, cmap='inferno', norm=norm)
        else: plot = ax.pcolormesh(binPlot, Yplot/1e3, hist2DArray_sig, cmap='inferno')
        ax.set_xlabel(r"Count")
        ax.set_ylabel(r"Readout time ($\mu$s)")
        if ifSubtractRef: ax.set_title("Histogram of counts - sig-ref")
        else: ax.set_title("Histogram of counts - sig")
        fig.colorbar(plot, orientation='vertical')

    return sig, ref, hist2DArray_sig, hist2DArray_ref, xplot, yplot, bins1DArray

def readDataLiveCounter(datafile, acqTimeMs=1, figsize=(4,4),
                     correctionFactor=1, ifLogColor=False, ifSubtractRef=False):
    readfile = np.loadtxt(datafile)
   
    iter = [xAxisAndResult[0] for xAxisAndResult in readfile]
    sig = [xAxisAndResult[1] for xAxisAndResult in readfile]

    iter = np.array(iter); sig = np.array(sig); 

    fig, ax = plt.subplots(figsize=figsize)
    ax.plot(iter*acqTimeMs/1e3*correctionFactor,sig)
    ax.set_xlabel(r"Time (s)")
    ax.set_ylabel(r"Count rate (c/s)")
    ax.set_title("Continuous count monitoring. Acq time = " + str(acqTimeMs) + " ms")

    return iter, sig

def readData(datafile, type=None, typeNorm=0, guess=None, ifPlot=True, ifPrint=True, ifSinusoid=False,
             ifFit=False, upperBounds=None, lowerBounds=None, endDataPoint=None, startDataPoint=None):
    readfile = np.loadtxt(datafile)
    x_s = [xAxisAndResult[0] for xAxisAndResult in readfile]
    ref = [xAxisAndResult[1] for xAxisAndResult in readfile]
    sig = [xAxisAndResult[2] for xAxisAndResult in readfile]
    sig = np.array(sig); ref = np.array(ref); sigOverRef = sig/ref; 

    x_s = np.array(x_s)
    if type == 'XY8': x_s = 8*x_s
    if type == 'T2E' or type == 'XY8': 
        x_s = x_s/1e3; 
        if endDataPoint != None: 
            x_s = x_s[0:endDataPoint]; ref = ref[0:endDataPoint]; sig = sig[0:endDataPoint]
        if startDataPoint != None: 
            x_s = x_s[startDataPoint:]; ref = ref[startDataPoint:]; sig = sig[startDataPoint:]
    if ifPlot:
        fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(12, 4))
        
        axs[0].plot(x_s, sig, 'o-', markersize=3, label='sig', color='C0')
        axs[0].plot(x_s, ref, 'o-', markersize=3, label='ref', color='C1')
        axs[0].legend(loc='best')
        axs[0].set_xlabel(r"$\tau$ (ns)")
        axs[0].set_ylabel("PL count rate (kc/s)")
        if type == 'T2E' or type == 'XY8': axs[0].set_xlabel(r"$\tau$ ($\mu$s)")
        axs[0].set_title(datafile.split('/')[4] + "/" + datafile.split('/')[5] + " - raw")
        
        if typeNorm == NO_MS_EQUALS_1:
            axs[1].plot(x_s, sig/ref, 'o-',markersize=3,  label='sig/ref', color='C0')
            axs[1].set_ylabel('sig/ref')
        elif typeNorm == THREE_PI_HALF_FINAL or typeNorm == Q_FINAL:
            axs[1].plot(x_s, (sig-ref)/(sig+ref), 'o-', markersize=3, label='(sig-ref)/(sig+ref)', color='C0')
            axs[1].set_ylabel('(sig-ref)/(sig+ref)')
        else:
            axs[1].plot(x_s, -(sig-ref)/(sig+ref), 'o-', markersize=3, label='(sig-ref)/(sig+ref)', color='C0')
            axs[1].set_ylabel('(sig-ref)/(sig+ref)')
        axs[1].legend(loc='best')
        axs[1].set_xlabel(r"$\tau$ (ns)")
        if type == 'T2E' or type == 'XY8': axs[1].set_xlabel(r"$\tau$ ($\mu$s)")
        axs[1].set_title(datafile.split('/')[4] + "/" + datafile.split('/')[5] + " - normalized")

    if type in ['Rabi', 'RabiDecay']and ifFit:
        fitFunc = fitSinusoid if type=='Rabi' else fitSinusoidDecay

        if typeNorm == NO_MS_EQUALS_1: y = sigOverRef
        else: y = np.abs((sig-ref)/(sig+ref))

        xfit, yfit, popt, perr = fitFunc(x_s, y, guess=guess)
        if ifPrint: print(popt)
    
        if ifPlot:
            axs[1].plot(xfit, yfit, color='C1')
            axs[1].set_title('$\pi$-pulse = %.2f $\pm$ %.2f ns' % (popt[1], perr[1]))
        
    if type in ['T2E'] and ifFit:
        if ifSinusoid: fitFunc = fitStrDecaySinusoid
        else: fitFunc = fitStrDecay

        if typeNorm == NO_MS_EQUALS_1: y = sigOverRef
        else: y = np.abs((sig-ref)/(sig+ref))

        xfit, yfit, popt, perr = fitFunc(x_s, y, guess=guess, upperBounds=upperBounds, lowerBounds=lowerBounds)
        if ifPrint: print(popt)
        if ifPlot:
            axs[1].plot(xfit, yfit, color='C1')
            axs[1].set_title(r'$T_{2E}$ = %.2f $\pm$ %.2f $\mu$s; $n$ = %.2f $\pm$ %.2f' % (popt[1], perr[1], popt[2], perr[2]))
    
    if type in ['ODMR'] and ifFit:
        fitFunc = fitLor

        if typeNorm == NO_MS_EQUALS_1: y = sigOverRef
        else: y = np.abs((sig-ref)/(sig+ref))

        xfit, yfit, popt, perr = fitFunc(x_s, y, guess=guess, upperBounds=upperBounds, lowerBounds=lowerBounds)
        if ifPrint: print(popt)
        if ifPlot:
            axs[1].plot(xfit, yfit, color='C1')
            axs[1].set_title("$f$ = %.2f $\pm$ %.2f MHz" % (popt[1]/1e6, perr[1]/1e6))
            axs[0].set_xlabel('f (GHz)')
            axs[1].set_xlabel('f (GHz)')

    if type in ['T1'] and ifFit:
        fitFunc = fitDecay

        if typeNorm == NO_MS_EQUALS_1: y = sigOverRef
        else: y = np.abs((sig-ref)/(sig+ref)); a = np.max(y)
        guess = (a,0,1e6)

        xfit, yfit, popt, perr = fitFunc(x_s, y, guess=guess)
        # if ifPrint: print(popt)
        if ifPlot:
            axs[1].plot(xfit, yfit, color='C1')
            s = "$T_{1}$=%.2f$\pm$%.2f ms" % (popt[2]/1e6, perr[2]/1e6)
            axs[1].set_title(s)
            axs[1].set_xscale('log')
            axs[0].set_xscale('log')
    
    
    plt.show()
    if not ifFit: popt=(1,1,1); perr=(1,1,1)    
    return sig, ref, popt, perr, x_s

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
def readDataNoPlotDual(datafile):
    readfile = np.loadtxt(datafile)
    
    x_s = [xAxisAndResult[0] for xAxisAndResult in readfile]
    ref = [xAxisAndResult[1] for xAxisAndResult in readfile]
    sig = [xAxisAndResult[3] for xAxisAndResult in readfile]
    sig2 = [xAxisAndResult[4] for xAxisAndResult in readfile]
    ref2 = [xAxisAndResult[2] for xAxisAndResult in readfile]
    
    return x_s, sig, ref, sig2, ref2

def readDataNoRef(datafile):
    readfile = np.loadtxt(datafile)
    # print(readfile)
    x_s = [xAxisAndResult[0] for xAxisAndResult in readfile]
    sig = [xAxisAndResult[1] for xAxisAndResult in readfile]
    
    return x_s, sig

def plotHistSweepTIon(sigs, tausArray, sweepWhat=None, ifPlot=1, ms=-1, power589 = 2, power532 = 1400, power635 = 9.5,
                t532 = 500e3, delay1 = 20e6, delay2 = 20, tsh = 100, delay3 = 600, ti=200, delay4 = 5e3, tr_ns = 250e6, finalDataFolder=None,
                gm=1.5,g0=2,nm=240,n0=90, ifLogPlot=0):
    ths = []; fids = []; pNVms = []; snrs = []; gms = []; g0s = []; nms = []; n0s = []; nMean0s = []; nMeanms = [];
    for i in range(len(tausArray)):
        if 'SCCPhStatSweepDelaySI' in finalDataFolder: 
            delay3=tausArray[i]
        elif sweepWhat =='tsh':
            tsh = tausArray[i]
        else:
            ti = tausArray[i]
        yPlot = sigs[i]
        acqTimeMs = tausArray[i]/1e6; binTimeMs = acqTimeMs
        binOverAcq = int(np.round(binTimeMs/acqTimeMs))

        reshaped_y = yPlot.reshape(-1, binOverAcq)
        summed_y = np.sum(reshaped_y, axis=1)

        nbins = int(np.max(summed_y) + 1)
        bins = np.arange(nbins)
        
        z = np.histogram(summed_y, bins=bins, density=True)
        p = z[0]; counts = bins[0:-1]; c = counts; pdata = p

        tr = tr_ns/1e9
        guess=(gm, g0, nm, n0, tr)

        start_idx = 0; end_idx = -1
        # if i == 5: start_idx = 0
        # p = p[start_idx:end_idx]; counts = counts[start_idx:end_idx]; counts=np.linspace(start_idx,counts[-1], len(counts))
        # counts0 = counts[0]; counts = counts-counts0

        fig, ax = plt.subplots(figsize=(5,3))
        ns = np.linspace(0,len(p)-1,len(p))
        pguess = ps(ns, gm, g0, nm, n0, tr)
        # ax.plot(ns,pguess, '--', label="guess")
        ax.plot(c, pdata, 'o', markersize=3, label="Data")
        xfit, yfit, popt, perr = fitBlinkTwoPois(c, pdata, guess=guess)
        ax.plot(xfit,yfit, label="Total")
        ax.set_xlabel(r"Count")
        ax.set_ylabel(r"Probability")

        (gm, g0, nm, n0, tr) = popt
        fitQual = np.sum(yfit)
        A0 = gm/(g0+gm); Am = g0/(g0+gm) 
        pfit0 = A0*p0s(c, *popt); pfit1 = Am*pms(c, *popt)
        pNV0 = np.sum(pfit0); pNVm = np.sum(pfit1)
        nMean0 = np.sum(c*pfit0)/pNV0; nMeanm = np.sum(c*pfit1)/pNVm
        ax.plot(xfit, pfit0, label="NV$^{0}$")
        ax.plot(xfit, pfit1, label="NV$^{-}$")

        fidel_arr = []; snr_arr = []
        for c in xfit:
            c = int(c-start_idx)
            p00_arr = pfit0[:c]; p10_arr = pfit0[c:-1]; p0 = np.sum(pfit0)  
            p01_arr = pfit1[:c]; p11_arr = pfit1[c:-1]; p1 = np.sum(pfit1)
            p00 = np.sum(p00_arr)/p0; p10 = np.sum(p10_arr)/p0; p01 = np.sum(p01_arr)/p1; p11 = np.sum(p11_arr)/p1
            fidel = 1 - 0.5*(p10 + p01) # same as fidel = 0.5*(p00 + p11)
            snr = (p00-p01)/np.sqrt(p00*p10 + p01*p11)
            if np.isnan(snr) or np.isinf(snr): snr = 0
            fidel_arr.append(fidel); snr_arr.append(snr)
        fidel_arr = np.array(fidel_arr); snr_arr = np.array(snr_arr)
        th_idx = np.argmax(fidel_arr); thres = th_idx + start_idx; 
        fidel = fidel_arr[th_idx]; snr = snr_arr[th_idx]
        sigma = np.sqrt(1+2/snr**2)

        s1 = r"$\tau$ = %.1f ms. $m_s$ = %.0f. " % (tr_ns/1e6, ms)
        s2 = "$g_{-0}$ = %.2f Hz. $g_{0-}$ = %.2f Hz. $\gamma_{-}$ = %.0f Hz. $\gamma_{0}$ = %.0f Hz" % (gm, g0, nm, n0)
        s3 = "$P_{589}$ = %.1f $\mu$W. $P_{532}$ = %.0f $\mu$W. $P_{635}$ = %.1f mW$. t_{532}$ = %.0f $\mu$s. $t_{in-MW}$ = %.0f $\mu$s" % (power589, power532, power635, t532/1e3, delay1/1e3)
        s4 = "$t_{MW-sh}$ = %.0f ns. $t_{sh}$ = %.0f ns. $t_{sh-i}$ = %.0f ns. $t_{i}$ = %.0f ns. $t_{i-r}$ = %.1f ms" % (delay2, tsh, delay3, ti, delay4/1e6)
        s5 = "Fit qual=%.3f. Thres=%.0f. Fidel=%.3f. $p_{NV^-}$=%.2f. SNR=%.2f. $\sigma_R$=%.2f"  % (fitQual, thres, fidel, pNVm, snr, sigma)
        s6 = "$n_{0, avg}$ = %.1f. $n_{-, avg}$ = %.1f" % (nMean0, nMeanm)
        ax.set_title(s1  + s2 + "\n" + s3 + "\n" + s4 + "\n" + s5 + "\n" + s6, fontsize=7);

        if ifLogPlot:
            ax.set_yscale('log')
        ax.set_ylim((1e-5,1.05*max(np.max(yfit),np.max(pdata))));
        ax.legend()
        ax.axvline(x=thres, ymin=0, ymax=1, color='black', linestyle='dashed', linewidth=1)
        if ifPlot:
            plt.show()
        else:
            plt.ioff(); plt.clf(); plt.close('all')
        ths.append(thres); fids.append(fidel); pNVms.append(pNVm); snrs.append(snr)
        gms.append(gm); g0s.append(g0); nms.append(nm); n0s.append(n0)
        nMeanms.append(nMeanm); nMean0s.append(nMean0)

        plt.tight_layout()

    ths = np.array(ths); fids = np.array(fids); pNVms = np.array(pNVms); snrs = np.array(snrs); tis = tausArray
    gms = np.array(gms); g0s = np.array(g0s); nms = np.array(nms); n0s = np.array(n0s)
    nMeanms = np.array(nMeanms); nMean0s = np.array(nMean0s)

    data={'tis':tis, 'ths':ths, 'fids':fids, 'pms':pNVms, 'snrs':snrs, 
          'gms':gms, 'g0s':g0s, 'nms':nms,   'n0s':n0s,   'nMeanms':nMeanms, 'nMean0s':nMean0s}
    df = pd.DataFrame(data)
    filename = finalDataFolder + "/fitParams_ms" + str(ms) + ".csv"
    df.to_csv(filename, index=False)

    
    return ths, fids, pNVms, snrs, gms, g0s, nms, n0s, nMeanms, nMean0s

def plotHistSweepTIonNoFit(sigs, tausArray, sweepWhat=None, ifPlot=1, ms=-1, power589 = 2, power532 = 1400, power635 = 9.5,
                t532 = 500e3, delay1 = 20e6, delay2 = 20, tsh = 100, delay3 = 600, ti=200, delay4 = 5e3, tr_ns = 250e6, finalDataFolder=None,
                gm=1.5,g0=2,nm=240,n0=90, ifLogPlot=0):
    for i in range(len(tausArray)):
        if 'SCCPhStatSweepDelaySI' in finalDataFolder: 
            delay3=tausArray[i]
        elif sweepWhat =='tsh':
            tsh = tausArray[i]
        else:
            ti = tausArray[i]
        yPlot = sigs[i]
        acqTimeMs = tausArray[i]/1e6; binTimeMs = acqTimeMs
        binOverAcq = int(np.round(binTimeMs/acqTimeMs))

        reshaped_y = yPlot.reshape(-1, binOverAcq)
        summed_y = np.sum(reshaped_y, axis=1)

        nbins = int(np.max(summed_y) + 1)
        bins = np.arange(nbins)
        
        z = np.histogram(summed_y, bins=bins, density=True)
        p = z[0]; counts = bins[0:-1]; c = counts; pdata = p

        tr = tr_ns/1e9
        guess=(gm, g0, nm, n0, tr)

        start_idx = 0; end_idx = -1
        # if i == 5: start_idx = 0
        # p = p[start_idx:end_idx]; counts = counts[start_idx:end_idx]; counts=np.linspace(start_idx,counts[-1], len(counts))
        # counts0 = counts[0]; counts = counts-counts0

        fig, ax = plt.subplots(figsize=(5,3))
        ns = np.linspace(0,len(p)-1,len(p))
        ax.plot(c, pdata, 'o', markersize=3, label="Data")
        ax.set_xlabel(r"Count")
        ax.set_ylabel(r"Probability")

        s1 = r"$\tau$ = %.1f ms. $m_s$ = %.0f. " % (tr_ns/1e6, ms)
        s3 = "$P_{589}$ = %.1f $\mu$W. $P_{532}$ = %.0f $\mu$W. $P_{635}$ = %.1f mW$. t_{532}$ = %.0f $\mu$s. $t_{in-MW}$ = %.0f $\mu$s" % (power589, power532, power635, t532/1e3, delay1/1e3)
        s4 = "$t_{MW-sh}$ = %.0f ns. $t_{sh}$ = %.0f ns. $t_{sh-i}$ = %.0f ns. $t_{i}$ = %.0f ns. $t_{i-r}$ = %.1f ms" % (delay2, tsh, delay3, ti, delay4/1e6)
        ax.set_title(s1  + "\n" + s3 + "\n" + s4, fontsize=7);

        if ifLogPlot:
            ax.set_yscale('log')
        
        if ifPlot:
            plt.show()
        else:
            plt.ioff(); plt.clf(); plt.close('all')
        
        plt.tight_layout()
    
def plotAnalysisSweepTIon(finalDataFolder, ifPlot=1, sweepWhat=None, power589 = 2, power532 = 1400, power635 = 9.5,
                t532 = 500e3, delay1 = 20e6, delay2 = 20, tsh = 100, delay3 = 600, ti=160, delay4 = 5e3, tr_ns = 250e6,
                ifYlim=0, ylim=None):
    ms=-1
    paramFile = finalDataFolder + "/fitParams_ms" + str(ms) + ".csv"
    df = pd.read_csv(paramFile)
    tis, ths, fids, pms, snrs, gms, g0s, nms, n0s, nMeanms, nMean0s = df['tis'].values, df['ths'].values, df['fids'].values, df['pms'].values, df['snrs'].values, df['gms'].values, df['g0s'].values, df['nms'].values, df['n0s'].values, df['nMeanms'].values, df['nMean0s'].values

    ms=0
    paramFile = finalDataFolder + "/fitParams_ms" + str(ms) + ".csv"
    df = pd.read_csv(paramFile)
    tis, thsref, fidsref, pmsref, snrsref, gmsref, g0sref, nmsref, n0sref, nMeanmsref, nMean0sref = df['tis'].values, df['ths'].values, df['fids'].values, df['pms'].values, df['snrs'].values, df['gms'].values, df['g0s'].values, df['nms'].values, df['n0s'].values, df['nMeanms'].values, df['nMean0s'].values

    b0 = pmsref; b1 = pms
    for i, bb1 in enumerate(b1):
        if b0[i]>bb1: b0[i] = 0.9999*bb1
    sigmaSCC = np.sqrt((b0+b1)*(2-b0-b1)/(b0-b1)**2)

    fig,axs = plt.subplots(4,2, figsize=(6,8))
    a00 = axs[0,0]; a01 = axs[0,1]; a10 = axs[1,0]; a11 = axs[1,1]
    a20 = axs[2,0]; a21 = axs[2,1]; a30 = axs[3,0]; a31 = axs[3,1]

    a00.plot(tis,ths,    'o-', markersize=3, label="$m_s$=-1")
    a00.plot(tis,thsref, 'o-', markersize=3, label="$m_s$=0")
    a00.set_ylabel("Threshold")
    a00.legend(fontsize=7)

    a01.plot(tis,pms,    'o-', markersize=3, label="$m_s$=-1")
    a01.plot(tis,pmsref, 'o-', markersize=3, label="$m_s$=0")
    a01.set_ylabel("$p_{NV^-}$")
    a01.legend(fontsize=7)

    a10.plot(tis,fids,    'o-', markersize=3, label="$m_s$=-1")
    a10.plot(tis,fidsref, 'o-', markersize=3, label="$m_s$=0")
    a10.set_ylabel("Fidelity")
    a10.legend(fontsize=7)

    #######################################################################
    n1s = nMeanms*pms + nMean0s*(1-pms)
    n0s = nMeanmsref*pmsref + nMean0sref*(1-pmsref)
    snr = (n1s-n0s)/np.sqrt(n1s+n0s)
    sigmaUnThres = np.sqrt(1 + 2/snr**2)

    a11.plot(tis,sigmaSCC,    'o-', markersize=3, label="Thresholded")
    a11.plot(tis,sigmaUnThres,'o-', markersize=3, label="Unthresholded")
    a11.set_ylabel("$\sigma^R_{SCC}$")
    a11.set_yscale('log')
    if ifYlim: a11.set_ylim(ylim)
    s1 = "$\sigma^R_{min, (un)thres}$=%.2f (%.2f) at $t_i$=%.0f (%.0f)" % (np.min(sigmaSCC), np.min(sigmaUnThres),tis[np.argmin(sigmaSCC)],tis[np.argmin(sigmaUnThres)])
    a11.set_title(s1, fontsize=8)
    a11.legend(fontsize=7)

    #######################################################################
    a20.plot(tis,gms,    'o-', markersize=3, label="$m_s$=-1")
    a20.plot(tis,gmsref, 'o-', markersize=3, label="$m_s$=0")
    a20.set_ylabel("$g_{-0}$ (Hz)")
    a20.legend(fontsize=7)

    a21.plot(tis,g0s,    'o-', markersize=3, label="$m_s$=-1")
    a21.plot(tis,g0sref, 'o-', markersize=3, label="$m_s$=0")
    a21.set_ylabel("$g_{0-}$ (Hz)")
    a21.legend(fontsize=7)

    a30.plot(tis,nms,    'o-', markersize=3, label="$m_s$=-1")
    a30.plot(tis,nmsref, 'o-', markersize=3, label="$m_s$=0")
    a30.set_ylabel("Bright: $\gamma_{-}$ (Hz)")
    a30.legend(fontsize=7)

    a31.plot(tis,n0s,    'o-', markersize=3, label="$m_s$=-1")
    a31.plot(tis,n0sref, 'o-', markersize=3, label="$m_s$=0")
    a31.set_ylabel("Dark: $\gamma_{0}$ (Hz)")
    a31.legend(fontsize=7)

    if 'SCCPhStatSweepDelaySI' in finalDataFolder: 
        lb = "$t_{sh-i}$ (ns)"
    elif sweepWhat == 'tsh':
        lb = "$t_{sh}$ (ns)"
    else:
        lb = "$t_i$ (ns)"

    a10.set_xlabel(lb)
    a11.set_xlabel(lb)
    a30.set_xlabel(lb)
    a31.set_xlabel(lb)

    s1 = "$P_{589}$ = %.1f $\mu$W. $P_{532}$ = %.0f $\mu$W. $P_{635}$ = %.1f mW$. t_{532}$ = %.0f $\mu$s. $t_{in-MW}$ = %.0f $\mu$s" % (power589, power532, power635, t532/1e3, delay1/1e3)
    if 'SCCPhStatSweepDelaySI' in finalDataFolder: 
        s2 = r"$t_{MW-sh}$ = %.0f ns. $t_{sh}$ = %.0f ns. $t_{i}$ = %.0f ns. $t_{i-r}$ = %.1f ms. $\tau$ = %.0f ms" % (delay2, tsh, ti, delay4/1e6, tr_ns/1e6)
    elif sweepWhat == 'tsh':
        s2 = r"$t_{MW-sh}$ = %.0f ns. $t_{sh-i}$ = %.0f ns. $t_{i}$ = %.0f ns. $t_{i-r}$ = %.1f ms. $\tau$ = %.1f ms" % (delay2, delay3, ti, delay4/1e6, tr_ns/1e6)
    else:
        s2 = r"$t_{MW-sh}$ = %.0f ns. $t_{sh}$ = %.0f ns. $t_{sh-i}$ = %.0f ns. $t_{i-r}$ = %.1f ms. $\tau$ = %.1f ms" % (delay2, tsh, delay3, delay4/1e6, tr_ns/1e6)
    plt.suptitle(s1  +  "\n" + s2, fontsize=9);

    plt.tight_layout()
    if ifPlot:
        plt.show()
    else:
        plt.ioff(); plt.clf(); plt.close('all')

    return tis, ths, fids, pms, snrs, gms, g0s, nms, n0s, nMeanms, nMean0s, thsref, fidsref, pmsref, snrsref, gmsref, g0sref, nmsref, n0sref, nMeanmsref, nMean0sref, sigmaSCC

def plotAnalysisT1SCC(finalDataFolder, ifPlot=1, power589 = 2, power532 = 1400, power635 = 9.5,
                t532 = 500e3, delay1 = 20e6, delay2 = 20, tsh = 100, delay3 = 600, ti=160, delay4 = 2e6, tr_ns = 250e6,
                ifYlim=0, ylim=None):
    datafile = finalDataFolder +'/T1SCCObject_sig_set.dat'
    sig, ref, hist2DArray_sig, hist2DArray_ref, xplot, yplot, bins1DArray = readDataFullData(
                                    datafile, num_of_bins=50, binwidth=1, plot_hist_every=100, ifPlot=False,
                                    ifDataSavedAsCountRate=False, ifLogColor=True, ifSubtractRef=False, ifPlotRef=False,)
    sigavg = np.average(sig,axis=1); refavg = np.average(ref,axis=1)

    ms=-1
    paramFile = finalDataFolder + "/fitParams_ms" + str(ms) + ".csv"
    df = pd.read_csv(paramFile)
    tis, ths, fids, pms, snrs, gms, g0s, nms, n0s, nMeanms, nMean0s = df['tis'].values, df['ths'].values, df['fids'].values, df['pms'].values, df['snrs'].values, df['gms'].values, df['g0s'].values, df['nms'].values, df['n0s'].values, df['nMeanms'].values, df['nMean0s'].values

    ms=0
    paramFile = finalDataFolder + "/fitParams_ms" + str(ms) + ".csv"
    df = pd.read_csv(paramFile)
    tis, thsref, fidsref, pmsref, snrsref, gmsref, g0sref, nmsref, n0sref, nMeanmsref, nMean0sref = df['tis'].values, df['ths'].values, df['fids'].values, df['pms'].values, df['snrs'].values, df['gms'].values, df['g0s'].values, df['nms'].values, df['n0s'].values, df['nMeanms'].values, df['nMean0s'].values

    b0 = pmsref; b1 = pms
    sigmaSCC = np.sqrt((b0+b1)*(2-b0-b1)/(b0-b1)**2)

    fig,axs = plt.subplots(5,2, figsize=(6,9.3))
    a00 = axs[0,0]; a01 = axs[0,1]; a10 = axs[1,0]; a11 = axs[1,1]
    a20 = axs[2,0]; a21 = axs[2,1]; a30 = axs[3,0]; a31 = axs[3,1]; a40 = axs[4,0]; a41 = axs[4,1]
    t1s = []; t1s_err = []

    ###################################################################################
    a00.plot(tis,ths,    'o-', markersize=3, label="$m_s$=-1")
    a00.plot(tis,thsref, 'o-', markersize=3, label="$m_s$=0")
    a00.set_ylabel("Threshold")
    a00.legend(fontsize=7)
    a00.set_xscale('log')

    ###################################################################################
    y = (pms-pmsref)/(pms+pmsref); y2 = (sigavg-refavg)/(sigavg+refavg); x = tis; a = np.max(y)
    guess = (a,0,1e6)
    xfit, yfit, popt, perr = fitDecay(x,y,guess=guess)
    xfit2, yfit2, popt2, perr2 = fitDecay(x,y2,guess=guess)

    # a01.plot(x,y,  'o-', linewidth=0.5, markersize=3, label="$C_{NV^-}$", color='C1')
    # a01.plot(xfit,yfit, color='C1')
    # s = "$T_{1,NV^-}$=%.2f$\pm$%.2f ms. $T_{1,tot}$=%.2f$\pm$%.2f ms" % (popt[2]/1e6, perr[2]/1e6, popt2[2]/1e6, perr2[2]/1e6)
    s = "$T_{1}$=%.2f$\pm$%.2f ms" % (popt2[2]/1e6, perr2[2]/1e6)
    a01.plot(x,y2, 'o-', linewidth=0.5, markersize=3, label="$C_{tot}$", color='C0')
    a01.plot(xfit2,yfit2, color='C0')
    a01.set_ylabel("Contrast")
    a01.legend(fontsize=7, loc='lower left')
    a01.set_xscale('log')
    a01.set_title(s, fontsize=7)
    t1s.append(popt[2]/1e6);  t1s_err.append(perr[2]/1e6)
    t1s.append(popt2[2]/1e6); t1s_err.append(perr2[2]/1e6)

    ##################################################################################
    a10.plot(tis,fids,    'o-', markersize=3, label="$m_s$=-1")
    a10.plot(tis,fidsref, 'o-', markersize=3, label="$m_s$=0")
    a10.set_ylabel("Fidelity")
    a10.legend(fontsize=7)
    a10.set_xscale('log')
    
    ##################################################################################
    y = sigmaSCC; x = tis
    a = np.max(y);    guess = (-a,a,1e6)
    xfit, yfit, popt, perr = fitDecay(x,y,guess=guess)
    s2 = "$T_{1,SNR}$=%.2f$\pm$%.2f ms" % (popt[2]/1e6, perr[2]/1e6)

    a11.plot(tis,sigmaSCC,'o-', linewidth=0.5, markersize=3, label=s2, color='C0')
    a11.plot(xfit, yfit, color='C0')
    a11.set_ylabel("$\sigma^R_{SCC}$")
    a11.set_xscale('log')
    if ifYlim: a11.set_ylim(ylim)
    s1 = "$\sigma^R_{SCC,min}$=%.2f at $t_{MW-sh}$=%.0f ns" % (np.min(sigmaSCC), tis[np.argmin(sigmaSCC)])
    a11.set_title(s1, fontsize=8)
    a11.legend(fontsize=7, loc='upper left')

    #################################################################################
    a20.plot(tis,gms,    'o-', markersize=3, label="$m_s$=-1")
    a20.plot(tis,gmsref, 'o-', markersize=3, label="$m_s$=0")
    a20.set_ylabel("$g_{-0}$ (Hz)")
    a20.legend(fontsize=7)
    a20.set_xscale('log')

    #################################################################################
    a21.plot(tis,g0s,    'o-', markersize=3, label="$m_s$=-1")
    a21.plot(tis,g0sref, 'o-', markersize=3, label="$m_s$=0")
    a21.set_ylabel("$g_{0-}$ (Hz)")
    a21.legend(fontsize=7)
    a21.set_xscale('log')

    #################################################################################
    a30.plot(tis,nms,    'o-', markersize=3, label="$m_s$=-1")
    a30.plot(tis,nmsref, 'o-', markersize=3, label="$m_s$=0")
    a30.set_ylabel("Bright: $\gamma_{-}$ (Hz)")
    a30.legend(fontsize=7)
    a30.set_xscale('log')

    #################################################################################
    a31.plot(tis,n0s,    'o-', markersize=3, label="$m_s$=-1")
    a31.plot(tis,n0sref, 'o-', markersize=3, label="$m_s$=0")
    a31.set_ylabel("Dark: $\gamma_{0}$ (Hz)")
    a31.legend(fontsize=7)
    a31.set_xscale('log')

    #################################################################################
    a40.plot(tis, pms,    'o-', markersize=3, label="$m_s$=-1")
    a40.plot(tis, pmsref, 'o-', markersize=3, label="$m_s$=0")
    a40.set_xscale('log')
    a40.legend(fontsize=7)
    a40.set_ylabel("$p_{NV^-}$")

    #################################################################################
    thresmin = np.min(ths); thresmax = np.max(ths); colors=['C0', 'C1']
    for iii, thres in enumerate(np.array((thresmin, thresmax))):
        (row, col) = np.shape(sig)

        sigsThd = np.zeros((row, col))
        for i in range(row):
            for j in range(col):
                if sig[i,j] > thres: sigsThd[i,j] = 1

        refsThd = np.zeros((row, col))
        for i in range(row):
            for j in range(col):
                if ref[i,j] > thres: refsThd[i,j] = 1
        sigThdAvg = np.average(sigsThd,axis=1); refThdAvg = np.average(refsThd,axis=1)
        c = (sigThdAvg-refThdAvg)/(sigThdAvg+refThdAvg)
        x = tis; a = np.max(c)
        guess = (a,0,1e6)
        xfit, yfit, popt, perr = fitDecay(x,c,guess=guess)
        s = "$T_{1}$=%.2f$\pm$%.2f ms. T=%.0f" % (popt[2]/1e6, perr[2]/1e6, thres)

        a41.plot(tis, c, 'o-', linewidth=0.5, markersize=3, color=colors[iii], label=s)
        a41.plot(xfit, yfit, color=colors[iii])    

        t1s.append(popt[2]/1e6); t1s_err.append(perr[2]/1e6)
        
    a41.set_title("Thresholded CSR", fontsize=8)
    a41.legend(fontsize=6, loc='lower left')
    a41.set_xscale('log')
    a41.set_ylabel("Contrast")
    ###################################################################################

    lb = "$t_{MW-sh}$ (ns)"
    a10.set_xlabel(lb)
    a11.set_xlabel(lb)
    a40.set_xlabel(lb)
    a41.set_xlabel(lb)

    s1 = "$P_{589}$ = %.1f $\mu$W. $P_{532}$ = %.1f mW. $P_{635}$ = %.1f mW$. t_{532}$ = %.0f $\mu$s. $t_{in-MW}$ = %.0f $\mu$s" % (power589, power532/1e3, power635, t532/1e3, delay1/1e3)
    s2 = r"$t_{sh}$ = %.0f ns. $t_{sh-i}$ = %.0f ns. $t_{i}$ = %.0f ns. $t_{i-r}$ = %.1f ms. $\tau$ = %.1f ms" % (tsh, delay3, ti, delay4/1e6, tr_ns/1e6)
    plt.suptitle(s1  +  "\n" + s2, fontsize=9);

    plt.tight_layout()
    if ifPlot:
        plt.show()
    else:
        plt.ioff(); plt.clf(); plt.close('all')

    return tis, ths, fids, pms, snrs, gms, g0s, nms, n0s, nMeanms, nMean0s, thsref, fidsref, pmsref, snrsref, gmsref, g0sref, nmsref, n0sref, nMeanmsref, nMean0sref, sigmaSCC, sigavg, refavg, t1s, t1s_err

def plotT1Simple(x, sig, ref, thresmax=2, ifPlot=1, power589 = 2, power532 = 1400, power635 = 9.5,
                t532 = 500e3, delay1 = 20e6, delay2 = 20, tsh = 100, delay3 = 600, ti=160, delay4 = 2e6, tr_ns = 250e6):
    
    fig, axs = plt.subplots(1,2, figsize=(6,3))
    a00 = axs[0]; a01 = axs[1]
    sigavg = np.average(sig,axis=1); refavg = np.average(ref,axis=1)
    t1s = []; t1s_err = []

    ###########################################################
    y = (sigavg-refavg)/(sigavg+refavg); a = np.max(y)
    guess = (a,0,1e6)
    xfit, yfit, popt, perr = fitDecay(x,y,guess=guess)
    s = "$T_{1}$=%.2f$\pm$%.2f ms" % (popt[2]/1e6, perr[2]/1e6)

    a00.plot(x,y, 'o-', linewidth=0.5, markersize=3, label="$C_{tot}$", color='C0')
    a00.plot(xfit,yfit, color='C0')
    a00.set_ylabel("Contrast")
    a00.legend(fontsize=7, loc='lower left')
    a00.set_xscale('log')
    a00.set_title(s, fontsize=7)
    a00.set_xlabel("$t_{MW-sh}$ (ns)")

    t1s.append(popt[2]/1e6); t1s_err.append(perr[2]/1e6)

    ###########################################################
    ratios = []
    for thres in range(thresmax+1):
        (row, col) = np.shape(sig)

        sigsThd = np.zeros((row, col))
        for i in range(row):
            for j in range(col):
                if sig[i,j] > thres: sigsThd[i,j] = 1

        refsThd = np.zeros((row, col))
        for i in range(row):
            for j in range(col):
                if ref[i,j] > thres: refsThd[i,j] = 1
        sigThdAvg = np.average(sigsThd,axis=1); refThdAvg = np.average(refsThd,axis=1)
        c = (sigThdAvg-refThdAvg)/(sigThdAvg+refThdAvg); a = np.max(c)
        guess = (a,0,1e6)
        xfit, yfit, popt, perr = fitDecay(x,c,guess=guess)
        ratios.append(perr[2]/popt[2])
    ratios = np.array(ratios); thres = np.argmin(ratios)
    ###########################################################

    (row, col) = np.shape(sig)
    sigsThd = np.zeros((row, col))
    for i in range(row):
        for j in range(col):
            if sig[i,j] > thres: sigsThd[i,j] = 1

    refsThd = np.zeros((row, col))
    for i in range(row):
        for j in range(col):
            if ref[i,j] > thres: refsThd[i,j] = 1
    sigThdAvg = np.average(sigsThd,axis=1); refThdAvg = np.average(refsThd,axis=1)
    c = (sigThdAvg-refThdAvg)/(sigThdAvg+refThdAvg); a = np.max(c)
    guess = (a,0,1e6)
    xfit, yfit, popt, perr = fitDecay(x,c,guess=guess)
    s = "$T_{1}$=%.2f$\pm$%.2f ms. Thres=%.0f" % (popt[2]/1e6, perr[2]/1e6, thres)
    
    a01.plot(x, c, 'o-', linewidth=0.5, markersize=3, label="Thresholded CSR")
    a01.plot(xfit, yfit, color='C0')    
        
    a01.set_title(s, fontsize=8)
    a01.legend(fontsize=6, loc='lower left')
    a01.set_xscale('log')
    a01.set_xlabel("$t_{MW-sh}$ (ns)")

    t1s.append(popt[2]/1e6); t1s_err.append(perr[2]/1e6)

    s1 = "$P_{589}$ = %.1f $\mu$W. $P_{532}$ = %.0f $\mu$W. $P_{635}$ = %.1f mW$. t_{532}$ = %.0f $\mu$s. $t_{in-MW}$ = %.0f $\mu$s" % (power589, power532, power635, t532/1e3, delay1/1e3)
    s2 = r"$t_{sh}$ = %.0f ns. $t_{sh-i}$ = %.0f ns. $t_{i}$ = %.0f ns. $t_{i-r}$ = %.1f ms. $\tau$ = %.1f ms" % (tsh, delay3, ti, delay4/1e6, tr_ns/1e6)
    plt.suptitle(s1  +  "\n" + s2, fontsize=9, y=0.935);

    plt.tight_layout()
    if ifPlot:
        plt.show()
    else:
        plt.ioff(); plt.clf(); plt.close('all')

    return t1s, t1s_err

def plotT2ESimple(x, sig, ref, xscale='linear', guess=(0.15,40e3,1,0), lowerBounds=(0, 0,0.1,-1), upperBounds=(np.inf,np.inf,4,1),
                thresmax=2, ifPlot=1, power589 = 2, power532 = 1400, power635 = 9.5,
                t532 = 500e3, delay1 = 20e6, delay2 = 20, tsh = 100, delay3 = 600, ti=160, delay4 = 2e6, tr_ns = 250e6):
    
    fig, axs = plt.subplots(1,2, figsize=(6,3))
    a00 = axs[0]; a01 = axs[1]
    sigavg = np.average(sig,axis=1); refavg = np.average(ref,axis=1)
    t2s = []; t2s_err = []; n2s = []; n2s_err = []

    ###########################################################
    y = -(sigavg-refavg)/(sigavg+refavg); a = np.max(y)

    # Find the best lower bound for n and upper bound for c
    ntry = np.linspace(2, 3.9, 20); ctry = np.linspace(0,0.2,11)
    perrs = np.zeros((len(ntry), len(ctry)))
    for i, nlow in enumerate(ntry):
        for j, chigh in enumerate(ctry):
            lowerBounds = (0,      0,      nlow, -1)
            upperBounds = (np.inf, np.inf, 4.01, chigh)
            xfit, yfit, popt, perr = fitStrDecay(x,y,guess=guess,upperBounds=upperBounds,lowerBounds=lowerBounds)
            perrs[i,j] = np.round(perr[1]/popt[1] + perr[2]/popt[2],3) # minimum relerr of T2 and n
    indices = np.unravel_index(np.argmin(perrs), perrs.shape)
    ngood = ntry[indices[0]]; cgood = ctry[indices[1]]

    lowerBounds = (0,      0,      ngood, -1)   
    upperBounds = (np.inf, np.inf, 4.0001, cgood)
    xfit, yfit, popt, perr = fitStrDecay(x,y,guess=guess,upperBounds=upperBounds,lowerBounds=lowerBounds)
    s = "$T_{2}$=%.2f$\pm$%.2f $\mu$s. $n$=%.1f" % (popt[1]/1e3, perr[1]/1e3, popt[2])

    a00.plot(x/1e3,y, 'o-', linewidth=0.5, markersize=3, label="$C_{tot}$", color='C0')
    a00.plot(xfit/1e3,yfit, color='C0')
    a00.set_ylabel("Contrast")
    a00.legend(fontsize=7, loc='lower left')
    a00.set_xscale(xscale)
    a00.set_title(s, fontsize=8)
    a00.set_xlabel("$\\tau$ ($\mu$s)")

    t2s.append(popt[1]/1e3); t2s_err.append(perr[1]/1e3); n2s.append(popt[2]); n2s_err.append(perr[2])

    ###########################################################
    # Use the same bounds as above, find the best threshold
    ratios = []
    for thres in range(thresmax+1):
        (row, col) = np.shape(sig)

        sigsThd = np.zeros((row, col))
        for i in range(row):
            for j in range(col):
                if sig[i,j] > thres: sigsThd[i,j] = 1
        refsThd = np.zeros((row, col))
        for i in range(row):
            for j in range(col):
                if ref[i,j] > thres: refsThd[i,j] = 1
        sigThdAvg = np.average(sigsThd,axis=1); refThdAvg = np.average(refsThd,axis=1)

        c = -(sigThdAvg-refThdAvg)/(sigThdAvg+refThdAvg)
        xfit, yfit, popt, perr = fitStrDecay(x,c,guess=guess,upperBounds=upperBounds,lowerBounds=lowerBounds)
        ratios.append(perr[1]/popt[1])
    ratios = np.array(ratios); thres = np.argmin(ratios)
    ###########################################################

    (row, col) = np.shape(sig)
    sigsThd = np.zeros((row, col))
    for i in range(row):
        for j in range(col):
            if sig[i,j] > thres: sigsThd[i,j] = 1
    refsThd = np.zeros((row, col))
    for i in range(row):
        for j in range(col):
            if ref[i,j] > thres: refsThd[i,j] = 1
    sigThdAvg = np.average(sigsThd,axis=1); refThdAvg = np.average(refsThd,axis=1)

    c = -(sigThdAvg-refThdAvg)/(sigThdAvg+refThdAvg)
    xfit, yfit, popt, perr = fitStrDecay(x,c,guess=guess,upperBounds=upperBounds,lowerBounds=lowerBounds)
    s = "$T_{2}$=%.2f$\pm$%.2f $\mu$s. $n$=%.1f. Thres=%.0f" % (popt[1]/1e3, perr[1]/1e3, popt[2], thres)
    
    a01.plot(x/1e3, c, 'o-', linewidth=0.5, markersize=3, label="Thresholded CSR")
    a01.plot(xfit/1e3, yfit, color='C0')    
        
    a01.set_title(s, fontsize=8)
    a01.legend(fontsize=6, loc='lower left')
    a01.set_xscale(xscale)
    a01.set_xlabel("$\\tau$ ($\mu$s)")

    t2s.append(popt[1]/1e3); t2s_err.append(perr[1]/1e3); n2s.append(popt[2]); n2s_err.append(perr[2])

    s1 = "$P_{589}$ = %.1f $\mu$W. $P_{532}$ = %.0f $\mu$W. $P_{635}$ = %.1f mW$. t_{532}$ = %.0f $\mu$s. $t_{in-MW}$ = %.0f $\mu$s" % (power589, power532, power635, t532/1e3, delay1/1e3)
    s2 = r"$t_{sh}$ = %.0f ns. $t_{sh-i}$ = %.0f ns. $t_{i}$ = %.0f ns. $t_{i-r}$ = %.1f ms. $\tau$ = %.1f ms" % (tsh, delay3, ti, delay4/1e6, tr_ns/1e6)
    plt.suptitle(s1  +  "\n" + s2, fontsize=9, y=0.935);

    plt.tight_layout()
    if ifPlot:
        plt.show()
    else:
        plt.ioff(); plt.clf(); plt.close('all')

    return t2s, t2s_err, n2s, n2s_err






























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

    mainFolder = 'C:/Users/lukin2dmaterials/data/2023-05-22/'
    for dataFolder in os.listdir(mainFolder):
        # print(dataFolder)
        if 'CalibReadoutPhotonStat' in dataFolder:
            idx = int(dataFolder[1:4])
            if idx >= 0:
                datafile = mainFolder + dataFolder +'/CalibReadoutPhotonStatObject_sig_set.dat'
                readDataFullData(datafile, type='CalibReadoutPhotonStat', typeNorm=1)

    # mainFolder = 'C:/Users/lukin2dmaterials/data/2022-12-19/'
    # for dataFolder in os.listdir(mainFolder):
    #     # print(dataFolder)
    #     if 'T2E' in dataFolder:
    #         idx = int(dataFolder[1:4])
    #         if idx >= 0:
    #             datafile = mainFolder + dataFolder +'/T2EObject_sig_set.dat'
    #             fig = readData(datafile, type='T2E', typeNorm=1)
    # plt.show()

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
    