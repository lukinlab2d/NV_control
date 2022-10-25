import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

# read dat file
datafile = 'C:/Users/lukin2dmaterials/data/2022-10-24/#029_Rabi_10-09-08/RabiObject_sig_set.dat'
readfile = np.loadtxt(datafile)
# print(readfile)
x_s = [xAxisAndResult[0] for xAxisAndResult in readfile]
sig = [xAxisAndResult[1] for xAxisAndResult in readfile]
ref = [xAxisAndResult[2] for xAxisAndResult in readfile]

fig,ax = plt.subplots()
ax.plot(x_s, sig, label='sig')
ax.plot(x_s, ref, label='ref')
ax.legend(loc='best')
ax.set_xlabel(r"$\tau$ (ns)")
plt.show()
