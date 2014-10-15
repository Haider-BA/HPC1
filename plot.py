import numpy as np
#import pylab as pl
import matplotlib.pyplot as plt 

array1=np.loadtxt('mkl.res')
array2=np.loadtxt('mkl_O2.res')
array3=np.loadtxt('mkl_O3.res')
x =range(20)
y1 =range(20)
y2 = range(20)
y3 = range(20)
for index in range(len(array1)):
	x[index] = array1[index][0]
	y1[index] = array1[index][1]
	y2[index] = array2[index][1]
	y3[index] = array3[index][1]
#pl.plot(x,y)
#pl.show()
plt.semilogx(x,y1,'r', label='MKL default')
plt.semilogx(x,y2,'g', label='MKL -O2')
plt.semilogx(x,y3,'b', label='MKL -O3')
plt.legend()
#plt.grid(True)
plt.ylabel('MFlop/s')
plt.xlabel('Vector Length')
plt.show()
