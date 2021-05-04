import numpy as np
import math
import cmath
from scipy import integrate
from pylab import plot,show,errorbar,linspace,array,zeros
from pylab import imshow,legend,xlabel,ylabel,colorbar,figure,title
import scipy.special as sc
import scipy.fft

a = float(input("sieze of canvas (in mm):"))*(10**(-3))
sigma=10**(-6)
N = 200
z = 4

M=np.arange(-21,22,2)
print(M)

#approximation of r
def r(x,y):
    return z + 0.5*(x**2+y**2)/z

#define the function of output
def I(x,y):
    I = 0
    for m in M:
        pi_m = m/100
        I = I + 2 / pi_m * math.exp(-( (x-2*pi_m / r(x,y) )**2 + y**2)/sigma)
    return I


g = zeros([N,N],complex)

x_1 = linspace (-a,a,N)    # recall a is the size of canvus
y_1 = linspace (-0.5*a,0.5*a,N)

# calculate the values before fourier tranform
i = 0
for x in x_1:
    j = 0
    for y in y_1:
        g[i][j] = I(x,y)
        j += 1
    i += 1
  
 
print(np.absolute(g))

imshow(np.real(g).T, origin="lower", extent = [-1000*a,1000*a, -500*a, 500*a],aspect = 1)
colorbar()
xlabel("x(mm)")
ylabel("y(mm)") 
title("z="+str(int(z)) + "m, sigma=10**(-6), ")
show()