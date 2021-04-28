import numpy as np
import math
import cmath
from scipy import integrate
from pylab import plot,show,errorbar,linspace,array,zeros
from pylab import imshow,legend,xlabel,ylabel,colorbar,figure,title
import scipy.special as sc
import scipy.fft

#input of variable
lambda_ = float(input("wavelength (in nm or 10e-9m):"))*(10**(-9))
z = float(input("distance of propergation (in m):"))
w = float(input("size of triangular aperture (in mm or 10e^-3):"))*(10**(-3))
a = float(input("sieze of canvas (in mm):"))*(10**(-3))
k = 2*math.pi/lambda_


#define a function h such that U-out is the convolution of U_in and h
def h(x,y):
    j = complex(0,1)
    return cmath.exp(j*k*z)  * cmath.exp((j*k*(x*x+y*y))/(2*z)) / (j*lambda_*z)


#define the function of U_in (1 if inside the area of apture, 0 otherwise)
def U(x,y):
    if y > math.sqrt(3)*w/4 or y < -math.sqrt(3)*w/4:
        return 0
    else:
        if y > math.sqrt(3)*x - math.sqrt(3)*w/4 and x > 0:
            return 1
        elif y > -math.sqrt(3)*x - math.sqrt(3)*w/4 and x < 0:
            return 1
        else:
            return 0

phase = np.vectorize(cmath.phase)       
N = 200
H = zeros([N,N],complex)   # for fouier transform of h
I = zeros([N,N],complex)   # for fouier transform of U_in
A = zeros([N,N],complex)   # for the product of I and H
x_1 = linspace (-a,a,N)    # recall a is the size of canvus
y_1 = linspace (-a,a,N)

# calculate the values before fourier tranform
i = 0
for x in x_1:
    j = 0
    for y in y_1:
        H[i][j] = h(x,y)
        I[i][j] = U(x,y)
        j += 1
    i += 1

# do fouier transform
H = scipy.fft.fft2(H)
I = scipy.fft.fft2(I)
A = np.multiply(H,I)

#intensity and phase after propagation
U_out = scipy.fft.ifft2(A)
U_phase = phase(U_out) 


#plot
imshow(np.absolute(U_out).T, origin="lower", extent = [-1000*a,1000*a, -1000*a, 1000*a],aspect = 1)
colorbar()
xlabel("x(mm)")
ylabel("y(mm)") 
title("wavelength:"+str(int(lambda_*(10**9))) +"nm" +".   "+ "distance:" + str(z) +"m.    "\
      +"size of aperture:" + str(w* 10 ** 3) +"mm")
show()

imshow(phase(U_out).T, origin="lower", extent = [-1000*a,1000*a, -1000*a, 1000*a],aspect = 1)
colorbar()
xlabel("x(mm)")
ylabel("y(mm)") 
title("Phase")
show()