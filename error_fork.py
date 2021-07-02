from setup import *
#constants and set up
Lx = 60*mm
Ly = 60*mm
Nx= 5000
Ny= 5000
lambda_ = 632*1e-9
k = 2*np.pi/lambda_
z=1
p=200*(10**-6)
l=7
sigma=0.05*(10**-2)
#output kspace size
kxmax = np.pi*Nx/(2*Lx)
kymax = np.pi*Ny/(2*Ly)

#output screen size
xmax = Nx*lambda_*z/(Lx*2)
ymax = Ny*lambda_*z/(Ly*2)



#denose data
def denoise(data):
    fft_data = fft2(data)
    fft_data[np.abs(fft_data) < 1e3] = 0
    good_data = absl(ifft2(fft_data))
    return good_data





#error
output = np.loadtxt("output_wave_fork.csv",dtype=np.complex_,delimiter=",",skiprows=0)
normal_output = normalization(output)


x_0 = linspace (-Lx,Lx,Nx)
y_0 = linspace (-Ly,Ly,Ny)
i = 0
for x in x_0:
    q = 0
    for y in y_0:
        if normal_output[q][i] < 0.0001:
            normal_output[q][i] = 0.00
        q += 1
    i += 1
    




output_expected = np.loadtxt("output_wave_expected.csv",dtype=np.complex_,delimiter=",",skiprows=0)

expected_phase = (phase(output_expected)+ 2*np.pi) % (2*np.pi)

error_phase = (phase(normal_output)+ 2*np.pi) % (2*np.pi) - expected_phase
error_intensity = normal_output - normalization(absl(output_expected))


plt.imshow(expected_phase , extent = [-xmax*1000,xmax*1000,-ymax*1000,ymax*1000], cmap ='gray', interpolation = "bilinear", aspect = 'auto')
colorbar()
xlabel("x(mm)")
ylabel("y(mm)")
title("expected_phase"+"wavelength:"+str(int(lambda_*(10**9))) +"nm" +".   "+ "distance:" + str(z) +"m.    "+"p:"+ str(p*10**6)+ " um"\
      +"   sigma:"+str(sigma*100)+"cm")
show()


plt.imshow( (phase(normal_output)+ 2*np.pi) % (2*np.pi) , extent = [-xmax*1000,xmax*1000,-ymax*1000,ymax*1000], cmap ='gray', interpolation = "bilinear", aspect = 'auto')
colorbar()
xlabel("x(mm)")
ylabel("y(mm)")
title("calculated output"+"wavelength:"+str(int(lambda_*(10**9))) +"nm" +".   "+ "distance:" + str(z) +"m.    "+"p:"+ str(p*10**6)+ " um"\
      +"   sigma:"+str(sigma*100)+"cm")
show()


plt.imshow( error_phase , extent = [-xmax*1000,xmax*1000,-ymax*1000,ymax*1000], cmap ='gray', interpolation = "bilinear", aspect = 'auto')
colorbar()
xlabel("x(mm)")
ylabel("y(mm)")
title("error_phase"+"wavelength:"+str(int(lambda_*(10**9))) +"nm" +".   "+ "distance:" + str(z) +"m.    "+"p:"+ str(p*10**6)+ " um"\
      +"   sigma:"+str(sigma*100)+"cm")
show()