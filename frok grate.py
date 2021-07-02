from setup import *




#grating function:
def fork_grate_(p,l,x,y):
    theta = atan(y/ x)
    if math.cos(2*math.pi*x/p + l*theta)>0:
        return cexp(np.pi/2*j)
    else:
        return cexp(-np.pi/2*j)
fork_grate = np.vectorize(fork_grate_)




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
sigma=0.1*(10**-2)

#output kspace size
kxmax = np.pi*Nx/(2*Lx)
kymax = np.pi*Ny/(2*Ly)

#output screen size
xmax = Nx*lambda_*z/(Lx*2)
ymax = Ny*lambda_*z/(Ly*2)



fork = zeros([Ny,Nx],complex)
x_0 = linspace (-Lx,Lx,Nx)
y_0 = linspace (-Ly,Ly,Ny)

i = 0
for x in x_0:
    q = 0
    for y in y_0:
        fork[q][i] = fork_grate_(p,l,x,y)
        q += 1
    i += 1
phase_fork = phase(fork)




input_wave = cexp(j*complex_(phase_fork)) 
output_wave = fftshift(fft2(input_wave))
np.savetxt('output_wave_fork.csv', output_wave, delimiter = ',')
output_intensity = absl(output_wave)
output_phase = phase(output_wave)

#expected output
sheet1 = Sheet(rangX = [-kxmax, kxmax] , rangY = [-kymax, kymax], Nx = Nx, Ny = Ny)
sheet1.sample_oam_output(np.arange(-9,10,2), z=z, sigma=sigma, p=p,l=l)
np.savetxt('output_wave_expected.csv', sheet1.Psi, delimiter = ',')






plt.imshow(phase_fork , extent = [-Lx*1000,Lx*1000,-Ly*1000,Ly*1000], cmap ='gray', interpolation = "bilinear", aspect = 'auto')
colorbar()
xlabel("x(mm)")
ylabel("y(mm)")
title("fork"+"wavelength:"+str(int(lambda_*(10**9))) +"nm" +".   "+ "distance:" + str(z) +"m.    "+"p:"+ str(p*10**6)+ " um"\
      +"   sigma:"+str(sigma*100)+"cm")
show()

plt.imshow(output_intensity , extent = [-xmax*1000,xmax*1000,-ymax*1000,ymax*1000], cmap ='gray', interpolation = "bilinear", aspect = 'auto')
colorbar()
xlabel("x(mm)")
ylabel("y(mm)")
title("output_intensity"+"wavelength:"+str(int(lambda_*(10**9))) +"nm" +".   "+ "distance:" + str(z) +"m.    "+"p:"+ str(p*10**6)+ " um"\
      +"   sigma:"+str(sigma*100)+"cm")
show()


plt.imshow(output_phase , extent = [-xmax*1000,xmax*1000,-ymax*1000,ymax*1000], cmap ='gray', interpolation = "bilinear", aspect = 'auto')
colorbar()
xlabel("x(mm)")
ylabel("y(mm)")
title("output_phase"+"wavelength:"+str(int(lambda_*(10**9))) +"nm" +".   "+ "distance:" + str(z) +"m.    "+"p:"+ str(p*10**6)+ " um"\
      +"   sigma:"+str(sigma*100)+"cm")
show()


plt.imshow(phase(sheet1.Psi) , extent = [-xmax*1000,xmax*1000,-ymax*1000,ymax*1000], cmap ='gray', interpolation = "bilinear", aspect = 'auto')
colorbar()
xlabel("x(mm)")
ylabel("y(mm)")
title("output_phase_expected"+"wavelength:"+str(int(lambda_*(10**9))) +"nm" +".   "+ "distance:" + str(z) +"m.    "+"p:"+ str(p*10**6)+ " um"\
      +"   sigma:"+str(sigma*100)+"cm")
show()



