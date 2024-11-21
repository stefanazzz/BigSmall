"""
S. Nielsen 28 / 11 / 2023

goood code version

1) Creates a random 2D distribution of uncorrelated values
2) Filters the distribution to produce a powerlaw w wavenumber dependance 
   i.e. introducing correlation and creating self-attine 'topography'
3) Option to create a bi-harmonic distrbution and a cut-off filter 
   mainly for purpose of testing that scaling and cutoff work as expected
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import sys, os
os.chdir('/home/stefan/Desktop/BigSmall/PythonScripts')
#%%

# Function to generate a 2D test image
def generate_test_image(size, dx):
    option = 'gaussian'
    if option == 'random':
        image = np.array(np.random.uniform(10, 110, size))
    elif option == 'gaussian':
        image = np.random.normal(500, 50, size)
    elif option == 'biharmonic':
        x = np.arange(0, size[1]) * dx
        y = np.arange(0, size[0]) * dx
        xx, yy = np.meshgrid(x, y)
        image = 4 + np.sin(2 * np.pi * 1.0 * xx) + np.cos(2 * np.pi * 1.0 * yy)  \
        + np.sin(2 * np.pi * 3.0 * xx) + np.cos(2 * np.pi * 3.0 * yy)
    #
    return image

#%% Function to apply 2D low-pass filter using FFT
def apply_2d_low_pass_filter(image, cutoff_frequency, dx):
    fft_result = np.fft.fft2(image)  #, norm = 'forward')
    u = np.fft.fftfreq(image.shape[0], d=dx)
    v = np.fft.fftfreq(image.shape[1], d=dx)
    uu, vv = np.meshgrid(u, v)
    option = 'powerlo'
    if option == 'powerlo':
        # avoid zero division and preserve mean value of original array:
        uu[0,0] = 1/ np.sqrt(2) 
        vv[0,0] = 1/ np.sqrt(2) 
        # exponent 1.2 Candela 2011 p 965 $P(k) \propto -2 -2 H_\tau$ 
        filter = np.sqrt(uu**2 + vv**2)**(-1.2) 
        fft_result_filtered = (fft_result * filter) 
        filtered_image = np.real((np.fft.ifft2((fft_result_filtered))))
        #mean_value = np.mean(image)
        #fft_result_filtered += mean_value
    elif option == 'low_pass':  
        filter = np.sqrt((uu**2 + vv**2)) <= cutoff_frequency
        fft_result_filtered = (fft_result * filter)
        filtered_image = np.real((np.fft.ifft2((fft_result_filtered))))
    # Perform inverse FFT to obtain the filtered image

    return filtered_image

######################################
def generate_and_filter():
    # Generate a test image
    test_image = generate_test_image(image_size, dx) 

    # Apply 2D low-pass filter
    filtered_image = apply_2d_low_pass_filter(test_image, 0.1, dx)
    # print(filtered_image)

    print('##########################')
    print('\n size', image_size[0])
    print('\n orig max min', np.max(test_image), np.min(test_image), 
          '\n mean', np.mean(test_image),'span',np.max(test_image) - np.min(test_image))
    print('\n filt max min', np.max(filtered_image), np.min(filtered_image),
          '\n mean', np.mean(filtered_image),'span',np.max(filtered_image) - np.min(filtered_image))
  
    return(filtered_image)

#%%###################################
# Set parametersc
#
# image_size = [1024,1024]
# subarrays = [[32,32],[64,64],[128,128],[256,256],[512,512],[1024,1024]]
#
# image_size = [4096,4096]
# subarrays = [[32,32],[64,64],[128,128],[256,256],[512,512],[1024,1024],[2048,2048],[4096,4096]]
#
image_size = [16384,16384]
# subarrays = [[32,32],[64,64],[128,128],[256,256],[512,512],[1024,1024],[2048,2048],[4096,4096],
#             [8192, 8192], [16384,16384]]
#
# image_size = [2048,2048]
# subarrays = [[32,32],[64,64],[128,128],[256,256],[512,512],[1024,1024],[2048,2048]]
#image_size = [8192,8192]
#subarrays = [[32,32],[64,64],[128,128],[256,256],[512,512],[1024,1024],[2048,2048],[4096,4096],[8192,8192]]


#sizes = [ [8192, 8192], [4096, 4096], [2048, 2048], [1024,1024], [512,512], [256,256], [128,128], [64,64], [32,32] ]
#subarrays = [ [8192, 8192], [4096, 4096], [2048, 2048], [1024,1024], [512,512], [256,256], [128,128], [64,64], [32,32] ]

siz = str(image_size[0])


dx = 1.0
filtered_image = generate_and_filter()

plt.figure(figsize=(4, 4))
plt.subplot(1, 1, 1)
plt.imshow(filtered_image, cmap='viridis')
plt.title('Filtered Array ($k^{-1.2}$)')

# plt.subplot(1, 3, 2)
# plt.imshow(np.log(np.abs(fftshift(fft2(random_array))) + 1), cmap='viridis')
# plt.title('Fourier Transform (log scale)')

# plt.subplot(1, 3, 3)
# plt.imshow(filtered_array, cmap='viridis')
# plt.title('Filtered Array ($k^{-1.2}$)')

# Save the filtered array to a file readable in Fortran
np.savetxt('filtered_array_fortran4.txt', filtered_array)

plt.show()

        
            
###
