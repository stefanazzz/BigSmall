"""
S. Nielsen 28 / 11 / 2023

goood code version

1) Creates a random 2D distribution of uncorrelated values
2) Filters the distribution to produce a powerlaw w wavenumber dependance 
   i.e. introducing correlation and creating self-attine 'topography'
3) Option to create a bi-harmonic distrbution and a cut-off filter 
   mainly for purpose of testing that scaling and cutoff work as expected
4) Create large random self-affine array and compute mean in each sub-array
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

#%%########################
'''
Choose on of TOTAL FAULT SIZE between 1024 or 2048, 4096, 8192 or 16384
Uncomment the corresponding image_size and subarrays
'''
###########
#
# SMALLEST FAULT CASE:
# image_size = [1024,1024]
# subarrays = [[32,32],[64,64],[128,128],[256,256],[512,512],[1024,1024]]
#
# SECOND SMALLEST:
# image_size = [2048,2048]
# subarrays = [[32,32],[64,64],[128,128],[256,256],[512,512],[1024,1024],[2048,2048]]
#
# etc. etc.
# image_size = [4096,4096]
# subarrays = [[32,32],[64,64],[128,128],[256,256],[512,512],[1024,1024],[2048,2048],[4096,4096]]
#
image_size = [8192,8192]
subarrays = [[32,32],[64,64],[128,128],[256,256],[512,512],[1024,1024],[2048,2048],[4096,4096],
             [8192, 8192]]
#
# image_size = [16384,16384]
# subarrays = [[32,32],[64,64],[128,128],[256,256],[512,512],[1024,1024],[2048,2048],[4096,4096],
#             [8192, 8192], [16384,16384]]
#

#%%
siz = str(image_size[0])

for it in range(6):
    f1 = open('tau_square_'+siz+'_iter_'+str(it)+'.csv', 'w') 
    f2 = open('Areas_stats_'+siz+'_iter_'+str(it)+'.csv', 'w')
    print('cycle:'+str(it))
    dx = 1.0
    filtered_image = generate_and_filter()
    #
    #sys.exit()
    for sub in subarrays:
        nsub = int(image_size[0]/sub[0])
        print('array:'+str(nsub))
        for i in range(nsub):
            x1 = i * sub[0]; x2 = (i+1) * sub[0]
            for j in range(nsub):
                y1 = j * sub[1]; y2 = (j+1) * sub[1]  
                #print('############ subarray ', sub[0],'sector ',x1,x2,y1,y2 )
                # f.close()
                tau_mean = np.mean (filtered_image[x1:x2, y1:y2])
                tau_sq = tau_mean**(2)
                tau_var = np.var(filtered_image[x1:x2, y1:y2])
                tau_std = np.std(filtered_image[x1:x2, y1:y2])
                f1.write("%s,%s\n"%(sub[0],tau_sq))
                f2.write("%s,%s,%s,%s\n"%(sub[0],tau_mean,tau_std,tau_var))
    f1.close()
    f2.close()   
#

        
            
###
