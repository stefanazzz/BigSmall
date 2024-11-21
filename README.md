2D dynamic rupture
------------------

2D_dynamic/

	fault2dPML.f: main code and subroutines 
 
	PMLpack.f: additional subroutines for absorbing boundaries
 
	arrays.f: array declaration for the fortran code
 
	fpar_model1: input parameter file for fortran code

	fpar_model2: alternate input parameter file
 
	Makefile: Makefile to compile fortran code
 
	RES_model1: output example of code
 
	RES_model2: output example of code
 
	README: instructions on use
 
	read_plot_dimensionless.py: plotting script for results

3D dynamic rupture
------------------

1) Homogeneous stress on fault:

3D_dynamic/Homogeneous/
	triffy.f: main fortran code and subroutines 
	triffy.dec: arrays declarations for fortran code
	fpar input parameters file for fortran code
	Makefile: Makefile to compile fortran code
	RES output of code: outoput example of code
	RES2 output of code: outoput example of code
	plot_area_beta.py: plotting script for results

2) Inhomogeneous stress on fault:

3D_dynamic/Inhomogeneous/
same file descriptors as above except for additional files:
	filtered_array_fortran3.txt: 
	inhomogeneous self-affine stress file as input for fotran code
	create_array.py:
	creates above file with self-affine stress array

Model to test scaling of critical breakout with self-affine ditribution:
------------------------------------------------------------------------

Self_affine/

Check_Areas_Separate.py: 
1) Creates a random 2D distribution of uncorrelated values
2) Filters the distribution to produce a powerlaw w wavenumber dependance 
   i.e. introducing correlation and creating self-attine 'topography'
3) Option to create a bi-harmonic distrbution and a cut-off filter 
   mainly for purpose of testing that scaling and cutoff work as expected
4) Computes mean and variance of the distribution for each sub-section of the
array, then saves it to a file

simplot_tau_square2_separate2.py:
Plots the outputs from Check_Areas.py
1) For different sub-fault sizes L, calculates F = mean(stress) x sqrt(L)
2) Finds max of F(L) 
3) Plots max F for each subarray in each maximum array size
4) Plots the max F in all subarrays as a function of maximum array size
	
