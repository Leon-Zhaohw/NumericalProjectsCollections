USAGE
Currently tested on ubuntu 14.04 LTS and ubuntu 16.04

1. Clone the code.
2. Install the related library (cmake, FFTW, glut, opengl, OpenMP, X11, libjpeg).
	cmake: [sudo apt-get install cmake]

	Install OpenMp before install FFTW, because we use fftw with omp support.
	
	FFTW install: For best performance, you need to install FFTW with multi-thread support. You can use the fftw in the repo.
	Download the fftw tar into a folder, unzip the tar ball.
	cd into fftw folder
	[./configure --enable-openmp] , you can also see the options by ./configure --help
	[make -j16]
	[sudo make install]

	x11 install: [sudo apt-get install libx11-dev]
	
	opengl and glut: If you don't have them [sudo apt-get install freeglut3 freeglut3-dev].

	OpenMP: seems most compliers have support for OpenMp
	
	libjpeg: [sudo apt-get install libjpeg-dev]

3. Compile the program.
	cd to the root folder of the repo, which contains a CMakeLists.txt
	[mkdir build]
	[cd build]
	[cmake ../]
	[make -j16]
	The binary then should be in ../bin

4. Run First example with diriclet bc and directable forward scattring.
	You need to first precompute a 3D tensor. To do this, at the root folder of the repo:
	
	[mkdir Tensor]
	
	[perl precompute_3Dtensor.pl]  This should precompute a all_dirichlet tensor with about 2700 basis functions
	It would take a little while to run. The tensor should be in the folder ./Tensor
	
	[mkdir preview]
	
	[perl fluid3D.perl] This would run the sphere passing through smoke example. After the program finishes, 
	the preview of the animation is saved in folder ./preview
	
	Enable forward scattering: To do this, change the parameter basisweightMultiplierC in fluid_moving_sphere.cfg
	to any value, this parameter corresponds to c in equation 22 of the paper. E.g. 0.0005
	To see the weights for each basis, after the program is run, a text file named "weight.txt" will be generated in the root folder, you can view what weight each basis get. The first column is |K|^2, and second column is basis weight.

	Precompute tensor with more basis functions: To do this, change the $basis_dim_des parameter in precompute_3Dtensor.pl. Detailed comments are in the precompute_3Dtensor.pl. To use the new tensor, change 
	basis_tensor_fname in the scene cfg file to point toward the tensor you computed.

5. Run second example with paddle_whell
	To run this, a Tensor with two neumann wall need to be generated.
	Go to  precompute_3Dtensor.pl and change basis_type to two_neumann_x
	Specify the basis_dim_des to 3000 (You can also use other number, but you need to change the basis_tensor_fname in
	fluid_3D_paddle.cfg)
	[perl precompute_3Dtensor.pl]
	Open fluid3D.perl, point the script to fluid_3D_paddle.cfg
	[perl fluid3D.perl]

6. Run the third smoke colliding example
	Go to fluid_two_phase_pulse.cfg and point basis_tensor_fname to a correct tensor.
	Open fluid3D.perl, point the script to fluid_two_phase_pulse.cfg
	[perl fluid3D.perl]

7. Out put the pbrt rendering files:
	In any of the cfg files, enable write_density_to_PBRT and point density_folder to a proper foler.
	If you want my rendering scripts, let me know.

8. If you have any questions regarding the code, email me at qiaodong@ucsb.edu