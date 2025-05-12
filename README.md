The two-dimensional Navier-Stokes equation in the stream function $\psi$ and vorticity $\omega$ formulation in the Fourier space is 
$$\hat{\psi} = -\hat{\omega}/k^2$$

$$\frac{\partial \hat{\omega}}{\partial t} = i (k_x \widehat{u_x \omega} + k_y \widehat{u_y \omega}) -(\nu k^2 + \mu) \hat{\omega} + \hat{f}$$,

where $\hat{.}$ is the Fourier transform, $(u_x, u_y) = (-\partial_y \psi , \partial_x \psi)$, and $\hat{f}$ is the curl of the forcing term in Fourier space. 
In two dimensions, the inverse cascade feeds energy into long wavelengths, which we avoid by including the Ekman friction $\mu$, in addition to the viscous 
dissipation $\nu$. We solve the spectral DNS in $512 \times 512$ numerical grids in a $2\pi$ periodic domain, and the following parameters

|Domain | $k_{f}$ | $F_{0}$ | $\nu$ | $\mu$ |
|--|--|--|--|--|
|512 $^{2}$ | 3 m $^{-1}$ | 0.1 ms $^{-2}$ | $5\times 10^{-6}$ m $^{2}$ s $^{-1}$ | 0.01 s $^{-1}$ |

We initialize the flow with a Taylor-Green vortex array

$$ u_x  = \frac{1}{2}\sin(x) \cos(y) \exp{(-2 \nu t)} $$

$$ u_y = - \frac{1}{2} \cos(x) \sin(y) \exp{(-2 \nu t)} $$

which gives the vorticity field $` \omega (x,y) = \sin x \sin y `$ , at $` t= 0 `$ . Since the flow in $` x `$ and $` y `$ directions are coupled, the forcing wavenumber along $` x `$ suffices to drive the flow in both directions in wavenumber space.

The Fourier and inverse Fourier transforms at each time step were carried out using the C library fftw3.
To eliminate aliasing errors, a third of the large wavenumber Fourier modes were set to zero at each time step. At each time step the fourier transform of forcing $` f = f_{0} k_{f} \cos (k_{f} x) `$ is added to the governing equation which drives the fluid. 

This at long enough time scales leads to a steady state turbulence, after which we introduce self propelled particles with uniformly random positions and orientations, and thereafter the particle and fluid evolves simultaneously. We calculate the fluid dynamics leading to the steady state turbulence only once, and use it as an initial condition (file initial.dat) for the DNS in subsequent simulations. The dynamics of particles happen in a continuous space, whereas fluid dynamics takes place on the numerical grid points separated by $` \Delta x = \Delta y = 2 \pi /512 `$ [see Figure below]. Particle dynamics at position $` (x,y) `$ require the velocity field $` \vec{U} `$ and the tensor field $` \nabla \vec{U} `$ at $` (x,y) `$. We use bilinear interpolation to get the value of $` \vec{U} = (u_x, u_y) `$ and $` \vec{\nabla U} `$ at a point $` (x,y) `$ lying in a grid specified by the indices $` (i,j) `$, $` (i+1,j) `$, $` (i,j+1) `$, $` (i+1,j+1) `$, which can be written in the matrix form to leading order in $` \Delta x, \Delta y `$ as given in https://en.wikipedia.org/wiki/Bilinear_interpolation

![bilinear](https://github.com/user-attachments/assets/44f3cbbe-7263-43d3-9078-f60435e9bfdd)

We treat other field components by similar interpolation. We evolve the particles using the equations (5) & (6) in R. Chajwa et al. https://arxiv.org/abs/2310.01829 and the flow fields using the Runge-Kutta-4 algorithm with the time step $` \Delta t = 10^{-3} `$ .


