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

which gives the vorticity field $\omega (x,y) = \sin x \sin y$, at $t= 0$. Since the flow in $x$ and $y$ directions are coupled, the forcing wavenumber along $x$ suffices to drive the flow in both directions in wavenumber space.

The Fourier and inverse Fourier transforms at each time step were carried out using the C library fftw3.
To eliminate aliasing errors, a third of the large wavenumber Fourier modes were set to zero at each time step \cite{Boffetta2012}. At each time step the fourier transform of forcing $f = f_{0} k_{f} \cos (k_{f} x) $ is added to \ref{eqn:520s} which drives the fluid. 
