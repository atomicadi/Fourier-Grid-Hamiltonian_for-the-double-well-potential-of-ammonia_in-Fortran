# 1D Fourier Grid Hamiltonian (FGH) for the double well potential of NH<sub>3</sub> in Fortran and Morse potential of H<sub>2</sub>
![image alt](https://github.com/atomicadi/Fourier-Grid-Hamiltonian_for-the-double-well-potential-of-ammonia_in-Fortran/blob/e337da37fbc176ae7765c56a2dc8750cec082dbe/Untitled.001.png)
Fourier Grid Hamiltonian (FGH) is a numerical technique used to solve the Schrodinger equation for bound states. Instead of using basis functions (like harmonic oscillators or plane waves), the method discretizes the coordinate space into a grid and constructs the Hamiltonian.
<p align="center">


$$
H =  T + V(x)  ...... (1)
$$


</p>

Where H = Hamiltonia, T = kinetic energy, and V(x) = Potential energy. In this method, the kinetic energy matrix is computed using Fourier transformation, which efficiently connects the momentum space and coordinate space representations and finally the form of the hamiltonian becomes,
<p align="center">


$$
H_{ij}^0 = \frac{2}{N} \sum_{l=1}^{n} (cos(\frac{l2\pi(i-j)}{N})\frac{\hbar^2}{2m}(lΔk)^2) + V(x_i)δ_{ij}  ...... (2)
$$


</p>

The delailed about the FGH is written in a pdf format and uploaded in this GitHub page (**FGH.pdf**).

Due to the rapid transformation of ammonia's equilibrium geometry (C<sub>3V</sub>) to its' mirror-image strucrure (which in invarient and same C<sub>3V</sub> point group) via ammonia's umbrella inversion mode, ammonia is trapped inside a symmetric double well potential (denoted as DWP in the rest of the README). This rapid transformation is occured through a planner sp<sup>2</sup> hybridiged geometry (D<sub>3h</sub>) which is responsiple for the potential barrier in this DWP. The form of the potential is,
<p align="center">


$$
H_{ij}^0 = \frac{2}{N} \sum_{l=1}^{n} (cos(\frac{l2\pi(i-j)}{N})\frac{\hbar^2}{2m}(lΔk)^2) + V(x_i)δ_{ij}  ...... (1)
$$


</p>

Where V<sub>max</sub> = Potential at the central maxima, ($\frac{a}{2}$ ± b) = Location of the two local minima, $\frac{a}{2}$ = Location of the central maxima.

The potential barrier in this DWP satisfies the conditions necessary for quantum tunneling, leading to transmission through the barrier. This, in turn, causes the splitting of the localized energy levels within each well of the DWP.

