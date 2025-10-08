# 1D Fourier Grid Hamiltonian (FGH) for the Morse potential of H<sub>2</sub> and double well potential of NH<sub>3</sub> in Fortran 
![image alt](https://github.com/atomicadi/Fourier-Grid-Hamiltonian_for-the-double-well-potential-of-ammonia_in-Fortran/blob/e337da37fbc176ae7765c56a2dc8750cec082dbe/Untitled.001.png)
Fourier Grid Hamiltonian (FGH) is a numerical technique used to solve the Schrodinger equation for bound states. Instead of using basis functions (like harmonic oscillators or plane waves), the method discretizes the coordinate space into a grid and constructs the Hamiltonian.
<p align="center">


$$
\hat{H} =  \hat{T} + \hat{V}(x)  ...... (1)
$$


</p>

Where $\hat{H}$ = Hamiltonia operator, $\hat{T}$ = Kinetic energy operator, and V(x) = Potential energy. In this method, the kinetic energy matrix is computed using Fourier transformation, which efficiently connects the momentum space and coordinate space representations and finally the form of the hamiltonian becomes,
<p align="center">


$$
\hat{H}_{ij}^0 = \frac{2}{N} \sum_{l=1}^{n} (cos(\frac{l2\pi(i-j)}{N})\frac{\hbar^2}{2m}(lΔk)^2) + V(x_i)δ_{ij}  ...... (3)
$$


</p>

The delails about the FGH is written in a pdf format and uploaded in this GitHub page (**FGH.pdf**).

# Example 1: Morse potential of H<sub>2</sub>
In this example the FGH is solved using the morse potential of H<sub>2</sub> in Fortran (**FGH_Morse.f90**). The form of the morse potential is,
<p align="center">


$$
V(r) = D_e (1 - e^{-β(r-r_0)})^2 ...... (2)
$$


</p>

Where, D<sub>e</sub> = The dissociation energy, r<sub>0</sub> = equilibrium bond distance, r = bond distance, and β = range of parameter related to the potential's curvature. For H<sub>2</sub> the chosen parameters are, D<sub>e</sub> = 0.1744 a.u., r<sub>0</sub> = 1.40201a.u., and β = 1.02764 a.u. Beside the morse parameters, the total distance (X<sub>min</sub>-X<sub>max</sub>) and the number of grid points are taken from 0-2 and 10 respectively for the FGH calculation.\
For more detalis: https://pubs.aip.org/aip/jcp/article/91/6/3571/221213/The-Fourier-grid-Hamiltonian-method-for-bound

# Example 2: Double well potential of NH<sub>3</sub>
In this example the FGH is solved using the double well potential of NH<sub>3</sub> in Fortran (**FGH_DWP.f90**). Due to the rapid transformation of ammonia's equilibrium geometry (C<sub>3V</sub>) to its' mirror-image strucrure (which in invarient and same C<sub>3V</sub> point group) via ammonia's umbrella inversion mode, ammonia is trapped inside a symmetric double well potential (denoted as DWP in the rest of the README). This rapid transformation is occured through a planner sp<sup>2</sup> hybridiged geometry (D<sub>3h</sub>) which is responsiple for the potential barrier in this DWP. The form of the potential is,
<p align="center">


$$
V(x) = \frac{V_{max}}{b^4}((x-\frac{a}{2})^2 - b^2)^2  ...... (1)
$$


</p>

Where V<sub>max</sub> = Potential at the central maxima, ($\frac{a}{2}$ ± b) = Location of the two local minima, $\frac{a}{2}$ = Location of the central maxima.

The potential barrier in this DWP satisfies the conditions necessary for quantum tunneling, leading to transmission through the barrier. This, in turn, causes the splitting of the localized energy levels within each well of the DWP.

