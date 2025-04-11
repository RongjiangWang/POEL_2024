FORTRAN code for simulating the diffusion and deformation process induced by pump tests in a layered poroelastic half-space.

Highlights:

(1) orthonormal propagator algorithm for numerical stability

(2) fully coupled diffusion and deformation system based on Biot’s poroelasticity theory

(3) wave-number-domain differential filter technique for fast convergence of numerical Hankel transform

For Windows user, the executable file is provided under folder "WindowsEXE". Linux user may compile the source codes with "gfortran" via a single command like, e.g.,

~>cd .../SourceCode

~>gfortran -o poel2024 *.f -O3

to get the excutable code poel2024.

After start the executable code, the program ask for an input file in the ASCII format. An example input file is provided under folder "InputFile". You may change the input data included in this file for your own applications.

References

Wang, R., and H.-J. Kümpel (2003), Poroelasticity: Efficient modelling of strongly coupled, slow deformation processes in a multilayered half-space, Geophysics, 68(2), 705-717.
