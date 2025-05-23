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

Silverii, F., F. Maccaferri, G. Richter, B. G. Cansado, R. Wang, S. Hainzl and T. Dahm (2021). Poroelastic model in a vertically sealed gas storage: a case study from cyclic injection/production in a carbonate aquifer. Geophysical Journal International. https://doi.org/10.1093/gji/ggab268.

Flóvenz, Ó. G., R. Wang, G. P. Hersir, T. Dahm, S. Hainzl, M. Vassileva, V. Drouin, S. Heimann, M. P. Isken, E. Á. Gudnason, K. Ágústsson, Th. Ágústsdóttir, J. Horalek, M. Motagh, Th. R. Walter, E. Rivalta, Ph. Jousset , Ch. M. Krawczyk and C. Milkereit (2022). Cyclic unrest in a geothermal field as a harbinger of the Fagradalsfjall eruption, Iceland. Nature Geoscience 15(5). doi:10.1038/s41561-022-00930-5.
