# Some-Developments-with-OpenMolcas

I'm developing analytic derivatives of CASPT2 and RASPT2 in OpenMolcas. Some ongoing and incomplete developments may be pushed here. Note that this repository is just a development snapshot. Some debug print may be shown, and the actual code is still messy (including this README). This code is based on OpenMolcas v23.10.

At present, analytic first-order derivatives (gradient and derivative coupling vectors) for single-state and all multistate variants ([X]MS, XDW, and RMS) CASPT2 and RASPT2 can be computed. Most functions have to be combined with the density-fitting or Cholesky decomposition approximation. Either real or imaginary level shift may be used. The IPEA shift can also be used, but note that CASPT2/RASPT2 with the IPEA shift is not invariant with respect to rotations among active orbitals. The performance is still poor, in particular with a large number of atomic orbitals. No symmetry constraints can be employed. MPI-like parallel calculations are now supported, I think.

2023 codes are already publicly available as a part of OpenMolcas

Reference:

- Nishimoto, Y. "Analytic Gradients for Restricted Active Space Second-order Perturbation Theory (RASPT2)" The Journal of Chemical Physics 2021, 154, 194103. DOI: 10.1063/5.0050074
- Nishimoto, Y.; Battaglia, S.; Lind, R. "Analytic First-Order Derivatives of (X)MS, XDW, and RMS Variants of the CASPT2 and RASPT2 Methods" Journal of Chemical Theory and Computation 2022, 18, 4269--4281. DOI: 10.1021/acs.jctc.2c00301
- Nishimoto, Y. "Analytic first-order derivatives of CASPT2 with IPEA shift" The Journal of Chemical Physics 2023, 158, 174112. DOI: 10.1063/5.0147611
- Nishimoto, Y. "Analytic First-Order Derivatives of CASPT2 combined with the Polarizable Continuum Model" ChemRxiv DOI: 10.26434/chemrxiv-2024-05mq8

***

Some history:

December 24, 2024: SA-CASSCF/PCM and parallel CASPT2/PCM gradient, hopefully

August 8, 2023: CASPT2-IPEA gradient (I hope I correctly merged)

July 31, 2023: Removed the wired loop for MS-CASPT2 variants (preliminary)

July 31, 2023: Now, based on v23.02

May 13, 2021: An old patch for RASPT2 gradient

March 22, 2021: State-specific density matrix (may work without the SADREF keyword)

March 19, 2021: Possibly RASPT2 (without the diagonal approximation)

March 11, 2021: Initial commit (RASPT2-D)

***

I think I will not make my previous NEVPT2 implementation in GAMESS-US publicly available for some reasons.
