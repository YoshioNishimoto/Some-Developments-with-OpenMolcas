# Some-Developments-with-OpenMolcas

I'm developing analytic derivatives of RASPT2 in OpenMolcas. I'm planning to merge the development in the main branch of OpenMolcas in the future, if possible, but for the moment, I will update here only. Note that this repository is just a development snapshot. Some debug print may be shown, and the actual code is still messy. At present, this development is based on OpenMolcas 19.11.

At present, the single-state CASPT2 and RASPT2 can be performed. The density-fitting may be used. Either real or imaginary level shift may be used. The IPEA shift is not yet; I think 99% OK but something is missing? At least, the SADREF keyword is required, if the IPEA is used. The performance is still poor, in particular with a large number of atomic orbitals (or impossible to complete). If the reference RASSCF averages more than one state, please add "SADREF" in &CASPT2. It uses the state-averaged density matrix for generating the one-electron Fock matrix. The energy is slightly different from that without the keyword. The gradient is implemented only with and without this keyword, but it is recommended that the keyword is used. No symmetry constraints can be employed.

As a byproduct, analytic gradients for state-averaged RASSCF can be computed. I implemented diagonal preconditioning for the active--active orbital rotation and fixed several things. However, no symmetry constrains can be employed, again.

Reference:

Nishimoto, Y. "Analytic Gradients for Restricted Active Space Second-order Perturbation Theory (RASPT2)" DOI: 10.26434/chemrxiv.14197457

***

Some history:

May 13, 2021: An old patch for RASPT2 gradient

March 22, 2021: State-specific density matrix (may work without the SADREF keyword)

March 19, 2021: Possibly RASPT2 (without the diagonal approximation)

March 11, 2021: Initial commit (RASPT2-D)

***

I think I will not make my previous NEVPT2 implementation in GAMESS-US publicly available for some reasons.
