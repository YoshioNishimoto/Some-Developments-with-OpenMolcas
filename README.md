# Some-Developments-with-OpenMolcas

I'm developing analytic derivatives of RASPT2 in OpenMolcas. I'm planning to merge the development in the main branch of OpenMolcas in the future, if possible, but for the moment, I will update here only. Note that this repository is just a development snapshot. Some debug print may be shown, and the actual code is still messy. At present, this development is based on OpenMolcas 19.11.

At present, the single-state CASPT2 and RASPT2 with the diagonal approximation (see the 1990 paper) can be performed. The density-fitting may be used. Either real or imaginary level shift may be used. The IPEA shift is not yet; I think 99% OK but something is missing. The performance is still poor, in particular with a large number of atomic orbitals (or impossible to complete). If the reference RASSCF averages more than one state, please add "SADREF" in &CASPT2. It uses the state-averaged density matrix for generating the one-electron Fock matrix. The energy is slightly different from that without the keyword, but the gradient is implemented only with this keyword (at least now). No symmetry constraints can be employed.

As a byproduct, analytic gradients for state-averaged RASSCF can be computed. I implemented diagonal preconditioning for the active--active orbital rotation and fixed several things. However, no symmetry constrains can be employed, again.

Reference:
Nishimoto, Y. "Analytic Gradients for Complete and Restricted Active Space Second-order Perturbation Theory within the Diagonal Approximation" DOI: 10.26434/chemrxiv.14197457

***

I think I will not make my previous NEVPT2 implementation in GAMESS-US publicly available for some reasons.
