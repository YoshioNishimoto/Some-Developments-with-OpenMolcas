# Some-Developments-with-OpenMolcas
Planning to implement analytic derivatives of RASPT2

I'm developing analytic derivatives of RASPT2 in OpenMolcas. I'm planning to merge these developments in the main branch of OpenMolcas in future, if possible, but for the moment, I will only update here. Note that this repository is just a development snapshot. Some debug print may be shown, and the actual code is still messy.

At present, single-state CASPT2 and RASPT2 with the diagonal approximation (see the 1990 paper) can be performed. The density-fitting may be used. Either real or imaginary level shift may be used. The IPEA shift is not yet; I think 99% OK but something is missing. The performance is still poor, in particular with a large number of atomic orbitals (or impossible to complete).

As a byproduct, analytic gradients for state-averaged RASSCF can be computed. I implemented diagonal preconditioning for the active--active orbital rotation and fixed several things. However, no symmetry constrains can be employed.

**********

I don't know how to put files. For the moment, I just put the zip file...
