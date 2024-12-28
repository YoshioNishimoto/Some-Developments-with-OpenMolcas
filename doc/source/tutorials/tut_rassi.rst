.. index::
   single: Program; RASSI
   single: RASSI
   single: Properties; With RASSI
   single: Properties; Expectation values
   single: Properties; Matrix elements
   single: Expectation values
   single: Matrix elements

.. _TUT\:sec\:rassi:

:program:`RASSI` --- A RAS State Interaction Program
====================================================

Program :program:`RASSI` (RAS State Interaction) computes matrix elements
of the Hamiltonian and other operators in a wave function basis, which
consists of individually optimized CI expansions from the :program:`RASSCF`
program. Also, it solves the Schrödinger equation within the space of
these wave functions. There are many possible applications for such type
of calculations. The first important consideration to have into account
is that :program:`RASSI` computes the interaction among RASSCF states
expanding the same set of configurations, that is,
having the same active space size and number of electrons.

The :program:`RASSI` program is routinely used to compute electronic
transition moments, as it is shown in the Advanced Examples in the
calculation of transition dipole moments for the
excited states of the thiophene molecule using CASSCF-type wave functions.
By default the program will compute the matrix elements and expectation values
of all the operators for which :program:`SEWARD` has computed the integrals
and has stored them in the :file:`ONEINT` file.

.. index::
   single: Non-orthogonal states

RASSCF (or CASSCF) individually optimized states are interacting and
non-orthogonal. It is imperative when the states involved have different
symmetry to transform the states to a common eigenstate basis in such
a way that the wave function remains unchanged. The State Interaction
calculation gives an unambiguous set of non-interacting and orthonormal
eigenstates to the projected Schrödinger equation and also the
overlaps between the original RASSCF wave functions and the eigenstates.
The analysis of the original states in terms of RASSI eigenstates is
very useful to identify spurious local minima and also to inspect the
wave functions obtained in different single-root RASSCF calculations,
which can be mixed and be of no help to compare the states.

Finally, the :program:`RASSI` program can be applied in situations when
there are two strongly interacting states and there are two very different
MCSCF solutions. This is a typical situation in transition metal chemistry
when there are many close states associated each one to a configuration
of the transition metal atom. It is also the case when there are two
close quasi-equivalent localized and delocalized solutions. :program:`RASSI`
can provide with a single set of orbitals able to represent, for instance,
avoided crossings. :program:`RASSI` will produce a
number of files containing the natural orbitals for
each one of the desired eigenstates to be used in subsequent calculations.

:program:`RASSI` requires as input files the :file:`ONEINT` and :file:`ORDINT`
integral files and the :file:`JOBIPH` files from the :program:`RASSCF` program
containing the states which are going to be computed. The :file:`JOBIPH` files
have to be named consecutively as :file:`JOB001`, :file:`JOB002`, etc.
The input for the :program:`RASSI` module has to contain at least
the definition of the number of states available in each of the input
:file:`JOBIPH` files. :numref:`block:rassi_input` lists the input file
for the :program:`RASSI` program in a calculation including two :file:`JOBIPH`
files (2 in the first line), the first one including three roots (3 in the first
line) and the second five roots (5 in the first line). Each one of the
following lines lists the number of these states within each :file:`JOBIPH` file.
Also in the input, keyword :kword:`NATOrb` indicates that three files
(named sequentially :file:`NAT001`, :file:`NAT002`, and :file:`NAT003`) will
be created for the three lowest eigenstates.

.. index::
   single: RASSI; Input

.. code-block:: none
   :caption: Sample input requesting the :program:`RASSI` module to calculate the matrix
             elements and expectation values for eight interacting RASSCF states
   :name: block:rassi_input

   &RASSI
   NROFjobiph= 2 3 5; 1 2 3; 1 2 3 4 5
   NATOrb= 3

.. index::
   single: RASSI; Output

:program:`RASSI` Output
-----------------------

The :program:`RASSI` section of the |molcas| output is basically divided
in three parts. Initially, the program prints the information about the
:file:`JOBIPH` files and input file, optionally prints the wave functions,
and checks that all the configuration spaces are the same in all the
input states. In second place :program:`RASSI` prints the expectation
values of the one-electron operators, the Hamiltonian matrix, the
overlap matrix, and the matrix elements of the one-electron operators,
all for the basis of input RASSCF states. The third part starts with
the eigenvectors and eigenvalues for the states computed in
the new eigenbasis, as well as the overlap of the computed eigenstates
with the input RASSCF states. After that, the expectation values and
matrix elements of the one-electron operators are repeated on the
basis of the new energy eigenstates. A final section informs about
the occupation numbers of the natural orbitals computed by
:program:`RASSI`, if any.

In the Advanced Examples a detailed example of how to interpret
the matrix elements output section for the thiophene molecule is
displayed. The rest of the output is self-explanatory. It has to be
remembered that to change the default origins for the one electron
operators (the dipole moment operator uses the nuclear charge
centroid and the higher order operators the center of the nuclear
mass) keyword :kword:`CENTer` in :program:`GATEWAY` must be used.
Also, if multipoles higher than order two are required, the
option :kword:`MULTipole` has to be used in :program:`GATEWAY`.

The program :program:`RASSI` can also be used to compute a spin--orbit Hamiltonian
for the input CASSCF wave functions as defined above. The keyword :kword:`AMFI`
has to be used in :program:`SEWARD` to ensure that the corresponding integrals
are available.

.. code-block:: none
   :caption: Sample input requesting the :program:`RASSI` module to calculate and diagonalize
             the spin--orbit Hamiltonian the ground and triplet excited state in water.
   :name: block:rassi_input1

   &RASSI
   NROFjobiph= 2 1 1; 1; 1
   Spinorbit
   Ejob

The first :file:`JOBMIX` file contains the wave function for the ground state and
the second file the :math:`^3B_2` state discussed above. The keyword :kword:`Ejob`
makes the :program:`RASSI` program use the CASPT2 energies which have been
written on the :file:`JOBMIX` files in the diagonal of the spin--orbit
Hamiltonian. The output of this calculation will give four spin--orbit states and
the corresponding transition properties, which can for example be used to
compute the radiative lifetime of the triplet state.

:program:`RASSI` --- Basic and Most Common Keywords
---------------------------------------------------

.. class:: keywordlist

:kword:`NROFjob`
  Number of input files, number of roots, and roots for each file

:kword:`EJOB`/:kword:`HDIAG`
  Read energies from input file / inline

:kword:`SPIN`
  Compute spin--orbit matrix elements for spin properties
