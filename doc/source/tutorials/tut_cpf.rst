.. index::
   single: Program; CPF
   single: CPF
   single: MCPF
   single: ACPF

.. _TUT\:sec\:cpf:

:program:`CPF` --- A Coupled-Pair Functional Program
====================================================

The :program:`CPF` program produces Single and Doubles Configuration
Interaction (SDCI), Coupled-Pair Functional (CPF), Modified Coupled-Pair
Functional (MCPF), and Averaged Coupled-Pair Functional (ACPF) wave
functions (see CPF section of the user's guide) from one
reference configuration. The difference between the :program:`MRCI` and
:program:`CPF` codes is that the former can handle Configuration
Interaction (CI) and Averaged Coupled-Pair Functional (ACPF) calculations
with more than one reference configuration. For a closed-shell reference
the wave function can be generated with the :program:`SCF` program. In
open-shell cases the :program:`RASSCF` has to be used.

The :kword:`TITLe` keyword behaviors in a similar fashion to the
other |molcas| modules. The :kword:`CPF` keyword requests an
Coupled-Pair Functional calculation.
This is the default and is mutually exclusive with keywords
:kword:`MCPF`, :kword:`ACPF`, and :kword:`SDCI` which request different
type of calculations. The input below lists the input files
for the :program:`guga` and :program:`cpf` programs to obtain the MCPF
energy for the lowest triplet state of :math:`B_2` symmetry in the water molecule.
The :program:`GUGA` module computes the coupling coefficients for a triplet
state of the appropriate symmetry and the :program:`CPF` module will
converge to the first excited triplet state. One orbital of the first
symmetry has been frozen in this case (core orbital) in the :program:`MOTRA`
step.

:program:`cpf` Output
---------------------

The :program:`cpf` section of the output lists the number of each type
of orbital in each symmetry including pre-frozen orbitals that were
frozen by the :program:`guga` module. After some information concerning the
total number of internal configurations used and storage data, it appears
the single reference configuration in the :program:`mrci` format: an empty
orbital is listed as "``0``" and a doubly occupied as "``3``". The
spin of a singly occupied orbital by "``1``" (spin up) or "``2``"
(spin down). The molecular orbitals are listed near the end of the output.

Sample input requested by the GUGA and CPF modules to calculate the MCPF energy for
the lowest :math:`B_1` triplet state of the water in :math:`C_{2v}` symmetry: ::

  &GUGA
  Title= H2O molecule. Triplet state.
  Electrons= 8; Spin= 3
  Inactive= 2 0 1 0; Active= 1 1 0 0
  CiAll= 2

  &CPF
  Title= MCPF of triplet state of C2v Water
  MCPF

There are four input files to the :program:`cpf` module; :file:`CIGUGA`
from :program:`GUGA`, :file:`TRAONE` and :file:`TRAINT` from
:program:`MOTRA` and :file:`ONEINT` from :program:`SEWARD`. The orbitals
are saved in :file:`CPFORB`.
