.. index::
   single: Reaction field
   single: Cavity
   single: Solvent

.. _TUT\:sec\:cavity:

Solvent models
==============

.. only:: html

  .. contents::
     :local:
     :backlinks: none

For isolated molecules of modest size the *ab initio* methods
have reached great accuracy at present both for ground and
excited states. Theoretical studies on isolated molecules, however,
may have limited value to bench chemists since most of the
actual chemistry takes place in a solvent. If solute--solvent interactions are
strong they may have a large impact on the electronic structure of a
system and then on its excitation spectrum, reactivity, and properties.
For these reasons, numerous models
have been developed to deal with solute--solvent interactions in *ab
initio* quantum chemical calculations. A microscopic
description of solvation effects can be obtained by a supermolecule
approach or by combining statistical mechanical simulation techniques
with quantum chemical methods.
Such methods, however, demand expensive computations. By contrast, at the
phenomenological level, the solvent can be regarded as a dielectric continuum,
and there are a number of approaches :cite:`Cossi:98,Cossi:01,Karlstroem:88,Serrano:97b,Tomasi:94`
based on the classical reaction field concept.

|molcas| can model the solvent
within the framework of SCF, RASSCF and CASPT2 programs, for the calculation of energies
and properties and also for geometry optimizations. The reaction field formalism
is based on a sharp partition of the system: the solute molecule (possibly
supplemented by some explicit solvent molecules) is placed in a cavity
surrounded by a polarizable dielectric.
The surrounding is characterized mainly by its dielectric constant and density:
an important parameter of the method is the size of the cavity;
the dielectric medium is polarized by the solute, and this polarization creates
a reaction field which perturbs the solute itself.

Two versions of the model are presently available: one is based on the Kirkwood model
:cite:`Karlstroem:88,Serrano:97b` and uses only spherical cavities; the other is
called PCM (polarizable continuum model) :cite:`Cossi:98,Cossi:01` and can
use cavities of general shape, modeled on the actual solute molecule. In the former
case, the reaction field is computed as a truncated multipolar expansion and added
as a perturbation to the one-electron Hamiltonian; in the latter case the reaction
field is expressed in terms of a collection of apparent charges (solvation charges)
spread on the cavity surface: the PCM reaction field perturbs both one- and
two-electron Hamiltonian operators. In both cases, the solvent effects can be
added to the Hamiltonian at any level of theory, including MRCI and CASPT2.

.. index::
   single: Kirkwood model

Kirkwood model
--------------

This version of the model only uses spherical cavities. In addition,
it includes Pauli repulsion due to the medium by introducing a repulsive
potential representing the exchange repulsion between the solute and the solvent.
This is done by defining a penalty function of Gaussian type, generating
the corresponding spherical well integrals, and adding them to the one-electron
Hamiltonian. When the repulsion potential is used, the size of the cavity should
be optimized for the ground state of the molecule (see below). If the repulsive
potential is not used and the cavity size is chosen to be smaller (molecular
size plus van der Waals radius as is the usual choice in the literature)
one must be aware of the consequences: larger solvent effects but also
an unknown presence of molecular charge outside the boundaries of the
cavity. This is not a consequence of the present model but it is a general
feature of cavity models :cite:`Serrano:97b`.

.. Comment: As long as this does not work correctly in the code we should not mention it here!

   Dielectric cavity models, in general, assume equilibrium between the
   electronic state of the solute and the reaction field. This condition is
   not fulfilled for an electronic excitation. Therefore, to compute excited
   states |molcas| introduces the time dependence into the :program:`RASSCF` program
   by partitioning the reaction field factor
   into two parts, a slow and a fast component. The slow
   field follows geometrical changes of the solute and is, in practice,
   determined by the properties of the ground state.
   The fast component can be considered as
   the instantaneous electronic polarization that follows the absorption of
   a photon and the strength of this fast component is proportional to the
   square of the refractive index, that is, the optical value of the
   dielectric constant :cite:`Serrano:97b`.

.. index::
   single: PCM

PCM
---

The cavities are defined as the envelope of spheres centered on solute
atoms or atomic groups (usually hydrogens are included in the same sphere
of the atoms they are bonded to). Two selection of radii are presently available,
i.e. Pauling radii, and the so-called UATM (united atom topological model) radii:
the latter is the default for PCM calculations; sphere radii can also be provided
by the user in the input file.
The solvation charges are placed
in the middle of small tiles (*tesserae*) drawn on the surface; the number of
solvation charges can be gauged by changing the average area of tesserae (keyword
AAre in :program:`SEWARD`).

The program prints some information related to the cavity, where one should
always check carefully the magnitude of sphere radii: the program adjusts them
automatically to the solute topology (each radius depends on
hybridization, bonds, etc.), and sometimes this causes some
problems (for instance, discontinuities could appear during the scan
of a potential energy surface): if this happens, it is preferable to provide
the desired radii in the input file, so that they will be kept at all
geometries.

When doing state-average RASSCF calculations, one has to specify which root is
to be used to generate the solvation charges: this means that the PCM reaction
field will be in equilibrium with a specific electronic state, while it
perturbs all the states included in the calculation.

In electronic transitions (e.g. photon absorption or emission) one has to include
non-equilibrium effects, due to the finite relaxation time of solvent molecules
following a sudden change in electronic distribution. This is done by partitioning the
reaction field in two components (fast and slow, the former always equilibrated,
the latter delayed), whose magnitude is determined by the static dielectric constant
and by a "fast" dielectric constant :cite:`Cossi:00` (for very fast processes, like
photon absorption, the fast constant is equal to the square of the refraction index).
To perform a non-equilibrium calculation, for example to study a ground-to-excited state
transition, one has to perform a regular calculation at equilibrium for the ground state,
followed by a calculation for the excited state specifying the keyword NONEQ in the
:program:`RASSCF` program. Failing to include the keyword NONEQ will cause the program
to compute equilibrium solvation also for the excited state, what would be appropriate
for an adiabatic, instead of a vertical, transition.

CASPT2 calculations can be performed as usual for isolated molecules, specifying
the keyword RFPERT. Geometry optimizations can be performed as usual: note that
the arrangement of solvation charges around the solute molecule is likely
to break the molecular symmetry. If the symmetry was explicitly requested in
:program:`SEWARD`, the system will keep it through the optimization even in the
presence of the solvent, otherwise the convergence could be more difficult, and
the final geometry could result of a lower symmetry.

.. index::
   single: Solvent
   single: Solvent; Ground state
   single: Kirkwood model

Calculation of solvent effects: Kirkwood model
----------------------------------------------

We begin by performing a CASSCF/CASPT2 reaction field calculation
on the ground state of a molecule.

.. compound::

  To use the Kirkwood model, the keyword ::

    REACtion field

  is needed; if no repulsive potential is going to be used the input
  simply consists in adding the
  appropriate data (dielectric constant of the medium, cavity size, and
  angular quantum number of the highest multipole moment of the charge
  distribution) into the :program:`SEWARD` input: ::

    &SEWARD &END
    ...
    ...
    RF-Input
    Reaction field
    80 8.0 4
    End of RF-Input
    ...
    ...
    End of Input

This will compute the reaction field at those levels. The dielectric
constant 80.0 correspond to water as solvent. The radius of the cavity
is 8.0 in atomic units. Finally 4 is the maximum angular moment
number used in the multipole expansion. The cavity origin is the
coordinate origin, thus the molecule must be placed accordingly.

If we want to include the reaction field (either PCM or Kirkwood
model)
at other levels of theory the keyword :kword:`RFPErt` must be
added to the :program:`MOTRA` or :program:`CASPT2` inputs.

We are, however, going to explain the more complicated situation where
a repulsive well potential has to be added to the model. In this case
it is convenient to optimize the size of the cavity, although
in so doing we obtain large cavity sizes and therefore
smaller solvent effects. More realistic results can be obtained if
additional and specific solvent molecules are added inside the cavity.

.. index::
   single: Repulsive well
   single: SEWARD; Well integrals

To define the well potential we have to add the keyword
:kword:`WELL Integrals` to the :program:`SEWARD` input to compute and add the
Pauli repulsion integrals to the bare Hamiltonian.

The requirements considered to build this potential
are that it shall reproduce solvation energies for spherical particles,
ions, and that it must be wide enough so that the electrons in the
excited state of the molecules are also confined to the cavity. Negative ions
have the property that their electrons are loosely bound and they are thus suited
for parametrizing the repulsive potential.
The final result of different calibration calculations
:cite:`Bernhardsson:96a,Serrano:97b` is a penalty function which includes four
Gaussians. If :math:`a` is
the radius of the cavity the Gaussians are placed at distances :math:`a+2.0`,
:math:`a+3.0`, :math:`a+5.0` and :math:`a+7.0` |a0| from the cavity's center with
exponents 5.0, 3.5, 2.0 and 1.4, respectively.

.. index::
   single: DMABN

As an example we will use the N,N-dimethylaminobenzonitrile (DMABN) molecule
(see :numref:`fig:dmabn`).
This is a well known system with large dipole moments both in ground and
excited states which suffer important effects due to the polar environment.

.. figure:: dmabn.*
   :name: fig:dmabn
   :width: 50%
   :align: center

   N,N-dimethylaminobenzonitrile (DMABN)

.. extractfile:: advanced/RF.input

  &SEWARD &END
  Title
  para-DMABN molecule. Cavity size: 10 au.
  Symmetry
   X XY
  Basis set
  N.ANO-S...3s2p1d.
  N1             0.0000000000        0.0000000000        4.7847613288
  N2             0.0000000000        0.0000000000       -8.1106617786
  End of basis
  Basis set
  C.ANO-S...3s2p1d.
  C1             0.0000000000        0.0000000000        2.1618352923
  C2             0.0000000000        2.2430930886        0.7747833630
  C3             0.0000000000        2.2317547910       -1.8500321252
  C4             0.0000000000        0.0000000000       -3.1917306021
  C5             0.0000000000        0.0000000000       -5.9242601761
  C6             0.0000000000        2.4377336900        6.0640991723
  End of basis
  Basis set
  H.ANO-S...2s.
  H1             0.0000000000        4.0043085530       -2.8534714086
  H2             0.0000000000        4.0326542950        1.7215314260
  H3             0.0000000000        2.1467175630        8.0879851846
  H4             1.5779129980        3.6622699270        5.5104123361
  End of basis

  RF-Input
  reaction field
  38.8 10.0 4
  End of RF-Input

  Well Int
  4
  1.0 5.0 12.0
  1.0 3.5 13.0
  1.0 2.0 15.0
  1.0 1.4 17.0
  End of Input

  &SCF &END
  TITLE
   DMABN molecule
  OCCUPIED
  20 2 12 5
  ITERATIONS
  50
  END OF INPUT

  &RASSCF &END
  TITLE
   p-DMABN
  SYMMETRY
      1
  SPIN
      1
  NACTEL
     10    0    0
  FROZEN
      8    0    3    0
  INACTIVE
     12    1    9    1
  RAS2
      0    2    0    7
  THRS
  1.0E-06,1.0E-03,1.0E-03
  ITER
  50,25
  LUMORB
  END OF INPUT

.. index::
   single: Cavity

In the :program:`SEWARD` input the :kword:`WELL Integrals` must include
first the number of Gaussians used (four), followed by the
coefficient and exponent of the Gaussian and the radius of
the cavity in the sequence explained above: first the most
compact Gaussian with the radius plus 2.0 |a0|, and so on
to the least compact Gaussian.
Here, we have defined a cavity size of 10 |a0|
(cavity centered at coordinate origin). The RASSCF program will read
the RCTFLD input, prepared this time for acetonitrile
(:math:`\epsilon = 38.8`), a cavity size of 10.0 |a0| (the same
as in the SEWARD input) and a multipole expansion up to the fourth order
which is considered sufficient :cite:`Serrano:97b`.
The active space includes the :math:`\pi` space over
the molecular plane, excluding the :math:`\pi` orbital of the :math:`\ce{CN}`
group which lies in the molecular plane.

We repeat the calculation for different cavity sizes in order
to find the radius which gives the lowest absolute energy at the CASSCF
level. The presence of the repulsive terms allows the cavity
radius to be computed by energy minimization. For the calculations
using different cavity sizes it is not necessary to repeat the
calculation of all the integrals, just those related to the
well potential. Therefore, the keyword :kword:`ONEOnly` can be
included in the SEWARD input. The :file:`ONEINT` file will be modified
and the :file:`ORDINT` file is kept the same for each molecular
geometry. The energies obtained are in :numref:`tab:cav1`.

.. table:: Ground state CASSCF energies for DMABN with different cavity sizes.
   :name: tab:cav1

   ============= ======================
   Radius (|a0|) CASSCF energies (|Eh|)
   ============= ======================
   no cav.       |-|\455.653242
   10.0          |-|\455.645550
   11.0          |-|\455.653486
   12.0          |-|\455.654483
   14.0          |-|\455.654369
   16.0          |-|\455.654063
   ============= ======================

Taking the gas-phase value (no cav.) as the reference, the CASSCF
energy obtained with a 10.0 |a0| cavity radius is higher. This is an
effect of the repulsive potential, meaning that the molecule is too close to
the boundaries. Therefore we discard this value and use the values
from 11.0 to 16.0 to make a simple second order fit and obtain a
minimum for the cavity radius at 13.8 |a0|.

.. index::
   single: Cavity

Once we have this value we also need to optimize the position of the
molecule in the cavity. Some parts of the molecule, especially those
with more negative charge, tend to move close to the
boundary. Remember than the sphere representing the cavity has
its origin in the cartesian coordinates origin. We use the radius of
13.8 |a0| and compute the CASSCF energy at different displacements
along the coordinate axis. Fortunately enough, this molecule has
:math:`C_{2v}` symmetry. That means that displacements along two of the
axis (:math:`x` and :math:`y`) are restricted by symmetry. Therefore it is
necessary to analyze only the displacements along the :math:`z` coordinate.
In a less symmetric molecule all
the displacements should be studied even including combination of the displacements.
The result may even be a three dimensional net, although no
great accuracy is really required. The results for DMABN in |Ctv| symmetry are compiled
in :numref:`tab:cav2`.

.. table:: Ground state CASSCF energies for different translations with respect to the initial position of of the DMABN molecule in a 13.8 |a0| cavity.
   :name: tab:cav2

   ========================= ======================
   Disp. in :math:`z` (|a0|) CASSCF energies (|Eh|)
   ========================= ======================
   +0.5                      -455.654325
   0.0                       -455.654400
   |-|\0.5                   -455.654456
   |-|\1.0                   -455.654486
   |-|\1.5                   -455.654465
   ========================= ======================

Fitting these values to a curve we obtain an optimal displacement of |-|\1.0 |a0|. We move the
molecule and reoptimize the cavity radius at the new position of the
molecule. The results are listed in :numref:`tab:cav3`.

.. table:: Ground state CASSCF energies for DMABN with different cavity sizes. The molecule position in the cavity has been optimized.
   :name: tab:cav3

   ============= ======================
   Radius (|a0|) CASSCF energies (|Eh|)
   ============= ======================
   11.8          |-|\455.653367
   12.8          |-|\455.654478
   13.8          |-|\455.654486
   ============= ======================

There is no significant change. The cavity radius is then selected as 13.8 |a0| and the
position of the molecule with respect to the cavity is kept as in
the last calculation. The calculation is carried out with the new values.
The SCF or RASSCF outputs will contain the information about the
contributions to the solvation energy. The CASSCF energy obtained will
include the reaction field effects and an analysis of the
contribution to the solvation energy for each value of the multipole
expansion: ::

        Reaction field specifications:
        ------------------------------

         Dielectric Constant :                  .388E+02
         Radius of Cavity(au):                  .138E+02
         Truncation after    :                 4

        Multipole analysis of the contributions to the dielectric solvation energy

        --------------------------------------
           l             dE
        --------------------------------------
           0            .0000000
           1           -.0013597
           2           -.0001255
           3           -.0000265
           4           -.0000013
        --------------------------------------

.. Note: contains a nbsp

.. Comment: As long as this does not work correctly in the code we should not mention it here!

  .. index::
     single: Excited states
     single: Solvent; Excited states

  Solvent effects on excited states
  ---------------------------------

  As was explained above, the calculation of the solvent effects
  on excited states must be done in a different way.
  The inclusion of the solvent effects on the ground state of one molecule leads to
  a relaxation both for nuclei and electrons. It is considered to be in equilibrium
  with the environment. However, when an absorption or emission occurs, a dynamic process
  starts. It is usually considered that the new electronic distribution of the final
  state polarizes the solvent, which affects again in a different way the electronic part
  of the molecule than in the ground state case. This is considered to be a fast motion,
  where the equilibrium is achieved by solvent and solute electrons. The situation is not
  the same for nuclei. Slow motions are included in this case, and the nuclei do not
  achieve the equilibrium with solvent. These slow motions and the related
  interactions remain as in the ground state. Therefore, the calculations must be
  carried out using this sequence of steps:

  * A normal reaction field calculation is done for the ground state. However, we
    have included in the RASSCF input the keyword

    .. index::
       single: RASSCF
       single: RASSCF; HFRctFld
       single: Reaction field

    ::

      HFrctfld
      1 1 1.796

    The first number, one, indicates that this is a ground (or initial) state calculation
    and the slow part or the reaction field is stored in the file :file:`COMFILE`.
    The second number is the number of the computed root, one as we are interested in the
    ground state at this stage. The third number is the square of the refractive index of
    the solvent, also called the dielectric constant at high frequencies or optical value
    (from an approximate expression of the Clausius--Mossoti formulae). The RF-Input section in the
    integral input
    must be included also, as in the previous case: ::

      &SEWARD &END
      ...
      ...
      RF-Input
      reaction field
      38.8 13.8 4
      end of RF-Input
      ...
      ...
      End of Input

  * .. compound::

      The RASSCF program will store in the communications file, :file:`COMFILE`,
      not only the full reaction field of the computed state (ground state in this case),
      but also the slow component of the reaction field.
      The partitioning of the reaction field is performed for each individual :math:`2^l`\-pole
      component. The slow part is calculated for the ground state wave function
      as the (:math:`c_l-c_l^\infty/c_l`) fraction of the reaction field,
      where :math:`c_l` is the reaction field factor

      .. math:: c_l = -\frac{l!(l+1)(\epsilon-1)}{(l+1)\epsilon+l}\frac{1}{a^{2l+1}},

      :math:`c_l^\infty` is the reaction field
      factor using :math:`\epsilon^\infty\approx n^2`, :math:`n` is the refraction
      index for some suitable frequency in the infrared or visible range,
      and :math:`M_l^m` is the component :math:`m` of a :math:`2^l`\-pole moment of order :math:`l`.
      This value, :math:`F_{l,\text{slow}}^m`, is stored in the communications
      file together with the corresponding
      self-energy used to create the cavity. These fields in the :file:`COMFILE` file
      will not be replaced until we perform another RASSCF calculation using the
      same :file:`COMFILE` and the first value in the keyword :kword:`HFRCtfld`
      set to one.

  * .. compound::

      The RASSCF input for the excited state must contain ::

        HFrctfld
        0 2 1.796

      where the zero indicates this is an excited state calculation, the two is the
      number of the computed root, and 1.796 the square of the refractive index of the solvent.
      The same RCTFLD namelist is used as before.

    The total reaction field of the excited state will be added to the :file:`COMFILE` file,
    which can be then used as a perturbation on other program for that particular excited state.
    The slow part of the reaction field of the initial state remains in the :file:`COMFILE`.

    .. index::
       single: Option; RFPert

  * Once the RASSCF calculation is finished one can use other programs such as MRCI
    or CASPT2 where the CASSCF reaction field is introduced as a perturbation by
    including the :kword:`RFPErt` in their inputs (in the MOTRA input for MRCI).

  Remember, however, that each time an SCF or RASSCF calculation is done the
  resulting total reaction field is stored in the communications file (:file:`COMFILE`).
  Therefore, one must control which file is being used. The following is an easy sequence
  for the calculations: ::

    &SEWARD &END
    Title
    para-DMABN molecule. Cavity size: 13.8 au.
    Symmetry
     X XY
    Basis set
    N.ANO-S...3s2p1d.
    N1             0.0000000000        0.0000000000        3.7847613288
    N2             0.0000000000        0.0000000000       -9.1106617786
    End of basis
    Basis set
    C.ANO-S...3s2p1d.
    C1             0.0000000000        0.0000000000        1.1618352923
    C2             0.0000000000        2.2430930886       -0.2252166370
    C3             0.0000000000        2.2317547910       -2.8500321252
    C4             0.0000000000        0.0000000000       -4.1917306021
    C5             0.0000000000        0.0000000000       -6.9242601761
    C6             0.0000000000        2.4377336900        5.0640991723
    End of basis
    Basis set
    H.ANO-S...2s.
    H1             0.0000000000        4.0043085530       -3.8534714086
    H2             0.0000000000        4.0326542950        0.7215314260
    H3             0.0000000000        2.1467175630        7.0879851846
    H4             1.5779129980        3.6622699270        4.5104123361
    End of basis

    RF-Input
    reaction field
    38.8 13.8 4
    End of RF-Input

    Well Int
    4
    1.0 5.0 15.8
    1.0 3.5 16.8
    1.0 2.0 18.8
    1.0 1.4 20.8
    End of Input

    &SCF &END
    TITLE
     DMABN molecule
    OCCUPIED
    20 2 12 5
    ITERATIONS
    50
    END OF INPUT

    &RASSCF &END
    TITLE
     p-DMABN
    SYMMETRY
        1
    SPIN
        1
    NACTEL
       10    0    0
    FROZEN
        8    0    3    0
    INACTIVE
       12    1    9    1
    RAS2
        0    2    0    7
    THRS
    1.0E-06,1.0E-03,1.0E-03
    ITER
    50,25
    LUMORB
    HFRCtfld
    1 1 1.796
    END OF INPUT

    &CASPT2 &END
    Title
     Solvent
    Maxit
    20
    Lroot
    1
    RFPert
    End of Input

    &RASSCF &END
    TITLE
     p-DMABN
    SYMMETRY
        1
    SPIN
        1
    NACTEL
       10    0    0
    FROZEN
        8    0    3    0
    INACTIVE
       12    1    9    1
    RAS2
        0    2    0    7
    THRS
    1.0E-06,1.0E-03,1.0E-03
    CIROot
    1 2
    2
    Iterations
    50,25
    LumOrb
    Hfrctfld
    0 2 1.796
    End of Input

    &CASPT2 &END
    Title
     Solvent
    Maxit
    20
    Lroot
    2
    RFPert
    End of Input

  This sequence of calculations is correct because we perform
  the CASPT2 calculations immediately after the RASSCF calculations
  for each root. They will not be correct if we first compute the
  two RASSCF states and then the two CASPT2 corrections, because
  in the :file:`COMFILE` file the reaction field stored would be
  that from the last RASSCF calculation. The safer way is to save
  explicitly the :file:`COMFILE` file each time.

  One additional problem that can occur in excited state calculations
  is that the :program:`RASSCF` program does not converge easily for single root
  calculations. One may have to do a state-average CASSCF calculation.
  The program can compute the reaction field in average calculations
  because at each iteration it takes the density matrix of the root
  we specify under the keyword :kword:`Hfrctfld` to compute the reaction
  field, although the optimized density matrix is the averaged one.
  As this is only an approximation, the recommended procedure is to
  increase as much as possible the weight of the computed state
  in the average procedure. For instance, for the :math:`2^1A_1` state
  in DMABN the RASSCF input could be: ::

    &RASSCF &END
    TITLE
     p-DMABN molecule. 21A1 averaged state.
    SYMMETRY
        1
    SPIN
        1
    NACTEL
       10    0    0
    FROZEN
        8    0    3    0
    INACTIVE
       12    1    9    1
    RAS2
        0    2    0    7
    CIROOT
    2 2
    1 2
    1 9
    LEVSHFT
    0.5
    ITER
    50,25
    CIMX
    25
    LUMORB
    HFRctFld
    0 2 1.801
    END OF INPUT

  The experimental and theoretical excitation energies, and computed
  dipole moments for the ground state (GS) and :math:`\pi\to\pi^*` :math:`2^1A_1`
  state of DMABN in different media are shown in :numref:`tab:dmabn`
  (see ref. :cite:`Serrano:97b`):

  .. index::
     single: Excited states; DMABN
     single: Excited states; Solvent effects

  .. table:: Excitation energies and dipole moments for the
             :math:`2^1A_1` state of DMABN in different media.
     :name: tab:dmabn

     ===================================== =========================== =========================== =========================== ===========================
     Solvent                               PT2 (eV)                    Exp (eV)                    :math:`\mu_{\text{GS}}` (D) :math:`\mu_{\text{ES}}` (D)
     ===================================== =========================== =========================== =========================== ===========================
     Gas phase                             4.41                        ---                         6.6                         14.2
     Cyclohexane   (:math:`\epsilon=4.30`) 4.38                        4.4                         6.7                         14.6
     Butylchloride (:math:`\epsilon=9.65`) 4.32                        4.3                         6.9                         15.1
     Acetonitrile  (:math:`\epsilon=38.8`) 4.31                        4.2                         7.0                         15.2
     ===================================== =========================== =========================== =========================== ===========================

  The large radius of the cavity and the empty space within the cavity
  lead to a clear underestimation of the solvent effects. In this case
  however the model accounts for the most important aspects of the
  interaction.

  We now consider a difficult case. We want to compute both :math:`\pi\pi^*` and
  :math:`n\pi^*` excited states in DMABN and we want to use different active
  spaces in each case. For the :math:`\pi\pi^*` states (:math:`2^1A_1` and :math:`1^1B_2`)
  we use a (0207) space (|ao|\ |at|\ |bt|\ |bo|) with 10 active electrons and for the :math:`n\pi^*`
  states (:math:`1^1A_2` and :math:`1^1B_1`)
  a (1207) active space with 12 active electrons. This means we need
  two different ground state calculations. In principle we could re-optimize
  the cavity size for the new ground state but this is going to be a very
  minor effect. The following is a Korn shell script designed to do the
  explained calculations.

  .. index::
     single: Shell script

  ::

    #!/bin/ksh
    ################################################################################
    export Project=DMABN
    export Solvent=ACN
    export HomeDir=/u/$LOGNAME/$Project
    export TempDir=/temp/$LOGNAME
    export WorkDir=$TempDir/$RANDOM
    mkdir $WorkDir
    cd $WorkDir
    ln -fs  $TempDir/$Project.OrdInt ORDINT
    ln -fs  $TempDir/$Project.OneInt.$Solvent ONEINT
    #------------------------------------------------------------------------------#
    # Compute integrals                                                            #
    #------------------------------------------------------------------------------#
    molcas run seward $HomeDir/seward.$Project.$Solvent.input
    #------------------------------------------------------------------------------#
    # Compute RHF-SCF wavefunction                                                 #
    #------------------------------------------------------------------------------#
    ln -fs  $TempDir/$Project.$Solvent.ScfOrb SCFORB
    molcas run scf $HomeDir/scf.$Project.input
    rm SCFORB
    #------------------------------------------------------------------------------#
    # Compute CASSCF state functions and CASPT2 energy corrections for the state   #
    #------------------------------------------------------------------------------#
    #
    #------------------------------------------------------------------------------#
    # Ground State CASSCF and CASPT2 calculations                                  #
    #------------------------------------------------------------------------------#
    Name_list='11A1_pipi 11A1_npi'
    for Name in $Name_list;
    do
      ln -fs  $TempDir/$Project.$Solvent.ScfOrb                INPORB
      ln -fs  $TempDir/$Project.$Name.$Solvent.JobIph          JOBIPH
      ln -fs  $TempDir/$Project.$Name.$Solvent.RasOrb          RASORB
      ln -fs  $TempDir/$Project.$Name.$Solvent.ComFile         COMFILE
      molcas run rasscf $HomeDir/rasscf.$Name.$Solvent.input
      molcas run caspt2 $HomeDir/caspt2.rfpert.1.input
      rm RASORB
      rm INPORB
      rm JOBIPH
      rm COMFILE
    done
    #------------------------------------------------------------------------------#
    # Excited States CASSCF calculations                                           #
    #------------------------------------------------------------------------------#
    Name_list='21A1 11B2'
    for Name in $Name_list;
    do
      #----------------------------------------------------------------#
      #  Changing the name of the COMFILE                              #
      #----------------------------------------------------------------#
      cp $Project.11A1_pipi.$Solvent.ComFile $Project.$Name.$Solvent.ComFile
      #----------------------------------------------------------------#
      ln -fs  $TempDir/$Project.11A1.$Solvent.RasOrb           INPORB
      ln -fs  $TempDir/$Project.$Name.$Solvent.JobIph          JOBIPH
      ln -fs  $TempDir/$Project.$Name.$Solvent.RasOrb          RASORB
      ln -fs  $TempDir/$Project.$Name.$Solvent.ComFile         COMFILE
      molcas run rasscf $HomeDir/rasscfs.$Name.$Solvent.input
      rm RASORB
      rm INPORB
      rm JOBIPH
      rm COMFILE
    done
    Name_list='11A2 11B1'
    for Name in $Name_list;
    do
      #----------------------------------------------------------------#
      #  Changing the name of the COMFILE                              #
      #----------------------------------------------------------------#
      cp $Project.11A1_npi.$Solvent.ComFile $Project.$Name.$Solvent.ComFile
      #----------------------------------------------------------------#
      ln -fs  $TempDir/$Project.11A1.$Solvent.RasOrb           INPORB
      ln -fs  $TempDir/$Project.$Name.$Solvent.JobIph          JOBIPH
      ln -fs  $TempDir/$Project.$Name.$Solvent.RasOrb          RASORB
      ln -fs  $TempDir/$Project.$Name.$Solvent.ComFile         COMFILE
      molcas run rasscf $HomeDir/rasscfs.$Name.$Solvent.input
      rm RASORB
      rm INPORB
      rm JOBIPH
      rm COMFILE
    done
    #------------------------------------------------------------------------------#
    # Excited States CASPT2 calculations                                           #
    #------------------------------------------------------------------------------#
    Name_list='21A1'
    for Name in $Name_list;
    do
      ln -fs  $TempDir/$Project.$Name.$Solvent.JobIph          JOBIPH
      ln -fs  $TempDir/$Project.$Name.$Solvent.ComFile        COMFILE
      molcas run caspt2 $HomeDir/caspt2.rfpert.2.input
      rm JOBIPH
      rm COMFILE
    done
    Name_list='11B2'
    for Name in $Name_list;
    do
      ln -fs  $TempDir/$Project.$Name.$Solvent.JobIph          JOBIPH
      ln -fs  $TempDir/$Project.$Name.$Solvent.ComFile        COMFILE
      molcas run caspt2 $HomeDir/caspt2.rfpert.1.input
      rm JOBIPH
      rm COMFILE
    done
    Name_list='11A2 11B1'
    for Name in $Name_list;
    do
      ln -fs  $TempDir/$Project.$Name.$Solvent.JobIph          JOBIPH
      ln -fs  $TempDir/$Project.$Name.$Solvent.ComFile        COMFILE
      molcas run caspt2 $HomeDir/caspt2.rfpert.1.input
      rm JOBIPH
      rm COMFILE
    done
    #------------------------------------------------------------------------------#
    cd $TempDir
    rm -r $WorkDir
    exit

.. index::
   single: SEWARD; OneOnly

Notice that the two-electron integral file is the same
independent of the cavity size, well integrals
used or translational movement of the molecule.
Therefore, if the well parameters are
changed only the one-electron integral file :file:`ONEINT` need to be
recomputed using the option :kword:`ONEOnly` in the SEWARD
input. In the previous script we have named the :file:`ONEINT`
as :file:`$Project.OneInt.$Solvent`, but not because it depends on
the solvent but on the cavity radius which should be different
for each solvent when the well potential is used.

.. index::
   single: Solvent
   single: PCM

Solvation effects in ground states. PCM model in formaldehyde
-------------------------------------------------------------

The reaction field parameters are added to the
:program:`SEWARD` program input through the keyword ::

  RF-Input

.. compound::

  To invoke the PCM model the keyword ::

    PCM-model

  is required. A possible input is ::

    RF-input
    PCM-model
    solvent
    acetone
    AAre
    0.2
    End of rf-input

  which requests a PCM calculation with acetone as solvent, with tesserae
  of average area 0.2 Å\ |2|. Note that the default parameters are
  solvent = water, average area 0.4 Å\ |2|; see the :program:`SEWARD`
  manual section for further PCM keywords. By default the PCM adds
  non-electrostatic terms (i.e. cavity formation energy, and dispersion
  and repulsion solute-solvent interactions) to the computed free-energy
  in solution.

A complete input for a ground state CASPT2 calculation on formaldehyde
(:math:`\ce{H2CO}`) in water is

.. extractfile:: advanced/CASPT2.RF.input

  &GATEWAY
  Title
  formaldehyde
  Coord
  4

  H      0.000000    0.924258   -1.100293    Angstrom
  H      0.000000   -0.924258   -1.100293    Angstrom
  C      0.000000    0.000000   -0.519589    Angstrom
  O      0.000000    0.000000    0.664765    Angstrom
  Basis set
   6-31G*
  Group
   X Y
  RF-input
  PCM-model
  solvent
  water
  end of rf-input
  End of input

  &SEWARD
  End of input

  &SCF
  Title
  formaldehyde
  Occupied
  5 1 2 0
  End of input

  &RASSCF
  Title
  formaldehyde
  nActEl
  4 0 0
  Inactive
  4 0 2 0
  Ras2
  1 2 0 0
  LumOrb
  End of input

  &CASPT2
  Frozen
  4  0  0  0
  RFPErt
  End of input

.. Originally written by Francesco Aquilante

Solvation effects in excited states. PCM model and acrolein
-----------------------------------------------------------

.. compound::

  In the PCM picture, the solvent reaction field is
  expressed in terms of a polarization charge density :math:`\sigma(\vec{s})` spread
  on the cavity surface, which, in the most recent version of the method,
  depends on the electrostatic potential
  :math:`V(\vec{s})` generated by the solute on the cavity according to

  .. _PCM:

  .. math:: \left[ \frac{\epsilon+1}{\epsilon-1} \hat{S}
            -\frac{1}{2\pi}\hat{S}\hat{D}^* \right]
            \sigma(\vec{s}) =
            \left[ - 1 + \frac{1}{2\pi} \hat{D}\right] V(\vec{s})

  where :math:`\epsilon` is the solvent dielectric constant and
  :math:`V(\vec{s})` is the (electronic+nuclear) solute potential at point
  :math:`\vec{s}` on the cavity surface.
  The :math:`\hat{S}` and :math:`\hat{D}^*` operators are related respectively to
  the electrostatic potential :math:`V^\sigma({\vec{s}})`
  and to the normal component of the
  electric field :math:`E_\perp^\sigma(\vec{s})`
  generated by the surface charge density :math:`\sigma(\vec{s})`.
  It is noteworthy that in this PCM formulation the polarization charge
  density :math:`\sigma(\vec{s})` is designed to take into account implicitly
  the effects of the fraction of solute electronic density lying outside the
  cavity.

In the computational practice, the surface charge distribution
:math:`\sigma(\vec{s})` is expressed in terms of a set of point charges `\vec{q}` placed
at the center of each surface tessera, so that operators are replaced by the
corresponding square matrices.
Once the solvation charges (:math:`\vec{q}`) have been determined,
they can be used to compute energies and properties in solution.

.. compound::

  The interaction energy between the solute and the solvation charges
  can be written

  .. math:: E_{\text{int}} = \mat{V}^{\text{T}} \vec{q} = \sum_i^{N_{\text{TS}}} V_i q_i

  where :math:`V_i` is the solute potential calculated at the representative point
  of tessera :math:`i`. The charges act as
  perturbations on the solute electron density :math:`\rho`: since the charges
  depend in turn on :math:`\rho` through the electrostatic potential, the solute
  density and the charges must be adjusted until self consistency.
  It can be shown :cite:`Tomasi:94` that for any SCF procedure including a
  perturbation linearly depending on the electron density,
  the quantity that is variationally minimized corresponds to a free energy
  (i.e. :math:`E_{\text{int}}` minus the work spent to polarize the dielectric and to create
  the charges).
  If :math:`E^0=E[\rho^0] + V_{\text{NN}}` is the solute energy in vacuo, the free energy
  minimized in solution is

  .. math:: G = E[\rho] + V_{\text{NN}} + \frac{1}{2} E_{\text{int}}

  where :math:`V_{\text{NN}}` is the solute nuclear repulsion energy, :math:`\rho^0` is the
  solute electronic density for the isolated molecule, and :math:`\rho` is the
  density perturbed by the solvent.

The inclusion of non-equilibrium solvation effects, like those
occurring during electronic excitations, is introduced in the model by
splitting the solvation charge on each surface element into
two components: :math:`q_{i,\text{f}}` is the charge due to electronic (fast) component
of solvent polarization, in equilibrium with the solute electronic density
upon excitations, and :math:`q_{i,\text{s}}`, the charge arising from the orientational
(slow) part, which is delayed when the solute undergoes a sudden transformation.

The photophysics and photochemistry of
acrolein are mainly controlled by the relative position of the
:math:`^1(n\to\pi^*)`, :math:`^3(n\to\pi^*)` and :math:`^3(\pi\to\pi^*)` states, which is,
in turn, very sensitive to the presence and the nature of the solvent.
We choose this molecule in order to show an example of how to
use the PCM model in a CASPT2 calculation of vertical excitation
energies.

The three states we want to compute are low-lying singlet
and triplet excited states of the *s-trans* isomer.
The :math:`\pi` space (4 :math:`\pi` MOs, 4 :math:`\pi`\-electrons)
with the inclusion of the lone-pair MO (:math:`n_y`) is a suitable choice
for the active space in this calculation.
For the calculation in aqueous solution, we need first to compute the CASPT2
energy of the ground state in presence of the solvent water.
This is done by including in the :program:`SEWARD` input for the corresponding gas-phase
calculation the section ::

  RF-input

  PCM-model
  solvent
   water
  DIELectric constant
   78.39
  CONDuctor version
  AARE
   0.4

  End of rf-input

If not specified, the default solvent is chosen to be water.
Some options are available. The value of the dielectric constant
can be changed for calculations at temperatures other than 298 K.
For calculations in polar solvents like water, the use of the conductor
model (C-PCM) is recommended.
This is an approximation that employs conductor rather than dielectric
boundary conditions. It works very well for polar solvents
(i.e. dielectric constant greater than about 5), and is based
on a simpler and more robust implementation. It can be useful also in cases when
the dielectric model shows some convergence problems.
Another parameter that can be varied in presence of convergency problem
is the average area of the tesserae of which the surface of the cavity is composed.
However, a lower value for this parameter may give poorer results.

Specific keywords are in general needed for the other modules to work with PCM, except for
the SCF. The keyword :kword:`NONEquilibrium` is
necessary when computing excited states energies in :program:`RASSCF`.
For a state specific calculation of the ground state CASSCF energy, the solvent effects
must be computed with an equilibrium solvation approach, so this keyword must be omitted.
None the less, the keyword :kword:`RFpert` must be included in the CASPT2 input
in order to add
the reaction field effects to the one-electron Hamiltonian as a constant perturbation. ::

  &RASSCF &END
  Title
  Acrolein GS + PCM
  Spin
   1
  Symmetry
   1
  nActEl
   6 0 0
  Frozen
   4 0
  Inactive
   8 0
  Ras2
   1 4
  LUMORB
  THRS
  1.0e-06 1.0e-04 1.0e-04
  ITERation
   100 100
  End of input

  &CASPT2 &END
  Title
   ground state + PCM
  RFpert
  End of Input

Information about the reaction field calculation employing
a PCM-model appear first in the SCF output ::

  Polarizable Continuum Model (PCM) activated
  Solvent:water
  Version: Conductor
  Average area for surface element on the cavity boundary: 0.4000 Angstrom2
  Minimum radius for added spheres: 0.2000 Angstrom



  Polarized Continuum Model Cavity
  ================================

   Nord Group  Hybr  Charge Alpha Radius          Bonded to
     1   O     sp2   0.00   1.20  1.590   C   [d]
     2   CH    sp2   0.00   1.20  1.815   O   [d]  C   [s]
     3   CH    sp2   0.00   1.20  1.815   C   [s]  C   [d]
     4   CH2   sp2   0.00   1.20  2.040   C   [d]
   ------------------------------------------------------------------------------

The following input is used for the CASPT2 calculation of the :math:`^3A''(n\to\pi^*)` state.
Provided that the same $WorkDir has been using, which contains all the files of of the
calculation done for the ground state, the excited state calculation is done
by using inputs for the :program:`RASSCF`
and the :program:`CASPT2` calculations: ::

  &RASSCF &END
  Title
  Acrolein n->pi* triplet state + PCM
  Spin
   3
  Symmetry
   2
  nActEl
   6 0 0
  Frozen
   4 0
  Inactive
   8 0
  Ras2
   1 4
  NONEquilibrium
  LUMORB
  ITERation
   100 100
  End of input

  &CASPT2 &END
  Title
   triplet state
  RFpert
  End of Input

Note the :program:`RASSCF` keyword NONEQ, requiring that the slow part of the reaction
field be frozen as in the ground state, while the fast part is
equilibrated to the new electronic distribution. In this case the fast
dielectric constant is the square of the refraction index, whose value
is tabulated for all the allowed solvents (anyway, it can be modified by
the user through the keyword :kword:`DIELectric` in :program:`SEWARD`).

The :program:`RASSCF` output include the line: ::

  Calculation type: non-equilibrium (slow component from JobOld)

   Reaction field from state:            1

This piece of information means that the program computes the solvent
effects on the energy of the :math:`^3A''(n\to\pi^*)`
by using a non-equilibrium approach.
The slow component of the solvent response is kept frozen in terms of the charges
that have been computed for the previous equilibrium calculation of
the ground state. The remaining part of the solvent response,
due to the fast charges, is instead computed self-consistently for the
state of interest (which is state 1 of the specified spatial and spin symmetry in this case).

The vertical excitations to the lowest valence states
in aqueous solution for *s-trans* acrolein are
listed in the :numref:`tab:acrol-water` and compared
with experimental data.
As expected by qualitative reasoning, the vertical excitation energy to the
:math:`^1A''(n\to\pi^*)` state exhibits a blue shift in water.
The value of the vertical transition energy computed with the inclusion of the
PCM reaction field is computed to be 3.96 eV at the CASPT2 level of theory.
The solvatochromic shift is thus of +0.33 eV.
Experimental data are available for the
excitation energy to the :math:`^1A''(n\to\pi^*)` state. The band shift in
going from isooctane to water is reported to be +0.24 eV which is in fair
agreement with the PCM result.

No experimental data are available for the excitation energies to the triplet
states of acrolein in aqueous solution. However it is of interest to see how
the ordering of these two states depends on solvent effects.
The opposing solvatochromic shifts produced by the solvent on these two electronic transitions
place the two triplet states closer in energy.
This result might suggest that a dynamical interconversion between
the :math:`n\pi^*` and :math:`\pi\pi^*` may occur more favorable in solution.

.. table:: Vertical excitation energies/eV (solvatochromic shifts)
           of *s-trans* acrolein in gas-phase and in aqueous solution.
   :name: tab:acrol-water

   ====================================== =================== =================== ===================
   State                                  Gas-phase           Water               Expt.\ [#a]_
   ====================================== =================== =================== ===================
   :math:`^1A''(n_y\to\pi^*)`             3.63                3.96 (+0.33)        3.94 (+0.24)\ [#b]_
   :math:`T_1` :math:`^3A''(n_y\to\pi^*)` 3.39                3.45 (+0.06)
   :math:`T_2` :math:`^3A'(\pi\to\pi^*)`  3.81                3.71 (\ |-|\0.10)
   ====================================== =================== =================== ===================

.. [#a] Ref. :cite:`Forbes:59`
.. [#b] Solvatochromic shifts derived by comparison of the
        absorption wave lengths in water and isooctane
