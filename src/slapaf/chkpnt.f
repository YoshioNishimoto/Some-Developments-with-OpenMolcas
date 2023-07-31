************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2018, Ignacio Fdez. Galvan                             *
************************************************************************
*                                                                      *
      Module Chkpnt
      Implicit None
      Character(Len=9) :: basename = 'SLAPAFCHK'
      Character(Len=12) :: filename
      Integer :: chkpnt_id, chkpnt_iter, chkpnt_hess
      Integer :: chkpnt_ener, chkpnt_coor, chkpnt_force, chkpnt_new
      Integer :: Iter_all

      Contains
*                                                                      *
************************************************************************
*                                                                      *
      Subroutine Chkpnt_open()
#ifdef _HDF5_
      Use Symmetry_Info, Only: nIrrep
      Use Slapaf_Info, Only: Coor
      Use Slapaf_Parameters, Only: IRC, iter
      Use mh5, Only: mh5_is_hdf5, mh5_open_file_rw, mh5_open_attr,
     &               mh5_open_dset, mh5_fetch_attr
      Character(Len=3) :: level
      Logical :: create
      Integer :: tmp

      Iter_all = Iter
      create = .True.
      Call GetEnvF('EMIL_InLoop',level)
      If ((level.eq.'0').or.(level.eq.'1')) level=''
      filename = basename//Trim(level)

      If (((Iter.gt.1).or.(IRC.eq.-1)).and.mh5_is_hdf5(filename)) Then
        create = .False.
        chkpnt_id = mh5_open_file_rw(filename)
        chkpnt_iter = mh5_open_attr(chkpnt_id, 'ITERATIONS')
        chkpnt_ener = mh5_open_dset(chkpnt_id, 'ENERGIES')
        chkpnt_coor = mh5_open_dset(chkpnt_id, 'COORDINATES')
        chkpnt_new = mh5_open_dset(chkpnt_id, 'CENTER_COORDINATES')
        chkpnt_force = mh5_open_dset(chkpnt_id, 'FORCES')
        chkpnt_hess = mh5_open_dset(chkpnt_id, 'HESSIAN')
        Call mh5_fetch_attr(chkpnt_id, 'NSYM', tmp)
        If (tmp.ne.nIrrep) create = .True.
        Call mh5_fetch_attr(chkpnt_id, 'NATOMS_UNIQUE', tmp)
        If (tmp.ne.SIZE(Coor,2)) create = .True.
        Call mh5_fetch_attr(chkpnt_id, 'ITERATIONS', tmp)
        If (IRC.eq.-1) Then
          If (Iter.ge.tmp) create = .True.
          Iter_all = tmp+1
        Else
          If (tmp.ge.Iter) create = .True.
        End If
      End If
      If (create) Then
        If ((Iter.gt.1).or.(IRC.eq.-1)) Then
          Call WarningMessage(2,
     &         'The HDF5 file does not exist or is inconsistent')
          Call AbEnd()
        Else
          Call Chkpnt_init()
        End If
      End If
#endif
      End Subroutine Chkpnt_open
*                                                                      *
************************************************************************
*                                                                      *
      Subroutine Chkpnt_init()
#ifdef _HDF5_
      Use Phase_Info
      Use Symmetry_Info, Only: nIrrep
      Use Slapaf_Info, Only: dMass, nStab, iCoSet, AtomLbl, Smmtrc,
     &                       Coor
      Use Slapaf_Parameters, Only: nDimBC, rMEP, MEP, dMEPStep
      Use mh5, Only: mh5_create_file, mh5_init_attr,
     &               mh5_create_attr_int, mh5_create_dset_str,
     &               mh5_create_dset_real, mh5_create_dset_int,
     &               mh5_put_dset, mh5_close_dset
#  include "Molcas.fh"
#  include "stdalloc.fh"
      Character :: lIrrep(24)
      Integer :: dsetid, mAtom, i, j, k
      Real*8, Allocatable :: charges(:)
      Integer, Allocatable :: desym(:,:), symdof(:,:)

      chkpnt_id = mh5_create_file(filename)

      Call mh5_init_attr(chkpnt_id, 'MOLCAS_MODULE', 'SLAPAF')

*     symmetry information
      Call mh5_init_attr(chkpnt_id, 'NSYM', nIrrep)
      Call Get_cArray('Irreps', lIrrep, 24)
      Call mh5_init_attr(chkpnt_id, 'IRREP_LABELS',
     &                   1, [nIrrep], lIrrep, 3)

      Call mh5_init_attr(chkpnt_id, 'NATOMS_UNIQUE', SIZE(Coor,2))

      Call mh5_init_attr(chkpnt_id, 'DOF', nDimBC)

*     atom labels
      dsetid = mh5_create_dset_str(chkpnt_id,
     &         'CENTER_LABELS', 1, [SIZE(Coor,2)], LENIN)
      Call mh5_init_attr(dsetid, 'DESCRIPTION',
     &     'Unique center labels arranged as one [NATOMS_UNIQUE] block')
      Call mh5_put_dset(dsetid, AtomLbl)
      Call mh5_close_dset(dsetid)

*     atom masses
      dsetid = mh5_create_dset_real(chkpnt_id,
     &         'CENTER_MASSES', 1, [SIZE(Coor,2)])
      Call mh5_init_attr(dsetid, 'DESCRIPTION',
     &    'Nuclear masses, stored as array of size [NATOMS_UNIQUE]')
      Call mh5_put_dset(dsetid, dMass)
      Call mh5_close_dset(dsetid)

*     atom charges
      dsetid = mh5_create_dset_real(chkpnt_id,
     &         'CENTER_CHARGES', 1, [SIZE(Coor,2)])
      Call mh5_init_attr(dsetid, 'DESCRIPTION',
     &    'Nuclear charges, stored as array of size [NATOMS_UNIQUE]')
      Call mma_allocate(charges, SIZE(Coor,2))
      Call Get_dArray('Nuclear Charge', charges, SIZE(Coor,2))
      Call mh5_put_dset(dsetid, charges)
      Call mma_deallocate(charges)
      Call mh5_close_dset(dsetid)

*     number of iterations
      chkpnt_iter = mh5_create_attr_int(chkpnt_id,
     &              'ITERATIONS')

*     atom coordinates (new iteration)
*       use the same dataset name as in run2hdf5
      chkpnt_new = mh5_create_dset_real(chkpnt_id,
     &             'CENTER_COORDINATES', 2, [3,SIZE(Coor,2)])
      Call mh5_init_attr(chkpnt_new, 'DESCRIPTION',
     &     'Atom coordinates for new iteration, matrix of size '//
     &     '[NATOMS_UNIQUE,3], stored with atom index varying slowest')

      If (nIrrep.gt.1) Then

        mAtom = 0
        Do i=1,SIZE(Coor,2)
          mAtom = mAtom+nIrrep/nStab(i)
        End Do
        Call mma_allocate(desym, 4, mAtom)
        Call mma_allocate(symdof, 2, nDimBC)
        mAtom = 0
        k = 0
        Do i=1,SIZE(Coor,2)
          Do j=0,nIrrep/nStab(i)-1
            mAtom = mAtom+1
            desym(1,mAtom) = i
            desym(2,mAtom) = iPhase(1,iCoSet(j,i))
            desym(3,mAtom) = iPhase(2,iCoSet(j,i))
            desym(4,mAtom) = iPhase(3,iCoSet(j,i))
          End Do
          Do j=1,3
            If (.Not.Smmtrc(j,i)) Cycle
            k = k+1
            symdof(1,k) = i
            symdof(2,k) = j
          End Do
        End Do

*     total number of atoms
        Call mh5_init_attr(chkpnt_id, 'NATOMS_ALL', mAtom)

*     desymmetrization factors
        dsetid = mh5_create_dset_int(chkpnt_id,
     &           'DESYM_FACTORS', 2, [4,mAtom])
        Call mh5_init_attr(dsetid, 'DESCRIPTION',
     &      'Factors for obtaining all coordinates, matrix of size '//
     &      '[NATOMS_ALL,4], each row contains the unique atom index'//
     &      'and the factors with which to multiply the x,y,z '//
     &      'coordinates')
        Call mh5_put_dset(dsetid, desym)
        Call mh5_close_dset(dsetid)
        Call mma_deallocate(desym)

*     symmetry-unique degrees of freedom (Cartesian indices)
        dsetid = mh5_create_dset_int(chkpnt_id,
     &           'DOF_INDICES', 2, [2, nDimBC])
        Call mh5_init_attr(dsetid, 'DESCRIPTION',
     &      'Indices of the Cartesian degrees of freedom, matrix of '//
     &      'size [DOF, 2], each row contains the atom index and the '//
     &      'Cartesian index (1=x, 2=y, 3=z)')
        Call mh5_put_dset(dsetid, symdof)
        Call mh5_close_dset(dsetid)
        Call mma_deallocate(symdof)
      End If

*     iteration data:

*     energies
      chkpnt_ener = mh5_create_dset_real(chkpnt_id,
     &              'ENERGIES', 1, [0], dyn=.True.)
      Call mh5_init_attr(chkpnt_ener, 'DESCRIPTION',
     &     'Energies for all iterations as a matrix of size '//
     &     '[ITERATIONS]')

*     atom coordinates
      chkpnt_coor = mh5_create_dset_real(chkpnt_id,
     &              'COORDINATES', 3, [3,SIZE(Coor,2),0], dyn=.True.)
      Call mh5_init_attr(chkpnt_coor, 'DESCRIPTION',
     &     'Atom coordinates, matrix of size [ITERATIONS,'//
     &     'NATOMS_UNIQUE,3], stored with iteration varying slowest, '//
     &     'then atom index')

*     Cartesian forces (F = -g)
      chkpnt_force = mh5_create_dset_real(chkpnt_id,
     &               'FORCES', 3, [3,SIZE(Coor,2),0], dyn=.True.)
      Call mh5_init_attr(chkpnt_force, 'DESCRIPTION',
     &     'Cartesian forces, matrix of size [ITERATIONS,'//
     &     'NATOMS_UNIQUE,3], stored with iteration varying slowest, '//
     &     'then atom index')

*     Cartesian Hessian
      chkpnt_hess = mh5_create_dset_real(chkpnt_id,
     &              'HESSIAN', 1, [nDimBC*(nDimBC+1)/2])
      Call mh5_init_attr(chkpnt_hess, 'DESCRIPTION',
     &     'Cartesian Hessian in triangular form, as a vector of '//
     &     'size [DOF*(DOF+1)/2]')

*     MEP/IRC information
      If (MEP.or.rMEP) Then
        Call mh5_init_attr(chkpnt_id, 'MEP_STEP', dMEPStep)
        Call mh5_init_attr(chkpnt_id, 'MEP_ITERATIONS', 0)

        dsetid = mh5_create_dset_int(chkpnt_id,
     &           'MEP_INDICES', 1, [0], dyn=.True.)
        Call mh5_init_attr(dsetid, 'DESCRIPTION',
     &       'Iteration number for each converged MEP step, as a '//
     &       'vector of size [MEP_ITERATIONS]')
      End If
#endif

      End Subroutine Chkpnt_init
*                                                                      *
************************************************************************
*                                                                      *
      Subroutine Chkpnt_update()
#ifdef _HDF5_
      Use Slapaf_Info, Only: Cx, Gx, Energy
      Use Slapaf_Parameters, Only: nDimBC, iter
      Use mh5, Only: mh5_put_attr, mh5_resize_dset,
     &               mh5_put_dset
#  include "WrkSpc.fh"
#  include "stdalloc.fh"
      Integer :: i, j
      Logical :: Found
      Real*8, Allocatable :: Hss_X(:)

      Call Qpg_dArray('Hss_X',Found,i)
      If (Found) Then
        If (i.ne.nDimBC**2) Then
          Call WarningMessage(2,'Hessian with wrong dimension')
          Call AbEnd()
        End If
        Call mma_allocate(Hss_X,i)
        Call Get_dArray('Hss_X',Hss_X,i)
        Do i=1,nDimBC
          Do j=1,nDimBC
            Hss_X(i*(i-1)/2+j) = Hss_X(nDimBC*(i-1)+j)
          End Do
        End Do
      End If

*     iterations
      Call mh5_put_attr(chkpnt_iter, Iter_all)
*     energies
      Call mh5_resize_dset(chkpnt_ener, [Iter_all])
      Call mh5_put_dset(chkpnt_ener,
     &     Energy(Iter:Iter), [1], [Iter_all-1])
*     coordinates
      Call mh5_resize_dset(chkpnt_coor, [3,SIZE(Cx,2),Iter_all])
      Call mh5_put_dset(chkpnt_coor,
     &     Cx(:,:,Iter), [3,SIZE(Cx,2),1], [0,0,Iter_all-1])
*     new coordinates
      Call mh5_put_dset(chkpnt_new,Cx(1,1,Iter+1))
*     forces
      Call mh5_resize_dset(chkpnt_force, [3,SIZE(Cx,2),Iter_all])
      Call mh5_put_dset(chkpnt_force,
     &     Gx(:,:,Iter), [3,SIZE(Cx,2),1], [0,0,Iter_all-1])
*     Hessian
      If (Found) Then
        Call mh5_put_dset(chkpnt_hess,Hss_X(1))
        Call mma_deallocate(Hss_X)
      End If
#endif
      End Subroutine Chkpnt_update
*                                                                      *
************************************************************************
*                                                                      *
      Subroutine Chkpnt_update_MEP(SaveMEP,IRCRestart)
#ifdef _HDF5_
      Use mh5, Only: mh5_init_attr, mh5_open_attr, mh5_get_attr,
     &               mh5_put_attr,  mh5_close_attr, mh5_open_dset,
     &               mh5_resize_dset, mh5_put_dset, mh5_close_dset
#endif
      Logical, Intent(In) :: SaveMEP, IRCRestart
#ifdef _HDF5_
      Integer :: attrid, dsetid, iMEP

      If (IRCRestart) Then
        Call mh5_init_attr(chkpnt_id, 'IRC_RESTART', Iter_all+1)
      End If
      If (SaveMEP) Then
        attrid = mh5_open_attr(chkpnt_id, 'MEP_ITERATIONS')
        Call mh5_get_attr(attrid, iMEP)
        iMEP = iMEP+1
        Call mh5_put_attr(attrid, iMEP)
        Call mh5_close_attr(attrid)
        dsetid = mh5_open_dset(chkpnt_id, 'MEP_INDICES')
        Call mh5_resize_dset(dsetid, [iMEP])
        Call mh5_put_dset(dsetid, [Iter_all], [1], [iMEP-1])
        Call mh5_close_dset(dsetid)
      End If
#else
      Return
#ifdef _WARNING_WORKAROUND_
      If (.False.) Then
        Call Unused_logical(SaveMEP)
        Call Unused_logical(IRCRestart)
      End If
#endif
#endif
      End Subroutine Chkpnt_update_MEP
*                                                                      *
************************************************************************
*                                                                      *
      Subroutine Chkpnt_close()
#ifdef _HDF5_
      Use mh5, Only: mh5_close_file
      Call mh5_close_file(chkpnt_id)
#endif
      End Subroutine Chkpnt_close
*                                                                      *
************************************************************************
*                                                                      *
      End Module Chkpnt
