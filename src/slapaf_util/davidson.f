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
* Copyright (C) 2014, Ignacio Fdez. Galvan                             *
************************************************************************
*  Davidson
*
*> @brief
*>   Compute the lowest \p k eigenvalues of a symmetric matrix.
*> @author Ignacio Fdez. Galv&aacute;n
*>
*> @details
*> Simple application of the Davidson procedure to obtain the lowest \p k eigenvalues
*> and corresponding eigenvectors of a symmetric matrix.
*> On input, \p Vec can contain an initial guess for the eigenvectors (from a previous
*> run with smaller \p k, for example), only the non-zero vectors are used.
*>
*> @param[in]     A   Symmetric matrix, in upper triangular packed format
*> @param[in]     n   Size of the matrix
*> @param[in]     k   Number of lowest eigenvalues to compute
*> @param[out]    Eig Lowest eigenvalues
*> @param[in,out] Vec Lowest eigenvectors
*> @param[out]    iRC Return code (0 if converged)
************************************************************************
      SUBROUTINE Davidson(A,n,k,Eig,Vec,iRC)
      IMPLICIT NONE
      INTEGER n,k,iRC
      REAL*8 A(n*(n+1)/2),Eig(k),Vec(n,k)
      REAL*8, DIMENSION(:), ALLOCATABLE :: Eig_old, EVec, Sub, Ab, Proj,
     &                                     EVal
      REAL*8 Aux,Thr,Thr2,Thr3,Conv,Alpha
      real*8 ddot_
      INTEGER mk,old_mk,mink,maxk,ig,info,nTmp,iter,maxiter
      INTEGER i,j,ii,jj
      LOGICAL Last,Augmented,Reduced
      external ddot_
      PARAMETER (Thr=1.0D-7, maxiter=300, Thr2=1.0D-16, Thr3=1.0D-16)
      Real*8 rDum(1)
      Real*8, Allocatable:: Vec2(:), Val(:), Tmp(:), Diag(:), TVec(:),
     &                      TAV(:), TRes(:)
      Integer, Allocatable:: Index(:)
*
#include "stdalloc.fh"
#include "real.fh"
*#define _DEBUGPRINT_

* Diagonal preconditioned residue (Davidson)
#define DAV_DPR 1
* Inverse-iteration generalized Davidson (Olsen et al.)
#define DAV_IIGD 2
* Generalized Jacobi-Davidson (Sleijpen et al.) (DO NOT USE IT)
#define DAV_GJD 3

#ifdef _DEBUGPRINT_
      CALL TriPrt('Initial matrix','',A,n)
#endif

*---- Initialize some parameters
*      mk   = subspace size (initially k)
*      maxk = maximum subspace size (25 if k=1)
*      mink = subspace size to reduce to when the maximum is exceeded (5 if k=1)
*
      IF (k.GT.n) THEN
        CALL SysAbendMsg('Davidson','Wrong k value.','')
      END IF
      mink=MIN(MAX(k+2,5),n)
      maxk=MIN(5*mink,n)
      mk=k
      iRC=0

*---- If all the eigenvalues are wanted, better solve the system directly
*      and return
*
      IF (mk.GE.n) THEN
        CALL mma_allocate(Val,n,Label='Val')
        CALL mma_allocate(Vec2,n*n,Label='Vec2')
        CALL FZero(Vec,n*n)
        DO j=1,n
          jj=(j-1)*n
          DO i=1,j
            Vec2(jj+i)=A(j*(j-1)/2+i)
          END DO
        END DO
        call dsyev_('V','U',n,Vec2,n,Val,rDum,-1,info)
        nTmp=INT(rDum(1))

        CALL mma_allocate(Tmp,nTmp,Label='Tmp')
        call dsyev_('V','U',n,Vec2,n,Val,Tmp,nTmp,info)
        CALL SortEig(Val,Vec2,n,n,1,.false.)
        CALL mma_deallocate(Tmp)

        call dcopy_(k,Val,1,Eig,1)
        call dcopy_(n*k,Vec2,1,Vec,1)

        CALL mma_deallocate(Val)
        CALL mma_deallocate(Vec2)
#ifdef _DEBUGPRINT_
        CALL RecPrt('Eigenvalues',' ',Eig,1,n)
        WRITE(6,*)
        WRITE(6,'(A)') 'Complete system solved'
#endif
        RETURN
      END IF

*---- Allocate matrices
*      Sub  = Vectors (columns) defining the subspace (maximum maxk vectors)
*      Ab   = A*b vectors (A * Sub)
*      EVal  = Set of computed eigenvalues (maximum maxk elements)
*      EVec  = Set of computed eigenvectors, in the subspace (maximum maxk*maxk)
*
      CALL mma_allocate(Eig_old,k,label="Eig_old")
      Call mma_allocate(Sub,n*maxk,Label='Sub')
      Call mma_allocate(Ab ,n*maxk,Label='Ab ')
      Call mma_allocate(Proj,maxk*maxk,Label='Proj')
      Call mma_allocate(EVal,maxk     ,Label='EVal')
      Call mma_allocate(EVec,maxk*maxk,Label='EVec')
      AB(:)   = Zero
      EVal(:) = Zero
      EVec(:) = Zero

*---- Build an index of sorted diagonal elements in A
*
      CALL mma_allocate(Index,n,Label='Index')
      DO i=1,n
        Index(i)=i
      END DO
      DO i=1,n
        ii=Index(i)
        Aux=A(ii*(ii+1)/2)
        ii=i
        DO j=i,n
          jj=Index(j)
          IF (A(jj*(jj+1)/2) .LT. Aux) THEN
            Aux=A(jj*(jj+1)/2)
            ii=j
          END IF
        END DO
        IF (ii .NE. i) THEN
          jj=Index(ii)
          Index(ii)=Index(i)
          Index(i)=jj
        END IF
      END DO

*---- Setup the initial subspace
*      Read the non-linear-dependent columns from the initial eigenvector matrix
*      Fill up to mk with selected base vectors from the initial matrix
*       (those corresponding to the lowest diagonal elements)
*      The rest is set to zero, just in case
*
      nTmp=0
      CALL mma_allocate(Tmp,n,Label='Tmp')
      DO i=1,k
        call dcopy_(n,Vec(1,i),1,Tmp,1)
        CALL Add_Vector(n,nTmp,Sub,Tmp,Thr3)
      END DO
*
      ii=0
      Tmp(:)=Zero
      DO WHILE ((nTmp .LT. mk) .AND. (ii .LT. n))
        ii=ii+1
        jj=Index(ii)
        Tmp(jj)=One
        CALL Add_Vector(n,nTmp,Sub,Tmp,Thr3)
        Tmp(jj)=Zero
      END DO
*     ig will be a global counter to loop across all n base vectors
      ig=ii
      CALL mma_deallocate(Tmp)
      CALL FZero(Sub(1+mk*n),(maxk-mk)*n)

*---- Iterative procedure starts here
*      mk     = subspace size at each iteration
*      old_mk = value of mk at the previous iteration
*
      Augmented=.FALSE.
      Reduced=.FALSE.
      Last=.FALSE.
      old_mk=0
      iter=0
      CALL mma_allocate(Diag,n,Label='Diag')
      CALL mma_allocate(TVec,n,Label='TVec')
      CALL mma_allocate(TAV,n,Label='TAV')
      CALL mma_allocate(TRes,n,Label='TRes')
      DO WHILE (.NOT. Last)
        iter=iter+1
        IF (iter .GT. 1) call dcopy_(k,Eig,1,Eig_old,1)
#ifdef _DEBUGPRINT_
        IF (.NOT. Reduced) THEN
          WRITE(6,'(A)') '---------------'
          WRITE(6,'(A,1X,I5)') 'Iteration',iter
        END IF
        CALL RecPrt('Orthonormalized subspace',' ',Sub,n,mk)
#endif

*----   Compute the matrix product
*        Ab = A * Sub
*        Only the new vectors since the last iterations need to be calculated
*
        CALL mma_allocate(Tmp,n,Label='Tmp')
        DO i=1,n
*         Reconstruct a row (or column) of the matrix A
          DO j=1,i
            Tmp(j)=A(i*(i-1)/2+j)
          END DO
          DO j=i+1,n
            Tmp(j)=A(j*(j-1)/2+i)
          END DO
*         Compute the i-th element of each new column
          DO j=old_mk,mk-1
            Ab(j*n+i)=DDot_(n,Tmp,1,Sub(1+j*n),1)
          END DO
        END DO
        CALL mma_deallocate(Tmp)

*----   Compute the matrix to diagonalize (symmetric)
*        Proj = Sub^t * Ab
*        Again, only the new rows/columns are needed
*
        IF (old_mk .EQ. 0) THEN
           CALL DGeMM_('T','N',
     &                 mk,mk,n,
     &                 One,Sub,n,
     &                     Ab,n,
     &                 Zero,Proj,maxk)
        ELSE
          DO i=0,mk-1
            DO j=MAX(old_mk,i),mk-1
              Proj(1+i*maxk+j)=DDot_(n,Sub(1+j*n),1,
     &                                     Ab(1+i*n),1)
            END DO
          END DO
        END IF

*----   Compute the eigenvalues of the projected matrix
*        Make sure the eigenpairs are sorted
*        If the subspace has been reduced, no need to compute new eigenpairs
*
        IF (.NOT. Reduced) THEN
#ifdef _DEBUGPRINT_
          WRITE(6,'(2X,A,1X,I5)') 'Solving for subspace size:',mk
#endif
          call dcopy_(maxk*maxk,Proj,1,EVec,1)
          call dsyev_('V','L',mk,EVec,maxk,EVal,rDum,-1,info)
          nTmp=INT(rDum(1))
          CALL mma_allocate(Tmp,nTmp,Label='Tmp')
          call dsyev_('V','L',mk,EVec,maxk,EVal,Tmp,nTmp,info)
          CALL mma_deallocate(Tmp)
          CALL SortEig(EVal,EVec,mk,maxk,1,.false.)
          call dcopy_(k,EVal,1,Eig,1)
#ifdef _DEBUGPRINT_
          CALL RecPrt('Current guess',' ',Eig,1,k)
#endif
        END IF
#ifdef _DEBUGPRINT_
        CALL RecPrt('Eigenvalues',' ',EVal,1,mk)
        CALL SubRecPrt('Subspace Eigenvectors',' ',EVec,maxk,mk,mk)
        WRITE(6,*)
#endif

*----   Check for convergence
*        Converge if the change in the eigenvalues is small
*         (but if a mink size has been reached)
*        Converge if the full system has been solved
*        Stop if the number of iterations exceeds the maximum
*        Stop if no new vectors to add are found
*
        IF (iter .GT. 1) THEN
          Conv=Zero
          DO i=1,k
            IF (ABS(Eig(i)) .GT. Thr2) THEN
              Conv=MAX(Conv,ABS((Eig(i)-Eig_old(i))/Eig(i)))
            ELSE
              Conv=MAX(Conv,ABS(Eig(i)-Eig_old(i)))
            END IF
          END DO
#ifdef _DEBUGPRINT_
          IF (Augmented)
     &      WRITE(6,'(2X,A,1X,G12.6)') 'Maximum relative eigenvalue '//
     &                                 'change:',Conv
#endif
        ELSE
          Conv=Ten*Thr
        END IF
        old_mk=mk
        IF (Augmented .AND. (Conv .LE. Thr) .AND. (mk .GE. mink)) THEN
#ifdef _DEBUGPRINT_
          WRITE(6,'(A)') 'Converged due to small change'
#endif
          Last=.TRUE.
        ELSE IF (mk .EQ. n) THEN
#ifdef _DEBUGPRINT_
          WRITE(6,'(A)') 'Complete system solved'
#endif
          Last=.TRUE.
        ELSE IF (iter .GE. maxiter) THEN
#ifdef _DEBUGPRINT_
          WRITE(6,'(A)') 'Not converged'
#endif
          Last=.TRUE.
          iRC=1

*----   Reduce the subspace size if it exceeds the maximum (maxk)
*        Sub' = Sub * Vec(1:mink)
*        Sub' should be orthonormal if Sub is orthonormal
*        (A reduction does not consume an iteration)
*       There is also a reduction if the process is stagnated
*
        ELSE IF ((MIN(mk+k,n) .GT. maxk) .OR. (iRC .EQ. 2)) THEN
          IF (iRC .EQ. 2) iRC=0
#ifdef _DEBUGPRINT_
          WRITE(6,'(2X,A,1X,I5)') 'Reducing search space to',mink
#endif
          CALL mma_allocate(Tmp,mink*n,Label='Tmp')
          CALL DGeMM_('N','N',n,mink,mk,One,Sub,n,
     &                                      EVec,maxk,
     &                                  Zero,Tmp,n)
          call dcopy_(mink*n,Tmp,1,Sub,1)
          CALL mma_deallocate(Tmp)

*----     To make sure Sub' is orthonormal, add the vectors one by one
*
          j=0
          i=0
          DO WHILE ((j .LT. mink) .AND. (i .LT. mk))
            i=i+1
            Call Add_Vector(n,j,Sub,Sub(1+(i-1)*n),Thr3)
          END DO

*----     j should be mink, but who knows...
*
#ifdef _DEBUGPRINT_
          IF (j .LT. mink) THEN
            WRITE(6,'(2X,A,1X,I5)') 'Fewer vectors found:',j
          END IF
#endif
          CALL FZero(Sub(1+j*n),(maxk-j)*n)
          CALL FZero(Ab(1+j*n),(maxk-j)*n)
          CALL FZero(EVec,maxk*maxk)
          DO i=0,j-1
            EVec(1+i*(maxk+1))=One
          END DO
          mk=j
          old_mk=0
          Augmented=.FALSE.
          Reduced=.TRUE.
          iter=iter-1

*----   Expand the subspace
*        For each eigenpair i of the first k,
*        check convergence for the residuals r:
*         r = Ab * Vec(i) - Val(i) * Sub * Vec(i)
*        Add a new vector, orthonormalized with the previous vectors,
*        computed from r and the eigenpair
*        (different possible variants)
*
        ELSE
          CALL mma_allocate(Tmp,n,Label='Tmp')
          Conv=Zero
          jj=0
          DO i=0,k-1
*           Vector in full space: Sub*Vec(i)
            Call dGeMV_('N',n,mk,One,Sub,n,
     &                              EVec(1+i*maxk),1,
     &                          Zero,TVec,1)
*           Product of matrix and vector: Ab*Vec(i)
            Call dGeMV_('N',n,mk,One,Ab,n,
     &                              EVec(1+i*maxk),1,
     &                          Zero,TAV,1)
*           Residual: (A-Val(i))*Vec(i) = Ab*Vec(i) - Val(i)*Sub*Vec(i)
            call dcopy_(n,TAV,1,TRes,1)
            call daxpy_(n,-EVal(1+i),TVec,1,TRes,1)
            Conv=MAX(Conv,DDot_(n,TRes,1,TRes,1))

*----       Scale vector, orthonormalize, and add to subspace
*
#define DAV_METH DAV_IIGD
#if DAV_METH == DAV_DPR
*           Diagonal matrix to scale the vectors: 1/(A(j,j)-Val(i))
            DO j=1,n
              Aux=A(j*(j+1)/2)-Eval(1+i)
              Diag(j)=One/SIGN(MAX(ABS(Aux),Thr2),Aux)
            END DO
*           scale
            DO j=1,n
              Tmp(j)=TRes(j)*Diag(j)
            END DO
#elif DAV_METH == DAV_IIGD
*           Diagonal matrix to scale the vectors: 1/(A(j,j)-Val(i))
            DO j=1,n
              Aux=A(j*(j+1)/2)-EVal(1+i)
              Diag(j)=One/SIGN(MAX(ABS(Aux),Thr2),Aux)
            END DO
*           scale
            DO j=1,n
              Tmp(j)=TRes(j)*Diag(j)
            END DO
            Alpha=Zero
            DO j=1,n
              Alpha=Alpha+Diag(j)*TVec(j)**2
            END DO
            Alpha=DDot_(n,TVec,1,Tmp,1)/Alpha
*           subtract
            DO j=1,n
              TVec(j)=TVec(j)*Diag(j)
            END DO
            call daxpy_(n,-Alpha,TVec,1,Tmp,1)
#elif DAV_METH == DAV_GJD
*   DO NOT USE THIS VARIANT!
* This is not practical as it stands, the equation should be
* solved only approximately, and it is not efficient to do this
* for each of the possibly many eigenpairs
            Block
            Real*8, Allocatable:: P1(:)
            CALL mma_allocate(P1,n*n,Label='P1')
*           project: (I-|v><v|) (A-e*I) (I-|v><v|) =
*                    A - e*I + (<v|P>+e)*|v><v| - |v><P| - |P><v|
*           e = Val(i); |P> = A|v>
            Aux=DDot_(n,TVec,1,TAV,1)
            DO kk=0,n-1
              ll=kk*(kk+1)/2
              DO ii=0,kk
                P1(1+ii*n+kk)=A(ll+ii+1)+
     &            (Aux+EVal(1+i))*TVec(1+ii)*TVec(1+kk)-
     &            TAV(1+ii)*TVec(1+kk)-
     &            TAV(1+kk)*TVec(1+ii)
              END DO
              P1(1+kk*n+kk)=P1(1+kk*n+kk)-EVal(1+i)
            END DO
            Block
              Integer iDum(1)
               iDum(1)=0
*              solve the equation
               CALL CG_Solver(n,n*n,P1,iDum,TRes,Tmp,info,5)
            End Block
#ifdef _DEBUGPRINT_
            WRITE(6,*) 'CG iterations',info
#endif
            CALL mma_deallocate(P1)
          End Block
#endif
            IF (mk+jj .LE. n-1) THEN
              jj=mk+jj
              CALL Add_Vector(n,jj,Sub,Tmp,Thr3)
              jj=jj-mk
            END IF
          END DO
#ifdef _DEBUGPRINT_
          WRITE(6,'(2X,A,1X,G12.6)') 'Maximum residual:',Conv
#endif
          IF ((Conv .LT. Thr3) .AND. (mk .GE. mink)) THEN
#ifdef _DEBUGPRINT_
            WRITE(6,'(A)') 'Converged due to small residual'
#endif
            Last=.TRUE.
          ELSE
            mk=MIN(mk+jj,n)

*----       If no new vector is found to add to the subspace, we are in trouble
*            -Try to find a non-linear-dependent base vector in the original matrix
*
            IF (jj .EQ. 0) THEN
#ifdef _DEBUGPRINT_
              WRITE(6,'(A)') 'Process stagnated'
#endif
              IF (mk .LT. maxk) THEN
                CALL FZero(Tmp,n)
                i=0
                DO WHILE ((jj .LT. 1) .AND. (i .LT. n))
                  i=i+1
                  ig=MOD(ig,n)+1
                  ii=Index(ig)
                  Tmp(ii)=One
                  jj=mk+jj
                  CALL Add_Vector(n,jj,Sub,Tmp,Thr3)
                  Tmp(ii)=Zero
                  jj=jj-mk
                END DO
                mk=MIN(mk+jj,n)
                IF (jj .GT. 0) Augmented=.TRUE.
              END IF
              IF (jj .EQ. 0) THEN
                Augmented=.FALSE.
                iRC=2
              ENDIF
            ELSE
              Augmented=.TRUE.
            END IF
          END IF
          CALL mma_deallocate(Tmp)
          Reduced=.FALSE.
        END IF
      END DO
      CALL mma_deallocate(Diag)
      CALL mma_deallocate(TVec)
      CALL mma_deallocate(TAV)
      CALL mma_deallocate(TRes)
      CALL mma_deallocate(Index)

*---- Store the current lowest k eigenvectors (in the full space)
*      Vec' = Sub * Vec(1:k)
*
      CALL DGeMM_('N','N',n,k,mk,One,Sub,n,EVec,maxk,
     &                           Zero,Vec,n)

      Call mma_deallocate(Sub)
      Call mma_deallocate(Ab)
      Call mma_deallocate(Proj)
      Call mma_deallocate(EVal)
      Call mma_deallocate(EVec)
      CALL mma_deallocate(Eig_old)

      END

************************************************************************
*  Add_Vector
*
*> @brief
*>   Add (or not) a vector to an orthonormal set
*> @author Ignacio Fdez. Galv&aacute;n
*>
*> @details
*> Adds a given vector to an existing set of orthonormal vectors.
*> The vector is orthogonalized against the existing set and, if the remainder
*> is large enough, it will be normalized and added to the set.
*> The input matrix with the initial set must have size at least \p n (\p m+1).
*> On output, \p m will be increased by 1 if the vector was added, and unchanged
*> otherwise.
*>
*> @param[in]     n   Size of the vectors
*> @param[in,out] m   Number of vectors in the subspace
*> @param[in,out] Sub Subspace of vectors
*> @param[in,out] Vec Vector to be added
*> @param[in]     Thr Threshold for linear dependence check
************************************************************************
      SUBROUTINE Add_Vector(n,m,Sub,Vec,Thr)
      IMPLICIT NONE
      INTEGER n,m,i
      REAL*8 Sub(n,m+1),Vec(n),Thr,Aux,DDot_
      external ddot_
#include "real.fh"

      DO i=1,m
        Aux=DDot_(n,Sub(1,i),1,Vec,1)
        call daxpy_(n,-Aux,Sub(1,i),1,Vec,1)
      END DO
      Aux=DDot_(n,Vec,1,Vec,1)
      IF (ABS(Aux) .GT. Thr) THEN
*----   Safety net: orthonormalize again before adding it
        call dscal_(n,One/SQRT(Aux),Vec,1)
        DO i=1,m
          Aux=DDot_(n,Sub(1,i),1,Vec,1)
          call daxpy_(n,-Aux,Sub(1,i),1,Vec,1)
        END DO
        m=m+1
        Aux=DDot_(n,Vec,1,Vec,1)
        call dscal_(n,One/SQRT(Aux),Vec,1)
        call dcopy_(n,Vec,1,Sub(1,m),1)
      END IF

      END
