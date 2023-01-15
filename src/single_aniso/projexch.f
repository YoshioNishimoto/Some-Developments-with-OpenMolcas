      Subroutine project_exchange_single (n1,n2,aniso_1,aniso_2,
     &                                    E,S,M,iPrint)
      Implicit None
      Integer, Parameter :: wp=selected_real_kind(p=15,r=307)
#include "stdalloc.fh"
      Integer, intent(in)          :: iPrint               ! print level
      Integer, intent(in)          :: n1, n2               ! size of the local pseudospins
      Real(kind=wp), intent(in)    :: E( n1*n2 )           ! spin-orbit energy
      Complex(kind=wp), intent(in) :: S( 3, n1*n2, n1*n2 ) ! spin momentum, total
      Complex(kind=wp), intent(in) :: M( 3, n1*n2, n1*n2 ) ! magnetic momentum, total
      Character(LEN=180), intent(in)   :: aniso_1              ! aniso file containing information about site 1
      Character(LEN=180), intent(in)   :: aniso_2              ! aniso file containing information about site 2
!     local variables
      Integer                       :: nExch, l, i, j, k, nLoc
      Integer                       :: neq(2), nexchange(2), lmax, nneq
      Integer                       :: nb, nb1, isite
      Integer                       :: intc(2), nind(2,2)

      Integer, allocatable          :: ibas(:,:)
      Real(kind=wp), allocatable    :: g(:), mg(:,:)
      Complex(kind=wp), allocatable :: S1(:,:,:), M1(:,:,:)
      Complex(kind=wp), allocatable :: S2(:,:,:), M2(:,:,:)
      Complex(kind=wp), allocatable :: Z(:,:),tmp(:,:)

      Integer                       :: i1,j1,i2,itmp,is1,js1,js2
      Integer                       :: is2,nb2,icoord(2),lp
      Integer                       :: k1,k2,q1,q2
      Character(len=1)              :: itype(2)
      Logical                       :: aniso_1_exists, aniso_2_exists
      Real(kind=wp), allocatable    :: g1(:),g2(:),mg1(:,:),mg2(:,:)
      Real(kind=wp), allocatable    :: eso_1(:), eso_2(:)
      Real(kind=wp)                 :: gtens_iso(3), D, EoverD, f1
      Real(kind=wp)                 :: rot(2,1,3,3), diff(n1*n2,n1*n2)
      Complex(kind=wp), allocatable :: mom1(:,:,:), mom2(:,:,:)
      Complex(kind=wp), allocatable :: som1(:,:,:), som2(:,:,:)
      Complex(kind=wp), allocatable :: lom1(:,:,:), lom2(:,:,:)
      Complex(kind=wp), allocatable :: mom1r(:,:,:), mom2r(:,:,:)
      Complex(kind=wp), allocatable :: som1r(:,:,:), som2r(:,:,:)
      Complex(kind=wp), allocatable :: Z1(:,:), Z2(:,:)

      Integer                       :: diff_min(n1*n2)
c      Real(kind=wp), allocatable    :: wout(:)
c      Complex(kind=wp), allocatable :: zout(:,:)

      ! exchange Hamiltonian
      Complex(kind=wp), allocatable :: M3(:,:,:),S3(:,:,:),Z12(:,:)
      Complex(kind=wp), allocatable :: MM(:,:,:,:),SM(:,:,:,:),ZM(:,:)
      ! exchange Hamiltonian:
      Real(kind=wp), allocatable    :: E1(:), rot1(:,:), rot2(:,:)
      Complex(kind=wp), allocatable :: H(:,:), H1(:,:,:,:), H2(:,:,:,:),
     &                                 H3(:,:,:,:)
      ! zero field splitting extracted from the general exchange matrix
      Complex(kind=wp), allocatable :: ZFS1(:,:), ZFS2(:,:)

      ! exchange parameters:
      Complex(kind=wp) :: JLin3(1:(n1-1),-(n1-1):(n1-1),
     &                          1:(n2-1),-(n2-1):(n2-1))
      Real(kind=wp) :: JLinC3(3,3)
      Logical :: dbg
      Integer :: norder
      External :: norder


      nExch=n1*n2 ! size of exchange matrix
      dbg=.true.
      ! sanity check:
      If ( nExch<=1 ) Then
          Call WarningMessage(1,'Project_exchange_single:: '//
     &                          'nExch is too small')
          Write(6,'(A,3(I5,1x))') 'n1, n2, nExch =',n1, n2, nExch
          Write(6,'(A         )') 'Calculation will continue, but'//
     &           'no projection of exchange interaction is done.'
          Return
      End If
      If(dbg) Write(6,*) 'PREX:  aniso_1=', trim(aniso_1)
      If(dbg) Write(6,*) 'PREX:  aniso_2=', trim(aniso_2)
      If(dbg) Write(6,*) 'PREX:  E() =', E(1:nExch)
!-----------------------------------------------------------------------
! allocate local variables:
      Call mma_allocate(g,3,'g main values')
      Call mma_allocate(mg,3,3,'g main axes')
      Call mma_allocate(s1,3,nExch,nExch,'rotated S')
      Call mma_allocate(m1,3,nExch,nExch,'rotated M')
      Call mma_allocate(s2,3,nExch,nExch,'pseudo S')
      Call mma_allocate(m2,3,nExch,nExch,'pseudo M')
      Call mma_allocate(Z,nExch,nExch,'rotated M')
      Call mma_allocate(ibas,nExch,2,'local basis')

      Call mma_allocate(g1,3,'local g main values, site 1')
      Call mma_allocate(g2,3,'local g main values, site 2')
      Call mma_allocate(mg1,3,3,'local g main axes, site 1')
      Call mma_allocate(mg2,3,3,'local g main axes, site 2')
!-----------------------------------------------------------------------
! 1. find total magnetic axis of the entire manifold
      Call dcopy_(3  ,[0.0_wp],0,g,1)
      Call dcopy_(3*3,[0.0_wp],0,mg,1)

      Call atens(M, n1*n2, g, mg, 2)

! 2. rotate the S and M to the main magnetic axes
      Call zcopy_(3*nExch*nExch,[(0.0_wp,0.0_wp)],0,s1,1)
      Call zcopy_(3*nExch*nExch,[(0.0_wp,0.0_wp)],0,m1,1)

      Call rotmom2(  S(1:3,1:nExch,1:nExch), nExch, mg,
     &              S1(1:3,1:nExch,1:nExch) )
      Call rotmom2(  M(1:3,1:nExch,1:nExch), nExch, mg,
     &              M1(1:3,1:nExch,1:nExch) )

!      Call atens(M1, n1*n2, g, mg, 2)

! 3. Find total exchange pseudospin:
      Call zcopy_(nExch*nExch,[(0.0_wp,0.0_wp)],0,Z,1)
      Call pseudospin(M1,nExch,Z,3,1,iPrint)

      Call mma_allocate(E1,n1*n2,'diagonal energies')
      Call dcopy_(n1*n2,[0.0_wp],0,E1,1)
      CALL rtrace(n1*n2,E,E1)
      ! re-write the CF matrix in J-pseudospin basis:
      ! energy units =  cm-1
      Call mma_allocate(H,n1*n2,n1*n2,'Hamiltonian')
      Call zcopy_(n1*n2*n1*n2,[(0.0_wp,0.0_wp)],0,H,1)
      Do i=1,n1*n2
        Do j=1,n1*n2
          Do k=1,n1*n2
            H(i,j)=H(i,j) + E1(k)*CONJG(Z(k,i))*Z(k,j)
          End Do
        End Do
      End Do

      Write(6,'(A)') 'Exchange Hamiltonian in pseudospin basis'
      Do i=1,n1*n2
        Write(6,'(12(2F7.4,1x))') (H(i,j),j=1,n1*n2)
      End Do




      Call mma_allocate(tmp,nExch,nExch,'scratch')
      Do L=1,3
         ! magnetic moment
         Call zcopy_(nExch*nExch,[(0.0_wp,0.0_wp)],0,TMP,1)
         Call zgemm_('C','N',nExch,nExch,nExch,
     &              (1.0_wp,0.0_wp),Z, nExch,
     &                              M1(L,:,:), nExch,
     &              (0.0_wp,0.0_wp),TMP, nExch )
         Call zcopy_(nExch*nExch,[(0.0_wp,0.0_wp)],0,M2(L,:,:),1)
         Call zgemm_('N','N',nExch,nExch,nExch,
     &              (1.0_wp,0.0_wp),TMP,nExch,
     &                                Z,nExch,
     &              (0.0_wp,0.0_wp), M2(L,:,:), nExch )

         ! spin moment
         Call zcopy_(nExch*nExch,[(0.0_wp,0.0_wp)],0,TMP,1)
         Call zgemm_('C','N',nExch,nExch,nExch,
     &              (1.0_wp,0.0_wp),Z,nExch,
     &                              S1(L,:,:), nExch,
     &              (0.0_wp,0.0_wp),TMP,nExch )
         Call zcopy_(nExch*nExch,[(0.0_wp,0.0_wp)],0,S2(L,:,:),1)
         Call zgemm_('N','N',nExch,nExch,nExch,
     &              (1.0_wp,0.0_wp),TMP,nExch,
     &                                Z,nExch,
     &              (0.0_wp,0.0_wp),S2(L,:,:),nExch )

      End Do ! L
      Call mma_deallocate(tmp)

      If (dbg) Then
         Write(6,*) 'Matrix M2'
         Do l=1,3
            Write(6,*)
            Write(6,'(A,i2)') 'Projection = ', l
            Write(6,*)
            Do i=1,nExch
c               Write(6,'(12(2F7.4,1x))') (M2(l,i,j),j=1,nExch)
            End Do
         End Do

c         Write(6,*) 'Matrix S2'
c         Do l=1,3
c            Write(6,*)
c            Write(6,'(A,i2)') 'Projection = ', l
c            Write(6,*)
c            Do i=1,nExch
c               Write(6,'(12(2F7.4,1x))') (S2(l,i,j),j=1,nExch)
c            End Do
c         End Do
      End If

      ! diagonal matrix M2(z,i,i) will be compared with the matrix we obtain
      ! from the coupling of two local momenta matrices, following their pseudospins
      ! We will attempt to assign states by comparing the diagonal elements.
      ! re-suffling may be neccecary;


! 4. Set up an indexing scheme for the states:
      nneq=2
      neq(1)=1
      neq(2)=1
      nexchange(1)=n1
      nexchange(2)=n2
      lmax=2
      nLoc=max(n1,n2)

      nind(:,:)=0
      intc(:)=0
      ibas(:,:)=0

      l=0
      Do i=1,nneq
        Do j=1,neq(i)
          l=l+1
          nind(l,1)=i
          nind(l,2)=j
        End Do
      End Do
      intc(1)=1
      If (lmax.gt.1) Then
        Do i=2,lmax
          isite=nind(i-1, 1)
          intc(i)=intc(i-1)*nexchange(isite)
        End Do
      End If
      Do nb=1,n1*n2
        nb1=nb-1
        Do i=1,lmax
          ibas(nb, lmax-i+1)= nb1 / intc(lmax-i+1)
          nb1=nb1-ibas(nb,lmax-i+1)*intc(lmax-i+1)
        End Do
      End Do

      If(dbg) Then
         Write(6,'(34x,A,1x,20i3)')  'site Nr.', (i,i=1,lmax)
         Do nb=1,n1*n2
           Write(6,'(A,i5,A,20i3)') 'COUPLING: basis set:  ibas(',nb,
     &                             ' ,isite) = ',(ibas(nb,i)+1,i=1,lmax)
         End Do
      End If ! dbg

! 5. Determine LOCAL matrices for the spin and orbital momenta
      Call mma_allocate(mom1,3,n1,n1,'local magnetic moment site 1')
      Call mma_allocate(mom1r,3,n1,n1,'local magnetic moment site 1')
      Call mma_allocate(som1r,3,n1,n1,'local magnetic moment site 1')
      Call mma_allocate(som1,3,n1,n1,'local spin moment site 1')
      Call mma_allocate(lom1,3,n1,n1,'local orbital moment site 1')
      Call mma_allocate(mom2,3,n2,n2,'local magnetic moment site 2')
      Call mma_allocate(mom2r,3,n2,n2,'local magnetic moment site 2')
      Call mma_allocate(som2r,3,n2,n2,'local magnetic moment site 2')
      Call mma_allocate(som2,3,n2,n2,'local spin moment site 2')
      Call mma_allocate(lom2,3,n2,n2,'local orbital moment site 2')
      Call mma_allocate(eso_1,n1,'local s-o states site 1')
      Call mma_allocate(eso_2,n2,'local s-o states site 2')
      Call mma_allocate(Z1,n1,n1,'local pseudospin site 1')
      Call mma_allocate(Z2,n2,n2,'local pseudospin site 2')

      Call zcopy_(3*n1*n1,[(0.0_wp,0.0_wp)],0,mom1,1)
      Call zcopy_(3*n1*n1,[(0.0_wp,0.0_wp)],0,mom1r,1)
      Call zcopy_(3*n1*n1,[(0.0_wp,0.0_wp)],0,som1r,1)
      Call zcopy_(3*n1*n1,[(0.0_wp,0.0_wp)],0,som1,1)
      Call zcopy_(3*n1*n1,[(0.0_wp,0.0_wp)],0,lom1,1)
      Call zcopy_(3*n2*n2,[(0.0_wp,0.0_wp)],0,mom2,1)
      Call zcopy_(3*n2*n2,[(0.0_wp,0.0_wp)],0,mom2r,1)
      Call zcopy_(3*n2*n2,[(0.0_wp,0.0_wp)],0,som2r,1)
      Call zcopy_(3*n2*n2,[(0.0_wp,0.0_wp)],0,som2,1)
      Call zcopy_(3*n2*n2,[(0.0_wp,0.0_wp)],0,lom2,1)
      Call dcopy_(n1,[0.0_wp],0,eso_1,1)
      Call dcopy_(n2,[0.0_wp],0,eso_2,1)

      ! Check if aniso files exist:
      aniso_1_exists=.false.
      aniso_2_exists=.false.
      Call f_inquire(aniso_1,aniso_1_exists)
      Call f_inquire(aniso_2,aniso_2_exists)
      ! fetch local data
      If (aniso_1_exists) Then
         itype(1)='A'
         Call read_aniso_old_exch( aniso_1,n1,eso_1,mom1,som1,lom1 )
      Else
         Write(6,'(A)') 'File '//trim(aniso_1)//'does not exist.'
         Write(6,'(A)') 'Generate the local data considerring '//
     &                  'the site 1 isotropic'
         gtens_iso(1)=2.0_wp
         gtens_iso(2)=2.0_wp
         gtens_iso(3)=2.0_wp
         D=0.0_wp
         EoverD=0.0_wp
         itype(1)='B'
         Call generate_isotrop_site( n1, itmp, n1, n1,
     &                             gtens_iso, D, EoverD,
     &                             eso_1,mom1,som1,lom1 )
c         Call rotmom2( som1, n2, mg, som1r )
c         Call rotmom2( lom1, n2, mg, lom1r )
c         Call rotmom2( mom1, n2, mg, mom1r )
      End If

      If (aniso_2_exists) Then
         itype(2)='A'
         Call read_aniso_old_exch( aniso_2,n2,eso_2,mom2,som2,lom2 )
      Else
         Write(6,'(A)') 'File '//trim(aniso_2)//' does not exist.'
         Write(6,'(A)') 'Generate the local data considerring '//
     &                  'the site 2 isotropic'
         gtens_iso(1)=2.0_wp
         gtens_iso(2)=2.0_wp
         gtens_iso(3)=2.0_wp
         D=0.0_wp
         EoverD=0.0_wp
         itype(2)='B'
         Call generate_isotrop_site( n2, itmp, n2, n2,
     &                             gtens_iso, D, EoverD,
     &                             eso_2,mom2,som2,lom2 )
         ! rotate to the magnetic axes of exchange set:
c         Call rotmom2( som2, n2, mg, som2r )
c         Call rotmom2( lom2, n2, mg, lom2r )
c         Call rotmom2( mom2, n2, mg, mom2r )

      End If
c---------------------------------------------------------------
      Call zcopy_(n1*n1,[(0.0_wp,0.0_wp)],0,Z1,1)
      Call zcopy_(n2*n2,[(0.0_wp,0.0_wp)],0,Z2,1)
      If(aniso_1_exists) Then
         g1=0.0_wp
         mg1=0.0_wp
         mom1r=(0.0_wp,0.0_wp)
         Call atens(mom1, n1, g1, mg1, 2)
         Call rotmom2(mom1, n1, mg1, mom1r )
         Call rotmom2(som1, n1, mg1, som1r )
         Call pseudospin(mom1r,n1,Z1,3,1,iPrint)
      Else
         Call zcopy_(3*n1*n1,mom1,1,mom1r,1)
         Call zcopy_(3*n1*n1,som1,1,som1r,1)
         Do i=1,n1
            Z1(i,i)=(1.0_wp,0.0_wp)
         End Do
      End If
      If(aniso_2_exists) Then
         g2=0.0_wp
         mg2=0.0_wp
         mom2r=(0.0_wp,0.0_wp)
         Call atens(mom2, n2, g2, mg2, 2)
         Call rotmom2(mom2, n2, mg2, mom2r )
         Call rotmom2(som2, n2, mg2, som2r )
         Call pseudospin(mom2r,n2,Z2,3,1,iPrint)
      Else
         Call zcopy_(3*n2*n2,mom2,1,mom2r,1)
         Call zcopy_(3*n2*n2,som2,1,som2r,1)
         Do i=1,n2
            Z2(i,i)=(1.0_wp,0.0_wp)
         End Do
      End If

      ! transform local momenta to their pseudospin basis
      Call mma_allocate(tmp,n1,n1,'scratch')
      Do L=1,3
         ! magnetic moment
         Call zcopy_(n1*n1,[(0.0_wp,0.0_wp)],0,TMP,1)
         Call zgemm_('C','N',n1,n1,n1,
     &              (1.0_wp,0.0_wp),Z1, n1,
     &                              mom1r(L,:,:), n1,
     &              (0.0_wp,0.0_wp),TMP, n1 )
         Call zcopy_(n1*n1,[(0.0_wp,0.0_wp)],0,mom1r(L,:,:),1)
         Call zgemm_('N','N',n1,n1,n1,
     &              (1.0_wp,0.0_wp),TMP,n1,
     &                                Z1,n1,
     &              (0.0_wp,0.0_wp), mom1r(L,:,:), n1 )

         ! spin moment
         Call zcopy_(n1*n1,[(0.0_wp,0.0_wp)],0,TMP,1)
         Call zgemm_('C','N',n1,n1,n1,
     &              (1.0_wp,0.0_wp),Z1, n1,
     &                              som1r(L,:,:), n1,
     &              (0.0_wp,0.0_wp),TMP, n1 )
         Call zcopy_(n1*n1,[(0.0_wp,0.0_wp)],0,som1r(L,:,:),1)
         Call zgemm_('N','N',n1,n1,n1,
     &              (1.0_wp,0.0_wp),TMP,n1,
     &                                Z1,n1,
     &              (0.0_wp,0.0_wp), som1r(L,:,:), n1 )

      End Do ! L
      Call mma_deallocate(tmp)


      ! transform local momenta to their pseudospin basis
      Call mma_allocate(tmp,n2,n2,'scratch')
      Do L=1,3
         ! magnetic moment
         Call zcopy_(n2*n2,[(0.0_wp,0.0_wp)],0,TMP,1)
         Call zgemm_('C','N',n2,n2,n2,
     &              (1.0_wp,0.0_wp),Z2, n2,
     &                              mom2r(L,:,:), n2,
     &              (0.0_wp,0.0_wp),TMP, n2 )
         Call zcopy_(n2*n2,[(0.0_wp,0.0_wp)],0,mom2r(L,:,:),1)
         Call zgemm_('N','N',n2,n2,n2,
     &              (1.0_wp,0.0_wp),TMP,n2,
     &                                Z2,n2,
     &              (0.0_wp,0.0_wp), mom2r(L,:,:), n2 )

         ! spin moment
         Call zcopy_(n2*n2,[(0.0_wp,0.0_wp)],0,TMP,1)
         Call zgemm_('C','N',n2,n2,n2,
     &              (1.0_wp,0.0_wp),Z2, n2,
     &                              som2r(L,:,:), n2,
     &              (0.0_wp,0.0_wp),TMP, n2 )
         Call zcopy_(n2*n2,[(0.0_wp,0.0_wp)],0,som2r(L,:,:),1)
         Call zgemm_('N','N',n2,n2,n2,
     &              (1.0_wp,0.0_wp),TMP,n2,
     &                                Z2,n2,
     &              (0.0_wp,0.0_wp), som2r(L,:,:), n2 )

      End Do ! L
      Call mma_deallocate(tmp)



      If (dbg) Then
         Write(6,*) 'Matrix mom1 site 1'
         Do l=1,3
            Write(6,*)
            Write(6,'(A,i2)') 'Projection = ', l
            Write(6,*)
            Do i=1,n1
               Write(6,'(12(2F7.4,1x))') (mom1r(l,i,j),j=1,n1)
            End Do
         End Do

         Write(6,*) 'Matrix mom2 site 2'
         Do l=1,3
            Write(6,*)
            Write(6,'(A,i2)') 'Projection = ', l
            Write(6,*)
            Do i=1,n2
               Write(6,'(12(2F7.4,1x))') (mom2r(l,i,j),j=1,n2)
            End Do
         End Do
      End If ! dbg
c---------------------------------------------------------------
!  build a test matrix of magnetic moment as direct product of local matrices
!  diagonalize and test the eigenvectors:

      Call mma_allocate(Z12,n1*n2,n1*n2,'Z12')
      Call mma_allocate(M3,3,n1*n2,n1*n2,'M3')
      Call mma_allocate(S3,3,n1*n2,n1*n2,'S3')
      Call mma_allocate(MM,2,3,nLoc,nLoc,'MM')
      Call mma_allocate(SM,2,3,nLoc,nLoc,'SM')
      Call zcopy_(  n1*n1*n2*n2,[(0.0_wp,0.0_wp)],0,Z12,1)
      Call zcopy_(3*n1*n1*n2*n2,[(0.0_wp,0.0_wp)],0,M3 ,1)
      Call zcopy_(3*n1*n1*n2*n2,[(0.0_wp,0.0_wp)],0,S3 ,1)
      Call zcopy_(2*3*nLoc*nLoc,[(0.0_wp,0.0_wp)],0,MM ,1)
      Call zcopy_(2*3*nLoc*nLoc,[(0.0_wp,0.0_wp)],0,SM ,1)
      rot(:,:,:,:)=0.0_wp
      rot(1,1,1,1)=1.0_wp
      rot(1,1,2,2)=1.0_wp
      rot(1,1,3,3)=1.0_wp
      rot(2,1,1,1)=1.0_wp
      rot(2,1,2,2)=1.0_wp
      rot(2,1,3,3)=1.0_wp

      Do i=1,n1
        Do j=1,n1
          Do l=1,3
            MM(1,l,i,j)=mom1r(l,i,j)
            SM(1,l,i,j)=som1r(l,i,j)
          End Do
        End Do
      End Do
      Do i=1,n2
        Do j=1,n2
          Do l=1,3
            MM(2,l,i,j)=mom2r(l,i,j)
            SM(2,l,i,j)=som2r(l,i,j)
          End Do
        End Do
      End Do

      Do L=1,3
        Do isite=1,lmax
          Do nb1=1,n1*n2
            Do lp=1,lmax
              icoord(lp)=ibas(nb1,lp)
            End Do
              i1=nind(isite,1)
              j1=nind(isite,2)
            is1=ibas(nb1,isite)+1

            Do js1=1,nexchange(i1)
              icoord(isite)=js1-1
              nb2=norder(icoord,intc,lmax)
              M3( l, nb1, nb2 ) = M3( l, nb1, nb2 ) +
     &              rot( i1, j1, l, 1 ) * MM( i1, 1, is1, js1 )
     &             +rot( i1, j1, l, 2 ) * MM( i1, 2, is1, js1 )
     &             +rot( i1, j1, l, 3 ) * MM( i1, 3, is1, js1 )
              S3( l, nb1, nb2 ) = S3( l, nb1, nb2 ) +
     &              rot( i1, j1, l, 1 ) * SM( i1, 1, is1, js1 )
     &             +rot( i1, j1, l, 2 ) * SM( i1, 2, is1, js1 )
     &             +rot( i1, j1, l, 3 ) * SM( i1, 3, is1, js1 )

            End Do  ! js1
          End Do  ! nb1
        End Do  ! isite
      End Do  ! L


      If (dbg) Then
         Write(6,*) 'Matrix M3'
         Do l=1,3
            Write(6,*)
            Write(6,'(A,i2)') 'Projection = ', l
            Write(6,*)
            Do i=1,nExch
               Write(6,'(12(2F7.4,1x))') (M3(l,i,j),j=1,nExch)
            End Do
         End Do

         Write(6,*) 'Matrix S3'
         Do l=1,3
            Write(6,*)
            Write(6,'(A,i2)') 'Projection = ', l
            Write(6,*)
            Do i=1,nExch
               Write(6,'(12(2F7.4,1x))') (S3(l,i,j),j=1,nExch)
            End Do
         End Do
      End If


      Write(6,*) 'Matrix M2 Z'
      Do i=1,nExch
      Write(6,'(12(2F7.4,1x))') (M2(3,i,j),j=1,nExch)
      End Do

      Write(6,*) 'Matrix M3 Z'
      Do i=1,nExch
      Write(6,'(12(2F7.4,1x))') (M3(3,i,j),j=1,nExch)
      End Do

      ! find the roatation matrix which makes M2 as close as possible to M3

      Call mma_allocate(ZM,nExch,nExch,'ZM')
      Call zcopy_(nExch*nExch,[(0.0_wp,0.0_wp)],0,ZM,1)

      diff=0.0_wp
      diff_min=0
      Do i=1,nExch
         Do j=1,nExch
             diff(i,j) = abs(M3(3,i,i) - M2(3,j,j))
         End Do
         diff_min(i)=MINLOC( diff(i,:),1)
      End Do

      Do i=1,nExch
        Write(6,'(A,i2,A,12(F8.4,1x),a,I3)') 'i= ',i,' j-> ',
     &    (diff(i,j),j=1,nExch), 'diff_min(i)=', diff_min(i)
      End Do

      Do i=1,nExch
        Do j=1,nExch
          if (diff_min(j)==i) ZM(i,j)=(1.0_wp,0.0_wp)
        End Do
      End Do

      Write(6,*) 'Matrix ZM'
      Do i=1,nExch
         Write(6,'(12(2F4.1,1x))') (ZM(i,j),j=1,nExch)
      End Do

      ! rotate magnetic and spin moments:

      Call mma_allocate(tmp,nExch,nExch,'scratch')
      Do L=1,3
         ! magnetic moment
         Call zcopy_(nExch*nExch,[(0.0_wp,0.0_wp)],0,TMP,1)
         Call zgemm_('C','N',nExch,nExch,nExch,
     &              (1.0_wp,0.0_wp),ZM, nExch,
     &                              M2(L,:,:), nExch,
     &              (0.0_wp,0.0_wp),TMP, nExch )
         Call zcopy_(nExch*nExch,[(0.0_wp,0.0_wp)],0,M2(L,:,:),1)
         Call zgemm_('N','N',nExch,nExch,nExch,
     &              (1.0_wp,0.0_wp),TMP,nExch,
     &                                Zm,nExch,
     &              (0.0_wp,0.0_wp), M2(L,:,:), nExch )

         ! spin moment
         Call zcopy_(nExch*nExch,[(0.0_wp,0.0_wp)],0,TMP,1)
         Call zgemm_('C','N',nExch,nExch,nExch,
     &              (1.0_wp,0.0_wp),ZM,nExch,
     &                              S2(L,:,:), nExch,
     &              (0.0_wp,0.0_wp),TMP,nExch )
         Call zcopy_(nExch*nExch,[(0.0_wp,0.0_wp)],0,S2(L,:,:),1)
         Call zgemm_('N','N',nExch,nExch,nExch,
     &              (1.0_wp,0.0_wp),TMP,nExch,
     &                                ZM,nExch,
     &              (0.0_wp,0.0_wp),S2(L,:,:),nExch )

      End Do ! L

      ! transform the exchange Hamiltonian accordingly:
      Call zcopy_(nExch*nExch,[(0.0_wp,0.0_wp)],0,TMP,1)
      Call zgemm_('C','N',nExch,nExch,nExch,
     &           (1.0_wp,0.0_wp),ZM,nExch,
     &                           H, nExch,
     &           (0.0_wp,0.0_wp),TMP,nExch )
      Call zcopy_(nExch*nExch,[(0.0_wp,0.0_wp)],0,H,1)
      Call zgemm_('N','N',nExch,nExch,nExch,
     &           (1.0_wp,0.0_wp),TMP,nExch,
     &                             ZM,nExch,
     &           (0.0_wp,0.0_wp),H,nExch )

      Call mma_deallocate(tmp)






      If (dbg) then
         Write(6,*) 'Matrix M2 (full ab initio):  X  FINAL'
         Do i=1,nExch
         Write(6,'(12(2F7.4,1x))') (M2(1,i,j),j=1,nExch)
         End Do

         Write(6,*) 'Matrix M3 (coupled basis): X  FINAL'
         Do i=1,nExch
         Write(6,'(12(2F7.4,1x))') (M3(1,i,j),j=1,nExch)
         End Do


         Write(6,*) 'Matrix M2 (full ab initio):  Y  FINAL'
         Do i=1,nExch
         Write(6,'(12(2F7.4,1x))') (M2(2,i,j),j=1,nExch)
         End Do

         Write(6,*) 'Matrix M3 (coupled basis): Y  FINAL'
         Do i=1,nExch
         Write(6,'(12(2F7.4,1x))') (M3(2,i,j),j=1,nExch)
         End Do



         Write(6,*) 'Matrix M2 (full ab initio):  Z  FINAL'
         Do i=1,nExch
         Write(6,'(12(2F7.4,1x))') (M2(3,i,j),j=1,nExch)
         End Do

         Write(6,*) 'Matrix M3 (coupled basis): Z  FINAL'
         Do i=1,nExch
         Write(6,'(12(2F7.4,1x))') (M3(3,i,j),j=1,nExch)
         End Do
      end if


      Write(6,'(A)') 'Exchange Hamiltonian in pseudospin basis --FIN'
      Do i=1,n1*n2
        Write(6,'(12(2F7.4,1x))') (H(i,j),j=1,n1*n2)
      End Do



      ! now we know the eigenstates ==  exact as in PA:
      ! proceed to extrat local M, S, L from the total exchange matrix
      Call zcopy_(3*n1*n1,[(0.0_wp,0.0_wp)],0,mom1,1)
      Call zcopy_(3*n1*n1,[(0.0_wp,0.0_wp)],0,som1,1)
      Call zcopy_(3*n1*n1,[(0.0_wp,0.0_wp)],0,lom1,1)
      Call zcopy_(3*n2*n2,[(0.0_wp,0.0_wp)],0,mom2,1)
      Call zcopy_(3*n2*n2,[(0.0_wp,0.0_wp)],0,som2,1)
      Call zcopy_(3*n2*n2,[(0.0_wp,0.0_wp)],0,lom2,1)

      f1=0.0_wp

      MM=(0.0_wp,0.0_wp)
      Do L=1,3
        Do isite=1,lmax
          Do nb1=1,n1*n2
            Do lp=1,lmax
              icoord(lp)=ibas(nb1,lp)
            End Do
              i1=nind(isite,1)
              j1=nind(isite,2)
            is1=ibas(nb1,isite)+1

            Do js1=1,nexchange(i1)
              icoord(isite)=js1-1
              nb2=norder(icoord,intc,lmax)

              if(i1==1) f1=1.0_wp/dble(n2)
              if(i1==2) f1=1.0_wp/dble(n1)

              MM(i1,l,is1,js1) = MM(i1,l,is1,js1) + f1*M3(l,nb1,nb2)
              SM(i1,l,is1,js1) = SM(i1,l,is1,js1) + f1*S3(l,nb1,nb2)

            End Do  ! js1
          End Do  ! nb1
        End Do  ! isite
      End Do  ! L

      Call prmom('extracted moment on site 1',MM(1,1:3,1:n1,1:n1),n1)
      Call prmom('extracted moment on site 2',MM(2,1:3,1:n2,1:n2),n2)

      Write(6,'(A)') 'g tensors for local pseudospins:'

      g1=0.0_wp
      g2=0.0_wp
      mg1=0.0_wp
      mg2=0.0_wp

      Call atens(MM(1,1:3,1:n1,1:n1), n1, g1, mg1, 2)
      Call atens(MM(2,1:3,1:n2,1:n2), n2, g2, mg2, 2)







!---------------------------------------------------------------
! re-write the exchange Hamiltonian as a four index array: H1(i1,j1,i2,j2)
      Call mma_allocate(H1,n1,n1,n2,n2,'H1')
      Call zcopy_(n1*n1*n2*n2,[(0.0_wp,0.0_wp)],0,H1,1)

      Do nb1 = 1,nExch
        Call icopy(lmax,0,0,icoord,1)
        Do i = 1,lmax
          icoord(i) = ibas(nb1,i)
        End Do

        i1  =   nind(1,1)
        i2  =   nind(2,1)
        is1 = icoord(1) + 1
        is2 = icoord(2) + 1
        Do js1 = 1, nexchange(i1)
          icoord(1) = js1-1
          Do js2 = 1, nexchange(i2)
            icoord(2) = js2-1
            nb2 = norder(icoord,intc,lmax)

            H1(is1,js1,is2,js2)=H(nb1,nb2)

            if(dbg) write(6,'(A,2I5,2x,4I2)')
     &          'nb1,nb2, is1,is2,js1,js2=', nb1,nb2, is1,is2,js1,js2
          End Do ! js2
        End Do ! js1
      End Do !nb1


!---------------------------------------------------------------
! 7. Compute the exchange contribution for the
!    LOCAL zero field splitting:
      Call mma_allocate(ZFS1,n1,n1,'ZFS1')
      Call mma_allocate(ZFS2,n2,n2,'ZFS2')
      Call zcopy_(n1*n1,[(0.0_wp,0.0_wp)],0,ZFS1,1)
      Call zcopy_(n2*n2,[(0.0_wp,0.0_wp)],0,ZFS2,1)

      If(n1>2) Then
         Do is1=1,n1
            Do js1=1,n1
               Do k=1,n2
                  ZFS1(is1,js1) = ZFS1(is1,js1) + H1(is1,js1,k,k)
               End Do
            End Do
         End Do
         If(dbg) Then
            Write(6,*)
            Write(6,'(A)') 'Local ZFS site 1'
            Do i=1,n1
               Write(6,'(12(2F7.4,1x))') (ZFS1(i,j),j=1,n1)
            End Do
         End If
      End If

      If(n2>2) Then
         Do is2=1,n2
            Do js2=1,n2
               Do k=1,n1
                  ZFS2(is2,js2) = ZFS2(is2,js2) + H1(k,k,is2,js2)
               End Do
            End Do
         End Do
         If(dbg) Then
            Write(6,*)
            Write(6,'(A)') 'Local ZFS site 2'
            Do i=1,n2
               Write(6,'(12(2F7.4,1x))') (ZFS2(i,j),j=1,n2)
            End Do
         End If
      End If


!---------------------------------------------------------------
! 8. Compute the exchange Hamiltonian, as difference between
!    TOTAL - ZFS(A) - ZFS(B)

      Call mma_allocate(H2,n1,n1,n2,n2,'H2')
      Call zcopy_(n1*n1*n2*n2,[(0.0_wp,0.0_wp)],0,H2,1)

      Do is1=1,n1
        Do js1=1,n1
            Do k=1,n2
c              H2(is1,js1,k,k) = H1(is1,js1,k,k) - ZFS1(is1,js1)
            End Do
        End Do
      End Do

      Do is2=1,n2
        Do js2=1,n2
            Do k=1,n1
c              H2(k,k,is2,js2) = H2(k,k,is2,js2) - ZFS2(is2,js2)
            End Do
        End Do
      End Do

c---------------------------------------------------------------
! 9. Project exchange matrix on products of ITO.

         Write(6,'(A)') 'Exchange Hamiltonian in pseudospin basis --FIN'
         Do i=1,n1*n2
           Write(6,'(12(2F7.3,1x))') (H(i,j),j=1,n1*n2)
         End Do

        ! transofrm the Hamiltonian:
          Call mma_allocate(H3,n1,n1,n2,n2,'H3')
          Call mma_allocate(rot1,3,3,'rot1')
          Call mma_allocate(rot2,3,3,'rot2')
          Call zcopy_(n1*n1*n2*n2,[(0.0_wp,0.0_wp)],0,H3,1)
          Call dcopy_(3*3,[0.0_wp],0,rot1,1)
          Call dcopy_(3*3,[0.0_wp],0,rot2,1)
          Do i=1,3
            rot1(i,i)=1.0_wp
            rot2(i,i)=1.0_wp
          End Do
c          iopt=2
c          Call transHam( n1, n2, rot1, rot2,
c     &                   MM(1,1:3,1:n1,1:n1), MM(2,1:3,1:n2,1:n2),
c     &                   itype(1), itype(2),
c     &                   H2(1:n1,1:n1, 1:n2,1:n2),
c     &                   H3(1:n1,1:n1, 1:n2,1:n2), iopt )
c
          JLin3( 1:(n1-1), -(n1-1):(n1-1),
     &           1:(n2-1), -(n2-1):(n2-1) )=(0.0_wp,0.0_wp)
          JLinC3(:,:)=0.0_wp

          Call JKQPar( n1, n2, H1(1:n1,1:n1,1:n2,1:n2),
     &                 JLin3( 1:(n1-1), -(n1-1):(n1-1),
     &                        1:(n2-1), -(n2-1):(n2-1) ) )

        ! print out the data:
        Write(6,'(A)')
        Write(6,'(10x,A)') 'Parameters of the ITOs:'
        Write(6,'( 5x,A)') 'with absolute values larger than:  0.5d-14 '
        Write(6,'(A)') '--- SITE 1 --|--- SITE 2 --|---------------'//
     &                 '----------------------|'
        Write(6,'(A)') ' rank | proj.| rank | proj.|    Lines  Exch'//
     &                 'ange  Interaction     |'
        Write(6,'(A)') '------|------|------|------|------ Real ---'//
     &                 '-------- Imag --------|'
        Do k1=1,n1-1,2
          Do q1=-k1,k1
            Do k2=1,n2-1,2
              Do q2=-k2,k2
                If(ABS(JLin3(k1,q1,k2,q2)) .gt. 0.5d-14 ) Then
                   Write(6,'(4(i4,2x,A),2(1x,E17.10),1x,A)')
     &                   k1,'|',q1,'|',k2,'|',q2,'|',
     &                   JLin3(k1,q1,k2,q2),'|'
                End If
              End Do
            End Do
          End Do
        End Do
        Write(6,'(A)') '------|------|------|------|---------------'//
     &                 '----------------------|'

        Call tensor2cart( JLin3( 1,-1:1, 1,-1:1), JLinC3 )

        Write(6,'(A)')
        Write(6,'(A)') 'Cartesian representation of the (rank-1)*'//
     &                 '(rank-1) exchange interaction: '
        Write(6,'(A)') 'Anisotropic exchange interaction:  J matrix:'
        Write(6,'(A)') 'LOCAL AXES:::'
        Write(6,'(A)') '     (  xx   xy  xz  )  '
        Write(6,'(A)') 'J =  (  yx   yy  yz  )  '
        Write(6,'(A)') '     (  zx   zy  zz  )  '
        Do i=1,3
          Write(6,'(3ES22.14)') (JLinC3(i,j),j=1,3)
        End Do

!-----------------------------------------------------------------------
      ! deallocate temporary memory
      Call mma_deallocate(g)
      Call mma_deallocate(mg)
      Call mma_deallocate(s1)
      Call mma_deallocate(m1)
      Call mma_deallocate(s2)
      Call mma_deallocate(m2)
      Call mma_deallocate(z)
      Call mma_deallocate(ibas)

      Call mma_deallocate(MM)
      Call mma_deallocate(SM)
      Call mma_deallocate(Z12)
      Call mma_deallocate(M3)
      Call mma_deallocate(S3)

      Call mma_deallocate(eso_1)
      Call mma_deallocate(eso_2)
      Call mma_deallocate(mom1)
      Call mma_deallocate(som1)
      Call mma_deallocate(lom1)
      Call mma_deallocate(mom2)
      Call mma_deallocate(som2)
      Call mma_deallocate(lom2)
      Call mma_deallocate(mom1r)
      Call mma_deallocate(mom2r)
      Call mma_deallocate(som1r)
      Call mma_deallocate(som2r)
      Call mma_deallocate(g1)
      Call mma_deallocate(g2)
      Call mma_deallocate(mg1)
      Call mma_deallocate(mg2)
      Call mma_deallocate(Z1)
      Call mma_deallocate(Z2)
      Call mma_deallocate(ZM)

      Call mma_deallocate(E1)
      Call mma_deallocate(H)
      Call mma_deallocate(H1)
      Call mma_deallocate(H2)
      Call mma_deallocate(rot1)
      Call mma_deallocate(rot2)
      Call mma_deallocate(H3)
      Call mma_deallocate(ZFS1)
      Call mma_deallocate(ZFS2)

      Return
      End subroutine project_exchange_single

