!DEC$ FREEFORM
!============================================================================
! COMMON BLOCK
!============================================================================
!
!----------------------------------------------------------------------------
INCLUDE 'GPARA.f'
!
! user varialbe field
double precision :: VARFLD(NVAROUT,NELOUT)
common VARFLD
!
! shape function matrix (interpolation matrix)
double precision, dimension(KGP,KNODE) :: KNN
common KNN
!
! integration point coordinates
double precision, dimension(KGP,KCORD) :: KGPCORD
common KGPCORD
!
! integration point weights
double precision, dimension(KGP) :: KWT
common KWT
!
! factor for quadrature rule
double precision :: KQUAD
common KQUAD
!
! identity tensor
double precision, dimension(KCORD,KCORD) :: KID
common KID
!
! flag parameter for UVARM routine using MPI
integer :: KFLAG
common KFLAG
!
!----------------------------------------------------------------------------