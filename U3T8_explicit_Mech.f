        !DEC$ FREEFORM
        !===============================================================================
        ! DEVELOPER: EMMA GRIFFITHS
        ! YEAR: 2018
        !===============================================================================
    
        ! include of ABAQUS subroutines
        !INCLUDE 'VEXTERNALDB.f'
        !INCLUDE 'UVARM.f'
        !INCLUDE 'SDVINI.f'
    
        module functions

            implicit none
            !---------------------------------------------------------------------------          
            public :: sgn,norm,trace,det,inv,dot,ddot,dya,cross,matvec,matmat
            
            contains
        
            !===========================================================================
            ! sign function of a scalar
            function sgn(sca) 
           
                !-----------------------------------------------------------------------  
                ! declaration
        
                double precision :: sca  
                double precision :: sgn    
                integer :: i
                !-----------------------------------------------------------------------

                if (sca .gt. 0.d0) then
                    sgn = 1.d0
                elseif (sca .lt. 0.d0) then
                    sgn = -1.d0
                elseif (sca == 0.d0) then
                    sgn = 0.d0
                end if
        
            end function sgn
            !===========================================================================
            
            !===========================================================================
            ! euclidean norm of a vector
            function norm(vec) 
           
                !-----------------------------------------------------------------------  
                ! declaration
        
                double precision, dimension(:) :: vec
                double precision :: norm
        
                integer :: i
                !-----------------------------------------------------------------------
        
                norm = 0.
                if (size(vec) == 3) then
                    norm = (vec(1)**2.d0+vec(2)**2.d0+vec(3)**2.d0)**(0.5d0)
                elseif (size(vec) == 2) then
                    norm = (vec(1)**2.d0+vec(2)**2.d0)**(0.5d0)
                else
                    do i=1,size(vec)
                        norm = norm + vec(i)**2.d0
                    end do
                    norm = sqrt(norm)
                end if
        
            end function norm
           !============================================================================
    
            !===========================================================================
            ! trace of a tensor
            function trace(mat) 
           
                !-----------------------------------------------------------------------  
                ! declaration
        
                double precision, dimension(:,:) :: mat  
                double precision :: trace    
                integer :: i
                !-----------------------------------------------------------------------
        
                trace = 0.d0
                if (size(mat,1) == size(mat,2)) then
                    do i=1,size(mat,1)
                        trace = trace + mat(i,i)
                    end do
                else
                    stop "Error in function `trace' - matrix is non-sqare!"
                end if
        
            end function trace
            !===========================================================================
    
            !===========================================================================
            ! determinant of a tensor
            function det(mat)
    
                !-----------------------------------------------------------------------    
                ! declaration
 
                double precision, dimension(:,:) :: mat
                double precision :: det
                !-----------------------------------------------------------------------
        
                det = 0.d0
                if (size(mat,1)==size(mat,2) .and. size(mat,1)==3) then
                    det = mat(1,1)*mat(2,2)*mat(3,3) &
                        + mat(1,2)*mat(2,3)*mat(3,1) &
                        + mat(1,3)*mat(2,1)*mat(3,2) &
                        - mat(1,3)*mat(2,2)*mat(3,1) &
                        - mat(1,2)*mat(2,1)*mat(3,3) &
                        - mat(1,1)*mat(2,3)*mat(3,2)
                else
                    stop "Error in function `det' - tensor is non-sqare!"
                end if
                    
            end function det
            !===========================================================================
    
            !===========================================================================
            ! inverse of a tensor
            function inv(mat)
           
                !-----------------------------------------------------------------------    
                ! declaration
        
                double precision, dimension(:,:) :: mat
                double precision, dimension(size(mat,1),size(mat,2)) :: inv
                !-----------------------------------------------------------------------
        
                inv = 0.d0
                if (size(mat,1)==size(mat,2) .and. size(mat,1)==3 ) then

                    inv(1,1) = mat(2,2)*mat(3,3)-mat(2,3)*mat(3,2)
                    inv(1,2) = mat(1,3)*mat(3,2)-mat(1,2)*mat(3,3)
                    inv(1,3) = mat(1,2)*mat(2,3)-mat(1,3)*mat(2,2)
                
                    inv(2,1) = mat(2,3)*mat(3,1)-mat(2,1)*mat(3,3)
                    inv(2,2) = mat(1,1)*mat(3,3)-mat(1,3)*mat(3,1)
                    inv(2,3) = mat(1,3)*mat(2,1)-mat(1,1)*mat(2,3)
                
                    inv(3,1) = mat(2,1)*mat(3,2)-mat(2,2)*mat(3,1)
                    inv(3,2) = mat(1,2)*mat(3,1)-mat(1,1)*mat(3,2)
                    inv(3,3) = mat(1,1)*mat(2,2)-mat(1,2)*mat(2,1)
                
                    inv=inv/det(mat)
                    
                elseif (size(mat,1)==size(mat,2) .and. size(mat,1)==2 ) then
                    
                    inv(1,1) =  mat(2,2)
                    inv(1,2) = -mat(1,2)
                    inv(2,1) = -mat(2,1)
                    inv(2,2) =  mat(1,1)
                    
                    inv=inv/det(mat)
                    
                elseif (size(mat,1)==size(mat,2) .and. size(mat,1)==1 ) then
                    
                    inv(1,1) = 1.d0/mat(1,1)
                    
                else
                    stop "Error in function `inv' - tensor is non-sqare or larger then 3x3!"
                end if

            end function inv
            !===========================================================================
    
            !===========================================================================
            ! scalar/dot product of two vectors of the same dimension
            function dot(vec1,vec2)
    
                !-----------------------------------------------------------------------   
                ! declaration
        
                double precision, dimension(:) :: vec1,vec2
                double precision :: dot
        
                integer :: i
                !-----------------------------------------------------------------------
        
                dot = 0.d0
                if (size(vec1)==size(vec2)) then  
                    do i=1,size(vec1)
                        dot = dot + vec1(i)*vec2(i)
                    end do
                else
                    stop "Error in function `dot' - vectors have not the same length!"
                end if
        
            end function dot
            !===========================================================================
            
            !============================================================================
            ! scalar/dot product of two vectors of the same dimension
            function cross(vec1,vec2)
    
                !-----------------------------------------------------------------------   
                ! declaration
        
                double precision, dimension(:) :: vec1,vec2
                double precision, dimension(size(vec1)) :: cross
        
                integer :: i
                !-----------------------------------------------------------------------
        
                cross = 0.d0
                if ((size(vec1)==size(vec2)) .and. (size(vec1)==3)) then  
                    cross(1) = vec1(2)*vec2(3)-vec1(3)*vec2(2)
                    cross(2) = vec1(3)*vec2(1)-vec1(1)*vec2(3)
                    cross(3) = vec1(1)*vec2(2)-vec1(2)*vec2(1)
                else
                    stop "error in `cross' - vector lengths are not the same or not given in 3D!"
                end if
        
            end function cross
           !============================================================================
    
            !===========================================================================
            ! scalar/dot product of two tensors of the same dimensions
            function ddot(mat1,mat2)
    
                !-----------------------------------------------------------------------   
                ! declaration
        
                double precision, dimension(:,:) :: mat1,mat2
                double precision :: ddot
        
                integer :: i,j
                !-----------------------------------------------------------------------
        
                ddot = 0.d0
                if (size(mat1,1)==size(mat2,1) .and. size(mat1,2)==size(mat2,2)) then  
                    do i=1,size(mat1,1)
                        do j=1,size(mat1,2)
                            ddot = ddot + mat1(i,j)*mat2(i,j)
                        end do
                    end do
                else
                    stop "Error in function `ddot' - tensor dimensions are not the same!"
                end if
        
            end function ddot
            !===========================================================================
    
            !===========================================================================
            ! dyadic prodcut of two vectors of the same length
            function dya(vec1,vec2)
           
                !-----------------------------------------------------------------------    
                ! declaration  
        
                integer :: i,j   
                double precision, dimension(:) :: vec1,vec2
                double precision, dimension(size(vec1),size(vec2)) :: dya
                !-----------------------------------------------------------------------

                dya = 0.d0
                if (size(vec1)==size(vec2)) then
        
                    do i=1,size(vec1)
                        do j=1,size(vec1)
            
                            dya(i,j) = vec1(i)*vec2(j)
            
                        end do
                    end do
        
                else
                    stop "Error in function `dya' - vector lengths are not the same!"
                end if
        
            end function dya
            !===========================================================================
        
            !===========================================================================
            ! matrix-vector operation
            function matvec(mat,vec)
           
                !-----------------------------------------------------------------------    
                ! declaration
        
                double precision, dimension(:,:) :: mat
                double precision, dimension(:) :: vec
                double precision, allocatable, dimension(:) :: matvec      
                integer :: i,j
                !-----------------------------------------------------------------------
        
                if (size(mat,2) == size(vec)) then
        
                    allocate(matvec(size(mat,1)))
                    matvec = 0.d0
                    do i=1,size(mat,1)
                        do j=1,size(mat,2)
                            matvec(i) = matvec(i) + mat(i,j)*vec(j)
                        end do
                    end do
            
                else 
                    stop "Dimension error in function `matvec' - dimension of vector must be consistent with the dimensions of the matrix"
                end if
         
            end function matvec
            !===========================================================================
    
            !===========================================================================
            ! matrix-matrix operation (only of square matrices)
            function matmat(mat1,mat2)
           
                !-----------------------------------------------------------------------    
                ! declaration
        
                double precision, dimension(:,:) :: mat1
                double precision, dimension(:,:) :: mat2        
                double precision, dimension(size(mat1,1),size(mat2,2)) :: matmat
                !-----------------------------------------------------------------------
        
                matmat=0.d0
                if (size(mat1,2) == size(mat2,1)) then

                    matmat(1,1) = mat1(1,1)*mat2(1,1)+mat1(1,2)*mat2(2,1)+mat1(1,3)*mat2(3,1)
                    matmat(1,2) = mat1(1,1)*mat2(1,2)+mat1(1,2)*mat2(2,2)+mat1(1,3)*mat2(3,2)
                    matmat(1,3) = mat1(1,1)*mat2(1,3)+mat1(1,2)*mat2(2,3)+mat1(1,3)*mat2(3,3)
                
                    matmat(2,1) = mat1(2,1)*mat2(1,1)+mat1(2,2)*mat2(2,1)+mat1(2,3)*mat2(3,1)
                    matmat(2,2) = mat1(2,1)*mat2(1,2)+mat1(2,2)*mat2(2,2)+mat1(2,3)*mat2(3,2)
                    matmat(2,3) = mat1(2,1)*mat2(1,3)+mat1(2,2)*mat2(2,3)+mat1(2,3)*mat2(3,3)
                
                    matmat(3,1) = mat1(3,1)*mat2(1,1)+mat1(3,2)*mat2(2,1)+mat1(3,3)*mat2(3,1)
                    matmat(3,2) = mat1(3,1)*mat2(1,2)+mat1(3,2)*mat2(2,2)+mat1(3,3)*mat2(3,2)
                    matmat(3,3) = mat1(3,1)*mat2(1,3)+mat1(3,2)*mat2(2,3)+mat1(3,3)*mat2(3,3)

                else 
                    stop "Dimension error in function `matmat' - matrix dimensions are not consistent or larger than 3x3"
                end if

            end function matmat
            !===========================================================================
        
        end module functions
        !===============================================================================
    

        module constitutive_relations

            implicit none
            !---------------------------------------------------------------------------          
            public :: stresses_lin_elastic
            
            contains
    
            !===========================================================================
            ! linear elasticity (Hooke's law)
            subroutine stresses_lin_elastic(SIG,P,S,F,Ee,pID,pGM,pLAM)
        
                use functions, only : trace,matmat,det,inv
       
                implicit none
                !---------------------------------------------------------------------------  
                ! declaration
         
                double precision, dimension(:,:), intent(inout) :: SIG,P,S
                double precision, dimension(:,:), intent(in) :: F,Ee,pID
                double precision, intent(in) :: pGM,pLAM
                !---------------------------------------------------------------------------

                S = 2.d0*pGM*Ee + pLAM*trace(Ee)*pID
            
                P = S
                SIG = S
            
                !---------------------------------------------------------------------------       
            
            end subroutine stresses_lin_elastic
            !===========================================================================
    
        end module constitutive_relations
        !===============================================================================
    
    
        !===============================================================================
        !-------------------------------------------------------------------------------
        ! USER SUBROUTINE - USER ELEMENT
        
        ! abaqus interface VUEL
        SUBROUTINE VUEL(nblock,rhs,amass,dtimeStable,svars,nsvars, &
                        energy, &
                        nnode,ndofel,props,nprops,jprops,njprops, &
                        coords,ncrd,u,du,v,a, &
                        jtype,jElem, &
                        time,period,dtimeCur,dtimePrev,kstep,kinc, &
                        lflags, &
                        dMassScaleFactor, &
                        predef,npredef, &
                        jdltyp, adlmag)
                        
            use functions
            use constitutive_relations
            
            implicit none ! use explicit=both and output_precision=full
            !include 'vaba_param.inc'
            ! if implicit none is not used and 'vaba_param.inc' is included - implicit declaration of variables:
            ! all variables for letters a-h, o-z is real and i-n are integers (i,j,k,l,m,n).
            !-------------------------------------------------------------------------------   
            ! declaration
            
            ! parameters - numbers
            !parameter ( zero=0.d0, half=0.5d0, one=1.d0, two=2.d0, four=4.d0, six=6.d0)
            double precision, parameter :: zero=0.d0
            double precision, parameter :: half=0.5d0
            double precision, parameter :: one=1.d0
            double precision, parameter :: two=2.d0
            double precision, parameter :: three=3.d0
            double precision, parameter :: four=4.d0
            double precision, parameter :: six=6.d0
            double precision, parameter :: eight=8.d0
            double precision, parameter :: Abcissa=SQRT(0.d0/3.d0)
            double precision, parameter :: factorStable=0.99d0
            
            ! parameters - problem specification
            !parameter ( iGP=8, iCORD=3, iNODE=8 )
            integer, parameter :: iGP=8
            integer, parameter :: iCORD=3
            integer, parameter :: iNODE=8
            
            !! parameters - material
            !double precision, parameter :: pEM=79000.d0
            !double precision, parameter :: pNU=0.44d0
            !double precision, parameter :: pRHO_Au=19.3d-06
            !double precision, parameter :: pGM=half*pEM/(one+pNU)
            !double precision, parameter :: pLAM = (pEM*pNU)/((one+pNU)*(one-two*pNU))

            ! Abaqus variables
            ! ================  
            
            ! predefined parameters - operational code keys
            integer, parameter :: jMassCalc = 1
            integer, parameter :: jIntForceAndDtStable = 2
            integer, parameter :: jExternForce = 3
            
            ! predefined parameters - flag indices
            integer, parameter :: iProcedure = 1
            integer, parameter :: iNlgeom = 2
            integer, parameter :: iOpCode = 3
            integer, parameter :: nFlags = 3
            
            !  predefined parameters - procedure flags
            integer, parameter :: jDynExplicit = 17
            
            ! predefined parameters - energy array indices
            integer, parameter :: iElPd = 1
            integer, parameter :: iElCd = 2
            integer, parameter :: iElIe = 3
            integer, parameter :: iElTs = 4
            integer, parameter :: iElDd = 5
            integer, parameter :: iElBv = 6
            integer, parameter :: iElDe = 7
            integer, parameter :: iElHe = 8
            integer, parameter :: iElKe = 9
            integer, parameter :: iElTh = 10
            integer, parameter :: iElDmd = 11
            integer, parameter :: iElDc = 12
            integer, parameter :: nElEnergy = 12

            ! predefined parameters - time indices
            integer, parameter :: iStepTime = 1
            integer, parameter :: iTotalTime = 2
            integer, parameter :: nTime = 2
            
            ! predefined parameters - predefined variables indices
            integer, parameter :: iPredValueNew = 1
            integer, parameter :: iPredValueOld = 2
            integer, parameter :: nPred = 2

            ! variables passed in for information
            integer, intent(in ) :: nblock                              ! number of user elements to be processed in this call to VUEL
            double precision, intent(in ) :: dtimeCur                   ! current time increment
            double precision, intent(in ) :: dtimePrev                  ! previous time increment
            double precision, intent(in ) :: period                     ! time period of the current step
            integer, intent(in ) :: ndofel                              ! number of dofs in the element
            integer, intent(in ) :: nsvars                    	        ! user-defined number of solution-dependent variables            
            integer, intent(in ) :: nprops                              ! user-defined number of real property values
            integer, intent(in ) :: njprops                             ! user-defined number of integer property values 
            integer, intent(in ) :: ncrd                                ! max. number of user-defined coordinates at any nodal point
            integer, intent(in ) :: nnode                               ! user-defined number of nodes on the element
            integer, intent(in ) :: jtype                               ! integer defining the element type
            integer, intent(in ) :: kstep 	                            ! current step number
            integer, intent(in ) :: kinc 	                            ! current increment number
            integer, intent(in ) :: npredef                             ! number of predefined field variables (incl. temperature)
            integer, intent(in ) :: jdltyp                              ! specifies the load type
            
            !dimension rhs( nblock,ndofel), amass(nblock,ndofel,ndofel), & 
            !            dtimeStable(nblock), & 
            !            svars(nblock,nsvars), energy(nblock,nElEnergy), & 
            !            props(nprops), jprops(njprops), & 
            !            jelem(nblock), time(nTime), lflags(nFlags), & 
            !            coords(nblock,nnode,ncrd), u(nblock,ndofel), & 
            !            du(nblock,ndofel), v(nblock,ndofel), a(nblock, ndofel), & 
            !            dMassScaleFactor(nblock), & 
            !            predef(nblock, nnode, npredef, nPred), adlmag(nblock) 
            
            double precision, intent(in ), dimension(nprops) :: props   ! real property values
            integer, intent(in ), dimension(njprops) :: jprops          ! integer property values
            double precision, dimension(nblock,nnode,ncrd) :: coords    ! block of original coordinates of nodes (reference configuration)
            integer, intent(in ), dimension(nblock) :: jElem 	        ! block of element numbers
            double precision, intent(in ), dimension(nblock) :: adlmag 	! block of total load magnitute of type jdltyp
            double precision, intent(in ), dimension(nblock,nnode,npredef,nPred) :: predef 	! block of values of predefined field variables (only uncoupled)
            integer, intent(in ), dimension(nFlags) :: lflags 	        ! array containing the flags of current solution procedure
            double precision, intent(in ), dimension(nblock) :: dMassScaleFactor            ! block containing mass scale factors (for each element)
            double precision, intent(in ), dimension(nTime) :: time     ! current value of step time
            double precision, intent(in ), dimension(nblock,ndofel) :: u    ! block of nodal solution variables (displacement, temperature, etc.)
            double precision, intent(in ), dimension(nblock,ndofel) :: du   ! block of incremental values (displacement, temperature, etc.)
            double precision, intent(in ), dimension(nblock,ndofel) :: v    ! block of time rates of nodal solution variables
            double precision, intent(in ), dimension(nblock,ndofel) :: a    ! block of acclerations of nodal solution variables
            
            ! variables to be defined (block-wise)
            double precision, intent(inout), dimension(nblock,ndofel) :: rhs 	        ! element contribution to right-hand-side vector
            double precision, intent(inout), dimension(nblock,ndofel,ndofel) :: amass  ! element contribution to mass matrix (symmetric)
            double precision, intent(inout), dimension(nblock) :: dtimeStable  	    ! scalar value, upper limit of time increment (Courant condition)
            double precision, intent(inout), dimension(nblock,nsvars) :: svars  	    ! solution-depending state variables
            double precision, intent(inout), dimension(nblock,nElEnergy) :: energy		! element energy quantities 
            
            ! user variables
            ! ==================
            
            double precision :: pGPCORD(iGP,iCORD)  ! integration point coordinates
            double precision :: pQUAD               ! factor for quadrature rule
            double precision :: pWT(iGP)            ! integration point weights
            double precision :: pNN(iGP,iNODE)      ! shape function matrix (interpolation matrix)
            double precision :: pID(iCORD,iCORD)    ! identity tensor
            
            double precision :: amass_row_sum       ! sum of array row
            double precision :: cd,volume,pd_min

            ! include others
            !INCLUDE 'COMMON.f'

            !===============================================================================

            double precision :: xi1,xi2,xi3         ! natural coordinates
            double precision :: X1(iNODE)           ! physical coordinates
            double precision :: X2(iNODE)           ! physical coordinates
            double precision :: X3(iNODE)           ! physical coordinates
            double precision :: dNdXi1(iNODE)       ! shape function derivatives
            double precision :: dNdXi2(iNODE)       ! shape function derivatives
            double precision :: dNdXi3(iNODE)       ! shape function derivatives
            double precision :: dNdX1(iGP,iNODE)    ! shape function derivatives
            double precision :: dNdX2(iGP,iNODE)    ! shape function derivatives
            double precision :: dNdX3(iGP,iNODE)    ! shape function derivatives
            double precision :: dX1dxi1             ! derivatives
            double precision :: dX1dxi2             ! derivatives
            double precision :: dX1dxi3             ! derivatives
            double precision :: dX2dxi1             ! derivatives
            double precision :: dX2dxi2             ! derivatives
            double precision :: dX2dxi3             ! derivatives
            double precision :: dX3dxi1             ! derivatives
            double precision :: dX3dxi2             ! derivatives
            double precision :: dX3dxi3             ! derivatives
            double precision :: JJ(iCORD,iCORD)     ! Jacobi matrix
            double precision :: detJ(iGP)           ! Jacobi-determinant (reference)
            double precision :: Ux(iNODE)           ! nodal displacement in x
            double precision :: Uy(iNODE)           ! nodal displacement in y
            double precision :: Uz(iNODE)           ! nodal displacement in z
            double precision :: H(iCORD,iCORD)      ! displacement gradient
            double precision :: F(iCORD,iCORD)      ! deformation gradient
            double precision :: Ee(iCORD,iCORD)     ! small /Gree-Lagrange strain tensor
            double precision :: SIG(iCORD,iCORD)    ! Cauchy stress tensor
            double precision :: S(iCORD,iCORD)      ! 2. PK
            double precision :: P(iCORD,iCORD)      ! 1. PK
            !double precision :: vMeq
            integer :: dofni(iCORD)        ! current dof
            integer :: dofnj(iCORD)        ! current dof
            double precision :: pRHO                 ! density
            double precision :: pEM
            double precision :: pNU
            double precision :: pGM
            double precision :: pLAM
                
            integer :: CordRange(iCORD) = (/(i, i=1,iCORD, 1)/)
            integer :: NODERange(iNODE) = (/(i, i=1,iNODE, 1)/)
            
            ! integer
            integer :: kblock,ip,nn,ni,nj,i

            !-------------------------------------------------------------------------------  

            ! element ingoing note =======================
            !write(*,*)" "
            !write(*,*)"VUEL in"
            !write(*,*)"kinc",kinc
            !write(*,*)"Abcissa",Abcissa
            !write(*,*)"ndofel",ndofel
            !write(*,*)"ncrd",ncrd
            !write(*,*)"kstep",kstep
            !write(*,*)"jtype",jtype
            !write(*,*)"nblock",nblock
            !write(*,*)"props",props
            !write(*,*)"block",jElem(1)," to ",jElem(nblock)
            !write(*,*)"lflags(iOpCode)",lflags(iOpCode)
            !write(*,*)"for mass calculation:",jMassCalc
            !write(*,*)"for internal force and stable increment calculation:",jIntForceAndDtStable
            !write(*,*)"(lflags(iOpCode).eq.jIntForceAndDtStable):",(lflags(iOpCode).eq.jIntForceAndDtStable)
            !write(*,*)" "
            ! ============================================
            
                        
            if (jtype.eq.38 .and. lflags(iProcedure).eq.jDynExplicit) then 
                
                pEM  = props(1)
                pNU  = props(2)
                pRHO = props(3)
                pGM  = half*pEM/(one+pNU)
                pLAM = (pEM*pNU)/((one+pNU)*(one-two*pNU))
                    
                ! integration point coordinates and weights --------------------------------
                if (iGP==8) then

                    pGPCORD(1,:) = (/ -Abcissa, -Abcissa, -Abcissa /)
                    pGPCORD(2,:) = (/  Abcissa, -Abcissa, -Abcissa /)
                    pGPCORD(3,:) = (/  Abcissa,  Abcissa, -Abcissa /)
                    pGPCORD(4,:) = (/ -Abcissa,  Abcissa, -Abcissa /)
                    pGPCORD(5,:) = (/ -Abcissa, -Abcissa,  Abcissa /)
                    pGPCORD(6,:) = (/  Abcissa, -Abcissa,  Abcissa /)
                    pGPCORD(7,:) = (/  Abcissa,  Abcissa,  Abcissa /)
                    pGPCORD(8,:) = (/ -Abcissa,  Abcissa,  Abcissa /)
                        
                    pWT(:) = (/ one, one, one, one, one, one, one, one /)
                        
                    pQUAD = one
                        
                else
        
                    stop "Error in computation of integration point coordinates. The number of IPs does not conform with the element type (4 node tetrahedral element with 4 integratin points)."
        
                end if
                !---------------------------------------------------------------------------
                    
                ! shape function for each ip -----------------------------------------------
                do ip=1,iGP
        
                    xi1=pGPCORD(ip,1)
                    xi2=pGPCORD(ip,2)
                    xi3=pGPCORD(ip,3)
        
                    if (iNODE==8) then ! cf. WRIGGERS - Nonlinear Finite Elemente Methods 2008 (p. 120)
        
                        pNN(ip,1) = one/eight*(1-xi1)*(1-xi2)*(1-xi3)
                        pNN(ip,2) = one/eight*(1+xi1)*(1-xi2)*(1-xi3)
                        pNN(ip,3) = one/eight*(1+xi1)*(1+xi2)*(1-xi3)
                        pNN(ip,4) = one/eight*(1-xi1)*(1+xi2)*(1-xi3)
                        pNN(ip,5) = one/eight*(1-xi1)*(1-xi2)*(1+xi3)
                        pNN(ip,6) = one/eight*(1+xi1)*(1-xi2)*(1+xi3)
                        pNN(ip,7) = one/eight*(1+xi1)*(1+xi2)*(1+xi3)
                        pNN(ip,8) = one/eight*(1-xi1)*(1+xi2)*(1+xi3)
                        
        
                    else
                        stop "Error in computation of shape functions. The number of nodes does not conform with the element type (4 node tetrahedral element)."
                    end if
        
                end do
                
                !---------------------------------------------------------------------------
                    
                ! identity matrix ----------------------------------------------------------
                pID=zero
                forall(i=1:iCORD) pID(i,i)=one
                !---------------------------------------------------------------------------
                    
                ! loop over element block
                do kblock=1,nblock ! ---------------------------------------------------
            
                    ! loop over all integration points (computation of FE variables)
                    do ip=1,iGP ! ------------------------------------------------------
            
                        ! get solution-dependent state variables (history variables)
            
                        ! natural coordinates of current ip
                        xi1 = pGPCORD(ip,1)
                        xi2 = pGPCORD(ip,2)
                        xi3 = pGPCORD(ip,3)
            
                        ! coordinate vectors of current element
                        X1 = coords(kblock,:,1)
                        X2 = coords(kblock,:,2)
                        x3 = coords(kblock,:,3)
        
                        ! --------------------------------------------------------------
                        if (iNODE==8) then
        
                            ! derivatives of shape functions with respect to natural coordinates                
                            dNdXi1(1) = -one/eight*(1-xi2)*(1-xi3)
                            dNdXi2(1) = -one/eight*(1-xi1)*(1-xi3)
                            dNdXi3(1) = -one/eight*(1-xi1)*(1-xi2)
            
                            dNdXi1(2) =  one/eight*(1-xi2)*(1-xi3)
                            dNdXi2(2) =  -one/eight*(1+xi1)*(1-xi3)
                            dNdXi3(2) =  -one/eight*(1+xi1)*(1-xi2)
            
                            dNdXi1(3) =  one/eight*(1+xi2)*(1-xi3)
                            dNdXi2(3) =  one/eight*(1+xi1)*(1-xi3)
                            dNdXi3(3) =  -one/eight*(1+xi1)*(1+xi2)
            
                            dNdXi1(4) =  -one/eight*(1+xi2)*(1-xi3)
                            dNdXi2(4) =  one/eight*(1-xi1)*(1-xi3)
                            dNdXi3(4) =  -one/eight*(1-xi1)*(1+xi2)
            
                            dNdXi1(5) =  -one/eight*(1-xi2)*(1+xi3)
                            dNdXi2(5) =  -one/eight*(1-xi1)*(1+xi3)
                            dNdXi3(5) =  one/eight*(1-xi1)*(1-xi2)
            
                            dNdXi1(6) =  one/eight*(1-xi2)*(1+xi3)
                            dNdXi2(6) =  -one/eight*(1+xi1)*(1+xi3)
                            dNdXi3(6) =  one/eight*(1+xi1)*(1-xi2)
            
                            dNdXi1(7) =  one/eight*(1+xi2)*(1+xi3)
                            dNdXi2(7) =  one/eight*(1+xi1)*(1+xi3)
                            dNdXi3(7) =  one/eight*(1+xi1)*(1+xi2)
            
                            dNdXi1(8) =  -one/eight*(1+xi2)*(1+xi3)
                            dNdXi2(8) =  one/eight*(1-xi1)*(1+xi3)
                            dNdXi3(8) =  one/eight*(1-xi1)*(1+xi2)
                         
                        else 
                            stop "Error in computation of shape function derivatives. The number of nodes does not conform with the element type (4 node tetrahedral element)."     
                        end if
        
                        ! derivatives of physical coordinates with respect to natural coordinates                
                        dX1dxi1=dot(X1,dNdXi1)
                        dX1dxi2=dot(X1,dNdXi2)
                        dX1dxi3=dot(X1,dNdXi3)
        
                        dX2dxi1=dot(X2,dNdXi1)
                        dX2dxi2=dot(X2,dNdXi2)
                        dX2dxi3=dot(X2,dNdXi3)
        
                        dX3dxi1=dot(X3,dNdXi1)
                        dX3dxi2=dot(X3,dNdXi2)
                        dX3dxi3=dot(X3,dNdXi3)
        
                        ! Jacobian determinant
                        JJ(1,:) = (/dX1dxi1, dX2dxi1, dX3dxi1/)
                        JJ(2,:) = (/dX1dxi2, dX2dxi2, dX3dxi2/)
                        JJ(3,:) = (/dX1dxi3, dX2dxi3, dX3dxi3/)
                        
                        detJ(ip) = det(JJ)

                        ! derivatives of shape functions with respect to physical coordinates  
                        do nn=1,iNODE
                            dNdX1(ip,nn) = one/detJ(ip)*( (dX2dxi2*dX3dxi3-dX3dxi2*dX2dxi3)*dNdXi1(nn) &
                                                        + (dX3dxi1*dX2dxi3-dX2dxi1*dX3dxi3)*dNdXi2(nn) &
                                                        + (dX2dxi1*dX3dxi2-dX3dxi1*dX2dxi2)*dNdXi3(nn) )
                            dNdX2(ip,nn) = one/detJ(ip)*( (dX3dxi2*dX1dxi3-dX1dxi2*dX3dxi3)*dNdXi1(nn) &
                                                        + (dX1dxi1*dX3dxi3-dX3dxi1*dX1dxi3)*dNdXi2(nn) &
                                                        + (dX3dxi1*dX1dxi2-dX1dxi1*dX3dxi2)*dNdXi3(nn) )
                            dNdX3(ip,nn) = one/detJ(ip)*( (dX1dxi2*dX2dxi3-dX2dxi2*dX1dxi3)*dNdXi1(nn) &
                                                        + (dX2dxi1*dX1dxi3-dX1dxi1*dX2dxi3)*dNdXi2(nn) &
                                                        + (dX1dxi1*dX2dxi2-dX2dxi1*dX1dxi2)*dNdXi3(nn) )
                        end do
                    end do
                                        
                    if (lflags(iOpCode).eq.jMassCalc) then ! compute mass matrix
                        !amass(kblock,1:ndofel,1:ndofel)=zero
                        ! loop over all integration points (computation of mass matrix)
                        do ip=1,iGP ! ------------------------------------------------------
    
                            ! summation over node_i
                            do ni=1,iNODE !-----------------------------loop-i--------------
                
                                ! current node dof
                                dofni(1:iCORD) = 1+iCORDTOTAL*((ni*1)-1)+(CordRange-1)
                                
                                ! summation over node_i
                                do nj=1,iNODE !-------------------------loop-j--------------
                
                                    ! current node dof
                                    dofnj(1:iCORD) = 1+iCORDTOTAL*((nj*1)-1)+(CordRange-1)
                                    
                                    ! regular mass matrix
                                    amass(kblock,dofni,dofnj) = amass(kblock,dofni,dofnj) &
                                                              + pQUAD*pWT(ip)*detJ(ip)*pRHO*matmat(pNN(ip,ni)*pID,pNN(ip,nj)*pID)
                                    
                                end do !--------------------------end-loop-j----------------
            
                            end do !------------------------------end-loop-i----------------
            
                        end do ! ----------------------------end-loop-ip--------------------
                        
                        ! mass lumbing
                        do i=1,ndofel
                            amass_row_sum = sum(amass(kblock,i,:))
                            amass(kblock,i,:) = zero
                            amass(kblock,i,i) = amass_row_sum
                        end do
                            
                            
                    elseif (lflags(iOpCode).eq.jIntForceAndDtStable) then !compute internal force + stable time increment    
                
                        rhs(kblock,1:ndofel)=zero
                        ! loop over all integration points (computation of residuum)
                        do ip=1,iGP ! ------------------------------------------------------
            
                            ! displacement gradient
                            Ux = u(kblock,1:iCORD*iNODE:iCORD)
                            Uy = u(kblock,2:iCORD*iNODE:iCORD)
                            Uz = u(kblock,3:iCORD*iNODE:iCORD)
                            H = zero
                            do ni=1,iNODE
                                H = H + dya( (/Ux(ni),Uy(ni),Uz(ni)/),(/dNdX1(ip,ni),dNdX2(ip,ni),dNdX3(ip,ni)/) )
                            end do
            
                            ! small strain tensor
                            Ee = half*(transpose(H) + H)
                            ! stress
                            call stresses_lin_elastic(SIG,P,S,F,Ee,pID,pGM,pLAM)
            
                            ! summation over node_i
                            do ni=1,iNODE !-----------------------------loop-i--------------
                
                                ! current node dof
                                dofni(1:iCORD) = 1+iCORDTOTAL*((ni*1)-1)+(CordRange-1)
                
                                ! internal residual force vector
                                rhs(kblock,dofni) = rhs(kblock,dofni) + pQUAD*pWT(ip)*detJ(ip)*matvec(P,(/dNdX1(ip,ni),dNdX2(ip,ni),dNdX3(ip,ni)/))
            
                            end do !------------------------------end-loop-i----------------
            
                        end do ! ----------------------------end-loop-ip--------------------
                        
                        cd = sqrt( (pEM*(one-pNU))/(pRHO*(one+pNU)*(one-two*pNU)) )
                        volume = detJ(1)+detJ(2)+detJ(3)+detJ(4)+detJ(5)+detJ(6)+detJ(7)+detJ(8)
                        pd_min = volume**(one/three)
                        dtimeStable(kblock) = factorStable*(pd_min/cd)
                        !dtimeStable(kblock) = 2.143102d-06
                        
                        !if (kblock==1) then
                        !    write(*,*)"dtimeStable(kblock)"
                        !    write(*,*)dtimeStable(kblock)
                        !    write(*,*)"rhs(kblock=1,:)"
                        !    write(*,*)rhs(kblock,:)
                        !end if
                                        
                    end if
                end do !----------------------------------------------------------------
                
            end if
                        
            ! ============================================
            ! element outgoing note
            !write(*,*)" "
            !write(*,*)"UEL out - ELEMNO",ELEMNO
            !write(*,*)" "
            ! ============================================

        return
        end subroutine VUEL
        !-------------------------------------------------------------------------------
        !===============================================================================
