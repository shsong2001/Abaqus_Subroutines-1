        !DEC$ FREEFORM
        !===============================================================================
        ! DEVELOPER: EDGAR HUSSER
        ! YEAR: 2017
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
            double precision, parameter :: factorStable=0.99d0
            
            ! parameters - problem specification
            !parameter ( iGP=1, iCORD=3, iNODE=4 )
            integer, parameter :: iGP=1
            integer, parameter :: iGPtri=3
            integer, parameter :: iCORD=3
            integer, parameter :: iNODE=4
            
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

            ! user variables- For triangular elements
            ! ==================
            
            double precision :: pGPCORDtri(iGPtri,iCORD-1)  ! integration point coordinates
            double precision :: pQUADtri               ! factor for quadrature rule
            double precision :: pWTtri(iGPtri)            ! integration point weights
            double precision :: pNNtri(iGPtri,iNODE-1)      ! shape function matrix (interpolation matrix)
            
            double precision :: amass_row_sum       ! sum of array row
            double precision :: cd,pd,pd_min
            double precision :: Xp(iCORD),Xa(iCORD),Xb(iCORD),Xc(iCORD)
            double precision :: Xba(iCORD),Xca(iCORD),Xpa(iCORD)
            double precision :: pNb(iCORD),pN(iCORD)

            ! include others
            !INCLUDE 'COMMON.f'

            !===============================================================================

            double precision :: xi1,xi2,xi3         ! natural coordinates
            double precision :: xi1tri,xi2tri,xi3tri         ! natural coordinates
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
            double precision :: DetjTri           ! Jacobi-determinant triangular faces (reference)
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
                
            double precision :: AB_edge(iCORD)           ! Edge between nodes 1 and 2
	    double precision :: BC_edge(iCORD)           ! Edge between nodes 2 and 3
	    double precision :: CD_edge(iCORD)           ! Edge between nodes 3 and 4
	    double precision :: AC_edge(iCORD)           ! Edge between nodes 1 and 3
	    double precision :: BD_edge(iCORD)           ! Edge between nodes 2 and 4

            double precision :: N1(iCORD)           ! Normal described by face given by vertex: 1,2,3
	    double precision :: N2(iCORD)           ! Normal described by face given by vertex: 2,3,4
	    double precision :: N3(iCORD)           ! Normal described by face given by vertex: 1,3,4
	    double precision :: N4(iCORD)           ! Normal described by face given by vertex: 1,2,4
	    double precision :: N_all(iCORD)           ! Normals described by different faces

	    double precision :: NODES(iCORD)           ! physical coordinates
            ! integer
            integer :: kblock,ip,nn,ni,nj,i,iptri,Filesize
            ! Allocatable arrays
	    integer, DIMENSION(:), ALLOCATABLE :: FrontEle 
	    integer, DIMENSION(:), ALLOCATABLE :: BackEle            
                        
            if (jtype.eq.34 .and. lflags(iProcedure).eq.jDynExplicit) then 
                
                pEM  = props(1)
                pNU  = props(2)
                pRHO = props(3)
                pGM  = half*pEM/(one+pNU)
                pLAM = (pEM*pNU)/((one+pNU)*(one-two*pNU))

		   
                ! integration point coordinates and weights --------------------------------
                if (iGP==1) then ! HUGHES - The Finite Element Method 1987 (p. 174)

                    pGPCORD(1,:) = (/ one/four, one/four, one/four /)
                        
                    pWT(1) = one
                        
                    pQUAD = one/six
		end if

                !---------------------------------------------------------------------------
                    
                ! shape function for each ip -----------------------------------------------
                do ip=1,iGP
        
                    xi1=pGPCORD(ip,1)
                    xi2=pGPCORD(ip,2)
                    xi3=pGPCORD(ip,3)
        
                    if (iNODE==4) then ! cf. WRIGGERS - Nonlinear Finite Elemente Methods 2008 (p. 120)
        
                        pNN(ip,1) = one-xi1-xi2-xi3
                        pNN(ip,2) = xi1
                        pNN(ip,3) = xi2
                        pNN(ip,4) = xi3
        
                    else
                        stop "Error in computation of shape functions. The number of nodes does not conform with the element type (4 node tetrahedral element)."
                    end if
        
                end do
                !---------------------------------------------------------------------------
                    
                ! identity matrix ----------------------------------------------------------
                pID=zero
                forall(i=1:iCORD) pID(i,i)=one
                !---------------------------------------------------------------------------
                
                if (lflags(iOpCode).eq.jMassCalc) then ! compute mass matrix
                    
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
                            if (iNODE==4) then
            
                                ! derivatives of shape functions with respect to natural coordinates                
                                dNdXi1(1) = -one
                                dNdXi2(1) = -one
                                dNdXi3(1) = -one
                
                                dNdXi1(2) =  one
                                dNdXi2(2) =  zero
                                dNdXi3(2) =  zero
                
                                dNdXi1(3) =  zero
                                dNdXi2(3) =  one
                                dNdXi3(3) =  zero
                
                                dNdXi1(4) =  zero
                                dNdXi2(4) =  zero
                                dNdXi3(4) =  one
                             
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
            
                            ! Jacobian determinant (detJ = 6V)
                            detJ(ip) = dX1dxi1*dX2dxi2*dX3dxi3 + dX2dxi1*dX3dxi2*dX1dxi3 + dX3dxi1*dX1dxi2*dX2dxi3 &
                                     - dX1dxi3*dX2dxi2*dX3dxi1 - dX2dxi3*dX3dxi2*dX1dxi1 - dX3dxi3*dX1dxi2*dX2dxi1
            
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
            
                        end do ! -----------------------------------------------------------
                        
                        !amass(kblock,1:ndofel,1:ndofel)=zero
                        ! loop over all integration points (computation of mass matrix)
                        do ip=1,iGP ! ------------------------------------------------------

                            ! summation over node_i
                            do ni=1,iNODE !-----------------------------loop-i--------------
                
                                ! current node dof
                                do i=1,iCORD
                                    dofni(i) = ni*iCORD-(iCORD-1)+(i-1)
                                end do
                                
                                ! summation over node_i
                                do nj=1,iNODE !-------------------------loop-j--------------
                
                                    ! current node dof
                                    do i=1,iCORD
                                        dofnj(i) = nj*iCORD-(iCORD-1)+(i-1)
                                    end do
                                    
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
                       

                    end do !----------------------------------------------------------------
                        
                elseif (lflags(iOpCode).eq.jIntForceAndDtStable) then !compute internal force + stable time increment
			open(unit=107, file='/home/cerecam/Desktop/MesoporousSilica/Short/BoundaryConditions/nodeSets/NumberPolymerEleSetFront.inp',status='old')!
			READ(107,*) Filesize
			Allocate ( FrontEle(Filesize) )
			close(107)

			open(unit=107, file='/home/cerecam/Desktop/MesoporousSilica/Short/BoundaryConditions/nodeSets/PolymerEleSetFront.csv',status='old')!
			READ(107,*) FrontEle
			close(107)

			open(unit=107, file='/home/cerecam/Desktop/MesoporousSilica/Short/BoundaryConditions/nodeSets/NumberPolymerEleSetBack.inp',status='old')!
			READ(107,*) Filesize
			Allocate ( BackEle(Filesize) )
			close(107)

			open(unit=107, file='/home/cerecam/Desktop/MesoporousSilica/Short/BoundaryConditions/nodeSets/PolymerEleSetBack.csv',status='old')!
			READ(107,*) BackEle
			close(107)

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
                            if (iNODE==4) then
            
                                ! derivatives of shape functions with respect to natural coordinates                
                                dNdXi1(1) = -one
                                dNdXi2(1) = -one
                                dNdXi3(1) = -one
                
                                dNdXi1(2) =  one
                                dNdXi2(2) =  zero
                                dNdXi3(2) =  zero
                
                                dNdXi1(3) =  zero
                                dNdXi2(3) =  one
                                dNdXi3(3) =  zero
                
                                dNdXi1(4) =  zero
                                dNdXi2(4) =  zero
                                dNdXi3(4) =  one
                             
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
            
                            ! Jacobian determinant (detJ = 6V)
                            detJ(ip) = dX1dxi1*dX2dxi2*dX3dxi3 + dX2dxi1*dX3dxi2*dX1dxi3 + dX3dxi1*dX1dxi2*dX2dxi3 &
                                     - dX1dxi3*dX2dxi2*dX3dxi1 - dX2dxi3*dX3dxi2*dX1dxi1 - dX3dxi3*dX1dxi2*dX2dxi1
            
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
            
                        end do ! -----------------------------------------------------------       
            
                        rhs(kblock,1:ndofel)=zero
			energy(kblock,iElIe)=zero
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
		            energy(kblock,iElIe)= energy(kblock,iElIe) + (detJ(ip)/6.0d0)*( ( pGM*ddot(Ee,Ee)+half*pLAM*trace(Ee)*trace(Ee)) )
                            ! von Mises equivalent stress
                            !vMeq = sqrt( half*( (SIG(1,1)- SIG(2,2))**(two) + (SIG(2,2)- SIG(3,3))**(two) + (SIG(3,3)- SIG(1,1))**(two) &
                            !     + six*(SIG(1,2)*SIG(1,2) + SIG(2,3)*SIG(2,3) + SIG(3,1)*SIG(3,1)) ) )
                
                            !VARFLD(VPGP*(ip-1)+1:VPGP*ip,JELEM) = (/ stress(1,1),stress(1,2),stress(1,3), &
                            !                                         stress(2,1),stress(2,2),stress(2,3), &
                            !                                         stress(3,1),stress(3,2),stress(3,3) /)
                            !VARFLD(VPGP*iGP+ip,JELEM) = vMeq
            
                            ! summation over node_i

                            do ni=1,iNODE !-----------------------------loop-i--------------
                		
                                ! current node dof
                                do i=1,iCORD
                                    dofni(i) = ni*iCORD-(iCORD-1)+(i-1)
                                end do
                
                                ! internal residual force vector
                                rhs(kblock,dofni) = rhs(kblock,dofni) + pQUAD*pWT(ip)*detJ(ip)*matvec(P,(/dNdX1(ip,ni),dNdX2(ip,ni),dNdX3(ip,ni)/))
            
                            end do !------------------------------end-loop-i----------------
            
                        end do ! ----------------------------end-loop-ip--------------------

			if (ANY(jElem(kblock).eq.FrontEle)) then
        	        	    pGPCORDtri(1,:) = (/ one/six, one/six /)
	                	    pGPCORDtri(2,:) = (/ two/three, one/six /)
	                	    pGPCORDtri(3,:) = (/ one/six, two/three /)
                	        
	                	    pWTtri = (/ one/six, one/six, one/six/)
        	
	                        do iptri=1,iGPtri
        	
        	            		xi1tri=pGPCORDtri(iptri,1)
        	            		xi2tri=pGPCORDtri(iptri,2)
        		      	
        	      			pNNtri(iptri,1) = xi1tri
        	              		pNNtri(iptri,2) = xi2tri
        	               		pNNtri(iptri,3) = 1.0d0-xi1tri-xi2tri
        	
        	        	end do
        	        	AB_edge = (/ (coords(kblock,1,1)-coords(kblock,2,1)) , (coords(kblock,1,2)-coords(kblock,2,2)), (coords(kblock,1,3)-coords(kblock,2,3)) /)	 
 				BC_edge = (/ (coords(kblock,2,1)-coords(kblock,3,1)) , (coords(kblock,2,2)-coords(kblock,3,2)), (coords(kblock,2,3)-coords(kblock,3,3)) /)
				CD_edge = (/ (coords(kblock,3,1)-coords(kblock,4,1)) , (coords(kblock,3,2)-coords(kblock,4,2)), (coords(kblock,3,3)-coords(kblock,4,3)) /)
				AC_edge = (/ (coords(kblock,1,1)-coords(kblock,3,1)) , (coords(kblock,1,2)-coords(kblock,3,2)), (coords(kblock,1,3)-coords(kblock,3,3)) /)
				BD_edge = (/ (coords(kblock,2,1)-coords(kblock,4,1)) , (coords(kblock,2,2)-coords(kblock,4,2)), (coords(kblock,2,3)-coords(kblock,4,3)) /)
	
				N1 = cross(AB_edge,BC_edge) ! Gives normal of the face defined by nodes (/1(A),2(B),3(C))
				N1 = N1/norm(N1)
				N2 = cross(BD_edge,CD_edge) ! Gives normal of the face defined by nodes (/2(B),3(C),4(D))
				N2 = N2/norm(N2)
				N3 = cross(AC_edge,CD_edge) ! Gives normal of the face defined by nodes (/1(A),3(C),4(D))
				N3 = N3/norm(N3)
				N4 = cross(AB_edge,BD_edge) ! Gives normal of the face defined by nodes (/1(A),2(B),4(D))
				N4 = N4/norm(N4)
	
				if (ABS(N1(3)).eq.1) then
					NODES = (/1, 2, 3/)
					DetjTri  =(coords(kblock,1,1)-coords(kblock,3,1))*(coords(kblock,2,2)-coords(kblock,3,2))& 
						-(coords(kblock,1,2) -coords(kblock,3,2))*(coords(kblock,2,1)-coords(kblock,3,1))
					if (DetjTri<0) then
						NODES = (/ 3,2,1/)
					end if 
!					write(*,*) "element with N1 normal along x-direction: ", jElem(kblock)
					write(*,*) "N1 normal : ", N1
!					write(*,*) "Jacobian : ", DetjTri
				elseif (ABS(N2(3)).eq.1) then
					NODES = (/2, 3, 4/)
					DetjTri  =(coords(kblock,2,1)-coords(kblock,4,1))*(coords(kblock,3,2)-coords(kblock,4,2))& 
						-(coords(kblock,2,2) -coords(kblock,4,2))*(coords(kblock,3,1)-coords(kblock,4,1))
					if (DetjTri<0) then
						NODES = (/ 4,3,2/)
					end if 
!					write(*,*) "element with N2 normal along x-direction: ", jElem(kblock)
					write(*,*) "N2 normal : ", N2
!					write(*,*) "Jacobian : ", DetjTri
				elseif (ABS(N3(3)).eq.1) then
					NODES = (/1, 3, 4/)
					DetjTri  =(coords(kblock,1,1)-coords(kblock,4,1))*(coords(kblock,3,2)-coords(kblock,4,2))& 
						-(coords(kblock,1,2) -coords(kblock,4,2))*(coords(kblock,3,1)-coords(kblock,4,1))
					if (DetjTri<0) then
						NODES = (/ 4,3,1/)
					end if 
!					write(*,*) "element with N3 normal along x-direction: ", jElem(kblock)
					write(*,*) "N3 normal : ", N3
!					write(*,*) "Jacobian : ", DetjTri
				elseif (ABS(N4(3)).eq.1) then
					NODES = (/1, 2, 4/)
					DetjTri  =(coords(kblock,1,1)-coords(kblock,4,1))*(coords(kblock,2,2)-coords(kblock,4,2))& 
						-(coords(kblock,1,2) -coords(kblock,4,2))*(coords(kblock,2,1)-coords(kblock,4,1))
					if (DetjTri<0) then
						NODES = (/ 4,2,1/)
					end if 
!					write(*,*) "element with N4 normal along x-direction: ", jElem(kblock)
					write(*,*) "N4 normal : ", N4
!					write(*,*) "Jacobian : ", DetjTri
				end if
! ------------------------------------ Application of flux vector -----------------------------------------------
				DO iptri=1,iGPtri
					do ni=1,size(NODES)
						nj = NODES(ni)
						do i=1,iCORD
        	                            		dofni(i) = nj*iCORD-(iCORD-1)+(i-1)
		                                end do
						RHS(kblock,dofni) = RHS(kblock,dofni) - pWTtri(iptri)*ABS(DetjTri)*pNNtri(iptri,ni)*(/0.d0,0.d0,0.1d0/)
		
					END DO ! ------------------------ ni-loop ------------------------
				END DO ! ------------------------ iptri-loop ------------------------
			end if ! ------------------------ Boundary jElem-loop ------------------------
	
                        cd = sqrt( (pEM*(one-pNU))/(pRHO*(one+pNU)*(one-two*pNU)) )
                        do i=1,iNODE
                            if (i==1) then
                                Xp = (/coords(kblock,1,1)+u(kblock,1),coords(kblock,1,2)+u(kblock,2),coords(kblock,1,3)+u(kblock,3)/)
                                Xa = (/coords(kblock,2,1)+u(kblock,4),coords(kblock,2,2)+u(kblock,5),coords(kblock,2,3)+u(kblock,6)/)
                                Xb = (/coords(kblock,3,1)+u(kblock,7),coords(kblock,3,2)+u(kblock,8),coords(kblock,3,3)+u(kblock,9)/)
                                Xc = (/coords(kblock,4,1)+u(kblock,10),coords(kblock,4,2)+u(kblock,11),coords(kblock,4,3)+u(kblock,12)/)
                            elseif (i==2) then
                                Xp = (/coords(kblock,2,1)+u(kblock,4),coords(kblock,2,2)+u(kblock,5),coords(kblock,2,3)+u(kblock,6)/)
                                Xa = (/coords(kblock,1,1)+u(kblock,1),coords(kblock,1,2)+u(kblock,2),coords(kblock,1,3)+u(kblock,3)/)
                                Xb = (/coords(kblock,3,1)+u(kblock,7),coords(kblock,3,2)+u(kblock,8),coords(kblock,3,3)+u(kblock,9)/)
                                Xc = (/coords(kblock,4,1)+u(kblock,10),coords(kblock,4,2)+u(kblock,11),coords(kblock,4,3)+u(kblock,12)/)
                            elseif (i==3) then
                                Xp = (/coords(kblock,3,1)+u(kblock,7),coords(kblock,3,2)+u(kblock,8),coords(kblock,3,3)+u(kblock,9)/)
                                Xa = (/coords(kblock,1,1)+u(kblock,1),coords(kblock,1,2)+u(kblock,2),coords(kblock,1,3)+u(kblock,3)/)
                                Xb = (/coords(kblock,2,1)+u(kblock,4),coords(kblock,2,2)+u(kblock,5),coords(kblock,2,3)+u(kblock,6)/)
                                Xc = (/coords(kblock,4,1)+u(kblock,10),coords(kblock,4,2)+u(kblock,11),coords(kblock,4,3)+u(kblock,12)/)
                            elseif (i==4) then
                                Xp = (/coords(kblock,4,1)+u(kblock,10),coords(kblock,4,2)+u(kblock,11),coords(kblock,4,3)+u(kblock,12)/)
                                Xa = (/coords(kblock,1,1)+u(kblock,1),coords(kblock,1,2)+u(kblock,2),coords(kblock,1,3)+u(kblock,3)/)
                                Xb = (/coords(kblock,2,1)+u(kblock,4),coords(kblock,2,2)+u(kblock,5),coords(kblock,2,3)+u(kblock,6)/)
                                Xc = (/coords(kblock,3,1)+u(kblock,7),coords(kblock,3,2)+u(kblock,8),coords(kblock,3,3)+u(kblock,9)/)
                            end if
                            Xba(1) = Xb(1)-Xa(1)
                            Xba(2) = Xb(2)-Xa(2)
                            Xba(3) = Xb(3)-Xa(3)
                            Xca(1) = Xc(1)-Xa(1)
                            Xca(2) = Xc(2)-Xa(2)
                            Xca(3) = Xc(3)-Xa(3)
                            Xpa(1) = Xp(1)-Xa(1)
                            Xpa(2) = Xp(2)-Xa(2)
                            Xpa(3) = Xp(3)-Xa(3)
                            pNb = cross(Xba,Xca)
                            pN = pNb/norm(pNb)
                            pd = abs(dot(Xpa,pN))
                            if (i==1) then
                                pd_min = pd
                            else
                                if ( pd .lt. pd_min ) then
                                    pd_min = pd
                                end if
                            end if
                        end do
                        dtimeStable(kblock) = factorStable*(pd_min/cd)
                        !dtimeStable(kblock) = 2.143102d-06
                        
                        !if (kblock==1) then
                        !    write(*,*)"dtimeStable(kblock)"
                        !    write(*,*)dtimeStable(kblock)
                        !    write(*,*)"rhs(kblock=1,:)"
                        !    write(*,*)rhs(kblock,:)
                        !end if
                        
                    end do !----------------------------------------------------------------
                    
                    
                    
                end if
                    
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
