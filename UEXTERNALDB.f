        !DEC$ FREEFORM
        !===============================================================================
        ! DEVELOPER: EDGAR HUSSER
        ! YEAR: 2017
        !===============================================================================

        !-------------------------------------------------------------------------------
        !USER SUBROUTINE - INITIAL CONDITIONS FOR SDV

        SUBROUTINE UEXTERNALDB(LOP,LRESTART,TIME,DTIME,KSTEP,KINC)

        implicit none
        !-------------------------------------------------------------------------------   
        ! declaration

        ! Abaqus variables
        ! ================================================================
        integer, intent(in ) :: LOP 	! flag parameter: 
                                        ! 0: start of the analysis
                                        ! 1: start of current increment
                                        ! 2: end of current increment
                                        ! 3: end of the analysis
                                        ! 4: beginning of restart analysis

        integer, intent(in ) :: LRESTART    ! 
                                            ! current step time/ total time
        integer, dimension(2), intent(in ) :: TIME 
        integer, intent(in ) :: DTIME       ! time increment
        integer, intent(in ) :: KSTEP       ! current step number
        integer, intent(in ) :: KINC        ! current increment number
        ! ================================================================


        ! global variables from COMMON block
        ! ================================================================
        INCLUDE 'COMMON.f'
        ! ================================================================


        ! internal variables
        ! ================================================================
        integer :: ip,i,j,alpha
        double precision :: xi1,xi2,xi3
        ! ================================================================


        if (LOP==0) then ! read user defined input file
     
            ! integration point coordinates and weights --------------------------------
            if (KGP==1) then ! HUGHES - The Finite Element Method 1987 (p. 174)

                KGPCORD(1,:) = (/ 1.d0/4.d0, 1.d0/4.d0, 1.d0/4.d0 /)
        
                KWT = (/ 1.d0 /)
                
                KQUAD = 1.d0/6.d0
        
            else
        
                stop "Error in computation of integration point coordinates. The number of IPs does not conform with the element type (4 node tetrahedral element with 4 integratin points)."
        
            end if
            !---------------------------------------------------------------------------
        

            ! shape function for each ip -----------------------------------------------
            do ip=1,KGP
        
                xi1=KGPCORD(ip,1)
                xi2=KGPCORD(ip,2)
                xi3=KGPCORD(ip,3)
        
                if (KNODE==4) then ! cf. WRIGGERS - Nonlinear Finite Elemente Methods 2008 (p. 120)
        
                    KNN(ip,1) = 1.d0-xi1-xi2-xi3
                    KNN(ip,2) = xi1
                    KNN(ip,3) = xi2
                    KNN(ip,4) = xi3
        
                else
                    stop "Error in computation of shape functions. The number of nodes does not conform with the element type (4 node tetrahedral element)."
                end if
        
            end do
            !---------------------------------------------------------------------------
    

            ! identity matrix ----------------------------------------------------------
            KID=0.d0
            forall(i=1:KCORD) KID(i,i)=1.d0
            !---------------------------------------------------------------------------
    
            
            ! initialize user output ---------------------------------------------------
            VARFLD = 0.d0
            !---------------------------------------------------------------------------   
    
        !elseif (LOP==1) then

        !elseif (LOP==2) then 
        
        end if

        RETURN
        END
        !-------------------------------------------------------------------------------
        !===============================================================================