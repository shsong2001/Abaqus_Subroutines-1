        !DEC$ FREEFORM
        !====================================================================
        ! Include file for material parameter
        !====================================================================
        !
        !--------------------------------------------------------------------
        ! number of dummy elements for user output
        integer, parameter :: NELOUT = 1
        ! number of variables per element for user output: (9+1+1+1)*KGP
        integer, parameter :: NVAROUT = 1 
        ! number of integration points per element
        integer, parameter :: KGP=1
        ! number of cartesian coordinates per node
        integer, parameter :: KCORD=3
        ! number of nodes per element
        integer, parameter :: KNODE=4
        ! number of slip systems (each system covers both dirrections)
        integer, parameter :: KSS = 1
        !
        ! ===================================================================
        ! material properties in kg,s,ym (MPa or N/mm2->N/ym2)
        !
        ! isotropic Young's modulus
        double precision, parameter :: KEM = 79000.d0 !/(1000.d0*1000.d0)   ! 130950 MPa
        ! Poisoon's ratio       
        double precision, parameter :: KNU = 0.44d0
        !
        ! Abaqus Neo-Hooke model
        !! parameter C10
        !double precision, parameter :: KGM = KEM/(4.d0*(1.d0+KNU))
        !! parameter D1
        !double precision, parameter :: KLAM = 6.d0*(1.d0-2.d0*KNU)/KEM
        !
        ! compressible Neo-Hookean
        ! shear modulus (Lame: mu)
        double precision, parameter :: KGM = 0.5d0*KEM/(1.d0+KNU)
        ! Lamé parameter lambda
        double precision, parameter :: KLAM = (KEM*KNU)/((1.d0+KNU)*(1.d0-2.d0*KNU))
        !
        ! --- plasticity parameters -----------------------------------------
        !
        ! saturation rate
        double precision, parameter :: KHSAT = 1000.d0
        ! rate-sensitivity exponent
        double precision, parameter :: KM = 20.d0    
        ! reference slip rate
        double precision, parameter :: KREFS = 0.001d0
        !
        ! ====================================================================
        ! variation value for numerical tangent (global)
        double precision, parameter :: KVNTG = 1.d-08
        ! abort criterium/tolerance for local iteration
        double precision, parameter :: KTOL = 1.d-09
        ! variation value for numerical tangent (local)
        double precision, parameter :: KVNTL = 1.d-07
        ! local iteration limit
        double precision, parameter :: KITER = 200
        !---------------------------------------------------------------------
