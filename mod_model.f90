module model 
    use mod_type 
    implicit none 
    private
    public:: sMODEL
    public:: sMODEL_temporal
    public:: read_sPAR
    public:: read_sINI
    public:: read_sINP

contains 

! read parameters from namelist 'para_namelist.nml'
subroutine read_sPAR(PAR)

    type(model_par):: PAR
    real(8):: Vd,KPOM,rMORT,Vsorp,rSorp,Ksorp,rDesorp
    real(8):: NUEx,CUEx,KNin,CNinp
    
    namelist /main_para/ Vd,KPOM,rMORT,Vsorp,rSorp,Ksorp,rDesorp, &
    &                   NUEx,CUEx,KNin,CNinp
    open(11,file='para_namelist.nml')
    read(11,nml=main_para)
    close(11)

    PAR%Vd      = Vd
    PAR%KPOM    = KPOM
    PAR%rMORT   = rMORT
    PAR%Vsorp   = Vsorp
    PAR%rSorp   = rSorp
    PAR%Ksorp   = Ksorp
    PAR%rDesorp = rDesorp
    PAR%NUEx    = NUEx
    PAR%CUEx    = CUEx
    PAR%KNin    = KNin
    PAR%CNinp   = CNinp
end subroutine read_sPAR

! read in initial conditions for pools and CN
subroutine read_sINI(sINI)
    type(model_init):: sINI
    sINI%CPOOL%MB=6.0
    sINI%CPOOL%POM=30.0
    sINI%CPOOL%MOM=40.0

    sINI%NPOOL%MB=1.0
    sINI%NPOOL%POM=1.0
    sINI%NPOOL%MOM=5.0

    sINI%CN%MB=6.0
    sINI%CN%POM=30.0 ! SAME AS CARBON INPUT CN RATIO
    sINI%CN%MOM=8.0

    sINI%MNPOOL%N=5.0
end subroutine read_sINI

! read in input
subroutine read_sINP(sINP)
    type(model_inp):: sINP
    sINP%carbon = 30.0d0/24.0d0 !mg C per hour
end subroutine read_sINP 

! the core model
subroutine sMODEL_temporal(PAR,INI,INP,OUT)
    type(model_par) :: PAR
    type(model_init):: INI 
    type(model_inp) :: INP ! forcing and carbon inputs
    type(model_out) :: OUT

! local variables
    real(8):: rTER  
    integer:: iday,ihour   

  associate(                &
      CFLUX => OUT%CFLUX,   &
      NFLUX => OUT%NFLUX,   &
      CPOOL => OUT%CPOOL,   &
      NPOOL => OUT%NPOOL,   &
      CN    => OUT%CN   ,   &
      MNPOOL => OUT%MNPOOL, &
      MNFLUX => OUT%MNFLUX, &
      rINP => INP%carbon    & 
  )

CPOOL = INI%CPOOL
NPOOL = INI%NPOOL
CN    = INI%CN
MNPOOL= INI%MNPOOL

if (CN%POM/=PAR%CNinp) then
    print*, 'make sure input CN equals POM cn'
else 
    print*, 'inputCN=POMCN'
end if 

!=====================================================!
! temperal running of model                           !
!=====================================================!
do iday=1,20
    do ihour=1,24
    ! ===================================================== !
    ! notes: 1. consider adding negative density dependence !
    ! ===================================================== !

    ! POM uptake
        CFLUX%POM_dec = PAR%Vd*CPOOL%MB*CPOOL%POM/(PAR%KPOM+CPOOL%POM)
        NFLUX%POM_dec = CFLUX%POM_dec/CN%POM

    ! microbial mortality
        CFLUX%MB_mortality = PAR%rMORT * CPOOL%MB !**1.5 (density dependence)
        NFLUX%MB_mortality = CFLUX%MB_mortality/CN%MB

    ! sorption of MB and POM by MOM
        CFLUX%TOT_sorp = PAR%Vsorp*(PAR%rSorp*CPOOL%POM+CFLUX%MB_mortality)/ &
        &               (PAR%Ksorp+(PAR%rSorp*CPOOL%POM+CFLUX%MB_mortality))
        CFLUX%MB_sorp = (CFLUX%MB_mortality/(CFLUX%MB_mortality+PAR%rSorp*CPOOL%POM)) &
        &               *CFLUX%TOT_sorp
        CFLUX%POM_sorp = (PAR%rSorp*CPOOL%POM/(CFLUX%MB_mortality+PAR%rSorp*CPOOL%POM)) &
        &               *CFLUX%TOT_sorp
        NFLUX%MB_sorp = CFLUX%MB_sorp/CN%MB
        NFLUX%POM_sorp = CFLUX%POM_sorp/CN%POM

    ! desorption of MOM
        CFLUX%MOM_desorp = PAR%rDesorp*CPOOL%MOM
        NFLUX%MOM_desorp = CFLUX%MOM_desorp/CN%MOM

    ! ===================================================== !
    ! notes: 2. N mineralization & immobilization: beg      !
    ! ===================================================== !
        rTER = CN%MB*(PAR%NUEx/PAR%CUEx) ! threshold element ratio

        if ((CFLUX%POM_dec*PAR%CUEx)/(NFLUX%POM_dec*PAR%NUEx).lt.rTER) then
            MNFLUX%Nmn = (NFLUX%POM_dec*PAR%NUEx-CFLUX%POM_dec*PAR%CUEx/CN%MB) &
            &           +NFLUX%POM_dec*(1-PAR%NUEx)
            MNFLUX%Nim = 0.d0
            MNFLUX%CO2_overflow = 0.d0
        else
            MNFLUX%Nmn = 0.d0
            MNFLUX%Nim = MNPOOL%N/(PAR%KNin+MNPOOL%N)
            MNFLUX%CO2_overflow = CFLUX%POM_dec*PAR%CUEx - CN%MB*(NFLUX%POM_dec*PAR%NUEx+MNFLUX%Nim)
            if (MNFLUX%Nim.gt.CFLUX%POM_dec*PAR%CUEx/CN%MB-NFLUX%POM_dec*PAR%NUEx) then
                MNFLUX%Nim = CFLUX%POM_dec*PAR%CUEx/CN%MB-NFLUX%POM_dec*PAR%NUEx
            end if 
        end if
    ! ===================================================== !
    ! notes: 2. N mineralization & immobilization: end      !
    ! ===================================================== !

    ! carbon balance
        CPOOL%POM = CPOOL%POM+rINP+(CFLUX%MB_mortality-CFLUX%MB_sorp) &
        &           +CFLUX%MOM_desorp-CFLUX%POM_dec-CFLUX%POM_sorp
        CPOOL%MOM = CPOOL%MOM+CFLUX%TOT_sorp-CFLUX%MOM_desorp
        CPOOL%MB = CPOOL%MB+CFLUX%POM_dec*PAR%CUEx-CFLUX%MB_mortality-MNFLUX%CO2_overflow

    ! nitrogen balance
        NPOOL%POM = NPOOL%POM+rINP/PAR%CNinp+(NFLUX%MB_mortality-NFLUX%MB_sorp) &
        &           +NFLUX%MOM_desorp-NFLUX%POM_dec-NFLUX%POM_sorp
        NPOOL%MOM = NPOOL%MOM+NFLUX%TOT_sorp-NFLUX%MOM_desorp
        NPOOL%MB  = NPOOL%MB+NFLUX%POM_dec*PAR%NUEx+MNFLUX%Nim-MNFLUX%Nmn-NFLUX%MB_mortality
        MNPOOL%N  = MNPOOL%N+MNFLUX%Nmn-MNFLUX%Nim

    ! update CN
        CN%MB = CPOOL%MB/NPOOL%MB
        CN%POM = CPOOL%POM/NPOOL%POM
        CN%MOM = CPOOL%MOM/NPOOL%MOM

    end do !hour
end do !day

end associate 

end subroutine sMODEL_temporal 

end module model 