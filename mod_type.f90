module mod_type 
    implicit none 

    type model_par 
        real(8):: Vd ! uptake rate of POM
        real(8):: KPOM ! half saturation constant for POM
        real(8):: rMORT ! microbial mortality rate
        real(8):: Vsorp ! sorption rate of MOM on dead microbes and POM
        real(8):: rSorp ! potential sorption rate of POM to MOM
        real(8):: Ksorp ! half saturation constant for sorption
        real(8):: rDesorp ! desorption rate of MOM
        real(8):: NUEx, CUEx ! 
        real(8):: KNin ! half saturation constant for mineral nitrogen
        real(8):: CNinp
    end type model_par 

    type model_flux
        real(8):: POM_dec 
        real(8):: MB_mortality
        real(8):: TOT_sorp ! total sorption of MB and POM
        real(8):: MB_sorp
        real(8):: POM_sorp
        real(8):: MOM_desorp
    end type model_flux

    type model_flux_mn
        real(8):: Nmn, Nim 
        real(8):: CO2_overflow
    end type model_flux_mn

    type model_pool
        real(8):: MB 
        real(8):: POM
        real(8):: MOM
    end type model_pool

    type model_pool_mn
        real(8):: N 
    end type model_pool_mn

    type model_init
        type(model_pool):: CPOOL
        type(model_pool):: NPOOL
        type(model_pool):: CN
        type(model_pool_mn):: MNPOOL
    end type model_init

    type model_inp
        real(8):: carbon 
        real(8):: ST 
        real(8):: SW
    end type model_inp

    type model_out
        type(model_flux):: CFLUX
        type(model_flux):: NFLUX
        type(model_pool):: CPOOL
        type(model_pool):: NPOOL
        type(model_pool):: CN
        type(model_pool_mn):: MNPOOL
        type(model_flux_mn):: MNFLUX 
    end type model_out
end module mod_type 