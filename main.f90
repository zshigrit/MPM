program main

use mod_type 
use model, only: sMODEL_temporal,read_sPAR,read_sINI,read_sINP 

implicit none 
    
type(model_par) :: sPAR 
type(model_init):: sINI
type(model_inp) :: sINP
type(model_out) :: sOUT

call read_sPAR(sPAR) ! read in parameters from a namelist
call read_sINI(sINI)
call read_sINP(sINP)

call sMODEL_temporal(sPAR,sINI,sINP,sOUT)

! print*, sINP!sINI!sPAR!sOUT%CPOOL
print*, sOUT%CPOOL

end program main 