!!****if* source/Simulation/SimulationMain/Sod/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer(IN) :: blockID) 
!!                       
!!
!!
!!
!! DESCRIPTION
!!
!!  Initializes fluid data (density, pressure, velocity, etc.) for
!!  a specified block.  This version sets up the Sod shock-tube
!!  problem.
!!
!!  Reference: Sod, G. A., 1978, J. Comp. Phys., 27, 1
!!
!! 
!! ARGUMENTS
!!
!!  blockID -           the number of the block to update
!!
!! PARAMETERS
!!
!!  sim_rhoLeft    Density in the left part of the grid
!!  sim_rhoRight   Density in the right part of the grid
!!  sim_pLeft      Pressure  in the left part of the grid
!!  sim_pRight     Pressure  in the righ part of the grid
!!  sim_uLeft      fluid velocity in the left part of the grid
!!  sim_uRight     fluid velocity in the right part of the grid
!!  sim_xangle     Angle made by diaphragm normal w/x-axis (deg)
!!  sim_ yangle    Angle made by diaphragm normal w/y-axis (deg)
!!  sim_posnR      Point of intersection between the shock plane and
!!                          the x-axis
!!  sim_cposx   Position of the clump in X
!!  sim_cposy   Position of the clump in Y
!!  sim_crad     Radius of the clump
!!  sim_crho     Density of the clump
!!  sim_cp         Pressure of the clump
!!
!!
!!***

subroutine Simulation_initBlock(blockID)

#include "constants.h"
#include "Flash.h"
#include "Eos.h"   

  use Simulation_data, ONLY: sim_posn, sim_xCos, sim_yCos, sim_zCos, &    
     &  sim_rhoLeft,  sim_tLeft, sim_uLeft, sim_rhoRight, sim_tRight, sim_uRight, &
     &  sim_smallX, sim_gamma, sim_smallP, &
     &  sim_cposx, sim_cposy, sim_crad, sim_crho, sim_ct

#ifdef FLASH_3T
  use Simulation_data, ONLY : &
       sim_pionLeft, sim_peleLeft, sim_pradLeft, &
       sim_pionRight, sim_peleRight, sim_pradRight
#endif
     
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
    Grid_getCellCoords, Grid_putPointData
  use Eos_interface, ONLY : Eos_wrapped, Eos

!------------------------------------------------------------------------------
! Ionization section
 use Simulation_speciesData, ONLY :  sim_specElementSymbol,sim_specElement,&
       sim_specSelected, sim_specNumElect, SPEC_NUM, SPEC_NIMAX, sim_specXfrac
  use Ionize_interface, ONLY : Ionize_equil
!------------------------------------------------------------------------------


  implicit none

!------------------------------------------------------------------------------
! Ionization variables

  integer :: id,nel,nion,inn 
  real             zel
  real,dimension(SPEC_NIMAX) :: delem
  real, dimension(NSPECIES) :: massFrac
!------------------------------------------------------------------------------



  ! compute the maximum length of a vector in each coordinate direction 
  ! (including guardcells)
  
  integer, intent(in) :: blockID
  

  integer :: i, j, k, n
  integer :: iMax, jMax, kMax
  


  real :: xx, yy,  zz, xxL, xxR

  real :: new_rho,new_t, distance, press
  
  real :: lPosn0, lPosn
  

  real,allocatable, dimension(:) ::xCenter,xLeft,xRight,yCoord,zCoord

  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer :: sizeX,sizeY,sizeZ
  integer, dimension(MDIM) :: axis

  
  real :: rhoZone, velxZone, velyZone, velzZone, presZone, & 
       eintZone, enerZone, ekinZone, tempZone

#ifdef FLASH_3T
  real :: peleZone, eeleZone
  real :: pionZone, eionZone
  real :: pradZone, eradZone
#endif

!For using Eos 
  real, dimension(EOS_NUM) :: eosData
  integer :: vecLen=1
  

 logical :: gcell = .true.

  
  ! dump some output to stdout listing the paramters
!!$  if (sim_meshMe == MASTER_PE) then
!!$     
!!$     
!!$1    format (1X, 1P, 4(A7, E13.7, :, 1X))
!!$2    format (1X, 1P, 2(A7, E13.7, 1X), A7, I13)
!!$     
!!$  endif
  
  
  ! get the integer index information for the current block
  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)
  
  sizeX = blkLimitsGC(HIGH,IAXIS)
  sizeY = blkLimitsGC(HIGH,JAXIS)
  sizeZ = blkLimitsGC(HIGH,KAXIS)
  allocate(xLeft(sizeX))
  allocate(xRight(sizeX))
  allocate(xCenter(sizeX))
  allocate(yCoord(sizeY))
  allocate(zCoord(sizeZ))
  xCenter = 0.0
  xLeft = 0.0
  xRight = 0.0
  yCoord = 0.0
  zCoord = 0.0

  if (NDIM == 3) call Grid_getCellCoords&
                      (KAXIS, blockId, CENTER,gcell, zCoord, sizeZ)
  if (NDIM >= 2) call Grid_getCellCoords&
                      (JAXIS, blockId, CENTER,gcell, yCoord, sizeY)

  call Grid_getCellCoords(IAXIS, blockId, LEFT_EDGE, gcell, xLeft, sizeX)
  call Grid_getCellCoords(IAXIS, blockId, CENTER, gcell, xCenter, sizeX)
  call Grid_getCellCoords(IAXIS, blockId, RIGHT_EDGE, gcell, xRight, sizeX)

!------------------------------------------------------------------------------

! Loop over cells in the block.  For each, compute the physical position of 
! its left and right edge and its center as well as its physical width.  
! Then decide which side of the initial discontinuity it is on and initialize 
! the hydro variables appropriately.


  do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     
     ! get the coordinates of the cell center in the z-direction
     zz = zCoord(k)
     
     ! Where along the x-axis the shock intersects the xz-plane at the current z.
     lPosn0 = sim_posn - zz*sim_zCos/sim_xCos
     
     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        
        ! get the coordinates of the cell center in the y-direction
        yy = yCoord(j)
        
        ! The position of the shock in the current yz-row.
        lPosn = lPosn0 - yy*sim_yCos/sim_xCos
        
        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
           
           ! get the cell center, left, and right positions in x
           xx  = xCenter(i)
           
           xxL = xLeft(i)
           xxR = xRight(i)

	   distance=sqrt(((xx-sim_cposx)**2) + ((yy-sim_cposy)**2)  )
	   ! if((xx >= sim_cposx-sim_crad).and.(xx<=sim_cposx+sim_crad).and.(yy >= sim_cposy-sim_crad).and.(yy<=sim_cposx+sim_crad)) then
           if (distance <= sim_crad) then 	  
		  tempZone=sim_ct
                  rhoZone = sim_crho
                  velxZone = 0
                  velyZone =0
                  velzZone = 0 

           ! initialize cells to the left of the initial shock.
           elseif (xxR <= lPosn) then
              tempZone = sim_tLeft


              rhoZone = sim_rhoLeft
              velxZone = sim_uLeft * sim_xCos
              velyZone = sim_uLeft * sim_yCos
              velzZone = sim_uLeft * sim_zCos 
              
              ! initialize cells which straddle the shock.  Treat them as though 1/2 of 
              ! the cell lay to the left and 1/2 lay to the right.
           elseif ((xxL < lPosn) .and. (xxR > lPosn)) then              

	      tempZone = 0.5 * (sim_tLeft+sim_tRight)

              
              rhoZone = 0.5 * (sim_rhoLeft+sim_rhoRight) 
              velxZone = 0.5 *(sim_uLeft+sim_uRight) * sim_xCos
              velyZone = 0.5 *(sim_uLeft+sim_uRight) * sim_yCos
              velzZone = 0.5 *(sim_uLeft+sim_uRight) * sim_zCos
              
              ! initialize cells to the right of the initial shock.
           else            
              
	      tempZone = sim_tRight

              rhoZone = sim_rhoRight
              velxZone = sim_uRight * sim_xCos
              velyZone = sim_uRight * sim_yCos
              velzZone = sim_uRight * sim_zCos
              
           endif


           axis(IAXIS) = i
           axis(JAXIS) = j
           axis(KAXIS) = k

           !put in value of default species
           if (NSPECIES > 0) then
              call Grid_putPointData(blockID, CENTER, SPECIES_BEGIN, EXTERIOR, &
                   axis, 1.0e0-(NSPECIES-1)*sim_smallX)


              !if there is only 1 species, this loop will not execute
              do n = SPECIES_BEGIN+1,SPECIES_END
                 call Grid_putPointData(blockID, CENTER, n, EXTERIOR, &
                      axis, sim_smallX)
              enddo
           end if

           ! Compute the gas energy and set the gamma-values needed for the equation of 
           ! state.
           ekinZone = 0.5 * (velxZone**2 + & 
                velyZone**2 + & 
                velzZone**2)
           
	
           !----------------------------------------------------------
	   !Assign equilibrium mass fractions

           !Copy in Electron & Hydrogen
           massFrac(1:2) = sim_specXfrac(1:2)

           id = 2
           do nel = 1,SPEC_NUM
              if (sim_specSelected(nel)) then
                 call Ionize_equil(tempZone,nel,nion,delem)
                 nion = sim_specNumElect(nel)+1
                 do inn = 1,nion
                    id = id+1
                    massFrac(id) = delem(inn)*sim_specXfrac(id)
                 enddo
              endif
           enddo

           !----------------------------------------------------------


           !----------------------------------------------------------
           ! get the pressure and internal energy corresponding 
           ! to the density and  temperature
           !-----------------------------------------------------------

           eosData(:) = 0.0
           eosData(EOS_DENS)=rhoZone
           eosData(EOS_TEMP)=tempZone
	   !Multigamma Eos requires mass fractions.
           call Eos(MODE_DENS_TEMP,vecLen,eosData, massFrac)

           presZone = eosData(EOS_PRES)

           eintZone = presZone / (sim_gamma-1.)
           eintZone = eintZone / rhoZone
           enerZone = eintZone + ekinZone
           enerZone = max(enerZone, sim_smallP)
           
           ! store the variables in the current zone via Grid put methods
           ! data is put stored one cell at a time with these calls to Grid_putData           


           call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, axis, rhoZone)
           call Grid_putPointData(blockId, CENTER, PRES_VAR, EXTERIOR, axis, presZone)
           call Grid_putPointData(blockId, CENTER, VELX_VAR, EXTERIOR, axis, velxZone)
           call Grid_putPointData(blockId, CENTER, VELY_VAR, EXTERIOR, axis, velyZone)
           call Grid_putPointData(blockId, CENTER, VELZ_VAR, EXTERIOR, axis, velzZone)
	   call Grid_putPointData(blockID, CENTER, SPECIES_BEGIN+n-1, EXTERIOR, &
                   axis, massFrac(n))
	   
#ifdef ENER_VAR
           call Grid_putPointData(blockId, CENTER, ENER_VAR, EXTERIOR, axis, enerZone)   
#endif
#ifdef EINT_VAR
           call Grid_putPointData(blockId, CENTER, EINT_VAR, EXTERIOR, axis, eintZone)   
#endif
#ifdef GAME_VAR          
           call Grid_putPointData(blockId, CENTER, GAME_VAR, EXTERIOR, axis, sim_gamma)
#endif
#ifdef GAMC_VAR
           call Grid_putPointData(blockId, CENTER, GAMC_VAR, EXTERIOR, axis, sim_gamma)
#endif
#ifdef TEMP_VAR
           call Grid_putPointData(blockId, CENTER, TEMP_VAR, EXTERIOR, axis, tempZone)
#endif

#ifdef FLASH_3T
           ! We must now compute the internal energy from the pressure
           ! for the ions, electrons, and radiation field:
           
           ! Electrons...
           eeleZone = peleZone / (sim_gamma - 1.0) / rhoZone
           eionZone = pionZone / (sim_gamma - 1.0) / rhoZone
           eradZone = 3.0 * pradZone / rhoZone
           
           call Grid_putPointData(blockId, CENTER, EELE_VAR, EXTERIOR, axis, eeleZone)
           call Grid_putPointData(blockId, CENTER, EION_VAR, EXTERIOR, axis, eionZone)
           call Grid_putPointData(blockId, CENTER, ERAD_VAR, EXTERIOR, axis, eradZone)
#endif
        enddo
     enddo
  enddo

! #ifdef EELE_VAR
!   call Eos_wrapped(MODE_DENS_EI_SCATTER,blkLimits,blockId)
! #endif

!   do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
!      do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
!         do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
!            axis(IAXIS) = i
!            axis(JAXIS) = j
!            axis(KAXIS) = k
! #ifdef ERAD_VAR
!            call Grid_putPointData(blockId, CENTER, ERAD_VAR, EXTERIOR, axis, 0.0  )   
! #endif
! #ifdef E3_VAR
!            call Grid_putPointData(blockId, CENTER, E3_VAR,   EXTERIOR, axis, 0.0  )   
! #endif

! #ifdef PRAD_VAR
!            call Grid_putPointData(blockId, CENTER, PRAD_VAR, EXTERIOR, axis, 0.0  )   
! #endif
! #ifdef TRAD_VAR
!            call Grid_putPointData(blockId, CENTER, TRAD_VAR, EXTERIOR, axis, 0.0  )   
! #endif
!         enddo
!      enddo
!   enddo

!! Cleanup!  Must deallocate arrays

  deallocate(xLeft)
  deallocate(xRight)
  deallocate(xCenter)
  deallocate(yCoord)
  deallocate(zCoord)

 
  return
end subroutine Simulation_initBlock
