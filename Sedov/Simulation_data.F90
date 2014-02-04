!!****if* source/Simulation/SimulationMain/Sedov/Simulation_data
!!
!! NAME
!!
!!  Simulation_data
!!
!! SYNOPSIS
!!
!!  use Simulation_data 
!!
!!  DESCRIPTION
!!
!!  Stores the local data for Simulation setup: Sedov
!!  
!! PARAMETERS
!!
!!  sim_pAmbient       Initial ambient pressure
!!  sim_rhoAmbient     Initial ambient density
!!  sim_expEnergy      Explosion energy (distributed over 2^dimen central zones)
!!  sim_rInit          Radial position of inner edge of grid (for 1D )
!!  sim_xctr           Explosion center coordinates
!!  sim_yctr           Explosion center coordinates
!!  sim_zctr           Explosion center coordinates
!!  sim_nsubzones      Number of `sub-zones' in cells for applying 1d profile
!!
!!
!!***

module Simulation_data

  implicit none
#include "constants.h"

  !! *** Runtime Parameters *** !!

  real, save    :: sim_pAmbient, sim_rhoAmbient, sim_expEnergy, sim_rInit
  real, save    :: sim_gamma, sim_xCenter, sim_yCenter, sim_zCenter
  real, save    :: sim_smallX, sim_smallRho, sim_smallP, sim_pi
  real, save    :: sim_xMin, sim_xMax, sim_yMin, sim_yMax, sim_zMin, sim_zMax
  integer, save :: sim_nSubZones
  integer, save :: sim_iFuel
  real, save    :: sim_tInitial

  !! *** Variables pertaining to this Simulation *** !!

  integer, parameter                  :: sim_nProfile = 1000
  real   , save                       :: sim_inSubZones, sim_inSubzm1
  real   , save                       :: sim_inszd
  real, dimension(sim_nProfile), save :: sim_rProf, sim_rhoProf, sim_pProf
  real, dimension(sim_nProfile), save :: sim_vProf
  real, save                          :: sim_drProf, sim_pExp, sim_vctr

integer, save :: sim_meshMe
end module Simulation_data
