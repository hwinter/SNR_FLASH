!!****if* source/Simulation/SimulationMain/Sedov/WriteParticleSubset/io_ptCreateSubset
!!
!! NAME
!!
!! io_ptCreateSubset
!!
!! SYNOPSIS
!!
!! io_ptCreateSubset(integer(IN)  :: numProcs,
!!                   integer(IN)  :: subsetIndex,
!!                   integer(IN)  :: numProperties,
!!                   integer(IN)  :: numParticles,
!!                   char(IN)     :: partLabels(OUTPUT_PROP_LENGTH,numParticles),
!!                   real(INOUT)  :: particlest(numProperties,numParticles),
!!                   real(OUT)    :: subsetSize(2),
!!                   char(OUT)    :: subsetName(MAX_STRING_LENGTH),
!!                   char(OUT)    :: subsetLabelName(MAX_STRING_LENGTH),
!!                   char(OUT)    :: subsetLabels(OUTPUT_PROP_LENGTH,numParticles),
!!                   logical(OUT) :: moreSubsetsRemain
!!
!! ARGUMENTS
!!
!!  numProcs          - Global number of processes.
!!  subsetIndex       - Current subset iteration.  This subroutine will
!!                      continually be called with incremented values of 
!!                      subsetIndex until the user sets moreSubsetsRemain to 
!!                      .false..
!!  numProperties     - Number of particle properties.
!!  numParticles      - Number of particles on myPE (local process ID).
!!  partLabels        - Text labels of particle properties.
!!  particlest        - Temporary particles array that is used to hold a
!!                      user-defined subset.  We reuse the same space for
!!                      different values of subsetIndex.
!!  subsetSize        - The size of user-defined subset.
!!  subsetName        - The particle dataset name in file.
!!  subsetLabelName   - The particle labels dataset name in file.
!!  subsetLabels      - The text labels of the particle properties in the
!!                      user-defined subset.
!!  moreSubsetsRemain - Whether there are more subsets after subset
!!                      "subsetIndex".
!! 
!! DESCRIPTION
!! 
!!  Creates custom particle subset(s) that are written to a particle file.
!!  These subsets do not get written to a checkpoint file.
!!
!! NOTES
!!  
!!  Visit will only read the particle dataset name "tracer particles"
!!  and the particle labels name "particle names".  You can replace these
!!  datasets with your custom data by setting the runtime parameters
!!  writeParticleSubset = .true. and writeParticleAll = .false..  Then
!!  return a subset from io_ptCreateSubset with names: 
!!  subsetName = "tracer particles" and subsetLabelName = "particle names"
!!
!!  The following example shows how to create a subset that replaces the
!!  standard particle datasets (and so can be visualized with Visit).  
!!  In this case the user particles array (particlest) is identical
!!  to the FLASH particles array (particles).  Customize as needed, e.g. 
!!  only include particles that meet a certain criterion.
!!
!!  if (subsetIndex == 0) then
!!     moreSubsetsRemain = .false.
!!     subsetName = "tracer particles"
!!     subsetLabelName = "particle names"
!!     subsetLabels(:) = partLabels(:)
!!     subsetSize(1) = numProperties
!!     subsetSize(2) = numParticles
!!     particlest(1:numProperties,1:numParticles) = &
!!          particles(1:numProperties,1:numParticles)
!!  end if
!!
!!***

subroutine io_ptCreateSubset ( subsetIndex, numProperties, &
     numParticles, partLabels, particlest, &
     subsetSize, subsetName, subsetLabelName, subsetLabels, &
     moreSubsetsRemain)
    
  use Particles_data, ONLY : particles
  implicit none

#include "Flash.h"
#include "constants.h"
    
  integer, intent(IN) ::  subsetIndex, numProperties, &
       numParticles
  character(len=OUTPUT_PROP_LENGTH), dimension(numProperties), &
       intent(IN) :: partLabels

  real, dimension(numProperties,numParticles), intent(INOUT) :: particlest

  integer, dimension(2), intent(OUT) :: subsetSize
  character(len=MAX_STRING_LENGTH), intent(OUT) :: subsetName, subsetLabelName
  character(len=OUTPUT_PROP_LENGTH), dimension(numProperties), &
       intent(OUT) :: subsetLabels
  logical, intent(OUT) :: moreSubsetsRemain


  integer :: i, matchingParticles


  !I will make 2 custom subsets. Check particle files using:  
  !
  !h5dump -d Tags_of_all_particles_labels \
  !-d Tags_of_all_particles \
  !-d Posx_greater_than_0.5_labels \
  !-d Posx_greater_than_0.5 \
  !sedov_2d_6lev_hdf5_part_0000

  if (subsetIndex == 0) then
     moreSubsetsRemain = .true.
     subsetName = "Tags_of_all_particles"
     subsetLabelName = "Tags_of_all_particles_labels"
     subsetLabels(1) = "tag"

     subsetSize(1) = 1  !We only want a single property.
     subsetSize(2) = numParticles
     do i = 1, numParticles
        particlest(1,i) = particles(TAG_PART_PROP,i)
     end do
  end if

  if (subsetIndex == 1) then
     moreSubsetsRemain = .false.
     subsetName = "Posx_greater_than_0.5"
     subsetLabelName = "Posx_greater_than_0.5_labels"
     subsetLabels(:) = partLabels(:)

     subsetSize(1) = numProperties  !We want all properties.
     matchingParticles = 0
     do i = 1, numParticles
        if (particles(POSX_PART_PROP,i) >= 0.5) then
           matchingParticles = matchingParticles + 1
           particlest(:,matchingParticles) = particles(:,i)
        end if
     end do
     subsetSize(2) = matchingParticles
  end if

end subroutine io_ptCreateSubset
