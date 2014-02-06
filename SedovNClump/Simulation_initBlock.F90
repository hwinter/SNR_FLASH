!!****if* source/Simulation/SimulationMain/Sedov/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!! 
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer :: blockId, 
!!                       
!!
!!
!! DESCRIPTION
!!
!!  Initializes fluid data (density, pressure, velocity, etc.) for
!!  a specified block.  This version sets up the Sedov spherical
!!  explosion problem.
!!
!!  References:  Sedov, L. I., 1959, Similarity and Dimensional Methods
!!                 in Mechanics (New York:  Academic)
!!
!!               Landau, L. D., & Lifshitz, E. M., 1987, Fluid Mechanics,
!!                 2d ed. (Oxford:  Pergamon)
!!
!! ARGUMENTS
!!
!!  blockId -        The number of the block to initialize
!!  
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

!!REORDER(4): solnData


subroutine Simulation_initBlock(blockId)

  use Simulation_data, ONLY: sim_xMax, sim_xMin, sim_yMax, sim_yMin, sim_zMax, sim_zMin, &
     &  sim_nProfile, sim_drProf, sim_rProf, sim_vProf, sim_pProf, sim_pExp, sim_rhoProf, &
     &  sim_tInitial, sim_gamma, sim_expEnergy, sim_pAmbient, sim_rhoAmbient, &
     &  sim_smallX, sim_smallRho, sim_smallP, sim_rInit, &
     &  sim_nSubZones, sim_xCenter, sim_yCenter, sim_zCenter, sim_inSubzm1, sim_inszd
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_releaseBlkPtr,&
    Grid_getCellCoords, Grid_putPointData
  
  implicit none

#include "constants.h"
#include "Flash.h"
  
  integer,intent(IN) ::  blockId
  
  
  integer  ::  i, j, k, n, jLo, jHi
  integer  ::  ii, jj, kk
  real     ::  distInv, xDist, yDist, zDist
  real     ::  sumRho, sumP, sumVX, sumVY, sumVZ
  real     ::  vel, diagonal
  real     ::  xx, dxx, yy, dyy, zz, dzz, frac
  real     ::  vx, vy, vz, p, rho, e, ek
  real     ::  dist
  logical  ::  validGeom
  integer :: istat

  real,allocatable,dimension(:) :: xCoord,yCoord,zCoord
  integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC
  integer :: sizeX,sizeY,sizeZ
  integer,dimension(MDIM) :: axis
  real, dimension(:,:,:,:),pointer :: solnData

  logical :: gcell = .true.

  !
  !  Construct the radial samples needed for the initialization.
  !
  diagonal = (sim_xMax-sim_xMin)**2
  diagonal = diagonal + K2D*(sim_yMax-sim_yMin)**2
  diagonal = diagonal + K3D*(sim_zMax-sim_zMin)**2
  diagonal = sqrt(diagonal)
  
  sim_drProf = diagonal / (sim_nProfile-1)
  
  do i = 1, sim_nProfile
     sim_rProf(i)   = (i-1) * sim_drProf
  enddo
  !
  !  If t>0, use the analytic Sedov solution to initialize the
  !  code.  Otherwise, just use a top-hat.
  !

  if (sim_tInitial .gt. 0.) then
     call set_analytic_sedov (sim_nProfile, sim_rProf, sim_rhoProf, sim_pProf, & 
          sim_vProf, sim_tInitial, sim_gamma, sim_expEnergy, & 
          sim_pAmbient, sim_rhoAmbient)
  else
     do i = 1, sim_nProfile
        sim_rhoProf(i) = sim_rhoAmbient
        sim_pProf(i)   = sim_pAmbient
        sim_vProf(i)   = 0.
        if (sim_rProf(i) .le. sim_rInit) sim_pProf(i) = sim_pExp
     enddo
     
  endif

  ! get the coordinate information for the current block from the database

  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)
  sizeX = blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1
  allocate(xCoord(sizeX),stat=istat); xCoord = 0.0
  sizeY = blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS) + 1
  allocate(yCoord(sizeY),stat=istat); yCoord = 0.0
  sizeZ = blkLimitsGC(HIGH,KAXIS) - blkLimitsGC(LOW,KAXIS) + 1
  allocate(zCoord(sizeZ),stat=istat); zCoord = 0.0

  if (NDIM == 3) call Grid_getCellCoords&
                      (KAXIS, blockId, CENTER, gcell, zCoord, sizeZ)
  if (NDIM >= 2) call Grid_getCellCoords&
                      (JAXIS, blockId, CENTER,gcell, yCoord, sizeY)
  call Grid_getCellCoords(IAXIS, blockId, CENTER, gcell, xCoord, sizeX)
  !
  !     For each cell
  !  
#ifdef FL_NON_PERMANENT_GUARDCELLS
  call Grid_getBlkPtr(blockId,solnData)
#endif
  do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
     ! Find a real difference between z's if problem is >= 3D
     if (NDIM > 2) then
        if (k .eq. 1) then
           dzz = zCoord(2) - zCoord(1) 
        else
           dzz = zCoord(k) - zCoord(k-1) 
        endif
     ! Otherwise this problem is <= 2D, so dzz is meaningless
     else
       dzz = 0.0
     endif
     zz = zCoord(k)
     
     do j = blkLimitsGC(LOW, JAXIS), blkLimitsGC(HIGH, JAXIS)
        ! Find a real difference between y's if problem is >= 2D
        if (NDIM > 1) then
           if (j .eq. 1) then
              dyy = yCoord(2) - yCoord(1) 
           else
              dyy = yCoord(j) - yCoord(j-1) 
           endif
        ! Otherwise this problem is <= 1D, so dyy is meaningless
        else
          dyy = 0.0
        endif
        yy = yCoord(j)
        

        do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH, IAXIS)
           xx = xCoord(i)
           if (i .eq. 1) then
              dxx = xCoord(2) - xCoord(1) 
           else
              dxx = xCoord(i) - xCoord(i-1) 
           endif
           
           sumRho = 0.
           sumP   = 0.
           sumVX  = 0.
           sumVY  = 0.
           sumVZ  = 0.
           
           !
           !       Break the cell into sim_nSubZones^NDIM sub-zones, and look up the
           !       appropriate quantities along the 1d profile for that subzone.  
           !
           !       Have the final values for the zone be equal to the average of
           !       the subzone values.
           ! 

           do kk = 0, (sim_nSubZones-1)*K3D
              zz    = zCoord(k) + (kk*sim_inSubzm1-.5)*dzz 
              zDist = (zz - sim_zCenter) * K3D
              
              do jj = 0, (sim_nSubZones-1)*K2D
                 yy    = yCoord(j) + (jj*sim_inSubzm1-.5)*dyy
                 yDist = (yy - sim_yCenter) * K2D
                 
                 do ii = 0, (sim_nSubZones-1)
                    xx    = xCoord(i) + (ii*sim_inSubzm1-.5)*dxx
                    xDist = xx - sim_xCenter
                    
                    dist    = sqrt( xDist**2 + yDist**2 + zDist**2 )
                    distInv = 1. / max( dist, 1.E-10 )
                    call sim_find (sim_rProf, sim_nProfile, dist, jLo)
                    !
                    !  a point at `dist' is frac-way between jLo and jHi.   We do a
                    !  linear interpolation of the quantities at jLo and jHi and sum those.
                    ! 
                    if (jLo .eq. 0) then
                       jLo = 1
                       jHi = 1
                       frac = 0.
                    else if (jLo .eq. sim_nProfile) then
                       jLo = sim_nProfile
                       jHi = sim_nProfile
                       frac = 0.
                    else
                       jHi = jLo + 1
                       frac = (dist - sim_rProf(jLo)) / & 
                            (sim_rProf(jHi)-sim_rProf(jLo))
                    endif
                    ! 
                    !   Now total these quantities.   Note that  v is a radial velocity; 
                    !   we multiply by the tangents of the appropriate angles to get
                    !   the projections in the x, y and z directions.
                    !
                    sumP = sumP +  & 
                         sim_pProf(jLo) + frac*(sim_pProf(jHi)  - sim_pProf(jLo))
                    
                    sumRho = sumRho + & 
                         sim_rhoProf(jLo) + frac*(sim_rhoProf(jHi)- sim_rhoProf(jLo))
                    
                    vel = sim_vProf(jLo) + frac*(sim_vProf(jHi)  - sim_vProf(jLo))
                    
                    sumVX  = sumVX  + vel*xDist*distInv
                    sumVY  = sumVY  + vel*yDist*distInv
                    sumVZ  = sumVZ  + vel*zDist*distInv
                    
                 enddo
              enddo
           enddo
           
           rho = max(sumRho * sim_inszd, sim_smallRho)
           p   = max(sumP   * sim_inszd, sim_smallP)
           vx  = sumVX  * sim_inszd
           vy  = sumVY  * sim_inszd
           vz  = sumVZ  * sim_inszd
           ek  = 0.5*(vx*vx + vy*vy + vz*vz)
           !
           !  assume gamma-law equation of state
           !
           e   = p/(sim_gamma-1.)
           e   = e/rho + ek
           e   = max (e, sim_smallP)
           
           axis(IAXIS)=i
           axis(JAXIS)=j
           axis(KAXIS)=k


#ifdef FL_NON_PERMANENT_GUARDCELLS
           if (NSPECIES > 0) then
              solnData(SPECIES_BEGIN,i,j,k)=1.0-(NSPECIES-1)*sim_smallX
              solnData(SPECIES_BEGIN+1:SPECIES_END,i,j,k)=sim_smallX
           end if
           solnData(DENS_VAR,i,j,k)=rho
           solnData(PRES_VAR,i,j,k)=p
           solnData(ENER_VAR,i,j,k)=e
           solnData(GAME_VAR,i,j,k)=sim_gamma
           solnData(GAMC_VAR,i,j,k)=sim_gamma
           solnData(VELX_VAR,i,j,k)=vx
           solnData(VELY_VAR,i,j,k)=vy
           solnData(VELZ_VAR,i,j,k)=vz
#else
           if (NSPECIES > 0) then
              ! putting in the value of the default species
              call Grid_putPointData(blockID, CENTER, SPECIES_BEGIN, EXTERIOR, &
                   axis, 1.0e0-(NSPECIES-1)*sim_smallX)


              !if there is only one species, this loop will not execute
              do n = SPECIES_BEGIN+1, SPECIES_END

                 call Grid_putPointData(blockID, CENTER, n, EXTERIOR, &
                      axis, sim_smallX)

              enddo
           end if


           call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, axis, rho)
           call Grid_putPointData(blockId, CENTER, PRES_VAR, EXTERIOR, axis, p)
           call Grid_putPointData(blockId, CENTER, ENER_VAR, EXTERIOR, axis, e)    
           call Grid_putPointData(blockId, CENTER, GAME_VAR, EXTERIOR, axis, sim_gamma)
           call Grid_putPointData(blockId, CENTER, GAMC_VAR, EXTERIOR, axis, sim_gamma)
           call Grid_putPointData(blockId, CENTER, VELX_VAR, EXTERIOR, axis, vx)
           call Grid_putPointData(blockId, CENTER, VELY_VAR, EXTERIOR, axis, vy)
           call Grid_putPointData(blockId, CENTER, VELZ_VAR, EXTERIOR, axis, vz)
#endif
        enddo
     enddo
  enddo
#ifdef FL_NON_PERMANENT_GUARDCELLS
  call Grid_releaseBlkPtr(blockID, solnData)
#endif
  deallocate(xCoord)
  deallocate(yCoord)
  deallocate(zCoord)
  return
end subroutine Simulation_initBlock




!  Subroutine:  set_analytic_sedov()

!  Description: Given a set of arrays to store radius, density, pressure, and
!               velocity, together with a time and a ratio of specific heats,
!               generate the analytical solution to the Sedov problem.

!               Currently this routine is limited to gamma=1.4 (I've hardwired
!               the dimensionless constant beta) and will complain if it gets
!               anything else.


subroutine set_analytic_sedov (N, r, rho, p, v, t, gamma, E, & 
     p_ambient, rho_ambient)

!==============================================================================

  implicit none

!  Arguments
  integer, intent(IN)     :: N
  real, intent(IN)        :: gamma, t, rho_ambient, E, p_ambient
  real, intent(IN), dimension(N)  :: r
  real, intent(OUT), dimension(N) :: rho, v, p

!  Local variables

  real    beta, R0, dr, xi, nu1, nu2, nu3, nu4, nu5, VV, G, Z, & 
       kappa, zeta, epsilon, c2sqr, k, gamp1, gam7, gamm1
  integer i

!==============================================================================

  if (gamma .ne. 1.4) then
     write (*,*) 'Warning!  Simulation_initBlock() found gamma<>1.4 and t>0.'
     write (*,*) '          Analytical initial conditions will be'
     write (*,*) '          wrong.  Assuming beta=1.033...'
  endif

  !               Compute dimensionless scaling constant and explosion radius.

  beta = 1.033
  R0   = beta * (E*t*t/rho_ambient)**0.2

  !               Compute exponents for self-similar solution.

  nu1 = - (13.*gamma*gamma - 7.*gamma + 12.) / & 
       ((3.*gamma - 1.) * (2.*gamma+1.))
  nu2 = 5. * (gamma - 1.) / (2.*gamma + 1.)
  nu3 = 3. / (2.*gamma + 1.)
  nu4 = - nu1 / (2. - gamma)
  nu5 = - 2. / (2. - gamma)

  !               Other useful combinations of gamma.

  gamp1 = gamma + 1.E0
  gamm1 = gamma - 1.E0
  gam7  = 7.E0 - gamma
  k     = gamp1 / gamm1

  !==============================================================================

  !               Generate the solution.

  do i = 1, N

     xi = r(i) / R0                ! Fraction of explosion radius.

     if (xi .le. 1.) then          ! Inside the explosion.

        ! Compute V(xi) using bisection.
        ! See Landau & Lifshitz for notation.
        call compute_sedov_v (xi, gamma, nu1, nu2, VV)

        G = k * (k*(gamma*VV-1.))**nu3 * & 
             (gamp1/gam7*(5.-(3.*gamma-1.)*VV))**nu4 * & 
             (k*(1.-VV))**nu5

        rho(i) = rho_ambient * G
        v(i)   = 2.*r(i)*VV / (5.*t)

        if (xi .le. 1.E-6) then     ! Use asymptotic r->0 solution.
           kappa = ( (0.5*gamp1/gamma)**2. * & 
                (gamp1/gam7* & 
                (5.-(3.*gamma-1.)/gamma))**nu1 )**(1./nu2) * & 
                gamm1/gamp1 / gamma

           epsilon = k**(nu5+1.) * (k*gamma*kappa)**nu3 * & 
                (gamp1/gam7*(3.*gamma-1.))**nu4 * & 
                ((2.*gamma+1)/gamma/(3.*gamma-1.)) * & 
                (gamm1/gamma)

           zeta = gamm1*gamm1/(2.*gamma*gamma*kappa)
           p(i) = rho_ambient/gamma * 0.16*(R0/t)**2 * epsilon*zeta
        else
           Z = gamma * gamm1 * (1.-VV) * VV**2 / (2.*(gamma*VV-1.))
           c2sqr = 0.16*(r(i)/t)**2 * Z
           p(i) = rho(i) * c2sqr / gamma
        endif

     else                          ! Outside the explosion.

        rho(i) = rho_ambient
        p(i)   = p_ambient
        v(i)   = 0.

     endif

  enddo

  return
end subroutine set_analytic_sedov



!  Subroutine:  compute_sedov_v()

!  Description: Compute the dimensionless velocity function V(xi) in the
!               Sedov problem using bisection.  See Landau & Lifshitz for
!               notation.

subroutine compute_sedov_v (xi, gamma, nu1, nu2, V)

  !==============================================================================

  implicit none

!  Arguments, LBR guessed intent on these
  real, intent(IN)  :: xi, gamma, nu1, nu2
  real, intent(OUT) :: V

! local variables
  real      VL, VR, xiL, xiR, Vmid, ximid, sedov_vfunc, tolerance, logxi
  integer   n_iter, n_iter_max
  parameter (n_iter_max = 500, tolerance = 1.E-6)

  !==============================================================================

  if (xi .le. 1.E-6) then         ! Use asymptotic xi->0 solution.

     !CD: The expression gamma*VV-1 on line 387 in set_analytic_sedov
     !(after this function call) should be 0.0 but generates a tiny value
     !(gdb) p gamma*VV-1.
     ! $1 = -1.1102230246251565e-16
     !This causes an FPE when it is raised to the power nu3.
     !When run in gdb it returns: 
     !"Cannot perform exponentiation: Numerical argument out of domain"
     !We can avoid the FPE by adding 1.E-8 to the original value of V:
     !V = 1./gamma
     V = 1./gamma + 1.E-8

  else                            ! Do bisection root search.


     logxi = alog(xi)             !was dlog but this fails on xlf compiler 
     VL = 1./gamma
     VR = 2./gamma
     xiL = sedov_vfunc(VL, gamma, nu1, nu2)
     xiR = sedov_vfunc(VR, gamma, nu1, nu2)
     n_iter = 1
10   Vmid = 0.5 * (VL + VR)
     ximid = sedov_vfunc(Vmid, gamma, nu1, nu2)
     if ((abs(ximid - logxi) .le. tolerance) .or. & 
          (n_iter .gt. n_iter_max)) goto 20
     n_iter = n_iter + 1
     if (ximid .gt. logxi) then
        VR = Vmid
     else
        VL = Vmid
     endif
     goto 10

#ifdef DEBUG_SIM
20   if (n_iter .gt. n_iter_max) & 
          write (*,*) 'compute_sedov_v:  did not reach ', & 
          'max precision for xi = ', xi
     V = Vmid
#else
20   V = Vmid
#endif

  endif

  return
end subroutine compute_sedov_v



!******************************************************************************

!  Function:    sedov_vfunc()

!  Description: Function to use in bisection root search (compute_sedov_v()).

real function sedov_vfunc (V, gamma, nu1, nu2)

  implicit none

! arguments, LBR guessed intent on these
  real, intent(IN)  :: V, gamma, nu1, nu2

!  local variables
  real gamp1, gamm1, gam7, k, xi, Vfpe


  !CD: The sub expression gamma*V-1 in xi expression causes a problem when 
  !V = 1./gamma.  Mathematically this is 0.0, but rounding gives a 
  !subexpression value of order -1.E-16.  This value is then passed
  !to the alog function (See alog(gamma*Vfpe-1.)) which gives -Infinity
  !for alog(0.0) and NaN for alog(-1.E-16).  Neither is good, so I 
  !guard against both cases by adding 1.E-8 to V.
  if ((gamma*V-1.) <= 0.0) then
     Vfpe = 1./gamma + 1.E-8
  else
     Vfpe = V
  end if
 

  gamp1 = gamma + 1.
  gamm1 = gamma - 1.
  gam7  = 7. - gamma
  k     = gamp1 / gamm1

  xi = nu1*alog(5.-(3.*gamma-1.)*Vfpe) + & 
       nu2*alog(gamma*Vfpe-1.) - & 
       nu1*alog(gam7/gamp1) - nu2*alog(gamm1/gamp1) - & 
       2.*alog(0.5*gamp1)
  sedov_vfunc = 0.2 * xi

  return
end function sedov_vfunc



!******************************************************************************

!  Routine:     sim_find()

!  Description: Given a monotonically increasing table x(N) and a test value
!               x0, return the index i of the largest table value less than
!               or equal to x0 (or 0 if x0 < x(1)).  Use binary search.

subroutine sim_find (x, N, x0, i)

  implicit none

! Arguments, LBR guessed intent on these
  integer, intent(IN) :: N
  integer, intent(OUT):: i
  real, intent(IN)    :: x(N), x0

! local variables
  integer  il, ir, im

  if (x0 .lt. x(1)) then

     i = 0

  elseif (x0 .gt. x(N)) then

     i = N

  else

     il = 1
     ir = N
10   if (ir .eq. il+1) goto 20
     im = (il + ir) / 2
     if (x(im) .gt. x0) then
        ir = im
     else
        il = im
     endif
     goto 10
20   i = il

  endif

  return
end subroutine sim_find
