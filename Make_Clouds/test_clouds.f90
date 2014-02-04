subroutine read_clouds(N,x_pos, y_pos, z_pos, cloud_rad)
    implicit none

    integer :: N, iii
    real*8 :: x_pos(N), y_pos(N), z_pos(N), cloud_rad(N)
    
    open (unit = 2, file = "clouds.txt")
    
    read(2,*) N

    
    print *,N
    
    
    DO iii=1, N
       PRINT *, iii
       read(2,*) x_pos(iii), y_pos(iii), z_pos(iii), cloud_rad(iii)
       
    ENDDO
    
    
    close(2)
    
    DO iii=1, N
       PRINT *, x_pos(iii), y_pos(iii), z_pos(iii), cloud_rad(iii)
    ENDDO
    return
  end subroutine read_clouds

integer function get_number_of_clouds()
  
  real*8, allocatable :: x_pos(:), y_pos(:), z_pos(:), cloud_rad(:)
  
  open (unit = 1, file = "clouds.txt")
  
  read(1,*) get_number_of_clouds
  close(1)
  return
end function get_number_of_clouds

integer function test_current_position(x,y,z,N, x_pos, y_pos, z_pos, cloud_rad)
   integer :: iii,N
   real*8 ::x,y,z,x_pos(N), y_pos(N), z_pos(N), cloud_rad(N), distance
   test_current_position=0
   
   do iii=1, N
      distance=sqrt(((x-x_pos(iii))**2) + ((y-y_pos(iii))**2)  )
      print *, distance, x, x_pos(iii), y, y_pos(iii), cloud_rad(iii)
      if (distance <= cloud_rad(iii)) then 
         test_current_position=1
         print *, 'PING'
      endif
   end do

   return

end function test_current_position
program test_clouds

  implicit none
  real*8, allocatable ::  x_pos(:), y_pos(:), z_pos(:), cloud_rad(:),xx,yy,zz
  integer :: N,get_number_of_clouds, test_current_position, test
  xx=15.378
  yy=8.94
  zz=1

  N=get_number_of_clouds()

  allocate(x_pos(N))
  allocate(y_pos(N))
  allocate(z_pos(N))
  allocate(cloud_rad(N))

  call read_clouds(N, x_pos, y_pos, z_pos, cloud_rad)

  test=test_current_position(xx,yy,zz,N, x_pos, y_pos, z_pos, cloud_rad)

  print *, test

end program test_clouds


