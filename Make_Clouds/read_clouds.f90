program read_clouds
  implicit none

    integer :: N, iii
    real*8, allocatable :: x_pos(:), y_pos(:), z_pos(:), cloud_rad(:)
    
    open (unit = 1, file = "clouds.txt")
    
    read(1,*) N


    print *,N
    
    allocate(x_pos(N))
    allocate(y_pos(N))
    allocate(z_pos(N))
    allocate(cloud_rad(N))
   
   
    
    DO iii=1, N
       PRINT *, iii
       read(1,*) x_pos(iii), y_pos(iii), z_pos(iii), cloud_rad(iii)

    ENDDO


    close(1)

     DO iii=1, N
        PRINT *, x_pos(iii), y_pos(iii), z_pos(iii), cloud_rad(iii)
     ENDDO
    
    
    
  end program read_clouds










