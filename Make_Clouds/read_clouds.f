program read_clouds
  implicit none

  character(1024) :: buffer
  integer :: pos
  real :: var1, var2,var3, var4

  OPEN(1,FILE='clouds.txt')
  read(1,"(A)") buffer
  pos = index(buffer, ",")
  var1 = buffer(1:pos-1)
  read(buffer(pos+1:), *) var2, var3, var4
  print *, var1, var2, var3, var4

end program read_clouds
