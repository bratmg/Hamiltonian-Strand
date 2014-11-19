program test_faces
  !
  integer, allocatable :: conn(:,:)
  integer, allocatable :: neig(:)
  integer, allocatable :: faces(:)
  !
  open(unit=1,file='con',status='unknown')
  read(1,*) ncells
  allocate(conn(3,ncells)
  read(1,*) conn
  close(1)
  !
  call find_edges(conn,neig,faces,ncells,3)
  !
end program test_faces
