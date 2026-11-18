program ProbandoCosas
    use FileTools
    use Derivatives
    implicit none
    character(len=*), parameter :: DataFile='PhiData'
    character(len=*), parameter :: OutMatrix='PhiMatrix.dat'
    double precision, allocatable, dimension(:) :: X,Y
    double precision, allocatable, dimension(:,:) :: Phi

    call MakePhiXY(X,Y,Phi,DataFile)

    call OutputMatrix(X,Y,Phi,OutMatrix)

    print*,size(Phi)

end program ProbandoCosas



 