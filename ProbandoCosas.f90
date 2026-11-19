program ProbandoCosas
    use FileTools
    use Derivatives
    implicit none
    character(len=*), parameter :: DataFile='PhiData'
    double precision, allocatable, dimension(:) :: X,Y
    double precision, allocatable, dimension(:,:) :: Phi,dxPhi,dyPhi

    call MakePhiXY(X,Y,Phi,DataFile)

    call OutputMatrix(X,Y,Phi,'PhiMatrix.dat')

    call partialx(X,Y,Phi,dxPhi) 

    call partialy(X,Y,Phi,dyPhi) 

    call OutputMatrix(X,Y,dxPhi,'dxPhiMatrix.dat')

    call OutputMatrix(X,Y,dyPhi,'dyPhiMatrix.dat')

end program ProbandoCosas



 
