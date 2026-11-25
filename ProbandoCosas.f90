program ProbandoCosas
    use FileTools
    use Derivatives
    implicit none
    character(len=*), parameter :: PhiFile='input/PhiData',GridFile='input/Gitterfuellung',EpsFile='input/EpsData'
    double precision, allocatable, dimension(:) :: X,Y,Xg,Yg
    double precision, allocatable, dimension(:,:) :: Phi,dxPhi,dyPhi,NablaPhi
    integer, allocatable, dimension(:,:) :: Material
    logical, allocatable, dimension(:,:) :: Metal,Deep
    integer :: iu

    call ReadPhiXY(X,Y,Phi,PhiFile)

    call ReadMaterial(Xg,Yg,Material,Metal,GridFile,EpsFile)

    call NodeDeepMetal(Metal,Deep)

    call partialx(X,Y,Phi,dxPhi,Deep)

    call partialy(X,Y,Phi,dyPhi,Deep) 

    NablaPhi = sqrt(dxPhi**2 + dyPhi**2)


    !OUTPUTS
        call WriteOutputMatrix(X,Y,Phi,'PhiMatrix.dat')
        call WriteOutputMatrix(X,Y,dxPhi,'dxPhiMatrix.dat')
        call WriteOutputMatrix(X,Y,dyPhi,'dyPhiMatrix.dat')
        
end program ProbandoCosas






 
