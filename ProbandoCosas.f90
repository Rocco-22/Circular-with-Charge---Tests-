program ProbandoCosas
    use FileTools
    use Derivatives
    implicit none
    character(len=*), parameter :: PhiFile='PhiData',GridFile='Gitterfuellung',EpsFile='EpsData'
    double precision, allocatable, dimension(:) :: X,Y,Xg,Yg
    double precision, allocatable, dimension(:,:) :: Phi,dxPhi,dyPhi
    integer, allocatable, dimension(:,:) :: Material
    logical, allocatable, dimension(:,:) :: Metal
    double precision, allocatable, dimension(:,:) :: RealMetal
    integer :: i,j

    call ReadPhiXY(X,Y,Phi,PhiFile)

    call ReadMaterial(Xg,Yg,Material,Metal,GridFile,EpsFile)

    

    call partialx(X,Y,Phi,dxPhi,Metal) 

    call partialy(X,Y,Phi,dyPhi,Metal) 

    call WriteOutputMatrix(X,Y,Phi,'PhiMatrix.dat')

    call WriteOutputMatrix(X,Y,dxPhi,'dxPhiMatrix.dat')

    call WriteOutputMatrix(X,Y,dyPhi,'dyPhiMatrix.dat')

    call WriteOutputMatrix(Xg,Yg,dble(Material),'MaterialMatrix.dat')

    allocate(RealMetal(LBOUND(Metal,1):UBOUND(Metal,1),LBOUND(Metal,2):UBOUND(Metal,2)))

    do i=LBOUND(RealMetal,1),UBOUND(RealMetal,1)
        do j=LBOUND(RealMetal,2),UBOUND(RealMetal,2)
            if (Metal(i,j)) then
                RealMetal(i,j) = 1
            else
                RealMetal(i,j) = 0 
            end if
        end do
    end do

    call WriteOutputMatrix(Xg,Yg,RealMetal,'MetalMatrix.dat')

    print   *,'Finished printing matrices to files: PhiMatrix.dat, dxPhiMatrix.dat, dyPhiMatrix.dat, MaterialMatrix.dat, MetalMatrix.dat'
    

end program ProbandoCosas






 
