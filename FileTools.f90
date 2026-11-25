module FileTools
    implicit none
contains
    subroutine ReadPhiXY(X,Y,Phi,PhiFile)
        character(len=*), intent(in) :: PhiFile
        double precision, allocatable, intent(out) :: X(:),Y(:),Phi(:,:)
        
        integer :: i,j,nx0,ny0,Nx,Ny
        integer :: iu


        open(newunit=iu,file=PhiFile)
        read(iu,*) nx0,Nx,ny0,Ny
        Nx = Nx-(nx0)
        Ny = Ny-(ny0)

        allocate(X(0:Nx), Y(0:Ny), Phi(0:Nx,0:Ny))

        read(iu,*)

        do i = 0,Nx
            read(iu,*) X(i)
        end do

        do j = 0,Ny
            read(iu,*) Y(j)
        end do

        do i = 1,6
            read(iu,*)
        end do

        do j = 0,Ny
            do i = 0,Nx
                read(iu,*) Phi(i,j)
            end do
        end do
        close(iu)
    end subroutine ReadPhiXY

    subroutine ReadMaterial(X,Y,Material,Metal,GridFile,EpsFile)
        character(len=*), intent(in) :: GridFile, EpsFile
        double precision, allocatable, intent(out) :: X(:), Y(:)
        integer, allocatable, intent(out) :: Material(:,:)
        logical, allocatable, intent(out) :: Metal(:,:)

        integer :: iu, nx, ny, nz, i, j, k, nrun, mat, nmat
        double precision :: eps, zk
        double precision :: Epsr(0:10)

        open(newunit=iu, file=GridFile, status='old')
        read(iu,*) nx, ny, nz
        read(iu,*)

        allocate(X(0:nx+1), Y(0:ny+1), Material(0:nx,0:ny), Metal(0:nx,0:ny))
        Material = 0

        do i = 0,nx+1
            read(iu,*) X(i)
        end do
        do j = 0,ny+1
            read(iu,*) Y(j)
        end do
        do k = 0,nz+1
            read(iu,*) zk        ! z planes, not used
        end do

        do k = 0,nz/2            ! the last plane read (nz/2) is the one kept
            do j = 0,ny
                i = 0
                do while (i <= nx)
                    read(iu,*) nrun, mat
                    if (i+nrun-1 > nx) error stop 'ReadMaterial: run exceeds row'
                    if (k == nz/2) Material(i:i+nrun-1,j) = mat
                    i = i + nrun
                end do
            end do
        end do
        close(iu)

        Epsr = 1.d0
        open(newunit=iu, file=EpsFile, status='old')
        read(iu,*) nmat
        do k = 1,nmat
            read(iu,*) mat, eps
            Epsr(mat) = eps
        end do
        close(iu)

        do i = 0,nx
            do j = 0,ny
                Metal(i,j) = Epsr(Material(i,j)) > 1.d10
            end do
        end do
    end subroutine ReadMaterial

    subroutine WriteOutputMatrix(X,Y,Matrix,OutMatrix)
        character(len=*), intent(in) :: OutMatrix
        double precision, dimension(0:), intent(in) :: X,Y
        double precision, dimension(0:,0:), intent(in) :: Matrix

        integer :: i
        integer :: iu

        open(newunit=iu,file=OutMatrix)
        write(iu,'(*(ES15.7))') real(size(X)), X         
        do i = LBOUND(Matrix,2), UBOUND(Matrix,2)
            write(iu,'(*(ES15.7))') Y(i), Matrix(:,i)
        end do
        close(iu)
    end subroutine WriteOutputMatrix


end module FileTools

