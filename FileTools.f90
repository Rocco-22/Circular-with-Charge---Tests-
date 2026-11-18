module FileTools
    implicit none
contains
    subroutine MakePhiXY(X,Y,Phi,DataFile)
        implicit none
        integer :: i,j,nx0,ny0,Nx,Ny
        character(len=*) :: DataFile
        double precision, allocatable, dimension(:) :: X,Y
        double precision, allocatable, dimension(:,:) :: Phi
        
        open(10,file=DataFile)
        read(10,*) nx0,Nx,ny0,Ny
        Nx = Nx-(nx0)
        Ny = Ny-(ny0)

        allocate(X(0:Nx), Y(0:Ny), Phi(0:Nx,0:Ny))

        read(10,*)

        do i = 0,Nx
            read(10,*) X(i)
        end do

        do j = 0,Ny
            read(10,*) Y(j)
        end do

        do i = 1,6
            read(10,*)
        end do

        do j = 0,Ny
            do i = 0,Nx
                read(10,*) Phi(i,j)
            end do
        end do
        close(10)
    end subroutine

    subroutine OutputMatrix(X,Y,Phi,OutMatrix)
        implicit none
        integer :: i
        character(len=*) :: OutMatrix
        double precision, allocatable, dimension(:) :: X,Y
        double precision, allocatable, dimension(:,:) :: Phi

        open(11, file=OutMatrix)
        write(11,'(*(ES15.7))') real(size(X)), X         
        do i = 0, size(Y)-1
            write(11,'(*(ES15.7))') Y(i), Phi(:,i)
        end do
        close(11)
    end subroutine

end module FileTools

