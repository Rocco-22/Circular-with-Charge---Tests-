module Derivatives
    implicit none
    
contains
    
    subroutine partialx(X,Y,Phi,dxPhi)
        implicit none
        integer :: i,j,Nx,Ny
        double precision :: h
        double precision, allocatable, dimension(:) :: X,Y
        double precision, allocatable, dimension(:,:) :: Phi,dxPhi
        
        Nx = size(X)
        Ny = size(Y)
        allocate(dxPhi(0:Nx-1,0:Ny-1))

        h = X(1)-X(0)
        do j = 0,Ny-1
            dxPhi(0,j) = (Phi(1,j)-Phi(0,j))/h
            do i = 1,Nx-2
                dxPhi(i,j) = (Phi(i+1,j)-Phi(i-1,j))/(2.d0*h)
            end do
            dxPhi(Nx-1,j) = (Phi(Nx-1,j)-Phi(Nx-2,j))/h
        end do

    end subroutine partialx

    subroutine partialy(X,Y,Phi,dyPhi)
        implicit none
        integer :: i,j,Nx,Ny
        double precision :: h
        double precision, allocatable, dimension(:) :: X,Y
        double precision, allocatable, dimension(:,:) :: Phi,dyPhi

        Nx = size(X)
        Ny = size(Y)
        allocate(dyPhi(0:Nx-1,0:Ny-1))

        h = Y(1)-Y(0)
        do i = 0,Nx-1
            dyPhi(i,0) = (Phi(i,1)-Phi(i,0))/h
            do j = 1,Ny-2
                dyPhi(i,j) = (Phi(i,j+1)-Phi(i,j-1))/(2.d0*h)
            end do
            dyPhi(i,Ny-1) = (Phi(i,Ny-1)-Phi(i,Ny-2))/h
        end do
        
    end subroutine partialy

end module Derivatives
