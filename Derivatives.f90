module Derivatives
    implicit none
    
contains

    subroutine partialx(X,Y,Phi,dxPhi,Metal)
        double precision, intent(in) :: X(0:),Y(0:),Phi(0:,0:)
        logical, intent(in) :: Metal(0:,0:)
        double precision, allocatable, intent(out) :: dxPhi(:,:)
        
        integer :: i,j,Nx,Ny
        double precision :: h

        Nx = size(X)
        Ny = size(Y)
        allocate(dxPhi(0:Nx-1,0:Ny-1))

        h = X(1)-X(0)
        do j = 0,Ny-1
            if (Metal(0,j).OR.Metal(1,j)) then
                dxPhi(0,j) = 0.d0
            else
                dxPhi(0,j) = (Phi(1,j)-Phi(0,j))/h
            end if
            do i = 1,Nx-2
                if (Metal(i-1,j).OR.Metal(i+1,j)) then
                    dxPhi(i,j) = 0.d0
                else
                    dxPhi(i,j) = (Phi(i+1,j)-Phi(i-1,j))/(2.d0*h)
                end if
            end do
            if (Metal(Nx-2,j).OR.Metal(Nx-1,j)) then
                dxPhi(Nx-1,j) = 0.d0
            else
                dxPhi(Nx-1,j) = (Phi(Nx-1,j)-Phi(Nx-2,j))/h
            end if
        end do

    end subroutine partialx

    subroutine partialy(X,Y,Phi,dyPhi,Metal)
        double precision, intent(in) :: X(0:),Y(0:),Phi(0:,0:)
        logical, intent(in) :: Metal(0:,0:)
        double precision, allocatable, intent(out) :: dyPhi(:,:)

        integer :: i,j,Nx,Ny
        double precision :: h
        
        Nx = size(X)
        Ny = size(Y)
        allocate(dyPhi(0:Nx-1,0:Ny-1))

        h = Y(1)-Y(0)
        do i = 0,Nx-1
            if (Metal(i,0).OR.Metal(i,1)) then
                dyPhi(i,0) = 0.d0
            else
                dyPhi(i,0) = (Phi(i,1)-Phi(i,0))/h
            end if
            do j = 1,Ny-2
                if (Metal(i,j-1).OR.Metal(i,j+1)) then
                    dyPhi(i,j) = 0.d0
                else
                    dyPhi(i,j) = (Phi(i,j+1)-Phi(i,j-1))/(2.d0*h)
                end if
            end do
            if (Metal(i,Ny-2).OR.Metal(i,Ny-1)) then
                dyPhi(i,Ny-1) = 0.d0
            else
                dyPhi(i,Ny-1) = (Phi(i,Ny-1)-Phi(i,Ny-2))/h
            end if
        end do
        
    end subroutine partialy

end module Derivatives