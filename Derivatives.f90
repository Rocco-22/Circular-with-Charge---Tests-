module Derivatives
    implicit none
    
contains

    subroutine partialx(X,Y,Phi,dxPhi,Deep)
        double precision, intent(in) :: X(0:),Y(0:),Phi(0:,0:)
        logical, intent(in) :: Deep(0:,0:)
        double precision, allocatable, intent(out) :: dxPhi(:,:)

        integer :: i,j,Nx,Ny
        double precision :: h
        logical :: lo,hi

        Nx = size(X); Ny = size(Y)
        allocate(dxPhi(0:Nx-1,0:Ny-1))
        dxPhi = 0.d0
        h = X(1)-X(0)
        do j = 0,Ny-1
            do i = 0,Nx-1
                if (Deep(i,j)) cycle
                lo = .false.; hi = .false.
                if (i > 0)    lo = .not. Deep(i-1,j)
                if (i < Nx-1) hi = .not. Deep(i+1,j)
                if (lo .and. hi) then
                    dxPhi(i,j) = (Phi(i+1,j)-Phi(i-1,j))/(2.d0*h)
                else if (hi) then
                    dxPhi(i,j) = (Phi(i+1,j)-Phi(i,j))/h
                else if (lo) then
                    dxPhi(i,j) = (Phi(i,j)-Phi(i-1,j))/h
                end if
            end do
        end do
    end subroutine partialx

    subroutine partialy(X,Y,Phi,dyPhi,Deep)
        double precision, intent(in) :: X(0:),Y(0:),Phi(0:,0:)
        logical, intent(in) :: Deep(0:,0:)
        double precision, allocatable, intent(out) :: dyPhi(:,:)

        integer :: i,j,Nx,Ny
        double precision :: h
        logical :: lo,hi
        
        Nx = size(X) 
        Ny = size(Y)
        allocate(dyPhi(0:Nx-1,0:Ny-1))

        h = Y(1)-Y(0)
        do i = 0,Nx-1
            do j = 0,Ny-1
                if (Deep(i,j)) cycle
                lo = .false.; hi = .false.
                if (j > 0)    lo = .not. Deep(i,j-1)
                if (j < Ny-1) hi = .not. Deep(i,j+1)
                if (lo .and. hi) then
                    dyPhi(i,j) = (Phi(i,j+1)-Phi(i,j-1))/(2.d0*h)
                else if (hi) then
                    dyPhi(i,j) = (Phi(i,j+1)-Phi(i,j))/h
                else if (lo) then
                    dyPhi(i,j) = (Phi(i,j)-Phi(i,j-1))/h
                end if
            end do
        end do
    end subroutine partialy
    ! Metal(0:nx,0:ny) are cells; Deep(0:nx,0:ny) is on the node grid (same as Phi).
    ! A node is "deep" when all 4 surrounding cells are metal: the solver zeroed it.
    ! Metal(0:nx,0:ny) are cells; Deep(0:nx,0:ny) is on the node grid (same as Phi).
    ! Copies the two rules of BuildAndSolve2DPoisson:
    !  1. A node is fixed (Variable = .false.) if one of its 4 surrounding cells is metal
    !     (dsoEpsA(3,..) < 1e-5), and every node on the last row/column (ix=nx, iy=ny)
    !     is never visited by the solver, so it stays fixed too.
    !  2. A fixed node whose 4 neighbouring nodes are all fixed is set to Phi = 0.
    subroutine NodeDeepMetal(Metal,Deep)
        logical, intent(in) :: Metal(0:,0:)
        logical, allocatable, intent(out) :: Deep(:,:)
        logical, allocatable :: Fixed(:,:)
        integer :: i,j,Nx,Ny

        Nx = ubound(Metal,1)
        Ny = ubound(Metal,2)

        ! Padded by one node on each side, so the neighbour lookups never go out of range.
        ! Outside the visited range 0..Nx-1, 0..Ny-1 the solver leaves Variable = .false.
        allocate(Fixed(-1:Nx+1,-1:Ny+1))
        Fixed = .true.
        do j = 0,Ny-1
            do i = 0,Nx-1
                Fixed(i,j) = CellIsMetal(i,j)   .or. CellIsMetal(i-1,j) .or. &
                            CellIsMetal(i-1,j-1) .or. CellIsMetal(i,j-1)
            end do
        end do

        allocate(Deep(0:Nx,0:Ny))
        do j = 0,Ny
            do i = 0,Nx
                Deep(i,j) = Fixed(i,j)   .and. &
                            Fixed(i-1,j) .and. Fixed(i+1,j) .and. &
                            Fixed(i,j-1) .and. Fixed(i,j+1)
            end do
        end do

        ! The cleanup loop never visits ix=nx or iy=ny, but Phi is 0 there anyway.
        Deep(Nx,:) = .true.
        Deep(:,Ny) = .true.

    contains

        ! Cells outside 0:Nx,0:Ny count as vacuum (Material = 0 in the solver).
        logical function CellIsMetal(ic,jc)
            integer, intent(in) :: ic,jc
            CellIsMetal = .false.
            if (ic >= 0 .and. ic <= Nx .and. jc >= 0 .and. jc <= Ny) then
                CellIsMetal = Metal(ic,jc)
            end if
        end function CellIsMetal

    end subroutine NodeDeepMetal

end module Derivatives