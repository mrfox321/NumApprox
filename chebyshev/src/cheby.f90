module cheby
implicit none

contains
!====================================================
subroutine lanczos_dense(ifc_s,n,omeg_new)


!use structures
!use mega_block
use sparsekit
implicit none

integer, intent(in) :: n
integer :: i,j
complex(8), intent(in) :: ifc_s(n,n)
complex(8) :: v_0(n),v_1(n),r(n)
real(8), intent(out) :: omeg_new
real(8) :: beta(1000), alpha(1000),omeg_old,eps,ratio
real(8), allocatable :: T(:,:),eig(:)
integer :: mm,clock
integer, allocatable :: seed(:)
real :: y

call random_seed(size=mm)
allocate(seed(mm))
CALL SYSTEM_CLOCK(COUNT=clock)
!seed = clock + 37 * (/ (i - 1, i = 1, mm) /)
seed  = 1367986020 + 37 * (/ (i - 1, i = 1, mm) /)
CALL RANDOM_SEED(PUT = seed)
deallocate(seed)

do i=1,n
    call random_number(y)
    v_0(i) = y
enddo

v_0 =  v_0/(L2_norm(v_0,n))
r = v_0

beta(1) = L2_norm(v_0,n)

omeg_old = 0d0
eps = 0.0001d0

do j=1,100
    v_1 = r/beta(j)
    r = matmul(ifc_s,v_1)
    r = r-beta(j)*v_0
    alpha(j) = (z_dot(v_1,r,n)) !v_1^dagger .dot. r
    r = r-alpha(j)*v_1
    beta(j+1) = L2_norm(r,n)
    v_0 = v_1
    if (mod(j,10).eq.0) then
        allocate(T(j,j))
        allocate(eig(j))
        T = 0d0
        do i=1,j
            if (i.ne.j) then
                T(i+1,i) = beta(i+1)
                T(i,i+1) = beta(i+1)
                T(i,i) = alpha(i)
            elseif (i.eq.j) then
                T(i,i) = alpha(i)
            endif
        enddo
        call diag_sym(j,T,eig)
        omeg_new = maxval(eig)
        ratio = abs((omeg_new-omeg_old)/(omeg_new))
        if (ratio.lt.eps) then
            exit
        else
            omeg_old = omeg_new
        endif
        deallocate(T)
        deallocate(eig)
    endif
enddo





end subroutine lanczos_dense
!====================================================
subroutine lanczos(ifc_s,ia,ja,n,nnz,omeg_new)


!use structures
!use mega_block
use sparsekit
implicit none

integer, intent(in) :: n,nnz
integer :: ia(nnz),ja(nnz),i,j
complex(8), intent(in) :: ifc_s(nnz)
complex(8) :: v_0(n),v_1(n),r(n)
real(8), intent(out) :: omeg_new
real(8) :: beta(1000), alpha(1000),omeg_old,eps,ratio
real(8), allocatable :: T(:,:),eig(:)
integer :: mm,clock
integer, allocatable :: seed(:)
real :: y

call random_seed(size=mm)
allocate(seed(mm))
CALL SYSTEM_CLOCK(COUNT=clock)
!seed = clock + 37 * (/ (i - 1, i = 1, mm) /)
seed  = 1367986020 + 37 * (/ (i - 1, i = 1, mm) /)
!print*,'the CLOCK VALUE is', clock
!print*,'the SEED IS',seed
CALL RANDOM_SEED(PUT = seed)
deallocate(seed)



v_0 =  v_0/(L2_norm(v_0,n))
r = v_0

beta(1) = L2_norm(v_0,n)

omeg_old = 0d0
eps = 0.0001d0

do j=1,100
    v_1 = r/beta(j)
    call amux(n,v_1,r,ifc_s,ja,ia)
    r = r-beta(j)*v_0
    alpha(j) = (z_dot(v_1,r,n)) !v_1^dagger .dot. r
    r = r-alpha(j)*v_1
    beta(j+1) = L2_norm(r,n)
    v_0 = v_1
    if (mod(j,10).eq.0) then
        allocate(T(j,j))
        allocate(eig(j))
        T = 0d0
        do i=1,j
            if (i.ne.j) then
                T(i+1,i) = beta(i+1)
                T(i,i+1) = beta(i+1)
                T(i,i) = alpha(i)
            elseif (i.eq.j) then
                T(i,i) = alpha(i)
            endif
        enddo
        call diag_sym(j,T,eig)
        omeg_new = maxval(eig)
        ratio = abs((omeg_new-omeg_old)/(omeg_new))
        if (ratio.lt.eps) then
            exit
        else
            omeg_old = omeg_new
        endif
        deallocate(T)
        deallocate(eig)
    endif
enddo

end subroutine lanczos
!====================================================
function L2_norm(a,n) result(norm)

implicit none
integer, intent(in) :: n
complex(8), intent(in) :: a(n)
integer :: i
real(8) :: norm


norm = 0d0
do i=1,n
    norm = norm + dreal(a(i))**2 + dimag(a(i))**2
enddo

norm = dsqrt(norm)

end function L2_norm
!====================================================
function z_dot(a,b,n) result(z)

implicit none

integer, intent(in) :: n
complex(8), intent(in) :: a(n), b(n)
complex(8) :: z
integer :: i


z = 0d0
do i=1,n
    z = z + dconjg(a(i))*b(i)
enddo

end function z_dot
!====================================================
subroutine rescale_ifc(ifc_s,nnz,n,ia,ja,omega_max)

implicit none

integer, intent(in) :: nnz,n
integer, intent(in) :: ia(nnz),ja(nnz)
real(8), intent(in) :: omega_max
integer :: i,j
complex(8) :: ifc_s(nnz)
real(8) :: eps,a,b


eps = 0.01d0
a = omega_max/(2d0-eps)
b = omega_max/2d0

do i=1,n
    do j=ia(i),(ia(i+1)-1)
        if (ja(j).eq.i) then
            ifc_s(j) = (ifc_s(j)-b)/a
        else
            ifc_s(j) = ifc_s(j)/a
        endif
    enddo
enddo

end subroutine rescale_ifc
!====================================================
subroutine jackson_kernel(mu,n)

implicit none

integer, intent(in) :: n
real(8) :: mu(0:n),g(0:n),pi
integer :: i

pi=4.D0*DATAN(1.D0)

do i=0,n
    g(i) = real(n-i+1)*cos((pi*i)/real(n+1))+sin((pi*i)/real(n+1))/(tan(pi/real(n+1)))
    g(i) = g(i)/(n+1)
    mu(i) = g(i)*mu(i)
enddo

end subroutine jackson_kernel
!====================================================
subroutine lorentz_kernel(mu,n,lam)

implicit none

integer, intent(in) :: n
real(8), intent(in) :: lam  !empirically between 3-5  (will stick to jackson...)
real(8) :: mu(0:n),g(0:n)
integer :: i

do i=0,n
    g(i) = sinh(lam*(1d0-real(i)/real(n)))/sinh(lam)
    mu(i) = g(i)*mu(i)
enddo

end subroutine lorentz_kernel
!====================================================
subroutine diag_sym(N,A,W)

implicit none
integer, intent(in) :: N
integer :: LDA
INTEGER          LWMAX
PARAMETER        ( LWMAX = 5000 )
integer :: INFO,LWORK

real(8) :: A( N, N ), W( N ), WORK( LWMAX )

LDA = N


LWORK = -1
CALL DSYEV( 'N', 'U', N, A, LDA, W, WORK, LWORK, INFO )
LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
CALL DSYEV( 'N', 'U', N, A, LDA, W, WORK, LWORK, INFO )

end subroutine diag_sym
!====================================================
subroutine dos_moment(ifc_s,ia,ja,nnz,n_mom,n_dim,mu)

implicit none


integer, intent(in) :: n_mom,n_dim,nnz
complex(8), intent(in) :: ifc_s(nnz)
integer, intent(in) :: ia(nnz),ja(nnz)
real(8), intent(out) :: mu(0:n_mom*2)
complex(8) :: vec_0(n_dim),vec_1(n_dim),vec_2(n_dim),vec_in(n_dim)
integer :: mm,clock,i,j
integer, allocatable :: seed(:)
real(8) :: y


call random_seed(size=mm)
allocate(seed(mm))
CALL SYSTEM_CLOCK(COUNT=clock)
!seed = clock + 37 * (/ (i - 1, i = 1, mm) /)
seed  = 1367986020 + 37 * (/ (i - 1, i = 1, mm) /)
CALL RANDOM_SEED(PUT = seed)
deallocate(seed)


do i=1,n_dim
    call random_number(y)
    y = sqrt(12d0)*(y-0.5d0)
    vec_in(i) = y
enddo


mu(0) = z_dot(vec_in,vec_in,n_dim)
vec_0 = vec_in
call amux(n_dim,vec_in,vec_1,ifc_s,ja,ia)
mu(1) = z_dot(vec_in,vec_1,n_dim)

do j=2,n_mom
    call amux(n_dim,vec_1,vec_2,ifc_s,ja,ia)
    vec_2 = 2d0*vec_2 - vec_0
    vec_0 = vec_1
    vec_1 = vec_2
    mu(j) = z_dot(vec_in,vec_1,n_dim)
    if (2*j.gt.n_mom) then
        mu(2*j-1) = 2d0*z_dot(vec_1,vec_0,n_dim) - mu(1)
        mu(2*j) = 2d0*z_dot(vec_1,vec_1,n_dim) - mu(0)
    endif
enddo

end subroutine dos_moment
!====================================================
subroutine cheby_calc(mu,n_mom,mesh)

implicit none

integer, intent(in) :: n_mom,mesh
real(8), intent(in) :: mu(0:n_mom)
real(8) :: omeg(mesh),cheb(mesh)

do i=1,mesh
    omeg(i) = real(i-1)*1.98d0/real(mesh-1) - 0.99d0
enddo

do i=1,mesh
    call cheby_sum(mu,n_mom,omeg(i),cheb(i))
enddo


end subroutine cheby_calc
!====================================================
subroutine cheby_sum(mu,n_mom,omeg,sum)

implicit none

integer, intent(in) :: n_mom
real(8), intent(in) :: mu(0:n_mom),omeg
real(8), intent(out) :: sum
integer :: i
real(8) :: pi,T(n_mom)



pi=4.D0*DATAN(1.D0)

T(0) = 1d0
T(1) = omeg

do i =2,n_mom
    T(i) = 2*omeg*T(i-1) - T(i-2)
enddo

sum  = T(0)*mu(0)

do i=1,n_mom
    sum = sum + 2d0*mu(i)*T(i)
enddo

sum = sum/(pi*sqrt(1d0-omeg**2))


end subroutine cheby_sum
!====================================================
end module cheby