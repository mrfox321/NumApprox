!===============================

module structures

implicit none

!===============================
type blocks
    integer :: c(8,3),t(8)            !n(atom,neighbor#,1:3) = -1:1:1 (x,y,z) block coord
    real(8) :: m(8)                   !n(atom,neighbor#,4)   = atom # in that block
end type blocks
!===============================

end module structures
!===============================


module mega_block

implicit none

contains
!===============================

subroutine blocked(n_b,nx,ny,nz,size_per,rough,r_np,cube)

use structures
implicit none
integer, intent(in) :: nx,ny,nz,size_per,rough,r_np
integer :: rl(8,3),n(2,46,3),tt
integer, intent(out) :: n_b(8,46,4)
integer :: i,j,k,l,a(3),b(3),c(3),cc,key(8)
type(blocks) :: cube(nx,ny,nz)


rl(1,1) = 0
rl(1,2) = 0
rl(1,3) = 0
rl(2,1) = 0
rl(2,2) = 2
rl(2,3) = 2
rl(3,1) = 2
rl(3,2) = 0
rl(3,3) = 2
rl(4,1) = 2
rl(4,2) = 2
rl(4,3) = 0
rl(5,1) = 3
rl(5,2) = 1
rl(5,3) = 1
rl(6,1) = 1
rl(6,2) = 3
rl(6,3) = 1
rl(7,1) = 1
rl(7,2) = 1
rl(7,3) = 3
rl(8,1) = 3
rl(8,2) = 3
rl(8,3) = 3

do i=1,8
    key(i) = rl(i,1)+3*rl(i,2)+9*rl(i,3)
enddo

do i=1,nx
    do j=1,ny
        do k=1,nz
            do l=1,8
                cube(i,j,k)%c(l,1) = rl(l,1) + (i-1)*4
                cube(i,j,k)%c(l,2) = rl(l,2) + (j-1)*4
                cube(i,j,k)%c(l,3) = rl(l,3) + (k-1)*4

                cube(i,j,k)%t(l) = mod(rl(l,1)+rl(l,2)+rl(l,3)+1,2)+1
            enddo
        enddo
    enddo
enddo

call read_neighbors(n)

!n_b(8,46,4)   !n(atom,neighbor#,1:3) = -1:1:1 (x,y,z) block coord
               !n(atom,neighbor#,4)   = atom # in that block

do i=1,8
    do j=1,46
        tt = cube(1,1,1)%t(i)
        a(:) = cube(1,1,1)%c(i,:)
        b(:) = a(:) + n(tt,j,:)
        n_b(i,j,1) = floor(real(b(1))/4)
        n_b(i,j,2) = floor(real(b(2))/4)
        n_b(i,j,3) = floor(real(b(3))/4)
        c(1) = modulo(b(1),4)
        c(2) = modulo(b(2),4)
        c(3) = modulo(b(3),4)
        cc = c(1)+3*c(2)+9*c(3)
        do k=1,8
            if (cc.eq.key(k)) then
                n_b(i,j,4) = k
            endif
        enddo
    enddo
enddo


call setup_mass(cube,size_per/2,nx,ny,nz,rough,r_np,n_b)

end subroutine blocked
!=======================================================================================
subroutine read_neighbors(n)

implicit none

integer :: i,j
integer, intent(out) :: n(2,46,3)

open(9999,file='neighbors.dat',status='old',action='read')

do i=1,46
read(9999,*)n(1,i,1),n(1,i,2),n(1,i,3)
enddo
close(9999)

do i=1,46
do j=1,3
n(2,i,j) = -n(1,i,j)
enddo
enddo

end subroutine read_neighbors
!======================================
subroutine setup_mass(cube,hper,nx,ny,nz,rough,r_np,n_b)

use structures
implicit none

real(8) :: m1,m2,m3,m4
integer :: i,j,k,l,m,cnt,aa,bb,cc
integer, intent(in) :: hper,rough,r_np,n_b(8,46,4),nx,ny,nz
type(blocks),intent(inout) :: cube(nx,ny,nz)
integer :: mm,clock
integer, allocatable :: seed(:)
real :: y
real(8) :: p
logical :: iss(4)



m1 = 69.723d0  !(gallium)
m2 = 26.9827d0 !(aluminum)
m3 = 74.9216d0  !(arsenic)
m4 = 167.259d0 !(erbium)


do i=1,nz
    do j=1,nx
        do k=1,ny
            do l=1,8
                if (mod(i-1,2*hper).lt.hper) then
                    if (cube(j,k,i)%t(l).eq.1) then
                        cube(j,k,i)%m(l) = m1  !gallium
                    elseif (cube(j,k,i)%t(l).eq.2) then
                        cube(j,k,i)%m(l) = m3 !arsenic
                    endif
                else
                    if (cube(j,k,i)%t(l).eq.1) then
                        cube(j,k,i)%m(l) = m2 !aluminum
                    elseif (cube(j,k,i)%t(l).eq.2) then
                        cube(j,k,i)%m(l) = m3 !arsenic
                    endif
                endif
            enddo
        enddo
    enddo
enddo

call random_seed(size=mm)
allocate(seed(mm))
CALL SYSTEM_CLOCK(COUNT=clock)
seed = clock + 37 * (/ (i - 1, i = 1, mm) /)
!seed  = 1367986020 + 37 * (/ (i - 1, i = 1, mm) /)
!print*,'the CLOCK VALUE is', clock
!print*,'the SEED IS',seed
CALL RANDOM_SEED(PUT = seed)
deallocate(seed)


p = 0.5d0
if (rough.eq.1) then
    do i=1,nz
        do j=1,nx
            do k=1,ny
                do l=1,8
                    if (cube(j,k,i)%t(l).eq.1) then
                        if (mod(i-1,2*hper).eq.0)then
                            call random_number(y)
                            if (y.lt.p) then
                                cube(j,k,i)%m(l) = m2
                            endif
                        elseif (mod(i-1,2*hper).eq.(hper-1)) then
                            call random_number(y)
                            if (y.lt.p) then   
                                cube(j,k,i)%m(l) = m2
                            endif
                        elseif (mod(i-1,2*hper).eq.(hper)) then
                            call random_number(y)
                            if (y.lt.p) then
                                cube(j,k,i)%m(l) = m1
                            endif
                        elseif (mod(i-1,2*hper).eq.(2*hper-1)) then
                            call random_number(y)
                            if (y.lt.p) then
                                cube(j,k,i)%m(l) = m1
                            endif
                        endif
                    endif
                enddo
            enddo
        enddo
    enddo
endif



call random_seed(size=mm)
allocate(seed(mm))
CALL SYSTEM_CLOCK(COUNT=clock)
seed = clock + 37 * (/ (i - 1, i = 1, mm) /)
!seed  = 1729137069 + 37 * (/ (i - 1, i = 1, mm) /)
!print*,'the CLOCK VALUE is', clock
!print*,'the SEED IS',seed
CALL RANDOM_SEED(PUT = seed)
deallocate(seed)


p = 0.004d0
!p = 0.0004d0
cnt = 0
if (r_np.eq.1) then
do i=1,nz
    do j=1,nx
        do k=1,ny
            do l=1,8
                iss(1) = (mod(i-1,2*hper).eq.0)
                iss(2) = (mod(i-1,2*hper).eq.(hper-1))
                iss(3) = (mod(i-1,2*hper).eq.(hper))
                iss(4) = (mod(i-1,2*hper).eq.(2*hper-1))
                if (iss(1).or.iss(2).or.iss(3).or.iss(4)) then
                    call random_number(y)
                    if (y.lt.p) then
                        cnt = cnt+1
                        print*,'nano',i,j,cnt
                        if (cube(j,k,i)%t(l).eq.1) then
                            cube(j,k,i)%m(l) = m4
                        endif
                        do m=1,46
                            aa = mod(n_b(l,m,1)+i,nz)+1
                            bb = mod(n_b(l,m,2)+j,nx)+1
                            cc = mod(n_b(l,m,3)+k,ny)+1
                            if (cube(bb,cc,aa)%t(n_b(l,m,4)).eq.1) then
                                cube(bb,cc,aa)%m(n_b(l,m,4)) = m4
                            endif
                        enddo
                    endif
                endif
            enddo
        enddo
    enddo
enddo
endif

!open(1234,file='mass.dat',status='unknown')
!do i=1,nx
!    do j=1,ny
!        do k=1,nz
!            do l=1,8
!                write(1234,*)cube(i,j,k)%c(l,1),cube(i,j,k)%c(l,2),cube(i,j,k)%c(l,3),cube(i,j,k)%m(l)
!            enddo
!        enddo
!    enddo
!enddo
!close(1234)

end subroutine setup_mass
!===========================================================================
subroutine setup_ifc(kpt,cube,nx,ny,nz,n_b,ifc_b)

use structures
implicit none

integer, intent(in) :: nx,ny,nz,n_b(8,46,4)
type(blocks), intent(in) :: cube(nx,ny,nz)
integer :: i,j,k,l,m,n,ref,num_1,num_2,cnt,n_vec(2,46,3),tt
real(8) :: ifc(2,46,3,3),sum_f(2,3)
complex(8), intent(out) :: ifc_b(-1:1,-1:1,-1:1,24,24)
complex(8) :: ph
real(8) :: pi,rk1,rk2,rk3,m_i,m_j
real(8), intent(in) :: kpt(3)


pi=4.D0*DATAN(1.D0)

call ifc_gaas(ifc)

call read_neighbors(n_vec)

ifc_b = 0d0

sum_f = 0d0
do i=1,2
    do j=1,46
        do k=1,3  !can set the threshold (eps) below which we dont include in the acoustic sum
            sum_f(i,1) = sum_f(i,1) + ifc(i,j,1,k) !need to also include this (eps) in the matrix below as well
            sum_f(i,2) = sum_f(i,2) + ifc(i,j,2,k)
            sum_f(i,3) = sum_f(i,3) + ifc(i,j,3,k)
        enddo
    enddo
    !print*,sum_f(i,:)
enddo

do i=1,8
    tt = cube(1,1,1)%t(i)
    do j=1,3
        ifc_b(0,0,0,3*(i-1)+j,3*(i-1)+j) = -sum_f(tt,j) !diagonal terms
    enddo

    do j=1,46
        ph = dcmplx(0,1)*(kpt(1)*n_vec(tt,j,1)+kpt(2)*n_vec(tt,j,2)+kpt(3)*n_vec(tt,j,3))/(4d0)
        do m=1,3
            do n=1,3
                ifc_b(n_b(i,j,1),n_b(i,j,2),n_b(i,j,3),3*(i-1)+m,3*(n_b(i,j,4)-1)+n) = ifc(tt,j,m,n)*zexp(ph)
            enddo
        enddo
    enddo
enddo


end subroutine setup_ifc
!======================================
subroutine convert_coord(cube,ifc_b,nx,ny,nz,ifc_c,coord)

use structures
use sparsekit
implicit none


integer, intent(in) :: nx,ny,nz
integer :: i,j,k,l,m,x,y,z,box(3),cnt
type(blocks), intent(in) :: cube(nx,ny,nz)
complex(8), intent(in) :: ifc_b(-1:1,-1:1,-1:1,24,24)
integer, allocatable :: coord(:,:)
integer :: iwk(24*nx*ny*nz+1)
complex(8), allocatable :: ifc_c(:)
real(8) :: m1,m2

cnt = 0
do i=-1,1
    do j=-1,1
        do k=-1,1
            do l=1,24
                do m=1,24
                    if (ifc_b(i,j,k,l,m).ne.0d0) then
                        cnt = cnt+1
                    endif
                enddo
            enddo
        enddo
    enddo
enddo

allocate(coord(cnt*nx*ny*nz,2))
allocate(ifc_c(cnt*nx*ny*nz))


coord = 0d0
ifc_c = 0d0

cnt = 0
do x=1,nx
    do y=1,ny
        do z=1,nz
            do i=-1,1
                do j=-1,1
                    do k=-1,1
                        do l=1,24
                            do m=1,24
                                if (ifc_b(i,j,k,l,m).ne.0d0) then
                                    cnt = cnt+1
                                    m1 = cube(x,y,z)%m((l+2)/3)
                                    box(1) = modulo(x-1+i,nx)+1
                                    box(2) = modulo(y-1+j,ny)+1
                                    box(3) = modulo(z-1+k,nz)+1
                                    m2 = cube(box(1),box(2),box(3))%m((m+2)/3)
                                    ifc_c(cnt) = ifc_b(i,j,k,l,m)/sqrt(m1*m2)
                                    coord(cnt,1) = l + 24*(x-1) + 24*nx*(y-1) + 24*nx*ny*(z-1)
                                    coord(cnt,2) = m + 24*(box(1)-1) + 24*nx*(box(2)-1) + 24*nx*ny*(box(3)-1)
                                endif
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo
enddo


!subroutine coicsr(n,nnz,job,a,ja,ia,iwk)

call coicsr(24*nx*ny*nz,size(ifc_c),1,ifc_c,coord(:,2),coord(:,1),iwk)

end subroutine convert_coord
!======================================
subroutine vec_neighbor(n_imp,atom_type,atom_number)

implicit none

integer, intent(in) :: n_imp(3),atom_type
integer :: n_vec(2,46,3),i
integer, intent(out) :: atom_number

call read_neighbors(n_vec)

atom_number = 0

do i=1,46
    if (n_imp(1).eq.n_vec(atom_type,i,1)) then
        if (n_imp(2).eq.n_vec(atom_type,i,2)) then 
            if (n_imp(3).eq.n_vec(atom_type,i,3)) then
                atom_number = i
            endif
        endif
    endif
enddo

end subroutine vec_neighbor
!======================================
subroutine ifc_gaas(ifc_ga)

implicit none

integer :: basis(3),lvec(3),cnt
integer :: i,j,k,l,xi,xj,ai,aj,m1,m2,m3,mm1,mm2,mm3,a_num
real(8) :: ifc
real(8) :: ifc_r(6,6,6,2,2,3,3)
real(8), intent(out) :: ifc_ga(2,46,3,3)

ifc_r = 0d0
ifc_ga = 0d0

open(9998,file='GaAs.ifc',status='old',action='read')
do i=1,7
    read(9998,*)
enddo

basis(1) = 1
basis(2) = 1
basis(3) = 1

do i=1,7812
    if (mod(i,6**3+1).eq.1) then
        read(9998,*) xi,xj,ai,aj
    else
        read(9998,*) m1,m2,m3,ifc
        ifc_r(m1,m2,m3,ai,aj,xi,xj) = ifc

        if (m1.gt.4) then
            mm1 = m1-6
        else   
            mm1 = m1
        endif
        
        if (m2.gt.4) then
            mm2 = m2-6
        else   
            mm2 = m2
        endif
        
        if (m3.gt.4) then
            mm3 = m3-6
        else   
            mm3 = m3
        endif
    
        lvec(1) = 2*(mm1-1)+0*(mm2-1)+2*(mm3-1) - (aj-ai)*basis(1)
        lvec(2) = 0*(mm1-1)+2*(mm2-1)+2*(mm3-1) - (aj-ai)*basis(2)
        lvec(3) = 2*(mm1-1)+2*(mm2-1)+0*(mm3-1) - (aj-ai)*basis(3)

        call vec_neighbor(lvec,aj,a_num)
        

        if (a_num.ne.0) then
            if ((xi.eq.1).or.(xj.eq.1)) then
                if (xi.eq.xj) then
                    ifc_ga(aj,a_num,xi,xj) = ifc
                else
                    ifc_ga(aj,a_num,xi,xj) = -ifc
                endif
            else
                ifc_ga(aj,a_num,xi,xj) = ifc
            endif            
        endif
    endif
enddo
close(9998)

end subroutine ifc_gaas
!======================================
end module mega_block












