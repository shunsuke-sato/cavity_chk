module global_variables
  implicit none
  real(8),parameter :: pi = 4d0*atan(1d0)
  complex(8),parameter :: zi = (0d0, 1d0)

  integer :: nx, nb
  real(8) :: length_x, dx
  real(8),allocatable :: xx_l(:)
  real(8),allocatable :: psi(:,:), eps(:)
  real(8),allocatable :: vpot(:)

  real(8),parameter :: gN0 = 0d0
  real(8),parameter :: gN1 = 4.d0/5.d0
  real(8),parameter :: gN2 = -1.d0/5.d0
  real(8),parameter :: gN3 = 4.d0/105.d0
  real(8),parameter :: gN4 = -1.d0/280.d0
  real(8),parameter :: cN0 = -205.d0/72.d0
  real(8),parameter :: cN1 = 8.d0/5.d0
  real(8),parameter :: cN2 = -1.d0/5.d0
  real(8),parameter :: cN3 = 8.d0/315.d0
  real(8),parameter :: cN4 = -1.d0/560.d0


end module global_variables
!-------------------------------------------------------------------------------
program main
  use global_variables
  implicit none

  call input
  call initialize

  call calc_groundstate

end program main
!-------------------------------------------------------------------------------
subroutine input
  use global_variables
  implicit none

  length_x = 40d0
  nx = 800
  nb = 4

  dx = length_x/nx

end subroutine input
!-------------------------------------------------------------------------------
subroutine initialize
  use global_variables
  implicit none
  integer :: ix

  allocate(xx_l(0:nx))
  allocate(psi(0:nx,nb),eps(nb))
  allocate(vpot(0:nx))


  xx_l = 0d0
  do ix = 0, nx
    xx_l(ix) = dx*ix -0.5d0*length_x
  end do

!potential
  do ix = 0, nx
    vpot(ix) = -0.5d0*xx_l(ix)**2 + 0.25d0*xx_l(ix)**4
  end do

end subroutine initialize
!-------------------------------------------------------------------------------
subroutine calc_groundstate
  use global_variables
  implicit none
  integer,parameter :: Ncg = 800
  real(8),allocatable :: xvec(:),pvec(:),rvec(:)
  real(8),allocatable :: hxvec(:),gvec(:),hpvec(:)
  real(8) :: xx,pp,xp,xhx,php,xhp,esp,esp_res,gg,gg0
  real(8) :: ss,lambda,alpha,beta,aa,bb,cc
  integer :: ib, ib_t,ix,iter_cg

  allocate(xvec(0:nx),pvec(0:nx),rvec(0:nx))
  allocate(hxvec(0:nx),gvec(0:nx),hpvec(0:nx))

  do ib = 1, nb

    call random_number(xvec)
    ss=sum(xvec(:)**2)*dx
    xvec = xvec/sqrt(ss)

!== Gramm-Schmidt start
    do ib_t = 1, ib-1
      ss = sum(psi(:,ib_t)*xvec(:))*dx
      xvec = xvec - ss*psi(:,ib_t)
    end do
    ss=sum(xvec(:)**2)*dx
    xvec = xvec/sqrt(ss)
!== Gramm-Schmidt end

    call hpsi(xvec,hxvec)
    xx = sum(xvec**2)*dx
    xhx = sum(xvec*hxvec)*dx
    lambda = xhx/xx
    do iter_cg = 1, Ncg
      gvec = 2d0*(hxvec-lambda*xvec)/xx
!== Gramm-Schmidt start
      do ib_t = 1, ib-1
        ss = sum(psi(:,ib_t)*gvec(:))*dx
        gvec = gvec - ss*psi(:,ib_t)
      end do
!== Gramm-Schmidt end

      gg0 = sum(gvec**2)*dx
      if(iter_cg == 1)then
        pvec = -gvec
      else
        beta = gg0/gg
        pvec = -gvec + beta*pvec
      end if
      gg=gg0
!== Gramm-Schmidt start
      do ib_t = 1, ib-1
        ss = sum(psi(:,ib_t)*pvec(:))*dx
        pvec = pvec - ss*psi(:,ib_t)
      end do
!== Gramm-Schmidt end

      call hpsi(pvec,hpvec)

      pp  = sum(pvec**2)*dx
      php = sum(pvec*hpvec)*dx
      xp  = sum(xvec*pvec)*dx
      xhp = sum(hxvec*pvec)*dx

      aa=php*xp-xhp*pp
      bb=php*xx-xhx*pp
      cc=xhp*xx-xhx*xp
      ss=bb**2-4d0*aa*cc
      if(ss > 0d0)then
        alpha=(-bb+sqrt(ss))/(2d0*aa)
      else
        exit
      end if
     
      xvec=xvec+alpha*pvec
!== Gramm-Schmidt start
      do ib_t = 1, ib-1
        ss = sum(psi(:,ib_t)*xvec(:))*dx
        xvec = xvec - ss*psi(:,ib_t)
      end do
!== Gramm-Schmidt end
      call hpsi(xvec, hxvec)
      xx = sum(xvec**2)*dx
      xhx = sum(xvec*hxvec)*dx
      lambda = xhx/xx
      write(*,"(2I7,2x,2e26.16e3)")ib,iter_cg,lambda,sum((hxvec-lambda*xvec)**2)*dx/xx

    end do

    xvec = xvec/sqrt(xx)
    call hpsi(xvec, hxvec)
    lambda = xhx/xx
    psi(:,ib) = xvec(:)
    eps(ib) = lambda

  end do


  open(20,file="gs_wf_1d.out")
  do ix = 0, nx
    write(20,"(999e26.16e3)")xx_l(ix),vpot(ix),psi(ix,:)
  end do
  close(20)

  open(20,file="gs_single_particle_energy.out")
  do ib = 1, nb
    write(20,"(I7,2x,999e26.16e3)")ib,eps(ib)
  end do
  close(20)


end subroutine calc_groundstate
!-------------------------------------------------------------------------------
subroutine hpsi(vec, hvec)
  use global_variables
  implicit none
  real(8),intent(in) :: vec(0:nx)
  real(8),intent(out) :: hvec(0:nx)
  real(8) :: c0,c1, c2,c3,c4
  integer :: ix

  c0 = -0.5d0*cN0/dx**2
  c1 = -0.5d0*cN1/dx**2
  c2 = -0.5d0*cN2/dx**2
  c3 = -0.5d0*cN3/dx**2
  c4 = -0.5d0*cN4/dx**2

  hvec = 0d0

  ix = 0
  hvec(ix) = c0*vec(ix) &
            +c1*(vec(ix+1)) &
            +c2*(vec(ix+2)) &
            +c3*(vec(ix+3)) &
            +c4*(vec(ix+4))

  ix = 1
  hvec(ix) = c0*vec(ix) &
            +c1*(vec(ix+1)+vec(ix-1)) &
            +c2*(vec(ix+2)) &
            +c3*(vec(ix+3)) &
            +c4*(vec(ix+4))

  ix = 2
  hvec(ix) = c0*vec(ix) &
            +c1*(vec(ix+1)+vec(ix-1)) &
            +c2*(vec(ix+2)+vec(ix-2)) &
            +c3*(vec(ix+3)) &
            +c4*(vec(ix+4))

  ix = 3
  hvec(ix) = c0*vec(ix) &
            +c1*(vec(ix+1)+vec(ix-1)) &
            +c2*(vec(ix+2)+vec(ix-2)) &
            +c3*(vec(ix+3)+vec(ix-3)) &
            +c4*(vec(ix+4))


  do ix = 0+4, nx-4
    hvec(ix) = c0*vec(ix) &
              +c1*(vec(ix+1)+vec(ix-1)) &
              +c2*(vec(ix+2)+vec(ix-2)) &
              +c3*(vec(ix+3)+vec(ix-3)) &
              +c4*(vec(ix+4)+vec(ix-4))
  end do

  ix = nx -3  
  hvec(ix) = c0*vec(ix) &
            +c1*(vec(ix+1)+vec(ix-1)) &
            +c2*(vec(ix+2)+vec(ix-2)) &
            +c3*(vec(ix+3)+vec(ix-3)) &
            +c4*(vec(ix-4))

  ix = nx -2
  hvec(ix) = c0*vec(ix) &
            +c1*(vec(ix+1)+vec(ix-1)) &
            +c2*(vec(ix+2)+vec(ix-2)) &
            +c3*(vec(ix-3)) &
            +c4*(vec(ix-4))

  ix = nx -1 
  hvec(ix) = c0*vec(ix) &
            +c1*(vec(ix+1)+vec(ix-1)) &
            +c2*(vec(ix-2)) &
            +c3*(vec(ix-3)) &
            +c4*(vec(ix-4))

  ix = nx
  hvec(ix) = c0*vec(ix) &
            +c1*(vec(ix-1)) &
            +c2*(vec(ix-2)) &
            +c3*(vec(ix-3)) &
            +c4*(vec(ix-4))

  hvec = hvec + vpot*vec
  
end subroutine hpsi
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
