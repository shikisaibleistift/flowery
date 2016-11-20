!
  subroutine calcpsmacdirect_ik_openmpi2
!
  use modulepara
  use moduleall
  use moduleu
  use modulev
  use modulew
  use moduledirectp_openmpi
  use modulevar
  use modulegeom
  use moduleflupr
  use modulecoef
  use moduletimemarching
  use modulemoni
  use moduleregion
  use moduleperiod
  use modulested
  use modulempi
  use moduleparamaster
  use modulempivar
!
  implicit none
!
  include "mpif.h"
!
  double precision,dimension(nio,nj,nk)   :: sppo
  double precision,dimension(nio,nj,nkkk) :: sppoo
  double precision,dimension(nio,nj,nkkk) :: ppo
  double precision                        :: pu,pv,pw
  double precision                        :: resmax,mbal,resmaxx,mbalt
  character(len=72)                       :: text
  integer                                 :: iflg  ,kflg !1:Neumann -1:Dirichlet
  parameter(                                 iflg=1,kflg=1)
!
!
!
! calculate coefficients 
!
  do 100 k = kst, ked
    do 100 j = jst, jed
      do 100 i = ist, ied
        sppoo(i,j,k)=                      &
         (up(i+1,j,k) - up(i,j,k))/sew(i)  &
        +(vp(i,j+1,k) - vp(i,j,k))/sns(j)  &
        +(wp(i,j,k+1) - wp(i,j,k))/stb(k)
100 continue
!
!
!
  call flipflop_mpi(            sppoo,iflg )
  call aggregate(               sppoo,kflg )
  call directpoisson_ik_openmpi(ppo  ,sppoo)
  call scatter(                 pp   ,ppo  )
!
!
!
  call jslipbc4j4(u,v,w,pp,2)
!
  if(myrank.eq.0)then
    write(*,*)'Property________',pp(imon,jmon,kmon)
  endif
!
!
!
  do 500 k = kst, ked
    do 500 j = jst, jed
      do 500 i = ist, ied
        p(i,j,k) = p(i,j,k)+densit*pp(i,j,k)/deltat
500 continue
!
  call sendrcvedgempi4(p ,ni,nj,nk)
  call sendrcvedgempi4(pp,ni,nj,nk)
!
  do 510 k = kst, ked
    do 510 j = jst, jed
      do 510 i = ist, ied
        pu=(pp(i-1,j,k)-pp(i,j,k))/sewu(i)
        pv=(pp(i,j-1,k)-pp(i,j,k))/snsv(j)
        pw=(pp(i,j,k-1)-pp(i,j,k))/stbw(k)
!
        u(i,j,k) = up(i,j,k) + pu
        v(i,j,k) = vp(i,j,k) + pv
        w(i,j,k) = wp(i,j,k) + pw
510 continue
!
  call sendrcvedgempi4(u,ni,nj,nk)
  call sendrcvedgempi4(v,ni,nj,nk)
  call sendrcvedgempi4(w,ni,nj,nk)
!
!
!
  resmax =0.d0
  mbal   =0.d0
  resmaxx=0.d0
  mbalt  =0.d0
!
  do k=kst,ked,1
    do j=jst,jed,1
      do i=ist,ied,1
        cellres(i,j,k)=                &
        abs(                           &
         (u(i+1,j,k)-u(i,j,k))/sew(i)  &
        +(v(i,j+1,k)-v(i,j,k))/sns(j)  &
        +(w(i,j,k+1)-w(i,j,k))/stb(k)  &
        )
!
        mbal=mbal+cellres(i,j,k)
        if(cellres(i,j,k).ge.resmax)then
          resmax=cellres(i,j,k)
          imax=i
          jmax=j
          kmax=k
        endif
!
      enddo
    enddo
  enddo
!
  call mpi_allreduce                           &
       (resmax,resmaxx,1,mpi_double_precision  &
       ,mpi_max,mpi_comm_world,ierr)
!
  call mpi_allreduce                       &
       (mbal,mbalt,1,mpi_double_precision  &
       ,mpi_sum,mpi_comm_world,ierr)
!
  if(myrank.eq.0)then
    write(*,*)"------------"
    write(*,*)"mass balance"
    write(*,*)mbalt/((nimaster-2*bff)*(njmaster-2*bff)*(nkmaster-2*bff))
    write(*,*)"max residual"
    write(*,*)resmaxx
  endif
!
  if(resmaxx.eq.resmax)write(*,*)"at ",myrank,imax,jmax,kmax
!
!
!
! check
!
  call integratempi4(pp,interu,ni,nj,nk,nimaster,njmaster,nkmaster)
  if(myrank.eq.0)then
    open(8,file="./output/pp.csv")
    do i=1,nimaster
      write(8,*)(interu(i,nj/2,k),k=1,nkmaster,1)
    enddo
    close(8)
  endif
!
!
!
  call integratempi4(u,interu,ni,nj,nk,nimaster,njmaster,nkmaster)
  if(myrank.eq.0)then
    open(8,file="./output/u.csv")
    do i=1,nimaster
      write(8,111)(interu(i,j,nkmaster/2),j=1,njmaster,1)
111 format(70(f15.10))
    enddo
    close(8)
  endif
!
!
!
  return
!
  end
!
