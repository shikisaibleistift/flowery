!
  subroutine propscsmmpi2
!
  use modulepara
  use moduleall
  use modulevar
  use modulegeom
  use moduleflupr
  use modulemoni
  use moduleprops
  use modulewalldpf
  use moduleperiod
  use moduleparamaster
  use modulempi
  use modulempivar
!
!
!
  implicit none
!
  double precision :: o5w2,o5s2
  double precision :: ccsm
!
!
!
  do 200 k=2, nkm1
    do 200 j=2, njm1
      do 200 i=2, nim1
        dudx = (u(i+1,j,k)-u(i,j,k))*mviu(i,j,k) / sew(i)

        dudy =                                                & 
        (                                                     &
         (u(i,j,k) + u(i+1,j,k) + u(i,j+1,k) + u(i+1,j+1,k))  &
        -(u(i,j,k) + u(i+1,j,k) + u(i,j-1,k) + u(i+1,j-1,k))  &
        )/sns(j)
!
        dudz =                                                &
        (                                                     &
         (u(i,j,k) + u(i+1,j,k) + u(i,j,k+1) + u(i+1,j,k+1))  &
        -(u(i,j,k) + u(i+1,j,k) + u(i,j,k-1) + u(i+1,j,k-1))  &
        )/stb(k)
!
        dvdx =                                                &
        (                                                     &
         (v(i,j,k) + v(i,j+1,k) + v(i+1,j,k) + v(i+1,j+1,k))  &
        -(v(i,j,k) + v(i,j+1,k) + v(i-1,j,k) + v(i-1,j+1,k))  &
        )/sew(i)
!
        dvdy = (v(i,j+1,k)-v(i,j,k))*mviv(i,j,k) / sns(j)

        dvdz =                                                &
        (                                                     &
         (v(i,j,k) + v(i,j+1,k) + v(i,j,k+1) + v(i,j+1,k+1))  &
        -(v(i,j,k) + v(i,j+1,k) + v(i,j,k-1) + v(i,j+1,k-1))  &
        )/stb(k)
!
        dwdx =                                                &
        (                                                     &
         (w(i,j,k) + w(i,j,k+1) + w(i+1,j,k) + w(i+1,j,k+1))  &
        -(w(i,j,k) + w(i,j,k+1) + w(i-1,j,k) + w(i-1,j,k+1))  &
        )/sew(i)
!
        dwdy =                                                &
        (                                                     &
         (w(i,j,k) + w(i,j,k+1) + w(i,j+1,k) + w(i,j+1,k+1))  &
        -(w(i,j,k) + w(i,j,k+1) + w(i,j-1,k) + w(i,j-1,k+1))  &
        )/sns(j)
!
        dwdz = (w(i,j,k+1) - w(i,j,k))*mviw(i,j,k) / stb(k)
!
!
!
        o5w2= 0.5d0*(                       &
        2.d0*(0.5d0*(dudy - dvdx))**2.d0 +  &
        2.d0*(0.5d0*(dudz - dwdx))**2.d0 +  &
        2.d0*(0.5d0*(dvdz - dwdy))**2.d0)
!
        o5s2= 0.5d0*(                       &
        (0.5d0*(dudx+dudx))**2.d0+          &
        (0.5d0*(dvdy+dvdy))**2.d0+          &
        (0.5d0*(dwdz+dwdz))**2.d0+          &
        2.d0*(0.5d0*(dudy + dvdx))**2.d0 +  &
        2.d0*(0.5d0*(dudz + dwdx))**2.d0 +  &
        2.d0*(0.5d0*(dvdz + dwdy))**2.d0)
!
!
!
        ccsm=                                     &
        (1.d0/20.d0*                              &
        (abs(((o5w2-o5s2)/(o5w2+o5s2)))**1.5d0))
!
        if(o5w2+o5s2.eq.0)ccsm=0.d0
!
        nut(i,j,k) =                             &
        ccsm*                                    &
        ((sew(i)*sns(j)*stb(k))**(2.d0/3.d0)) *  &
        sqrt(                                    &
         2.d0*dudx**2.d0    +                    &
         2.d0*dvdy**2.d0    +                    &
         2.d0*dwdz**2.d0    +                    &
        (dudy + dvdx)**2.d0 +                    &
        (dudz + dwdx)**2.d0 +                    &
        (dvdz + dwdy)**2.d0                      &
        )
!
200 continue
!
!
!
  do 100 i = 2, nim1
    do 100 k = 2, nkm1
      do 100 j = 2, njm1
        vistl(i,j,k) = nut(i,j,k) * densit
        vis(i,j,k)   = vistl(i,j,k)+viscos
100 continue
!
  call sendrcvedgempi4(vis,ni,nj,nk)
!
!
!
  return
!
  end

