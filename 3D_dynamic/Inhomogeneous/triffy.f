c*************************
c PROGRAM TriffyBloc
c A QuarterSpace with SYMMETRY on three sides
c and absorbing boundaries on the three others.
c for simulation of an infinite medium
c *************************
c Stefan Nielsen September 2000
c
c PILOT
c
       include 'triffy.dec'            !declaration file
       open (unit=17,file='check',form='formatted')
       it=0
       call zero                    !initialize variables & tools
       do while(it.le.ntime)  !start time iteration

        it = it+1
        call stress                !calculate the stresses
c       call imposed_source
        call rupture
        call save               !visualize and/or store variables 
        call extrapolate           !extrapolate in time
        call bound                 !dummy boundaries
       end do                      !end of time iteration
       close (12)
       close (26)
       stop 
       end
c **********************************
c SUBROUTINES:
c **********************************
       subroutine zero
       include 'triffy.dec'
       real co,si,xpo,ypo,gam
c FIXED PARAMETERS:
       nos1=20
       savefields=1        ! 1 if want to save 
       savetraces=1        ! 1 if want to save 
       rn=9./8.            ! interpolation coefficients
       rnn=-1./24.         !  "  "
       pi=acos(-1.)
       is1=0               ! markers initialization
       is2=0               ! 
       nxs=nx/2
       nys=ny/2
       nzs=nz/2
c MEDIUM PARAMETERS - NO INPUT FILE:
       do k=1,nz
       do j=1,ny
       do i=1,nx
        mu(i,j,k)= 1.
        lam(i,j,k)=1.
        rho(i,j,k)=1.
       enddo
       enddo
       enddo
c Fault parameters- NO INPUT FILE:
      nzf=3
      nxf1=3 
      nxf2=nx-5
      nyf1=3
      nyf2=ny-5 
c     peak=1.
c     mu_s=0.
c     mu_d=0.
c     delta=2.
c     vmax=0.
c     dissipation=0.00
c READ INPUT PARAMETERS:
      open (11,file='./fpar',form='formatted',status='old')
      read (11,'(a)')
      read (11,*) dx    
      read (11,*) dt
      read (11,*) ntime
      read (11,'(a)')
      read (11,*) tw  
      read (11,*) vforce
      read (11,*) peak 
      read (11,*) ssini
      read (11,*) delta
      read (11,*) vmax
      read (11,*) dissipation
      read (11,*) RadiusInit
      read (11,*) mu_s 
      read (11,*) mu_d 
      close(11)
c
      call system('cp fpar RES/.')
      open(unit=10, file="filtered_array_fortran3.txt", action='read')
      read(10, *) strini
      close(10)
c     strini = strini/100.
c     convert vforce in fraction of shear wave velocity
      vforce=vforce*sqrt(mu(1,1,1)/rho(1,1,1))
      do j=nyf1,ny
      do i=nxf1,nx
       rdis(i,j)=0.
       broken(i,j)=0
       incrack(i,j)=1
       strini(i,j)=1.7*strini(i,j)
c      initial asperity of RadiusInit
       xpo=float(i-nxf1)*dx;ypo=float(j-nyf1)*dx
       raa=sqrt(xpo**2+ypo**2)
       if (ypo.eq.0) then
          co=1.
       else
          co=xpo/ypo
       endif
       if (xpo.eq.0) then
          si=1.
       else
          si=ypo/xpo
       endif
       gam=atan(si/co)
       open (unit=12,file='source',form='formatted')
       open (unit=26, file='iter',form='formatted')
       open (unit=38, file='./RES/frate',form='formatted')
       open (unit=37, file='./RES/ruptures',form='formatted')
       co=(cos(gam));si=(sin(gam))
       xp=float(i-nxf1);yp=float(j-nyf1);
       elrad=dx*sqrt((xp/sqrt(3.0))**2+yp**2)
       if ( (raa. lt. sqrt(co)*RadiusInit)
     & .or.(raa . lt.  sqrt(si)*RadiusInit) ) then
       !if (((nxf1-i)**2+(nyf1-j)**2)*dx.lt.(RadiusInit)**2) then
       !if (elrad.le.RadiusInit) then
            peak_31(i,j)=(s31(i,j,nzf)+strini(i,j))*.99
            broken(i,j)=1
            rdis(i,j)=0.
            nbr=nbr+1
            write (37,*) i,j,it,nbr
        !strini(i,j)=peak*1.05;else;strini(i,j)=ssini
        !peak_31(i,j)=(s31(i,j,nzf)+strini(i,j))*.99
        !broken(i,j)=1
       endif
      enddo
      enddo
      call flush(37)
c the dirac functions remove forth order in layers nzf+1 and nzf-1 
c form some terms, eventually add non zero dissipation:
      do k=1,nz
       dirac(k)=0.
      enddo
      dirac(nzf)=1. 
c EFFECTIVE MEDIUM PARAMETERS: computed in iterations to save mem.
c SPONGE BOUNDARIES:
       spong=.false.
       call init_sponge
c OUTPUT FILES, OPEN:
       return
       end
c************************************
       subroutine stress
       include 'triffy.dec'      
       real fsa !,rdis
       fsa=3.
c Sij = lam*dij*Ekk + 2*mu*Eij
c 4th order interpolation:
       do k=3,nz-1
       do j=3,ny-2
       do i=3,nx-1
        s11(i,j,k)=  s11(i,j,k) + (dt/dx) * (
     &  (lam(i,j,k)+2*mu(i,j,k)) * (
     &   rn*(v1(i,j,k)-v1(i-1,j,k))+rnn*(v1(i+1,j,k)-v1(i-2,j,k)) )
     &  +lam(i,j,k) * (
     &   rn*(v2(i,j+1,k)-v2(i,j,k))+rnn*(v2(i,j+2,k)-v2(i,j-1,k))+ 
     &   rn*(v3(i,j,k)-v3(i,j,k-1))+rnn*(v3(i,j,k+1)-v3(i,j,k-2))))
     &  +dissipation*(dirac(k)+dirac(k-1))*(v1(i,j,k)-v1(i-1,j,k))

        s22(i,j,k)= s22(i,j,k) + (dt/dx) * (
     &  (lam(i,j,k)+2*mu(i,j,k)) * (
     &   rn*(v2(i,j+1,k)-v2(i,j,k))+rnn*(v2(i,j+2,k)-v2(i,j-1,k)))
     &  +lam(i,j,k) * (
     &   rn*(v1(i,j,k)-v1(i-1,j,k))+rnn*(v1(i+1,j,k)-v1(i-2,j,k))+ 
     &   rn*(v3(i,j,k)-v3(i,j,k-1))+rnn*(v3(i,j,k+1)-v3(i,j,k-2))))
     &  +dissipation*(dirac(k)+dirac(k-1))*(v2(i,j+1,k)-v2(i,j,k))

        s33(i,j,k)= s33(i,j,k) + (dt/dx) * (
     &  (lam(i,j,k)+2*mu(i,j,k)) * (
     &   rn*(v3(i,j,k)-v3(i,j,k-1))+rnn*(v3(i,j,k+1)-v3(i,j,k-2)))
     &  +lam(i,j,k) * (
     &   rn*(v1(i,j,k)-v1(i-1,j,k))+rnn*(v1(i+1,j,k)-v1(i-2,j,k))+ 
     &   rn*(v2(i,j+1,k)-v2(i,j,k))+rnn*(v2(i,j+2,k)-v2(i,j-1,k))))
        enddo
        enddo
        enddo
c
        do k=3,nz-1  ! no der
        do j=3,ny-1
        do i=3,nx-2
        s12(i,j,k)=  s12(i,j,k) +
     &   (4./(1/mu(i,j,k)+ 1/mu(i+1,j,k)+ 
     &   1/mu(i+1,j-1,k)+ 1/mu(i,j-1,k)))*(dt/dx)*(
     &   rn*(v1(i,j,k)-v1(i,j-1,k))+rnn*(v1(i,j+1,k)-v1(i,j-2,k))+
     &   rn*(v2(i+1,j,k)-v2(i,j,k))+rnn*(v2(i+2,j,k)-v2(i-1,j,k)))
     &  +dissipation*(dirac(k)+dirac(k-1))*
     &   ((v1(i,j,k)-v1(i,j-1,k))+(v2(i+1,j,k)-v2(i,j,k)))
        enddo
        enddo
        enddo
c
        do k=3,nz-2
        do j=3,ny-1
        do i=3,nx-2 !no der
        s23(i,j,k)=  s23(i,j,k) + mu(i,j,k)*(dt/dx)*(
     &   rn*(v3(i,j,k)-v3(i,j-1,k))+rnn*(v3(i,j+1,k)-v3(i,j-2,k))+
     &   rn*(v2(i,j,k+1)-v2(i,j,k))+rnn*(v2(i,j,k+2)-v2(i,j,k-1)))
        enddo
        enddo
        enddo
c
c the dirac functions remove forth order in layers nzf+1 and nzf-1 :
        do k=3,nz-2
        do j=3,ny-2 !no der
        do i=3,nx-2
        s31(i,j,k)=s31(i,j,k)+mu(i,j,k)*(dt/dx)*(
     &   (1.-dirac(k-1))*(1.-dirac(k+1))*
     &   (rn*(v1(i,j,k+1)-v1(i,j,k))+rnn*(v1(i,j,k+2)-v1(i,j,k-1)))
     &   +(dirac(k-1)+dirac(k+1))*(v1(i,j,k+1)-v1(i,j,k))
     &   +rn*(v3(i+1,j,k)-v3(i,j,k))+rnn*(v3(i+2,j,k)-v3(i-1,j,k)))
        enddo
        enddo
        enddo
c second order interpolation:
c                  kk, 12, 23, 31
       call stressi(nx,nx-1,nx,nx-1) ! i=nx
c      call stressi(2,  1,  2,  1)   ! i=1 ! apply SYMMETRY instead:
       do k=1,nz
        do j=1,ny
         s11(nxf1  ,j,k)=-s11(nxf1+1,j,k)
         s11(nxf1-1,j,k)=-s11(nxf1+2,j,k)
         s22(nxf1  ,j,k)=-s22(nxf1+1,j,k)
         s22(nxf1-1,j,k)=-s22(nxf1+2,j,k)
         s33(nxf1  ,j,k)=-s33(nxf1+1,j,k)
         s33(nxf1-1,j,k)=-s33(nxf1+2,j,k)
         s12(nxf1-1,j,k)= s12(nxf1+1,j,k)
         s12(nxf1-2,j,k)= s12(nxf1+2,j,k)
         s23(nxf1  ,j,k)=-s23(nxf1+1,j,k)
         s23(nxf1-1,j,k)=-s23(nxf1+2,j,k)
         s31(nxf1-1,j,k)= s31(nxf1+1,j,k)
         s31(nxf1-2,j,k)= s31(nxf1+2,j,k)
        enddo
       enddo
       call stressj(ny-1,ny,ny,ny-1) ! j=ny
c      call stressj(1,  2,  2,  1)   ! j=1! apply SYMMETRY instead:
       do k=1,nz
        do i=1,nx
         s11(i,nyf1-1,k)= s11(i,nyf1+1,k)
         s11(i,nyf1-2,k)= s11(i,nyf1+2,k)
         s22(i,nyf1-1,k)= s22(i,nyf1+1,k)
         s22(i,nyf1-2,k)= s22(i,nyf1+2,k)
         s33(i,nyf1-1,k)= s33(i,nyf1+1,k)
         s33(i,nyf1-2,k)= s33(i,nyf1+2,k)
         s12(i,nyf1  ,k)=-s12(i,nyf1+1,k)
         s12(i,nyf1-1,k)=-s12(i,nyf1+2,k)
         s23(i,nyf1  ,k)=-s23(i,nyf1+1,k)
         s23(i,nyf1-1,k)=-s23(i,nyf1+2,k)
         s31(i,nyf1-1,k)= s31(i,nyf1+1,k)
         s31(i,nyf1-2,k)= s31(i,nyf1+2,k)
        enddo
       enddo
       call stressk(nz,nz,nz-1,nz-1) ! k=nz
c      call stressk(2,  2,  1,  1)   ! k=1! apply SYMMETRY instead:
       do j=1,ny
        do i=1,nx
         s11(i,j,nzf  )=-s11(i,j,nzf+1)
         s11(i,j,nzf-1)=-s11(i,j,nzf+2)
         s22(i,j,nzf  )=-s22(i,j,nzf+1)
         s22(i,j,nzf-1)=-s22(i,j,nzf+2)
         s33(i,j,nzf  )=-s33(i,j,nzf+1)
         s33(i,j,nzf-1)=-s33(i,j,nzf+2)
         s12(i,j,nzf  )=-s12(i,j,nzf+1)
         s12(i,j,nzf-1)=-s12(i,j,nzf+2)
         s23(i,j,nzf-1)= s23(i,j,nzf+1)
         s23(i,j,nzf-2)= s23(i,j,nzf+2)
         s31(i,j,nzf-1)= s31(i,j,nzf+1)
         s31(i,j,nzf-2)= s31(i,j,nzf+2)
        enddo
       enddo
       return
       end
c *******************************************
       subroutine stressi(i11,i12,i23,i31)
       include 'triffy.dec'
       i=i11
       do k=3,nz
       do j=3,ny-1
        s11(i,j,k)=  s11(i,j,k) + (dt/dx) * (
     &  (lam(i,j,k)+2*mu(i,j,k))*(v1(i,j,k)-v1(i-1,j,k))+
     &  lam(i,j,k)*(
     &  v2(i,j+1,k)-v2(i,j,k)+v3(i,j,k)-v3(i,j,k-1) ) )
        s22(i,j,k)= s22(i,j,k) + (dt/dx) * (
     &  (lam(i,j,k)+2*mu(i,j,k))*(v2(i,j+1,k)-v2(i,j,k))+
     &  lam(i,j,k)*(
     &  v1(i,j,k)-v1(i-1,j,k)+v3(i,j,k)-v3(i,j,k-1)))
        s33(i,j,k)= s33(i,j,k) + (dt/dx) * (
     &  (lam(i,j,k)+2*mu(i,j,k))*(v3(i,j,k)-v3(i,j,k-1))+
     &  lam(i,j,k)*(
     &  v1(i,j,k)-v1(i-1,j,k)+v2(i,j+1,k)-v2(i,j,k)))
       enddo
       enddo
       i=i12
       do k=3,nz
       do j=3,ny
        s12(i,j,k)=  s12(i,j,k) +
     &  (4./(1/mu(i,j,k)+ 1/mu(i+1,j,k)+ 
     &  1/mu(i+1,j-1,k)+ 1/mu(i,j-1,k)))*(dt/dx)*(
     &  v1(i,j,k)-v1(i,j-1,k)+v2(i+1,j,k)-v2(i,j,k))
       enddo
       enddo
       i=i31
       do k=3,nz-1
       do j=3,ny
        s31(i,j,k)=  s31(i,j,k) + mu(i,j,k)*(dt/dx)*(
     &  v1(i,j,k+1)-v1(i,j,k)+v3(i+1,j,k)-v3(i,j,k))
       enddo
       enddo
       return
       end
c *******************************************
       subroutine stressk(i11,i12,i23,i31)
       include 'triffy.dec'
       k=i11
       do j=3,ny-1
       do i=3,nx
        s11(i,j,k)=  s11(i,j,k) + (dt/dx) * (
     &  (lam(i,j,k)+2*mu(i,j,k))*(v1(i,j,k)-v1(i-1,j,k))+
     &  lam(i,j,k) * (
     &  v2(i,j+1,k)-v2(i,j,k)+v3(i,j,k)-v3(i,j,k-1)))
        s22(i,j,k)= s22(i,j,k) + (dt/dx) * (
     &  (lam(i,j,k)+2*mu(i,j,k))*(v2(i,j+1,k)-v2(i,j,k))+
     &  lam(i,j,k) * (
     &  v1(i,j,k)-v1(i-1,j,k)+v3(i,j,k)-v3(i,j,k-1)))
        s33(i,j,k)= s33(i,j,k) + (dt/dx) * (
     &  (lam(i,j,k)+2*mu(i,j,k))*(v3(i,j,k)-v3(i,j,k-1))+
     &  lam(i,j,k) * (
     &  v1(i,j,k)-v1(i-1,j,k)+v2(i,j+1,k)-v2(i,j,k)))
       enddo
       enddo
       k=i23
       do j=3,ny
       do i=3,nx
        s23(i,j,k)=  s23(i,j,k) + mu(i,j,k)*(dt/dx)*(
     &   v3(i,j,k)-v3(i,j-1,k)+v2(i,j,k+1)-v2(i,j,k))
       enddo
       enddo
       k=i31
       do j=3,ny
       do i=3,nx-1
        s31(i,j,k)=  s31(i,j,k) + mu(i,j,k)*(dt/dx)*(
     &  v1(i,j,k+1)-v1(i,j,k)+v3(i+1,j,k)-v3(i,j,k))
       enddo
       enddo
       return
       end
c *******************************************
       subroutine stressj(i11,i12,i23,i31)
       include 'triffy.dec'
       j=i11
       do k=3,nz
       do i=3,nx
        s11(i,j,k)=  s11(i,j,k) + (dt/dx) * (
     &  (lam(i,j,k)+2*mu(i,j,k))*(v1(i,j,k)-v1(i-1,j,k))+
     &  lam(i,j,k) * (
     &  v2(i,j+1,k)-v2(i,j,k)+v3(i,j,k)-v3(i,j,k-1)))
        s22(i,j,k)= s22(i,j,k) + (dt/dx) * (
     &  (lam(i,j,k)+2*mu(i,j,k))*(v2(i,j+1,k)-v2(i,j,k))+
     &  lam(i,j,k) * (
     &  v1(i,j,k)-v1(i-1,j,k)+v3(i,j,k)-v3(i,j,k-1)))
        s33(i,j,k)= s33(i,j,k) + (dt/dx) * (
     &  (lam(i,j,k)+2*mu(i,j,k))*(v3(i,j,k)-v3(i,j,k-1))+
     &  lam(i,j,k) * (
     &  v1(i,j,k)-v1(i-1,j,k)+v2(i,j+1,k)-v2(i,j,k)))
       enddo
       enddo
       j=i12
       do k=3,nz
       do i=3,nx-1
        s12(i,j,k)=  s12(i,j,k) +
     &  (4./(1/mu(i,j,k)+ 1/mu(i+1,j,k)+ 
     &  1/mu(i+1,j-1,k)+ 1/mu(i,j-1,k)))*(dt/dx)*(
     &   v1(i,j,k)-v1(i,j-1,k)+v2(i+1,j,k)-v2(i,j,k))
       enddo
       enddo
       j=i23
       do k=3,nz-1
       do i=3,nx
        s23(i,j,k)=  s23(i,j,k) + mu(i,j,k)*(dt/dx)*(
     &   v3(i,j,k)-v3(i,j-1,k)+v2(i,j,k+1)-v2(i,j,k))
       enddo
       enddo
       return
       end
c *******************************************
       subroutine bound
       include 'triffy.dec'
c paraxial velo. in "first" layer:
c                    pos ,shift
         call paraxi( nx , -1 )
         call paraxj( ny , -1 ) 
         call paraxk( nz , -1 )
         call parax_diag
c        call paraxi( 1  ,  1 ) ! apply SYMMETRY instead:
         do k=1,nz
          do j=1,ny
           v1(nxf1-1,j,k)= v1(nxf1+1,j,k)
           v1(nxf1-2,j,k)= v1(nxf1+2,j,k)
           v2(nxf1  ,j,k)=-v2(nxf1+1,j,k)
           v2(nxf1-1,j,k)=-v2(nxf1+2,j,k)
           v3(nxf1  ,j,k)=-v3(nxf1+1,j,k)
           v3(nxf1-1,j,k)=-v3(nxf1+2,j,k)
          enddo
         enddo
c        call paraxj( 1  ,  1 ) ! apply SYMMETRY instead:
         do k=1,nz
          do i=1,ny
           v1(i,nyf1-1,k)= v1(i,nyf1+1,k)
           v1(i,nyf1-2,k)= v1(i,nyf1+2,k)
           v2(i,nyf1  ,k)=-v2(i,nyf1+1,k)
           v2(i,nyf1-1,k)=-v2(i,nyf1+2,k)
           v3(i,nyf1-1,k)= v3(i,nyf1+1,k)
           v3(i,nyf1-2,k)= v3(i,nyf1+2,k)
          enddo
         enddo
c        call paraxk( 1  ,  1 ) ! apply SYMMETRY instead:
         do i=1,nx
          do j=1,ny
           v1(i,j,nzf  )=-v1(i,j,nzf+1)
           v1(i,j,nzf-1)=-v1(i,j,nzf+2)
           v2(i,j,nzf  )=-v2(i,j,nzf+1)
           v2(i,j,nzf-1)=-v2(i,j,nzf+2)
           v3(i,j,nzf-1)= v3(i,j,nzf+1)
           v3(i,j,nzf-2)= v3(i,j,nzf+2)
          enddo
         enddo
c diagonal paraxial velo. at "pathological" edges:
       return
       end
c **********************************
       subroutine extrapolate
       include 'triffy.dec'
       integer iu,ju,ku
c               23........
c             . .      . .  
c           .   .    .   . 
c         .     .  .     . 
c       .       ..       . STAGGERING
c     .        ..        .   OF THE 
c   .        .  3.......31    GRID
c 2.......12  .        .    
c .        ..        .           3
c .       ..       .            /
c .     .  .     .             .---> 1
c .   .    .   .               |
c . .      . .                 |
c kk.......1                   2
c                    
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      3
c      ^        |  3   31 |  3   31 |  3   31 |
c      |        |         |         |       : |
c      o-->1    | kk    1 | kk    1 | kk    1 |  NZF+1 
c...............|.........|.........|..............
c    fault>--------3---31----3---31----3---31---
c               |         |         |       : |  NZF 
c               | kk    1 | kk    1 | kk    1 |
c...............|.........|.........|.................
c               |      31 |      31 |      31 |
c               |         |         |       : |  NZF-1 
c               | kk    1 | kk    1 | kk    1 |
c...............|.........|.........|.................
c               |      31 |      31 |      31 |
c               |         |         |       : |  NZF-2 
c               | kk    1 | kk    1 | kk    1 |
c...............|.........|.........|.................
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      3
c      ^        | 23      | 23      | 23      |
c      |        |         |         |  :      | NZF+1
c      o-->1    |  2   12 |  2   12 |  2   12 |
c...............|.........|.........|..............
c    fault>-------23--------23--------23--------
c               |         |         |  :      | NZF
c               |  2   12 |  2   12 |  2   12 |
c...............|.........|.........|.................
c               | 23      | 23      | 23      |
c               |         |         |  :      | NZF-1
c               |  2   12 |  2   12 |  2   12 |
c...............|.........|.........|.................
c               | 23      | 23      | 23      |
c               |         |         |  :      | NZF-2
c               |  2   12 |  2   12 |  2   12 |
c...............|.........|.........|.................
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                              
c EXTRAPOL SPEED:
c  4th order interpolation velo. in"third" layer:
      do k=3,nz-2
      do j=3,ny-2
      do i=3,nx-2
      v1(i,j,k)=qf(i,j,k)*( v1(i,j,k)+
     &(dt/(dx*2./((1/rho(i,j,k)+1/rho(i+1,j,k)))))*(
     &rn*(s11(i+1,j,k)-s11(i,j,k))+rnn*(s11(i+2,j,k)-s11(i-1,j,k))+
     &rn*(s12(i,j+1,k)-s12(i,j,k))+rnn*(s12(i,j+2,k)-s12(i,j-1,k))+
     &rn*(s31(i,j,k)-s31(i,j,k-1))+rnn*(s31(i,j,k+1)-s11(i,j,k-2))))
      v2(i,j,k)=qf(i,j,k)*( v2(i,j,k)+
     &(dt/(dx*2./((1/rho(i,j,k)+1/rho(i,j-1,k)))))*(
     &rn*(s22(i,j,k)-s22(i,j-1,k))+rnn*(s22(i,j+1,k)-s22(i,j-2,k))+
     &rn*(s23(i,j,k)-s23(i,j,k-1))+rnn*(s23(i,j,k+1)-s23(i,j,k-2))+
     &rn*(s12(i,j,k)-s12(i-1,j,k))+rnn*(s12(i+1,j,k)-s12(i-2,j,k))))
      v3(i,j,k)=qf(i,j,k)*( v3(i,j,k)+
     &(dt/(dx*2./((1/rho(i,j,k)+1/rho(i,j,k+1)))))*(
     &rn*(s33(i,j,k+1)-s33(i,j,k))+rnn*(s33(i,j,k+2)-s33(i,j,k-1))+
     &rn*(s23(i,j+1,k)-s23(i,j,k))+rnn*(s23(i,j+2,k)-s23(i,j-1,k))+
     &rn*(s31(i,j,k)-s31(i-1,j,k))+rnn*(s31(i+1,j,k)-s31(i-2,j,k))))
      enddo
      enddo
      enddo
c 2nd order interpolation velo. in "second" layer:
      call extrapoli(nx-1)
      call extrapolj(ny-1)
      call extrapolk(nz-1)
c     call extrapoli(2) Subst. by SYMMETRY
c     call extrapolj(2) Subst. by SYMMETRY
c     call extrapolk(2) Subst. by SYMMETRY
c
      return         
      end
ccccccccccccccccccccccccccccccccccccccccc
      subroutine extrapoli(ii)
      include 'triffy.dec'
      i=ii
      do k=3,nz-1
       do j=3,ny-1
        v1(i,j,k)=qf(i,j,k)*( v1(i,j,k)+
     &  (dt/(dx*2./((1/rho(i,j,k)+1/rho(i+1,j,k)))))*(
     &  s11(i+1,j,k)-s11(i,j,k)+s12(i,j+1,k)-s12(i,j,k)+
     &  s31(i,j,k)-s31(i,j,k-1)                            ) )
        v2(i,j,k)=qf(i,j,k)*( v2(i,j,k)+
     &  (dt/(dx*2./((1/rho(i,j,k)+1/rho(i,j-1,k)))))*(
     &  s22(i,j,k)-s22(i,j-1,k)+s23(i,j,k)-s23(i,j,k-1)+
     &  s12(i,j,k)-s12(i-1,j,k)                            ) )
        v3(i,j,k)=qf(i,j,k)*( v3(i,j,k)+
     &  (dt/(dx*2./((1/rho(i,j,k)+1/rho(i,j,k+1)))))*(
     &  s33(i,j,k+1)-s33(i,j,k)+s23(i,j+1,k)-s23(i,j,k)+
     &  s31(i,j,k)-s31(i-1,j,k)                            ) )
       enddo
      enddo
      return
      end
ccccccccccccccccccccccccccccccccccccccccc
       subroutine extrapolk(ii)
       include 'triffy.dec'
       k=ii
       do j=3,ny-1
        do i=3,nx-1
         v3(i,j,k)=qf(i,j,k)*( v3(i,j,k)+
     &   (dt/(dx*2./((1/rho(i,j,k)+1/rho(i,j,k+1)))))*(
     &   s33(i,j,k+1)-s33(i,j,k)+s23(i,j+1,k)-s23(i,j,k)+
     &   s31(i,j,k)-s31(i-1,j,k)                            ) )
         v1(i,j,k)=qf(i,j,k)*( v1(i,j,k)+
     &   (dt/(dx*2./((1/rho(i,j,k)+1/rho(i+1,j,k)))))*(
     &   s11(i+1,j,k)-s11(i,j,k)+s12(i,j+1,k)-s12(i,j,k)+
     &   s31(i,j,k)-s31(i,j,k-1)                            ) )
         v2(i,j,k)=qf(i,j,k)*( v2(i,j,k)+
     &   (dt/(dx*2./((1/rho(i,j,k)+1/rho(i,j-1,k)))))*(
     &   s22(i,j,k)-s22(i,j-1,k)+s23(i,j,k)-s23(i,j,k-1)+
     &   s12(i,j,k)-s12(i-1,j,k)                            ) )
        enddo
       enddo
       return
       end
ccccccccccccccccccccccccccccccccccccccccc
       subroutine extrapolj(ii)
       include 'triffy.dec'
       j=ii
       do k=3,nz-1
        do i=3,nx-1
         v1(i,j,k)=qf(i,j,k)*( v1(i,j,k)+
     &   (dt/(dx*2./((1/rho(i,j,k)+1/rho(i+1,j,k)))))*(
     &   s11(i+1,j,k)-s11(i,j,k)+s12(i,j+1,k)-s12(i,j,k)+
     &   s31(i,j,k)-s31(i,j,k-1)                            ) )
         v3(i,j,k)=qf(i,j,k)*( v3(i,j,k)+
     &   (dt/(dx*2./((1/rho(i,j,k)+1/rho(i,j,k+1)))))*(
     &   s33(i,j,k+1)-s33(i,j,k)+s23(i,j+1,k)-s23(i,j,k)+
     &   s31(i,j,k)-s31(i-1,j,k)                            ) )
         v2(i,j,k)=qf(i,j,k)*( v2(i,j,k)+
     &   (dt/(dx*2./((1/rho(i,j,k)+1/rho(i,j-1,k)))))*(
     &   s22(i,j,k)-s22(i,j-1,k)+s23(i,j,k)-s23(i,j,k-1)+
     &   s12(i,j,k)-s12(i-1,j,k)                            ) )
        enddo
       enddo
       return
       end
ccccccccccccccccccccccccccccccccccccccccc
       subroutine save
       include 'triffy.dec'
       real bluff(nx/5,ny/5)
       real temp(nx,ny)
       character*15 namp
       character*15 namt
       character*15 nams
       real*8 r8
       integer*4 it4
       integer*4 ip4
       integer kmi
c      here files trace* are replenished with vertical
c      displacement at the surface. v2 can be replaced by v1
c      or any combination of the 2 components. The convention
c      here is that the velocity is the interpolated value at
c      node s11/s22, so that all velocity components are
c      known at the same spot.
c
       is1=is1+1
       is2=is2+1
cccccc
c%%% save fields for snapshots:
       if (savefields.eq.1.and.is1.ge.nos1) then
        backspace (26)
        write (26,*) it
        call flush(26)
         print*,it
         is1=0
         numf=numf+1
         write (namp,'(a,I5.5)') './RES/fldx',it
         write (nams,'(a,I5.5)') './RES/fldy',it
         write (namt,'(a,I5.5)') './RES/fldz',it
cc field1:
         open (51,file=namp,form='formatted')
         open (52,file=nams,form='formatted')
         open (53,file=namt,form='formatted')
         do k=3,nz
          write(51,79) (v1(nx/2,j,k),j=1,ny)
         enddo
         do i=3,nx
          write(52,79) (v1(i,ny/2,k),k=1,nz)
         enddo
         do j=3,ny
          write(53,79) (v1(i,j,nzf),i=1,nx)
         enddo
         close (51)
         close (52)
         close (53)
ccccc
       endif
78     continue
79     format(1000G17.9)
       return
       end 
ccccccccccccccccccccccccccccccccccccccccccc
       subroutine init_sponge
       include 'triffy.dec'
       integer TT
       real sig
       TT=10  !Thickness of sponge
       sig=sqrt( -float(TT)**2 / log( .95) )
c
       do k=1,nz
        do j=1,ny
         do i=1,nx
          qf(i,j,k)=1.
         enddo
        enddo
       enddo
c
       if (spong) then 
c
       do i=1,nx
        do k=1,nz
         do j=ny-TT,ny !(bottom)
c          qf(i,j,k)=exp(sig*float(ny-j-TT))
c          qf(i,j,k)=1.-.2*float(j-ny+TT)/float(TT)
           qf(i,j,k)=qf(i,j,k)*exp(-((float(ny-j-TT)/sig)**2))
         enddo
         do j=1,TT !(top)
c          qf(i,j,k)=exp(sig*float(j-TT))
c          qf(i,j,k)=.8-.2*float(j)/float(TT)
           qf(i,j,k)=qf(i,j,k)*exp(-((float(j-TT)/sig)**2))
         enddo
        enddo
       enddo
       do j=1,ny
        do k=1,nz
         do i=1,TT     !(left side)
c          qf(i,j,k)=exp(sig*float(i-TT))
c          qf(i,j,k)=.8-.2*float(i)/float(TT)
           qf(i,j,k)=qf(i,j,k)*exp(-((float(i-TT)/sig)**2))
         enddo
         do i=nx-TT,nx !(right side)
c          qf(i,j)=exp(sig*float(nx-i-TT))
c          qf(i,j)=1.-.2*float(i-nx+TT)/float(TT)
           qf(i,j,k)=qf(i,j,k)*exp(-((float(nx-i-TT)/sig)**2))
         enddo
        enddo
       enddo
       do j=1,ny
        do i=1,nx
         do k=1,TT     !(front side)
c          qf(i,j,k)=exp(sig*float(k-TT))
c          qf(i,j,k)=.8-.2*float(k)/float(TT)
           qf(i,j,k)=qf(i,j,k)*exp(-((float(k-TT)/sig)**2))
         enddo
         do k=nz-TT,nz !(back side)
c          qf(i,j,k)=exp(sig*float(nz-k-TT))
c          qf(i,j,k)=1.-.2*float(k-nz+TT)/float(TT)
           qf(i,j,k)=qf(i,j,k)*exp(-((float(nz-k-TT)/sig)**2))
         enddo
        enddo
       enddo
c
       endif
c
       return
       end
ccccccccccccccccccccccccccccccccccccccccccccc
       subroutine paraxi(ipos,ish)
       include 'triffy.dec'
       real Vp,Vs,vtemp(nx,ny,nz)
       integer ish,ipos,isi
       i=ipos
       do k=3,nz-1
       do j=3,ny-1
        Vp=sqrt((2*mu(i,j,k)+lam(i,j,k))/rho(i,j,k))
        Vs=sqrt(mu(i,j,k)/rho(i,j,k))
        vtemp(i,j,k)=v1(i,j,k)-
     &   (v1(i,j,k)-v1(i+ish,j,k))*Vp*dt/dx
     &   +(v2(i+0,j+1,k)-v2(i+0,j,k))*(Vp-Vs)*dt/dx
     &   +(v3(i+0,j,k)-v3(i+0,j,k-1))*(Vp-Vs)*dt/dx
       enddo
       do j=3,ny-1
        v2(i,j,k)=v2(i,j,k)-
     &  (v2(i,j,k)-v2(i+ish,j,k))*Vs*dt/dx 
c    &  +(v1(i+ish,j,k)-v1(i+ish,j-1,k))*(Vp-Vs)*dt/dx 
       enddo
       do j=3,ny-1
        v3(i,j,k)=v3(i,j,k)-
     &  (v3(i,j,k)-v3(i+ish,j,k))*Vs*dt/dx 
c    &  +(v1(i+ish,j,k+1)-v1(i+ish,j,k))*(Vp-Vs)*dt/dx 
       enddo
       do j=3,ny-1
        v1(i,j,k)=vtemp(i,j,k)
       enddo
       enddo
       return
       end
ccccccccccccccccccccccccccccccccccccccc
       subroutine paraxk(ipos,ish)
       include 'triffy.dec'
       real Vp,Vs,vtemp(nx,ny,nz)
       integer ish,ipos,isi
       k=ipos
       do i=3,nx-1
       do j=3,ny-1
        Vp=sqrt((2*mu(i,j,k)+lam(i,j,k))/rho(i,j,k))
        Vs=sqrt(mu(i,j,k)/rho(i,j,k))
        vtemp(i,j,k)=v3(i,j,k)-
     &   (v3(i,j,k)-v3(i,j,k+ish))*Vp*dt/dx
     &   +(v2(i,j+1,k+0)-v2(i,j,k+0))*(Vp-Vs)*dt/dx
     &   +(v1(i,j,k+0)-v1(i-1,j,k+0))*(Vp-Vs)*dt/dx
       enddo
       do j=3,ny-1
        v2(i,j,k)=v2(i,j,k)-
     &  (v2(i,j,k)-v2(i,j,k+ish))*Vs*dt/dx
c    &  +(v3(i,j,k+ish)-v3(i,j-1,k+ish))*(Vp-Vs)*dt/dx 
       enddo
       do j=3,ny-1
        v1(i,j,k)=v1(i,j,k)-
     &  (v1(i,j,k)-v1(i,j,k+ish))*Vs*dt/dx 
c    &  +(v3(i+1,j,k+ish)-v3(i,j,k+ish))*(Vp-Vs)*dt/dx 
       enddo
       do j=3,ny-1
        v3(i,j,k)=vtemp(i,j,k)
       enddo
       enddo
       return
       end
ccccccccccccccccccccccccccccccccccccccc
       subroutine paraxj(ipos,ish)
       include 'triffy.dec'
       real Vp,Vs,vtemp(nx,ny,nz)
       integer ish,ipos,isi
       j=ipos
       do i=3,nx-1
       do k=3,nz-1
        Vp=sqrt((2*mu(i,j,k)+lam(i,j,k))/rho(i,j,k))
        Vs=sqrt(mu(i,j,k)/rho(i,j,k))
        vtemp(i,j,k)=v2(i,j,k)-
     &   (v2(i,j,k)-v2(i,j+ish,k))*Vp*dt/dx
     &   +(v3(i,j-1  ,k)-v3(i,j-1,k-1))*(Vp-Vs)*dt/dx
     &   +(v1(i,j-1  ,k)-v1(i-1,j-1,k))*(Vp-Vs)*dt/dx
       enddo
       do k=3,nz-1
        v1(i,j,k)=v1(i,j,k)-
     &  (v1(i,j,k)-v1(i,j+ish,k))*Vs*dt/dx
c    &  +(v2(i+1,j+ish,k)-v2(i,j+ish,k))*(Vp-Vs)*dt/dx 
       enddo
       do k=3,nz-1
        v3(i,j,k)=v3(i,j,k)-
     &  (v3(i,j,k)-v3(i,j+ish,k))*Vs*dt/dx 
c    &  +(v2(i,j+ish,k+1)-v2(i,j+ish,k))*(Vp-Vs)*dt/dx
       enddo
       do k=3,nz-1
        v2(i,j,k)=vtemp(i,j,k)
       enddo
       enddo
       return
       end
ccccccccccccccccccccccccccccccccccccccc
       subroutine parax_diag
       include 'triffy.dec'
       real V,dtx
c diagonal paraxial approx for 9 out of 12 edges
c where velocity is in the inside of the corner
c      sqrt(2.)/2. diagonal incidence
c      dxx=sqrt(2.)*dx diagonal step
c      * radiation at 45^o
c      * diagonal step of sqrt(dx^2+dx^2)
c      * mean of P and S influence
c             => V=V/4

       dtx=dt/(sqrt(2.)*dx)
c edges xy:
c      i=1 ; j=1
c      do k=1,nz
c      intermediate S/P velocity:
c      V=(sqrt((2*mu(i,j,k)+lam(i,j,k))/rho(i,j,k))
c    &  +sqrt(mu(i,j,k)/rho(i,j,k)))/2.
c       v1(i,j,k)=v1(i,j,k)-(v1(i,j,k)-v1(i+1,j+1,k))*V*dtx
c       v2(i,j,k)=v2(i,j,k)-(v2(i,j,k)-v2(i+1,j+1,k))*V*dtx
c       v3(i,j,k)=v3(i,j,k)-(v3(i,j,k)-v3(i+1,j+1,k))*V*dtx
c      enddo
c      i=1 ; j=ny
c      do k=1,nz
c      V=(sqrt((2*mu(i,j,k)+lam(i,j,k))/rho(i,j,k))
c    &  +sqrt(mu(i,j,k)/rho(i,j,k)))/2.
c       v1(i,j,k)=v1(i,j,k)-(v1(i,j,k)-v1(i+1,j-1,k))*V*dtx
c       v2(i,j,k)=v2(i,j,k)-(v2(i,j,k)-v2(i+1,j-1,k))*V*dtx
c       v3(i,j,k)=v3(i,j,k)-(v3(i,j,k)-v3(i+1,j-1,k))*V*dtx
c      enddo
c      i=nx ; j=1   
c      do k=1,nz
c      V=(sqrt((2*mu(i,j,k)+lam(i,j,k))/rho(i,j,k))
c    &  +sqrt(mu(i,j,k)/rho(i,j,k)))/2.
c       v1(i,j,k)=v1(i,j,k)-(v1(i,j,k)-v1(i-1,j+1,k))*V*dtx
c       v2(i,j,k)=v2(i,j,k)-(v2(i,j,k)-v2(i-1,j+1,k))*V*dtx
c       v3(i,j,k)=v3(i,j,k)-(v3(i,j,k)-v3(i-1,j+1,k))*V*dtx
c      enddo
       i=nx ; j=ny   
       do k=3,nz
       V=(sqrt((2*mu(i,j,k)+lam(i,j,k))/rho(i,j,k))
     &  +sqrt(mu(i,j,k)/rho(i,j,k)))/2.
        v1(i,j,k)=v1(i,j,k)-(v1(i,j,k)-v1(i-1,j-1,k))*V*dtx
        v2(i,j,k)=v2(i,j,k)-(v2(i,j,k)-v2(i-1,j-1,k))*V*dtx
        v3(i,j,k)=v3(i,j,k)-(v3(i,j,k)-v3(i-1,j-1,k))*V*dtx
       enddo
c Edges yz:
c      j=1 ; k=1 
c      do i=1,nx
c       V=(sqrt((2*mu(i,j,k)+lam(i,j,k))/rho(i,j,k))
c    &   +sqrt(mu(i,j,k)/rho(i,j,k)))/2.
c       v1(i,j,k)=v1(i,j,k)-(v1(i,j,k)-v1(i,j+1,k+1))*V*dtx
c       v2(i,j,k)=v2(i,j,k)-(v2(i,j,k)-v2(i,j+1,k+1))*V*dtx
c       v3(i,j,k)=v3(i,j,k)-(v3(i,j,k)-v3(i,j+1,k+1))*V*dtx
c      enddo
       j=ny ; k=nz   
       do i=3,nx
        V=(sqrt((2*mu(i,j,k)+lam(i,j,k))/rho(i,j,k))
     &   +sqrt(mu(i,j,k)/rho(i,j,k)))/2.
        v1(i,j,k)=v1(i,j,k)-(v1(i,j,k)-v1(i,j-1,k-1))*V*dtx
        v2(i,j,k)=v2(i,j,k)-(v2(i,j,k)-v2(i,j-1,k-1))*V*dtx
        v3(i,j,k)=v3(i,j,k)-(v3(i,j,k)-v3(i,j-1,k-1))*V*dtx
       enddo
c      j=1 ; k=nz 
c      do i=1,nx
c       V=(sqrt((2*mu(i,j,k)+lam(i,j,k))/rho(i,j,k))
c    &   +sqrt(mu(i,j,k)/rho(i,j,k)))/2.
c       v1(i,j,k)=v1(i,j,k)-(v1(i,j,k)-v1(i,j+1,k-1))*V*dtx
c       v2(i,j,k)=v2(i,j,k)-(v2(i,j,k)-v2(i,j+1,k-1))*V*dtx
c       v3(i,j,k)=v3(i,j,k)-(v3(i,j,k)-v3(i,j+1,k-1))*V*dtx
c      enddo
c      j=ny ; k=1
c      do i=1,nx
c       V=(sqrt((2*mu(i,j,k)+lam(i,j,k))/rho(i,j,k))
c    &   +sqrt(mu(i,j,k)/rho(i,j,k)))/2.
c       v1(i,j,k)=v1(i,j,k)-(v1(i,j,k)-v1(i,j-1,k+1))*V*dtx
c       v2(i,j,k)=v2(i,j,k)-(v2(i,j,k)-v2(i,j-1,k+1))*V*dtx
c       v3(i,j,k)=v3(i,j,k)-(v3(i,j,k)-v3(i,j-1,k+1))*V*dtx
c      enddo
c Edges zx:
c      i=1 ; k=nz
c      do j=1,ny
c       V=(sqrt((2*mu(i,j,k)+lam(i,j,k))/rho(i,j,k))
c    &   +sqrt(mu(i,j,k)/rho(i,j,k)))/2.
c       v1(i,j,k)=v1(i,j,k)-(v1(i,j,k)-v1(i+1,j,k-1))*V*dtx
c       v2(i,j,k)=v2(i,j,k)-(v2(i,j,k)-v2(i+1,j,k-1))*V*dtx
c       v3(i,j,k)=v3(i,j,k)-(v3(i,j,k)-v3(i+1,j,k-1))*V*dtx
c      enddo
       i=nx ; k=nz
       do j=3,ny
        V=(sqrt((2*mu(i,j,k)+lam(i,j,k))/rho(i,j,k))
     &   +sqrt(mu(i,j,k)/rho(i,j,k)))/2.
        v1(i,j,k)=v1(i,j,k)-(v1(i,j,k)-v1(i-1,j,k-1))*V*dtx
        v2(i,j,k)=v2(i,j,k)-(v2(i,j,k)-v2(i-1,j,k-1))*V*dtx
        v3(i,j,k)=v3(i,j,k)-(v3(i,j,k)-v3(i-1,j,k-1))*V*dtx
       enddo
c      i=nx ; k=1
c      do j=1,ny
c       V=(sqrt((2*mu(i,j,k)+lam(i,j,k))/rho(i,j,k))
c    &   +sqrt(mu(i,j,k)/rho(i,j,k)))/2.
c       v1(i,j,k)=v1(i,j,k)-(v1(i,j,k)-v1(i-1,j,k+1))*V*dtx
c       v2(i,j,k)=v2(i,j,k)-(v2(i,j,k)-v2(i-1,j,k+1))*V*dtx
c       v3(i,j,k)=v3(i,j,k)-(v3(i,j,k)-v3(i-1,j,k+1))*V*dtx
c      enddo
c      i=1 ; k=1
c      do j=1,ny
c       V=(sqrt((2*mu(i,j,k)+lam(i,j,k))/rho(i,j,k))
c    &   +sqrt(mu(i,j,k)/rho(i,j,k)))/2.
c       v1(i,j,k)=v1(i,j,k)-(v1(i,j,k)-v1(i+1,j,k+1))*V*dtx
c       v2(i,j,k)=v2(i,j,k)-(v2(i,j,k)-v2(i+1,j,k+1))*V*dtx
c       v3(i,j,k)=v3(i,j,k)-(v3(i,j,k)-v3(i+1,j,k+1))*V*dtx
c      enddo
c
       return
       end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rupture
      include 'triffy.dec'
      real friction,dd,dv,rns,rad,elrad,ftw,xp,yp
      k=nzf
ccc   check for rupture:
       do i=nxf1,nxf2
        do j=nyf1,nyf2
          rns=( s33(i,j,k+1)+s33(i,j,k) ) / 2.
          rdis(i,j)=v1(i,j,nzf+1)-v1(i,j,nzf)
          xp=float(i-nxf1);yp=float(j-nyf1);
          elrad=dx*sqrt((xp/sqrt(3.0))**2+yp**2)
          rad=dx*sqrt(xp**2+yp**2)
          if(
     &       (
     &        (s31(i,j,nzf)+strini(i,j).ge.peak+rns*mu_s)!.or.
     &       !(elrad.le.(float(it)*dt*vforce+RadiusInit))
     &       ).and.
     &        broken(i,j).eq.0.and.incrack(i,j).eq.1
     &      ) then
            peak_31(i,j)=(s31(i,j,nzf)+strini(i,j))*.99
            broken(i,j)=1
            rdis(i,j)=0.
            nbr=nbr+1
            write (37,*) i,j,it,nbr
          endif
ccc update friction:
          if (broken(i,j).eq.1) then
            if(delta.eq.0.)then
              dd=0.
            else
c             dd=delta/(delta+slip(i,j))
              dd=1. - slip(i,j)/delta
              if (dd.lt.0.) dd=0.
            endif
            if (vmax.eq.0.) then
              dv=0.
            else
              dv=vmax/(vmax+rdis(i,j))
            endif
            !if (dv.gt.dd) dd=dv
ccc friction with normal stress dependence:
c           friction=mu_d*(rns+peak/mu_s)*(1.-dd)+dd*peak_31(i,j) 
c           friction=dd*peak_31(i,j)+(1.-dd)*(mu_d*rns)
            friction=dd*peak_31(i,j) 
            !ftw=peak_31(i,j)*(1.-(float(it)*dt-rad/vforce))/tw
            !if (ftw.lt.0.) ftw=0.
            !if(friction.gt.ftw) friction=ftw
c           friction=peak_31(i,j)*exp(-(float(it)/50.)**2)
c           rfa=sqrt(1./(rho(i,j,k)*mu(i,j,k)))
c           rdis(i,j)=2.*rfa*(s31(i,j,k)+strini(i,j)-friction)
            s31(i,j,k)=friction-strini(i,j)
ccc  healing criterion:
c           if (rdis(i,j).lt.0.) then
c             broken(i,j)=0
c             rdis(i,j)=0.
c             nbr=nbr-1
c             strini(i,j)=0.
c             write(39,*) i,j,it
c           endif
          else
            rdis(i,j)=0.
          endif
          slip(i,j)=slip(i,j)+rdis(i,j)*dt
        enddo
       enddo
       write (38,'(1000G17.9)') (rdis(i,3),i=3,nxf2)
       call flush (38)
       call flush (39)
       call flush (37)
       return
       end
cccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine imposed_source
      include 'triffy.dec'
c STRESS SOURCE:
      rf=0.
      if (float(it)*dt.le.6*tw) then
        rf=100000.*exp(-((float(it)*dt-2.*tw)/tw)**2)
        s11(nxs,nys,nzs)=s11(nxs,nys,nzs)+rf
        s22(nxs,nys,nzs)=s22(nxs,nys,nzs)+rf
        s33(nxs,nys,nzs)=s33(nxs,nys,nzs)+rf
        write (12,*) it,nxs,nys,nzs,rf
        call flush(12)
      endif
      return
      end
ccccccccccccccccccccccccccccccccc
