c ******************
c PROGRAM fault2dPML
c ******************
c 2-d rectangle bounded by PML with fault line
c Gaetano Festa - Stefan Nielsen - January 2002
c
       include 'arrays.f'            !declaration file
c Timing definition
!      integer (kind=8) it1,it2,irtc_rate,rate
!      external irtc
!      rate=irtc_rate()
!      call irtc(it1)
c
       open (unit=25, file='iter',form='formatted')
       write (25,*) 'initializing...'
       call flush (25)
       call zero
       it=0
       do while(it.le.ntime)  !start time iteration
         it = it+1
         call stress                !calculate the stresses
         call rupture
c        call source
         call extrapolate           !extrapolate in time
         call save                  !visualize and/or store variables 
       end do                       !end of time iteration
       close (25)
!      call irtc(it2)
!      sec=float((it2-it1))/float(rate)
!      print *,'Time elapsed',sec
       end
c **********************************
c SUBROUTINES:
c **********************************
       subroutine zero
       include 'arrays.f'
c GENERAL PARAMETERS:
       rn=9./8.            ! interpolation coefficients
       rnn=-1./24.         !  "  "
       pi=acos(-1.)
       is1=0               ! markers initialization
       is2=0               ! 
       l1=3.               ! width of space smoot. funct. (in dx)
c MEDIUM PARAMETERS - INPUT FILE:
c here read model parameters if desired:
c (elastic parmeters and density and quality)
c
       open(11,file='fpar',form='formatted',status='old')
       read(11,'(a)')
       read(11,*) dx
       read(11,*) dt
       read(11,*) ntime
       read(11,*) peak
       read(11,*) ssini
       read(11,'(a)') FricType
       read(11,*) delta, vcha, vmax, tau_f
       read(11,*) rvr,rvh
       read(11,*) tau_s
       read(11,*) InitAspLen
       read(11,*) mu_s
       read(11,*) mu_d
       read(11,*) nxs, nys
       read(11,*) nos1
       read(11,*) nos2
       read(11,*) diss
       close(11)
       call fillmodel
c
       open (unit=31,file='./RES/trace',form='formatted')
c fault parameters init:
       do i=1,nx
         strini(i)=ssini
         broken(i)=0
         slip(i)=0.
         rdis(i)=0.
         peak_12(i)=0.
       enddo
       do j=1,ny
        fault(j)=0.
       enddo
       fault(nys)=1.
       print*,'nx,ny,nxs,nys,InitAspLen',nx,ny,nxs,nys,InitAspLen
       np1=nxs-int(InitAspLen/2);np2=nxs+int(InitAspLen/2)+1
       epoints(1)=np1;epoints(2)=np2
       newfail=0
       waittime=0
       do i=nxs-int(InitAspLen/2),nxs+int(InitAspLen/2)+1
         broken(i)=1
         peak_12(i)=(s12(nys,i)+strini(i))*.99
         write (37,*) i,j,it,nbr
       enddo
       call flush(37)
c      print*, 'files read, initializing stuff...'
c EFFECTIVE MEDIUM PARAMETERS:
c this is meant to increase accuracy in heterogeneous media:
c!$omp parallel do  private(i,j)
        do i=1,nx
          do j=1,ny
            muu(j,i)=mu(j,i)
            ro1(j,i)=rho(j,i)
            ro2(j,i)=rho(j,i)
           enddo
        enddo
c!$omp parallel do private(i,j)
        do i=1,nx-1
         do j=2,ny
          muu(j,i)=4./(1/mu(j,i)+ 1/mu(j,i+1)+ 
     &             1/mu(j-1,i+1)+ 1/mu(j-1,i))  
          ro1(j,i)=2./((1/rho(j,i)+1/rho(j,i+1)))
          ro2(j,i)=2./((1/rho(j,i)+1/rho(j-1,i)))
         enddo
        enddo
c
c PML elastic parameters
c 
        call ParElPML
        CALL Initdistance
c
c!$omp parallel do private(i,j)
      do i=1,nx
       do j=1,ny
        v1(j,i)=0.
        v2(j,i)=0.
        s11(j,i)=0.
        s22(j,i)=0.
        s12(j,i)=0.
       enddo
      enddo
      open (unit=37, file='./RES/ruptures',form='formatted')

c
c  PML ARRAYS INITIALIZING
c
      call IniziaPML
c
      print*,'starting computation...'
      return
      end
c************************************
       subroutine stress
       include 'arrays.f'      
c Sij = lam*dij*Ekk + 2*mu*Eij
c 4th order interpolation:
c
       call strekko4
       call streijo4

       call drive_PML_stress
 
       return
       end 
c**************************************
       subroutine drive_PML_stress
       include 'arrays.f'
c Driving PML boundary condition
c
c   Left boundary
       call strekko2PMLleft
       call streijo2PMLleft
c   Right boundary
       call strekko2PMLright
       call streijo2PMLright
c   Bottom boundary
       call strekko2PMLbottom
       call streijo2PMLbottom
c   Top boundary
       call strekko2PMLtop
       call streijo2PMLtop
c   Bottom left boundary
       call strekko2PMLbottomleft
       call streijo2PMLbottomleft
c   Bottom right boundary
       call strekko2PMLbottomright
       call streijo2PMLbottomright
c   Top left boundary
       call strekko2PMLtopleft
       call streijo2PMLtopleft
c   Top right boundary
       call strekko2PMLtopright
       call streijo2PMLtopright
c   Paste stress between different blocks
       call PasteStressPML
c
       return
       end
ccccccccccccccccccccccccccccccccccc
       subroutine extrapolate
       include 'arrays.f'
c
c STAGGERING:         .---> 1
c                     | 
c       |   X | |     | 
c             |       V 2
c       O   - | O        
c       ------|----   X=s12
c       |   X | |     O=s11&s22
c             |       -=v1
c                     |=v2
c EXTRAPOL SPEED:
c  4th order interpolation:
!$omp parallel do private(i,j)
       do i=3,nx-2
       do j=3,ny-2
         v1(j,i)= v1(j,i)+(dt/(dx*ro1(j,i)))*(
     &    (rn*(s11(j,i+1)-s11(j,i))+rnn*(s11(j,i+2)-s11(j,i-1)))
     &   +(rn*(s12(j+1,i)-s12(j,i))+rnn*(s12(j+2,i)-s12(j-1,i))) )
         v2(j,i)= v2(j,i)+(dt/(dx*ro2(j,i)))*(
     &    (rn*(s22(j,i)-s22(j-1,i))+rnn*(s22(j+1,i)-s22(j-2,i)))
     &   +(rn*(s12(j,i)-s12(j,i-1))+rnn*(s12(j,i+1)-s12(j,i-2))) )
       enddo
       enddo

       call drive_PML_velocity
       return
       end 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine drive_PML_velocity
       include 'arrays.f'
c
c DRIVING PML VELOCITY ROUTINE  
c
c  LEFT-RIGHT
c 
         call velo2PMLleft
         call velo2PMLright
c
c BOTTOM-TOP
c
         call velo2PMLbottom
         call velo2PMLtop
c
c  BOTTOM-TOP LEFT-RIGHT
c
         call velo2PMLbottomleft
         call velo2PMLbottomright
         call velo2PMLtopleft
         call velo2PMLtopright
c
c  PASTE VELOCITY BETWEEN DIFFERENT BLOCKS
c
       call PasteVeloPML
c
       return
       end
ccccccccccccccccccccccccccccccccccccccccc
       subroutine save
       include 'arrays.f'
       character*15 nam1
       character*15 nam2
       character*15 nam3
       character*15 nam4
       character*15 nam5
       real*8 r8,rp8,rt8
       real buff(nx,ny)
       integer*4 it4
       integer*4 ip4
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
       backspace (25)
       write (25,*) 'current iteration ',it
       call flush(25)
c%%% save trace:
       if (is2.ge.nos2) then
         is2=0
         write(31,79) (slip(i),i=1,nx)
         call flush(31)
       endif
c%%% save fields for snapshots:
      if (is1.ge.nos1) then
       is1=0
       numf=numf+1
       write (nam1,'(a,I5.5)') './RES/fvex',it
       write (nam2,'(a,I5.5)') './RES/fvey',it
       write (nam3,'(a,I5.5)') './RES/fs11',it
       write (nam4,'(a,I5.5)') './RES/fs22',it
       write (nam5,'(a,I5.5)') './RES/fs12',it
       open  (51,file=nam1,form='formatted')
       open  (52,file=nam2,form='formatted')
       open  (53,file=nam3,form='formatted')
       open  (54,file=nam4,form='formatted')
       open  (55,file=nam5,form='formatted')
       do i=1,nx
        write(51,79)(v1(j,i),j=1,ny)
        write(52,79)(v2(j,i),j=1,ny)
        write(53,79)(s11(j,i),j=1,ny)
        write(54,79)(s22(j,i),j=1,ny)
        write(55,79)(s12(j,i),j=1,ny)
       enddo
       do jar=51,55
         close (jar)
       enddo
78    endif
ccccc
79    format(1000G17.9)
89    format(G17.9)
      return
      end

cccccccccccccccccccccccccccccccccccccccc
       subroutine source
       include 'arrays.f'
c STRESS SOURCE: (this is an explosive source, gaussian in time)
       rf=0.
       rf=(1./tau_s)*exp(-((float(it)*dt-2.*tau_s)/tau_s)**2)
       s22(nys,nxs)=s22(nys,nxs)+rf
       s11(nys,nxs)=s11(nys,nxs)+rf
c      write(26,*) rf
c      call flush(26)
       return
       end
cccccccccccccccccccccccccccccccccccccccc
       subroutine fillmodel
       include 'arrays.f'
       integer jr,nyb
       real rmu,rlam,rrho
       parameter (nyb = ny+2*npm-4)
       character*80 fmodel
c      rmu=25798500000.
c      rlam=52851500000.
c      rrho=2600.
c      vs=3220. ; vp=5577.
       rmu= 28000000000.
       rlam=rmu
       rrho=2700.
c  end of declarations       
       do i=1,npm+1
c FILL TOP-LEFT BORDER:
        do j=1,npm+1
         muTL(j,i)=rmu
         lamTL(j,i)=rlam
         rhoTL(j,i)=rrho
        enddo
c FILL LEFT BORDER:
        do j=1,ny
         jr = j+npm-2
         muL(j,i)= rmu
         lamL(j,i)=rlam
         rhoL(j,i)=rrho
        enddo
c FILL BOTTOM LEFT BORDER:
        do j=1,npm+1
         jr = j+npm+ny-5
         muBL(j,i)=rmu
         lamBL(j,i)=rlam
         rhoBL(j,i)=rrho
        enddo
       enddo
c CENTRAL 3 BACKSPACE (overlap left/center)
       do i=1,nx
c FILL TOP BORDER:
        do j=1,npm+1
         muT(j,i)=rmu
         lamT(j,i)=rlam
         rhoT(j,i)=rrho
        enddo
c FILL CENTRAL BORDER:
        do j=1,ny
         jr = j+npm-2
         mu(j,i)=rmu
         lam(j,i)=rlam
         rho(j,i)=rrho
        enddo
c FILL BOTTOM BORDER:
        do j=1,npm+1
         jr = j+npm+ny-5
         muB(j,i)=rmu
         lamB(j,i)=rlam
         rhoB(j,i)=rrho
        enddo
       enddo
c RIGHT 3 BACKSPACE (overlap center/right)
       do i=1,npm+1
c FILL TOP-RIGHT BORDER:
        do j=1,npm+1
         muTR(j,i)=rmu
         lamTR(j,i)=rlam
         rhoTR(j,i)=rrho
        enddo
c FILL RIGHT BORDER:
        do j=1,ny
         jr = j+npm-2
         muR(j,i)=rmu
         lamR(j,i)=rlam
         rhoR(j,i)=rrho
        enddo
c FILL BOTTOM RIGHT BORDER:
        do j=1,npm+1
         jr = j+npm+ny-5
         muBR(j,i)=rmu
         lamBR(j,i)=rlam
         rhoBR(j,i)=rrho
        enddo
       enddo
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccc
       subroutine strekko4
       include 'arrays.f'
!$omp parallel do private(i,j)
       do i=3,nx-2
         rbr=float(broken(i))
       do j=3,ny-2
         s11(j,i)=  s11(j,i) + (dt/dx) * (
     &    (lam(j,i)+2*mu(j,i)) * 
     &    (rn*(v1(j,i)-v1(j,i-1))+rnn*(v1(j,i+1)-v1(j,i-2)))
     &    +lam(j,i) * 
     &    (rn*(v2(j+1,i)-v2(j,i))+rnn*(v2(j+2,i)-v2(j-1,i))) )
         s22(j,i)=s22(j,i)+ (dt/dx)*(
     &    (lam(j,i)+2*mu(j,i))*
     &    (rn*(v2(j+1,i)-v2(j,i))+rnn*(v2(j+2,i)-v2(j-1,i)))
     &    + lam(j,i)*(
     &    (rn*(v1(j,i)-v1(j,i-1))+rnn*(v1(j,i+1)-v1(j,i-2)))
     &               ))
!    &    (1.-rbr*fault(j)-rbr*fault(j-1))* 
!    &    (rn*(v1(j,i)-v1(j,i-1))+rnn*(v1(j,i+1)-v1(j,i-2)))
!    &    +(rbr*fault(j)+rbr*fault(j+1))*(v1(j,i)-v1(j,i-1))
!    &               ))
      enddo
      enddo
      return 
      end
ccccccccccccccccccccccccccc
       subroutine streijo4
       include 'arrays.f'
!$omp parallel do private(i,j)
       do i=3,nx-2
        rbr=float(broken(i))
       do j=3,ny-2 !fault(j) impose order 2 across fault
        s12(j,i)=  s12(j,i) + muu(j,i)*(dt/dx)*(
     &   (1.-rbr*fault(j-1)-rbr*fault(j+1)-rbr*fault(j))*
     &   (rn*(v1(j,i)-v1(j-1,i))+rnn*(v1(j+1,i)-v1(j-2,i)))+
     &   (rbr*fault(j-1)+rbr*fault(j+1)+rbr*fault(j))
     &   *(v1(j,i)-v1(j-1,i)) +
     &   (rn*(v2(j,i+1)-v2(j,i))+rnn*(v2(j,i+2)-v2(j,i-1)))
     &   ) 
     &   +diss*(rbr*fault(j)+rbr*fault(j+1)+rbr*fault(j-1))
     &   *(v2(j,i+1)-v2(j,i))
       enddo
       enddo
       return
       end
cccccccccccccccccccccccccccccccccccccccc
      subroutine rupture
      include 'arrays.f'
      real friction,dd,dv,rns
      real rtim,rpos
      j=nys
ccc fail new points if not spontaneous rupture yet:
       if (newfail .lt. 2) then
         waittime=waittime+1
         if (waittime > 500) then
           waittime=0
c
           i=epoints(1)
           broken(i)=1
           nbr=nbr+1
           peak_12(i)=(s12(j,i)+strini(i))*.99
           write (37,*) i,j,it,nbr
           epoints(1)=i-1
c
           i=epoints(2)
           broken(i)=1
           nbr=nbr+1
           peak_12(i)=(s12(j,i)+strini(i))*.99
           write (37,*) i,j,it,nbr
           epoints(2)=i+1
         endif
       endif
ccc check for spontaneous rupture:

       do i=1,nx
        strini(i)=strini(i)+vcha
        rns=( s22(j,i)+s22(j-1,i) ) / 2.
        rdis(i)=v1(j,i)-v1(j-1,i)
        rpos=float(i-nxs)*dx
        rtim=float(it)*dt
        if( 
     &   (s12(j,i)+strini(i).ge.peak+rns*mu_s).and.
     &   (broken(i).eq.0) ) then
c       if(rvr*rtim.gt.rpos.and.rvh*(rtim-tau_s).lt.rpos.
c    &    and.broken(i).eq.0.and.rpos.gt.0.) then
          peak_12(i)=(s12(j,i)+strini(i))*.99
          broken(i)=1
          itim(i)=it
          rdis(i)=0.
          nbr=nbr+1
          newfail=newfail+1
          if (nbr .gt. 2) then 
           vcha=0.0
          endif
          write (37,*) i,j,it,nbr
        endif
ccc update friction:
        if (broken(i).eq.1) then
         if (FricType.eq.'srw'.or.
     &       FricType.eq.''.or.
     &       FricType.eq.'sw') then
            if(delta.eq.0.)then
              dd=0.
            else
c             dd=delta/(delta+slip(i))
              dd=1.-slip(i)/delta
              if (dd.lt.0.) dd=0.
              if (dd.gt.1.) dd=1.
            endif
            if (FricType.eq.'srw') then
             if (vcar.eq.0.) then
               dv=0.
             else
               dv=vcar/(vcar+rdis(i))
             endif
             if (dv.gt.dd) dd=dv
            endif
         elseif (FricType.eq.'time') then
            dd=1.-float(it-itim(i))*dt/tau_f
            if (dd.lt.0.) dd=0.
            if (dd.gt.1.) dd=1.
         endif
ccc friction law enforcement! :
         friction=dd*peak_12(i)
         s12(j,i)=friction-strini(i)
ccc  healing criterion:
cc        if (rdis(i).lt.0.) then
c         if (rvh*(rtim-tau_s).gt.rpos) then
c            broken(i)=0
c            rdis(i)=0.
c            nbr=nbr-1
cc           strini(i)=0.
c         endif
        else
            rdis(i)=0.
        endif
        slip(i)=slip(i)+rdis(i)*dt
       enddo
c
       call flush(37)
       return
       end
cccccccccccccccccccccccccccccccccccccccccccccccc
