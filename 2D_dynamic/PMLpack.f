cccccccccccccccccccccccccccccccccccccccccc
       subroutine strekko2PMLleft
       include 'arrays.f'
 
!$omp parallel do private(i,j,daf,c0)
       do i=2,npm
         do j=3,ny-2
         daf=DiXvL(j,i)
         c0=1+daf
         soL11(j,i)=(1-daf)/c0*soL11(j,i) + (dt/dx)/c0 *
     &   (lamL(j,i)+2*muL(j,i)) *
     &     (v1L(j,i)-v1L(j,i-1))
         spL11(j,i)=spL11(j,i) + (dt/dx) *
     &   lamL(j,i) *(
     &     rn*(v2L(j+1,i)-v2L(j,i))+rnn*(v2L(j+2,i)-v2L(j-1,i)) )
         sL11(j,i)=spL11(j,i)+soL11(j,i)
         spL22(j,i)=spL22(j,i)+ (dt/dx)*
     &         (lamL(j,i)+2*muL(j,i))*
     &   (rn*(v2L(j+1,i)-v2L(j,i))+rnn*(v2L(j+2,i)-v2L(j-1,i)))
         soL22(j,i)=(1-daf)/c0*soL22(j,i) + (dt/dx)/c0
     &          * lamL(j,i)*
     &   (v1L(j,i)-v1L(j,i-1))
         sL22(j,i)=spL22(j,i)+soL22(j,i)
        enddo
       enddo

       return
       end
ccccccccccccccccccccccccccccccccccccccc
       subroutine streijo2PMLleft
       include 'arrays.f'
 
!$omp parallel do private(i,j,daf,c0)
       do i=1,npm
        do j=3,ny-2
 
        daf=DiXoL(j,i)
        c0=1+daf
        soL12(j,i)=(1-daf)/c0*soL12(j,i) + (dt/dx)/c0
     &  *muuL(j,i) *
     &    (v2L(j,i+1)-v2L(j,i))
        spL12(j,i)=spL12(j,i) + (dt/dx) *
     &  muuL(j,i) * (
     &    rn*(v1L(j,i)-v1L(j-1,i))+rnn*(v1L(j+1,i)-v1L(j-2,i)) )
        sL12(j,i)=spL12(j,i)+soL12(j,i)
        enddo
       enddo
       return
       end
ccccccccccccccccccccccccccccccccccccccc
       subroutine strekko2PMLright
       include 'arrays.f'
!$omp parallel do private(i,j,daf,c0)
       do i=2,npm+1
        do j=3,ny-2
         daf=DiXvR(j,i)
         c0=1+daf
         soR11(j,i)=(1-daf)/c0*soR11(j,i) + (dt/dx)/c0 *
     &   (lamR(j,i)+2*muR(j,i)) *
     &     (v1R(j,i)-v1R(j,i-1))
         spR11(j,i)=spR11(j,i) + (dt/dx) *
     &   lamR(j,i) *(
     &     rn*(v2R(j+1,i)-v2R(j,i))+rnn*(v2R(j+2,i)-v2R(j-1,i)) )
         sR11(j,i)=spR11(j,i)+soR11(j,i)
         spR22(j,i)=spR22(j,i)+ (dt/dx)*
     &         (lamR(j,i)+2*muR(j,i))*
     &   (rn*(v2R(j+1,i)-v2R(j,i))+rnn*(v2R(j+2,i)-v2R(j-1,i)))
         soR22(j,i)=(1-daf)/c0*soR22(j,i) + (dt/dx)/c0
     &          * lamR(j,i)*
     &   (v1R(j,i)-v1R(j,i-1))
         sR22(j,i)=spR22(j,i)+soR22(j,i)
        enddo
       enddo
       return
       end
ccccccccccccccccccccccccccccccccccccccccc
       subroutine streijo2PMLright
       include 'arrays.f'
 
!$omp parallel do private(i,j,daf,c0)
       do i=2,npm
        do j=3,ny-2

        daf=DiXoR(j,i)
        c0=1+daf
        soR12(j,i)=(1-daf)/c0*soR12(j,i) + (dt/dx)/c0
     &  *muuR(j,i) *
     &    (v2R(j,i+1)-v2R(j,i))
        spR12(j,i)=spR12(j,i) + (dt/dx) *
     &  muuR(j,i) * (
     &    rn*(v1R(j,i)-v1R(j-1,i))+rnn*(v1R(j+1,i)-v1R(j-2,i)) )
        sR12(j,i)=spR12(j,i)+soR12(j,i)

        enddo
       enddo
       return
       end
ccccccccccccccccccccccccccccccccccccccc
       subroutine strekko2PMLbottom
       include 'arrays.f'

!$omp parallel do private(i,j,daf,c0)
       do i=3,nx-2
        do j=2,npm
 
        daf=DiYoB(j,i)
        c0=1.+daf
        spB11(j,i)=(1.-daf)/c0*spB11(j,i) + (dt/dx)/c0 *
     &  lamB(j,i) *
     &    (v2B(j+1,i)-v2B(j,i))
        soB11(j,i)= soB11(j,i) +(dt/dx) *
     &  (lamB(j,i)+2*muB(j,i)) * (
     &    rn*(v1B(j,i)-v1B(j,i-1))+rnn*(v1B(j,i+1)-v1B(j,i-2)) )
        sB11(j,i)=spB11(j,i)+soB11(j,i)
        spB22(j,i)=(1-daf)/c0*spB22(j,i) + (dt/dx)/c0 *
     &        (lamB(j,i)+2*muB(j,i))*
     &  (v2B(j+1,i)-v2B(j,i))
        soB22(j,i)=soB22(j,i)+ (dt/dx)
     &         * lamB(j,i)*
     &  (rn*(v1B(j,i)-v1B(j,i-1))+rnn*(v1B(j,i+1)-v1B(j,i-2)))
        sB22(j,i)=spB22(j,i)+soB22(j,i)

        enddo
       enddo

       return
       end
ccccccccccccccccccccccccccccccccccccccc
       subroutine streijo2PMLbottom
       include 'arrays.f'
 
!$omp parallel do private(i,j,daf,c0)
       do i=3,nx-2
        do j=2,npm+1

        daf=DiYvB(j,i)
        c0=1.+daf
        spB12(j,i)=(1.-daf)/c0*spB12(j,i) + (dt/dx)/c0 *
     &  muuB(j,i) *
     &    (v1B(j,i)-v1B(j-1,i))
        soB12(j,i)= soB12(j,i) +(dt/dx) *
     &  muuB(j,i) * (
     &    rn*(v2B(j,i+1)-v2B(j,i))+rnn*(v2B(j,i+2)-v2B(j,i-1)) )
        sB12(j,i)=spB12(j,i)+soB12(j,i)

        enddo
       enddo
       return
       end
ccccccccccccccccccccccccccccccccccccccc
       subroutine strekko2PMLtop
       include 'arrays.f'

!$omp parallel do private(i,j,daf,c0)
       do i=3,nx-2
        do j=1,npm
 
        daf=DiYoT(j,i)
        c0=1.+daf
        spT11(j,i)=(1.-daf)/c0*spT11(j,i) + (dt/dx)/c0 *
     &  lamT(j,i) *
     &    (v2T(j+1,i)-v2T(j,i))
        soT11(j,i)= soT11(j,i) +(dt/dx) *
     &  (lamT(j,i)+2*muT(j,i)) * (
     &    rn*(v1T(j,i)-v1T(j,i-1))+rnn*(v1T(j,i+1)-v1T(j,i-2)) )
        sT11(j,i)=spT11(j,i)+soT11(j,i)
        spT22(j,i)=(1-daf)/c0*spT22(j,i) + (dt/dx)/c0 *
     &        (lamT(j,i)+2*muT(j,i))*
     &  (v2T(j+1,i)-v2T(j,i))
        soT22(j,i)=soT22(j,i)+ (dt/dx)
     &         * lamT(j,i)*
     &  (rn*(v1T(j,i)-v1T(j,i-1))+rnn*(v1T(j,i+1)-v1T(j,i-2)))
        sT22(j,i)=spT22(j,i)+soT22(j,i)
        enddo
       enddo
       return
       end
c*******************************************************       
       subroutine streijo2PMLtop
       include 'arrays.f'
 
!$omp parallel do private(i,j,daf,c0)
       do i=3,nx-2
        do j=2,npm

        daf=DiYvT(j,i)
        c0=1.+daf
        spT12(j,i)=(1.-daf)/c0*spT12(j,i) + (dt/dx)/c0 *
     &  muuT(j,i) *
     &    (v1T(j,i)-v1T(j-1,i))
        soT12(j,i)= soT12(j,i) +(dt/dx) *
     &  muuT(j,i) * (
     &    rn*(v2T(j,i+1)-v2T(j,i))+rnn*(v2T(j,i+2)-v2T(j,i-1)) )
        sT12(j,i)=spT12(j,i)+soT12(j,i)

        enddo
       enddo
       return
       end
ccccccccccccccccccccccccccccccccccccccc
       subroutine strekko2PMLbottomleft
       include 'arrays.f'
!$omp parallel do private(i,j,dafx,dafy,c0x,c0y)
       do i=2,npm
        do j=2,npm

        dafy=DiYoBL(j,i)
        dafx=DiXvBL(j,i)
        c0x=1.+dafx
        c0y=1.+dafy
        spBL11(j,i)=(1.-dafy)/c0y*spBL11(j,i) + (dt/dx)/c0y *
     &  lamBL(j,i) *
     &    (v2BL(j+1,i)-v2BL(j,i))
        soBL11(j,i)=(1.-dafx)/c0x*soBL11(j,i) + (dt/dx)/c0x *
     &  (lamBL(j,i)+2*muBL(j,i)) *
     &    (v1BL(j,i)-v1BL(j,i-1))
        sBL11(j,i)=soBL11(j,i)+spBL11(j,i)
        spBL22(j,i)=(1-dafy)/c0y*spBL22(j,i) + (dt/dx)/c0y *
     &        (lamBL(j,i)+2*muBL(j,i))*
     &  (v2BL(j+1,i)-v2BL(j,i))
        soBL22(j,i)=(1-dafx)/c0x*soBL22(j,i) + (dt/dx)/c0x
     &         * lamBL(j,i)*
     &  (v1BL(j,i)-v1BL(j,i-1))
        sBL22(j,i)=soBL22(j,i)+spBL22(j,i)
       
        enddo
       enddo
       return
       end
ccccccccccccccccccccccccccccccccccccccc
       subroutine streijo2PMLbottomleft
       include 'arrays.f'
!$omp parallel do private(i,j,dafx,dafy,c0x,c0y)
       do i=1,npm
        do j=2,npm+1

        dafy=DiYvBL(j,i)
        dafx=DiXoBL(j,i)
        c0x=1.+dafx
        c0y=1.+dafy
        spBL12(j,i)=(1.-dafy)/c0y*spBL12(j,i) + (dt/dx)/c0y *
     &  muuBL(j,i) *
     &    (v1BL(j,i)-v1BL(j-1,i))
        soBL12(j,i)=(1-dafx)/c0x*soBL12(j,i) + (dt/dx)/c0x
     &  *muuBL(j,i) *
     &    (v2BL(j,i+1)-v2BL(j,i))
        sBL12(j,i)=soBL12(j,i)+spBL12(j,i)
        enddo
       enddo
       return
       end
cccccccccccccccccccccccccccccccccccccccc
       subroutine strekko2PMLbottomright
       include 'arrays.f'
!$omp parallel do private(i,j,dafx,dafy,c0x,c0y) 
       do i=2,npm+1
       do j=2,npm
 
        dafy=DiYoBR(j,i)
        dafx=DiXvBR(j,i)
        c0x=1.+dafx
        c0y=1.+dafy
        spBR11(j,i)=(1.-dafy)/c0y*spBR11(j,i) + (dt/dx)/c0y *
     &  lamBR(j,i) *
     &    (v2BR(j+1,i)-v2BR(j,i))
        soBR11(j,i)=(1.-dafx)/c0x*soBR11(j,i) + (dt/dx)/c0x *
     &  (lamBR(j,i)+2*muBR(j,i)) *
     &    (v1BR(j,i)-v1BR(j,i-1))
        sBR11(j,i)=soBR11(j,i)+spBR11(j,i)
        spBR22(j,i)=(1-dafy)/c0y*spBR22(j,i) + (dt/dx)/c0y *
     &        (lamBR(j,i)+2*muBR(j,i))*
     &  (v2BR(j+1,i)-v2BR(j,i))
        soBR22(j,i)=(1-dafx)/c0x*soBR22(j,i) + (dt/dx)/c0x
     &         * lamBR(j,i)*
     &  (v1BR(j,i)-v1BR(j,i-1))
        sBR22(j,i)=soBR22(j,i)+spBR22(j,i)

        enddo
       enddo
       return
       end
ccccccccccccccccccccccccccccccccccccccc
       subroutine streijo2PMLbottomright
       include 'arrays.f'
!$omp parallel do private(i,j,dafx,dafy,c0x,c0y) 
       do i=2,npm
        do j=2,npm+1

        dafy=DiYvBR(j,i)
        dafx=DiXoBR(j,i)
        c0x=1.+dafx
        c0y=1.+dafy
        spBR12(j,i)=(1.-dafy)/c0y*spBR12(j,i) + (dt/dx)/c0y *
     &  muuBR(j,i) *
     &    (v1BR(j,i)-v1BR(j-1,i))
        soBR12(j,i)=(1-dafx)/c0x*soBR12(j,i) + (dt/dx)/c0x
     &  *muuBR(j,i) *
     &    (v2BR(j,i+1)-v2BR(j,i))
        sBR12(j,i)=soBR12(j,i)+spBR12(j,i)
      
        enddo
       enddo
       return
       end
ccccccccccccccccccccccccccccccccccccccc
       subroutine strekko2PMLtopleft
       include 'arrays.f'
!$omp parallel do private(i,j,dafx,dafy,c0x,c0y)
       do i=2,npm
        do j=1,npm

        dafy=DiYoTL(j,i)
        dafx=DiXvTL(j,i)
        c0x=1.+dafx
        c0y=1.+dafy
        spTL11(j,i)=(1.-dafy)/c0y*spTL11(j,i) + (dt/dx)/c0y *
     &  lamTL(j,i) *
     &    (v2TL(j+1,i)-v2TL(j,i))
        soTL11(j,i)=(1.-dafx)/c0x*soTL11(j,i) + (dt/dx)/c0x *
     &  (lamTL(j,i)+2*muTL(j,i)) *
     &    (v1TL(j,i)-v1TL(j,i-1))
        sTL11(j,i)=soTL11(j,i)+spTL11(j,i)
        spTL22(j,i)=(1-dafy)/c0y*spTL22(j,i) + (dt/dx)/c0y *
     &        (lamTL(j,i)+2*muTL(j,i))*
     &  (v2TL(j+1,i)-v2TL(j,i))
        soTL22(j,i)=(1-dafx)/c0x*soTL22(j,i) + (dt/dx)/c0x
     &         * lamTL(j,i)*
     &  (v1TL(j,i)-v1TL(j,i-1))
        sTL22(j,i)=soTL22(j,i)+spTL22(j,i)
       
        enddo
       enddo
       return
       end
ccccccccccccccccccccccccccccccccccccccc
       subroutine streijo2PMLtopleft
       include 'arrays.f'
!$omp parallel do private(i,j,dafx,dafy,c0x,c0y)
       do i=1,npm
        do j=2,npm

        dafy=DiYvTL(j,i)
        dafx=DiXoTL(j,i)
        c0x=1.+dafx
        c0y=1.+dafy
        spTL12(j,i)=(1.-dafy)/c0y*spTL12(j,i) + (dt/dx)/c0y *
     &  muuTL(j,i) *
     &    (v1TL(j,i)-v1TL(j-1,i))
        soTL12(j,i)=(1-dafx)/c0x*soTL12(j,i) + (dt/dx)/c0x
     &  *muuTL(j,i) *
     &    (v2TL(j,i+1)-v2TL(j,i))
        sTL12(j,i)=soTL12(j,i)+spTL12(j,i)
        enddo
       enddo
       return
       end
cccccccccccccccccccccccccccccccccccccccc
       subroutine strekko2PMLtopright
       include 'arrays.f'
!$omp parallel do private(i,j,dafx,dafy,c0x,c0y) 
       do i=2,npm+1
       do j=1,npm
 
        dafy=DiYoTR(j,i)
        dafx=DiXvTR(j,i)
        c0x=1.+dafx
        c0y=1.+dafy
        spTR11(j,i)=(1.-dafy)/c0y*spTR11(j,i) + (dt/dx)/c0y *
     &  lamTR(j,i) *
     &    (v2TR(j+1,i)-v2TR(j,i))
        soTR11(j,i)=(1.-dafx)/c0x*soTR11(j,i) + (dt/dx)/c0x *
     &  (lamTR(j,i)+2*muTR(j,i)) *
     &    (v1TR(j,i)-v1TR(j,i-1))
        sTR11(j,i)=soTR11(j,i)+spTR11(j,i)
        spTR22(j,i)=(1-dafy)/c0y*spTR22(j,i) + (dt/dx)/c0y *
     &        (lamTR(j,i)+2*muTR(j,i))*
     &  (v2TR(j+1,i)-v2TR(j,i))
        soTR22(j,i)=(1-dafx)/c0x*soTR22(j,i) + (dt/dx)/c0x
     &         * lamTR(j,i)*
     &  (v1TR(j,i)-v1TR(j,i-1))
        sTR22(j,i)=soTR22(j,i)+spTR22(j,i)

        enddo
       enddo
       return
       end
ccccccccccccccccccccccccccccccccccccccc
       subroutine streijo2PMLtopright
       include 'arrays.f'
!$omp parallel do private(i,j,dafx,dafy,c0x,c0y) 
       do i=2,npm
        do j=2,npm

        dafy=DiYvTR(j,i)
        dafx=DiXoTR(j,i)
        c0x=1.+dafx
        c0y=1.+dafy
        spTR12(j,i)=(1.-dafy)/c0y*spTR12(j,i) + (dt/dx)/c0y *
     &  muuTR(j,i) *
     &    (v1TR(j,i)-v1TR(j-1,i))
        soTR12(j,i)=(1-dafx)/c0x*soTR12(j,i) + (dt/dx)/c0x
     &  *muuTR(j,i) *
     &    (v2TR(j,i+1)-v2TR(j,i))
        sTR12(j,i)=soTR12(j,i)+spTR12(j,i)
      
        enddo
       enddo
       return
       end
ccccccccccccccccccccccccccc
       subroutine velo2PMLleft
       include 'arrays.f'
!$omp parallel do private(i,j,c0,daf)
       do i=2,npm
        do j=3,ny-2
        daf=DiXoL(j,i)
        c0=1.+daf
         v1oL(j,i)=(1.-daf)/c0*v1oL(j,i)+(dt/(dx*ro1L(j,i)*c0))*
     &    (sL11(j,i+1)-sL11(j,i))
         v1pL(j,i)=v1pL(j,i)+(dt/(dx*ro1L(j,i)))*
     &   (rn*(sL12(j+1,i)-sL12(j,i))+rnn*(sL12(j+2,i)-sL12(j-1,i)))
         v1L(j,i)=v1pL(j,i)+v1oL(j,i)
        daf=DiXvL(j,i)
        c0=1.+daf
         v2oL(j,i)=(1.-daf)/c0*v2oL(j,i)+(dt/(dx*ro2L(j,i)*c0))*
     &   (sL12(j,i)-sL12(j,i-1))
         v2pL(j,i)= v2pL(j,i)+(dt/(dx*ro2L(j,i)))*(
     &    rn*(sL22(j,i)-sL22(j-1,i))+rnn*(sL22(j+1,i)-sL22(j-2,i)))
         v2L(j,i)=v2oL(j,i)+v2pL(j,i)
        enddo
        enddo
       return
       end
ccccccccccccccccccccccccccc
       subroutine velo2PMLright
       include 'arrays.f'
!$omp parallel do private(i,j,daf,c0)
        do i=2,npm
        do j=3,ny-2
        daf=DiXoR(j,i)
        c0=1.+daf
         v1oR(j,i)=(1.-daf)/c0*v1oR(j,i)+(dt/(dx*ro1R(j,i)*c0))*
     &    (sR11(j,i+1)-sR11(j,i))
         v1pR(j,i)=v1pR(j,i)+(dt/(dx*ro1R(j,i)))*
     &   (rn*(sR12(j+1,i)-sR12(j,i))+rnn*(sR12(j+2,i)-sR12(j-1,i)))
         v1R(j,i)=v1pR(j,i)+v1oR(j,i)
        daf=DiXvR(j,i)
        c0=1.+daf
         v2oR(j,i)=(1.-daf)/c0*v2oR(j,i)+(dt/(dx*ro2R(j,i)*c0))*
     &   (sR12(j,i)-sR12(j,i-1))
         v2pR(j,i)= v2pR(j,i)+(dt/(dx*ro2R(j,i)))*(
     &    rn*(sR22(j,i)-sR22(j-1,i))+rnn*(sR22(j+1,i)-sR22(j-2,i)))
         v2R(j,i)=v2oR(j,i)+v2pR(j,i)
        enddo
        enddo

       return
       end
ccccccccccccccccccccccccc
       subroutine velo2PMLbottom
       include 'arrays.f'
!$omp parallel do private(i,j,c0,daf)
       do i=3,nx-2
        do j=2,npm
       daf=DiYoB(j,i)
       c0=1.+daf
       v1oB(j,i)=v1oB(j,i)+(dt/(dx*ro1B(j,i)))*
     &    (rn*(sB11(j,i+1)-sB11(j,i))+rnn*(sB11(j,i+2)-sB11(j,i-1)))
         v1pB(j,i)=(1.-daf)/c0*v1pB(j,i)+(dt/(dx*ro1B(j,i)*c0))*
     &   (sB12(j+1,i)-sB12(j,i))
         v1B(j,i)=v1pB(j,i)+v1oB(j,i)
       daf=DiYvB(j,i)
       c0=1.+daf
         v2oB(j,i)=v2oB(j,i)+(dt/(dx*ro2B(j,i)))*
     &   (rn*(sB12(j,i)-sB12(j,i-1))+rnn*(sB12(j,i+1)-sB12(j,i-2)))
         v2pB(j,i)=(1.-daf)/c0*v2pB(j,i)+(dt/(dx*ro2B(j,i)*c0))*(
     &    sB22(j,i)-sB22(j-1,i))
         v2B(j,i)=v2oB(j,i)+v2pB(j,i)
        enddo
       enddo
       return
       end
ccccccccccccccccccccccc
       subroutine velo2PMLtop
       include 'arrays.f'
!$omp parallel do private(i,j,c0,daf)
       do i=3,nx-2
        do j=2,npm
       daf=DiYoT(j,i)
       c0=1.+daf
       v1oT(j,i)=v1oT(j,i)+(dt/(dx*ro1T(j,i)))*
     &    (rn*(sT11(j,i+1)-sT11(j,i))+rnn*(sT11(j,i+2)-sT11(j,i-1)))
         v1pT(j,i)=(1.-daf)/c0*v1pT(j,i)+(dt/(dx*ro1T(j,i)*c0))*
     &   (sT12(j+1,i)-sT12(j,i))
         v1T(j,i)=v1pT(j,i)+v1oT(j,i)
       daf=DiYvT(j,i)
       c0=1.+daf
         v2oT(j,i)=v2oT(j,i)+(dt/(dx*ro2T(j,i)))*
     &   (rn*(sT12(j,i)-sT12(j,i-1))+rnn*(sT12(j,i+1)-sT12(j,i-2)))
         v2pT(j,i)=(1.-daf)/c0*v2pT(j,i)+(dt/(dx*ro2T(j,i)*c0))*(
     &    sT22(j,i)-sT22(j-1,i))
         v2T(j,i)=v2oT(j,i)+v2pT(j,i)
        enddo
       enddo
       return
       end
cccccccccccccccccccccccccccccccccccccccccccc
       subroutine velo2PMLbottomleft
       include 'arrays.f'
!$omp parallel do private(i,j,c0x,c0y,dafx,dafy)
       do i=2,npm
        do j=2,npm

       dafy=DiYoBL(j,i)
       dafx=DiXoBL(j,i)
       c0x=1.+dafx
       c0y=1.+dafy
       v1oBL(j,i)=(1.-dafx)/c0x*v1oBL(j,i)+(dt/(dx*ro1BL(j,i)*c0x))*
     &    (sBL11(j,i+1)-sBL11(j,i))
       v1pBL(j,i)=(1.-dafy)/c0y*v1pBL(j,i)+(dt/(dx*ro1BL(j,i)*c0y))*
     &   (sBL12(j+1,i)-sBL12(j,i))
       v1BL(j,i)=v1oBL(j,i)+v1pBL(j,i)
       dafy=DiYvBL(j,i)
       dafx=DiXvBL(j,i)
       c0x=1.+dafx
       c0y=1.+dafy
       v2oBL(j,i)=(1.-dafx)/c0x*v2oBL(j,i)+(dt/(dx*ro2BL(j,i)*c0x))*
     &   (sBL12(j,i)-sBL12(j,i-1))
       v2pBL(j,i)=(1.-dafy)/c0y*v2pBL(j,i)+(dt/(dx*ro2BL(j,i)*c0y))*
     &    (sBL22(j,i)-sBL22(j-1,i))
       v2BL(j,i)=v2oBL(j,i)+v2pBL(j,i)
        enddo
       enddo
      return
      end
ccccccccccccccccccccccc
       subroutine velo2PMLbottomright
       include 'arrays.f'
!$omp parallel do private(i,j,c0x,c0y,dafx,dafy) 
       do i=2,npm
        do j=2,npm
       dafy=DiYoBR(j,i)
       dafx=DiXoBR(j,i)
       c0x=1.+dafx
       c0y=1.+dafy
       v1oBR(j,i)=(1.-dafx)/c0x*v1oBR(j,i)+(dt/(dx*ro1BR(j,i)*c0x))*
     &    (sBR11(j,i+1)-sBR11(j,i))
       v1pBR(j,i)=(1.-dafy)/c0y*v1pBR(j,i)+(dt/(dx*ro1BR(j,i)*c0y))*
     &   (sBR12(j+1,i)-sBR12(j,i))
       v1BR(j,i)=v1oBR(j,i)+v1pBR(j,i)
       dafy=DiYvBR(j,i)
       dafx=DiXvBR(j,i)
       c0x=1.+dafx
       c0y=1.+dafy
       v2oBR(j,i)=(1.-dafx)/c0x*v2oBR(j,i)+(dt/(dx*ro2BR(j,i)*c0x))*
     &   (sBR12(j,i)-sBR12(j,i-1))
       v2pBR(j,i)=(1.-dafy)/c0y*v2pBR(j,i)+(dt/(dx*ro2BR(j,i)*c0y))*(
     &    sBR22(j,i)-sBR22(j-1,i))
       v2BR(j,i)=v2oBR(j,i)+v2pBR(j,i)
        enddo
       enddo
      return
      end
cccccccccccccccccccccccccccccccccccccccccccc
       subroutine velo2PMLtopleft
       include 'arrays.f'
!$omp parallel do private(i,j,c0x,c0y,dafx,dafy)
       do i=2,npm
        do j=2,npm

       dafy=DiYoTL(j,i)
       dafx=DiXoTL(j,i)
       c0x=1.+dafx
       c0y=1.+dafy
       v1oTL(j,i)=(1.-dafx)/c0x*v1oTL(j,i)+(dt/(dx*ro1TL(j,i)*c0x))*
     &    (sTL11(j,i+1)-sTL11(j,i))
       v1pTL(j,i)=(1.-dafy)/c0y*v1pTL(j,i)+(dt/(dx*ro1TL(j,i)*c0y))*
     &   (sTL12(j+1,i)-sTL12(j,i))
       v1TL(j,i)=v1oTL(j,i)+v1pTL(j,i)
       dafy=DiYvTL(j,i)
       dafx=DiXvTL(j,i)
       c0x=1.+dafx
       c0y=1.+dafy
       v2oTL(j,i)=(1.-dafx)/c0x*v2oTL(j,i)+(dt/(dx*ro2TL(j,i)*c0x))*
     &   (sTL12(j,i)-sTL12(j,i-1))
       v2pTL(j,i)=(1.-dafy)/c0y*v2pTL(j,i)+(dt/(dx*ro2TL(j,i)*c0y))*
     &    (sTL22(j,i)-sTL22(j-1,i))
       v2TL(j,i)=v2oTL(j,i)+v2pTL(j,i)
        enddo
       enddo
      return
      end
ccccccccccccccccccccccc
       subroutine velo2PMLtopright
       include 'arrays.f'
!$omp parallel do private(i,j,c0x,c0y,dafx,dafy) 
       do i=2,npm
        do j=2,npm
       dafy=DiYoTR(j,i)
       dafx=DiXoTR(j,i)
       c0x=1.+dafx
       c0y=1.+dafy
       v1oTR(j,i)=(1.-dafx)/c0x*v1oTR(j,i)+(dt/(dx*ro1TR(j,i)*c0x))*
     &    (sTR11(j,i+1)-sTR11(j,i))
       v1pTR(j,i)=(1.-dafy)/c0y*v1pTR(j,i)+(dt/(dx*ro1TR(j,i)*c0y))*
     &   (sTR12(j+1,i)-sTR12(j,i))
       v1TR(j,i)=v1oTR(j,i)+v1pTR(j,i)
       dafy=DiYvTR(j,i)
       dafx=DiXvTR(j,i)
       c0x=1.+dafx
       c0y=1.+dafy
       v2oTR(j,i)=(1.-dafx)/c0x*v2oTR(j,i)+(dt/(dx*ro2TR(j,i)*c0x))*
     &   (sTR12(j,i)-sTR12(j,i-1))
       v2pTR(j,i)=(1.-dafy)/c0y*v2pTR(j,i)+(dt/(dx*ro2TR(j,i)*c0y))*(
     &    sTR22(j,i)-sTR22(j-1,i))
       v2TR(j,i)=v2oTR(j,i)+v2pTR(j,i)
        enddo
       enddo
      return
      end
     
ccccccccccccccccccccccccccccccccccccc
      subroutine PasteStressPML
      include 'arrays.f'
c
c   File paste to the border
c LEFT PASTE
c
      do j=2,ny-2
       s11(j,1)=sL11(j,npm-1)
       s11(j,2)=sL11(j,npm)
       sL11(j,npm+1)=s11(j,3)
       s12(j,1)=sL12(j,npm-1)
       s12(j,2)=sL12(j,npm)
       sL12(j,npm+1)=s12(j,3)
       s22(j,1)=sL22(j,npm-1)
       s22(j,2)=sL22(j,npm)
       sL22(j,npm+1)=s22(j,3)
c
c RIGHT PASTE
c   
       s11(j,nx-1)=sR11(j,2)
       s11(j,nx)=sR11(j,3)
       sR11(j,1)=s11(j,nx-2)
       s12(j,nx-1)=sR12(j,2)
       s12(j,nx)=sR12(j,3)
       sR12(j,1)=s12(j,nx-2)
       s22(j,nx-1)=sR22(j,2)
       s22(j,nx)=sR22(j,3)
       sR22(j,1)=s22(j,nx-2)
      enddo
c
c BOTTOM PASTE
c
      do  i=3,nx-2
       s11(ny-1,i)=sB11(2,i)
       s11(ny,i)=sB11(3,i)
       sB11(1,i)=s11(ny-2,i)
       s12(ny-1,i)=sB12(2,i)
       s12(ny,i)=sB12(3,i)
       sB12(1,i)=s12(ny-2,i)
       s22(ny-1,i)=sB22(2,i)
       s22(ny,i)=sB22(3,i)
       sB22(1,i)=s22(ny-2,i)
      enddo
c
c TOP PASTE
c
      do  i=3,nx-2
       s11(1,i)=sT11(npm-1,i)
       s11(2,i)=sT11(npm,i)
       sT11(npm+1,i)=s11(3,i)
       s12(1,i)=sT12(npm-1,i)
       s12(2,i)=sT12(npm,i)
       sT12(npm+1,i)=s12(3,i)
       s22(1,i)=sT22(npm-1,i)
       s22(2,i)=sT22(npm,i)
       sT22(npm+1,i)=s22(3,i)
      enddo
c
c LEFT/RIGHT - BOTTOM LEFT/RIGHT PASTE
c
      do i=1,npm
       sBL11(1,i)=sL11(ny-2,i)
       sL11(ny-1,i)=sBL11(2,i)
       sL11(ny,i)=sBL11(3,i)
       sBL12(1,i)=sL12(ny-2,i)
       sL12(ny-1,i)=sBL12(2,i)
       sL12(ny,i)=sBL12(3,i)
       sBL22(1,i)=sL22(ny-2,i)
       sL22(ny-1,i)=sBL22(2,i)
       sL22(ny,i)=sBL22(3,i)
      enddo 
      do i=2,npm+1
       sBR11(1,i)=sR11(ny-2,i)
       sR11(ny-1,i)=sBR11(2,i)
       sR11(ny,i)=sBR11(3,i)
       sBR12(1,i)=sR12(ny-2,i)
       sR12(ny-1,i)=sBR12(2,i)
       sR12(ny,i)=sBR12(3,i)
       sBR22(1,i)=sR22(ny-2,i)
       sR22(ny-1,i)=sBR22(2,i)
       sR22(ny,i)=sBR22(3,i)
      enddo
c
c LEFT/RIGHT - TOP LEFT/RIGHT PASTE
c
      do i=1,npm
       sTL11(npm+1,i)=sL11(3,i)
       sL11(1,i)=sTL11(npm-1,i)
       sL11(2,i)=sTL11(npm,i)
       sTL12(npm+1,i)=sL12(3,i)
       sL12(1,i)=sTL12(npm-1,i)
       sL12(2,i)=sTL12(npm,i)
       sTL22(npm+1,i)=sL22(3,i)
       sL22(1,i)=sTL22(npm-1,i)
       sL22(2,i)=sTL22(npm,i)
      enddo 
      do i=2,npm+1
       sTR11(npm+1,i)=sR11(3,i)
       sR11(1,i)=sTR11(npm-1,i)
       sR11(2,i)=sTR11(npm,i)
       sTR12(npm+1,i)=sR12(3,i)
       sR12(1,i)=sTR12(npm-1,i)
       sR12(2,i)=sTR12(npm,i)
       sTR22(npm+1,i)=sR22(3,i)
       sR22(1,i)=sTR22(npm-1,i)
       sR22(2,i)=sTR22(npm,i)
      enddo
c
c BOTTOM -BOTTOM LEFT/RIGHT PASTE
c       
      do j=2,npm+1
       sBL11(j,npm+1)=sB11(j,3)
       sB11(j,1)=sBL11(j,npm-1)
       sB11(j,2)=sBL11(j,npm)
       sBL12(j,npm+1)=sB12(j,3)
       sB12(j,1)=sBL12(j,npm-1)
       sB12(j,2)=sBL12(j,npm)
       sBL22(j,npm+1)=sB22(j,3)
       sB22(j,1)=sBL22(j,npm-1)
       sB22(j,2)=sBL22(j,npm)
c       
       sBR11(j,1)=sB11(j,nx-2)
       sB11(j,nx-1)=sBR11(j,2)
       sB11(j,nx)=sBR11(j,3)
       sBR12(j,1)=sB12(j,nx-2)
       sB12(j,nx-1)=sBR12(j,2)
       sB12(j,nx)=sBR12(j,3)
       sBR22(j,1)=sB22(j,nx-2)
       sB22(j,nx-1)=sBR22(j,2)
       sB22(j,nx)=sBR22(j,3)
      enddo
c
c TOP -TOP LEFT/RIGHT PASTE
c       
      do j=1,npm
       sTL11(j,npm+1)=sT11(j,3)
       sT11(j,1)=sTL11(j,npm-1)
       sT11(j,2)=sTL11(j,npm)
       sTL12(j,npm+1)=sT12(j,3)
       sT12(j,1)=sTL12(j,npm-1)
       sT12(j,2)=sTL12(j,npm)
       sTL22(j,npm+1)=sT22(j,3)
       sT22(j,1)=sTL22(j,npm-1)
       sT22(j,2)=sTL22(j,npm)
c       
       sTR11(j,1)=sT11(j,nx-2)
       sT11(j,nx-1)=sTR11(j,2)
       sT11(j,nx)=sTR11(j,3)
       sTR12(j,1)=sT12(j,nx-2)
       sT12(j,nx-1)=sTR12(j,2)
       sT12(j,nx)=sTR12(j,3)
       sTR22(j,1)=sT22(j,nx-2)
       sT22(j,nx-1)=sTR22(j,2)
       sT22(j,nx)=sTR22(j,3)
      enddo
      return
      end
cccccccccccccccccccccccccccccccccccccc
      subroutine PasteVeloPML
      include 'arrays.f'

c LEFT-RIGHT PASTE      
       do j=3,ny-2
        v1(j,1)=v1L(j,npm-1)
        v1(j,2)=v1L(j,npm)
        v1L(j,npm+1)=v1(j,3)
        v2(j,1)=v2L(j,npm-1)
        v2(j,2)=v2L(j,npm)
        v2L(j,npm+1)=v2(j,3)
c         
        v1(j,nx-1)=v1R(j,2)
        v1(j,nx)=v1R(j,3)
        v1R(j,1)=v1(j,nx-2)
        v2(j,nx-1)=v2R(j,2)
        v2(j,nx)=v2R(j,3)
        v2R(j,1)=v2(j,nx-2)
       enddo
c BOTTOM-TOP PASTE       
       do i=3,nx-2
        v1(ny-1,i)=v1B(2,i)
        v1(ny,i)=v1B(3,i)
        v1B(1,i)=v1(ny-2,i)
        v2(ny-1,i)=v2B(2,i)
        v2(ny,i)=v2B(3,i)
        v2B(1,i)=v2(ny-2,i)
c        
        v1(1,i)=v1T(npm-1,i)
        v1(2,i)=v1T(npm,i)
        v1T(npm+1,i)=v1(3,i)
        v2(1,i)=v2T(npm-1,i)
        v2(2,i)=v2T(npm,i)
        v2T(npm+1,i)=v2(3,i)
       enddo 
c LEFT/RIGHT - BOTTOM LEFT/RIGHT PASTE       
       do i=1,npm
        v1L(ny-1,i)=v1BL(2,i)
        v1L(ny,i)=v1BL(3,i)
        v1BL(1,i)=v1L(ny-2,i)
        v2L(ny-1,i)=v2BL(2,i)
        v2L(ny,i)=v2BL(3,i)
        v2BL(1,i)=v2L(ny-2,i)
       enddo 
       do i=2,npm+1
        v1R(ny-1,i)=v1BR(2,i)
        v1R(ny,i)=v1BR(3,i)
        v1BR(1,i)=v1R(ny-2,i)
        v2R(ny-1,i)=v2BR(2,i)
        v2R(ny,i)=v2BR(3,i)
        v2BR(1,i)=v2R(ny-2,i)
       enddo
c LEFT/RIGHT - TOP LEFT/RIGHT PASTE
      do i=1,npm
       v1TL(npm+1,i)=v1L(3,i)
       v1L(1,i)=v1TL(npm-1,i)
       v1L(2,i)=v1TL(npm,i)
       v2TL(npm+1,i)=v2L(3,i)
       v2L(1,i)=v2TL(npm-1,i)
       v2L(2,i)=v2TL(npm,i)
      enddo 
      do i=2,npm+1
       v1TR(npm+1,i)=v1R(3,i)
       v1R(1,i)=v1TR(npm-1,i)
       v1R(2,i)=v1TR(npm,i)
       v2TR(npm+1,i)=v2R(3,i)
       v2R(1,i)=v2TR(npm-1,i)
       v2R(2,i)=v2TR(npm,i)
      enddo 
c BOTTOM - BOTTOM LEFT/RIGHT PASTE      
       do j=2,npm+1
        v1BL(j,npm+1)=v1B(j,3)
        v1B(j,1)=v1BL(j,npm-1)
        v1B(j,2)=v1BL(j,npm)
        v2BL(j,npm+1)=v2B(j,3)
        v2B(j,1)=v2BL(j,npm-1)
        v2B(j,2)=v2BL(j,npm)
c        
        v1BR(j,1)=v1B(j,nx-2)
        v1B(j,nx-1)=v1BR(j,2)
        v1B(j,nx)=v1BR(j,3)
        v2BR(j,1)=v2B(j,nx-2)
        v2B(j,nx-1)=v2BR(j,2)
        v2B(j,nx)=v2BR(j,3)
       enddo
c TOP - TOP LEFT/RIGHT PASTE      
       do j=1,npm
        v1TL(j,npm+1)=v1T(j,3)
        v1T(j,1)=v1TL(j,npm-1)
        v1T(j,2)=v1TL(j,npm)
        v2TL(j,npm+1)=v2T(j,3)
        v2T(j,1)=v2TL(j,npm-1)
        v2T(j,2)=v2TL(j,npm)
c        
        v1TR(j,1)=v1T(j,nx-2)
        v1T(j,nx-1)=v1TR(j,2)
        v1T(j,nx)=v1TR(j,3)
        v2TR(j,1)=v2T(j,nx-2)
        v2T(j,nx-1)=v2TR(j,2)
        v2T(j,nx)=v2TR(j,3)
       enddo 
       return
       end
ccccccccccccccccccccccccccccccccccccccccccccc
       subroutine IniziaPML
       include 'arrays.f'

       do i=1,npm+1
        do j=1,ny
         soR11(j,i)=0
         spR11(j,i)=0
         sR11(j,i)=0
         soR12(j,i)=0
         spR12(j,i)=0
         sR12(j,i)=0
         soR22(j,i)=0
         spR22(j,i)=0
         sR22(j,i)=0
         v1oR(j,i)=0
         v1pR(j,i)=0
         v1R(j,i)=0
         v2oR(j,i)=0
         v2pR(j,i)=0
         v2R(j,i)=0
         soL11(j,i)=0
         spL11(j,i)=0
         sL11(j,i)=0
         soL12(j,i)=0
         spL12(j,i)=0
         sL12(j,i)=0
         soL22(j,i)=0
         spL22(j,i)=0
         sL22(j,i)=0
         v1oL(j,i)=0
         v1pL(j,i)=0
         v1L(j,i)=0
         v2oL(j,i)=0
         v2pL(j,i)=0
         v2L(j,i)=0
        enddo
       enddo
       do i=1,nx
        do j=1,npm+1
         soB11(j,i)=0
         spB11(j,i)=0
         sB11(j,i)=0
         soB12(j,i)=0
         spB12(j,i)=0
         sB12(j,i)=0
         soB22(j,i)=0
         spB22(j,i)=0
         sB22(j,i)=0
         v1oB(j,i)=0
         v1pB(j,i)=0
         v1B(j,i)=0
         v2oB(j,i)=0
         v2pB(j,i)=0
         v2B(j,i)=0
         soT11(j,i)=0
         spT11(j,i)=0
         sT11(j,i)=0
         soT12(j,i)=0
         spT12(j,i)=0
         sT12(j,i)=0
         soT22(j,i)=0
         spT22(j,i)=0
         sT22(j,i)=0
         v1oT(j,i)=0
         v1pT(j,i)=0
         v1T(j,i)=0
         v2oT(j,i)=0
         v2pT(j,i)=0
         v2T(j,i)=0
        enddo
       enddo
       do i=1,npm+1
        do j=1,npm+1
         soBR11(j,i)=0
         spBR11(j,i)=0
         sBR11(j,i)=0
         soBR12(j,i)=0
         spBR12(j,i)=0
         sBR12(j,i)=0
         soBR22(j,i)=0
         spBR22(j,i)=0
         sBR22(j,i)=0
         v1oBR(j,i)=0
         v1pBR(j,i)=0
         v1BR(j,i)=0
         v2oBR(j,i)=0
         v2pBR(j,i)=0
         v2BR(j,i)=0
         soBL11(j,i)=0
         spBL11(j,i)=0
         sBL11(j,i)=0
         soBL12(j,i)=0
         spBL12(j,i)=0
         sBL12(j,i)=0
         soBL22(j,i)=0
         spBL22(j,i)=0
         sBL22(j,i)=0
         v1oBL(j,i)=0
         v1pBL(j,i)=0
         v1BL(j,i)=0
         v2oBL(j,i)=0
         v2pBL(j,i)=0
         v2BL(j,i)=0
         soTR11(j,i)=0
         spTR11(j,i)=0
         sTR11(j,i)=0
         soTR12(j,i)=0
         spTR12(j,i)=0
         sTR12(j,i)=0
         soTR22(j,i)=0
         spTR22(j,i)=0
         sTR22(j,i)=0
         v1oTR(j,i)=0
         v1pTR(j,i)=0
         v1TR(j,i)=0
         v2oTR(j,i)=0
         v2pTR(j,i)=0
         v2TR(j,i)=0
         soTL11(j,i)=0
         spTL11(j,i)=0
         sTL11(j,i)=0
         soTL12(j,i)=0
         spTL12(j,i)=0
         sTL12(j,i)=0
         soTL22(j,i)=0
         spTL22(j,i)=0
         sTL22(j,i)=0
         v1oTL(j,i)=0
         v1pTL(j,i)=0
         v1TL(j,i)=0
         v2oTL(j,i)=0
         v2pTL(j,i)=0
         v2TL(j,i)=0
        enddo
       enddo
       return
       end
ccccccccccccccccccccccccccccccccccca
       subroutine ParElPML
       include 'arrays.f'
c Calculate effective parameters
c
      do i=1,npm+1
       do j=1,ny
        muuL(j,i)=muL(j,i)
        ro1L(j,i)=rhoL(j,i)
        ro2L(j,i)=rhoL(j,i)
        muuR(j,i)=muR(j,i)
        ro1R(j,i)=rhoR(j,i)
        ro2R(j,i)=rhoR(j,i)
       enddo
      enddo
     
      do i=1,npm
       do j=2,ny
        muuL(j,i)=4./(1/muL(j,i)+ 1/muL(j,i+1)+
     &             1/muL(j-1,i+1)+ 1/muL(j-1,i))
        ro1L(j,i)=2./((1/rhoL(j,i)+1/rhoL(j,i+1)))
        ro2L(j,i)=2./((1/rhoL(j,i)+1/rhoL(j-1,i)))
c        
        muuR(j,i)=4./(1/muR(j,i)+ 1/muR(j,i+1)+
     &             1/muR(j-1,i+1)+ 1/muR(j-1,i))
        ro1R(j,i)=2./((1/rhoR(j,i)+1/rhoR(j,i+1)))
        ro2R(j,i)=2./((1/rhoR(j,i)+1/rhoR(j-1,i)))
       enddo
      enddo
c
      do i=1,nx
       do j=1,npm+1
        muuB(j,i)=muB(j,i)
        ro1B(j,i)=rhoB(j,i)
        ro2B(j,i)=rhoB(j,i)
        muuT(j,i)=muT(j,i)
        ro1T(j,i)=rhoT(j,i)
        ro2T(j,i)=rhoT(j,i)
       enddo
      enddo
c
      do i=1,nx-1
       do j=2,npm+1
        muuB(j,i)=4./(1/muB(j,i)+ 1/muB(j,i+1)+
     &             1/muB(j-1,i+1)+ 1/muB(j-1,i))
        ro1B(j,i)=2./((1/rhoB(j,i)+1/rhoB(j,i+1)))
        ro2B(j,i)=2./((1/rhoB(j,i)+1/rhoB(j-1,i)))
        muuT(j,i)=4./(1/muT(j,i)+ 1/muT(j,i+1)+
     &             1/muT(j-1,i+1)+ 1/muT(j-1,i))
        ro1T(j,i)=2./((1/rhoT(j,i)+1/rhoT(j,i+1)))
        ro2T(j,i)=2./((1/rhoT(j,i)+1/rhoT(j-1,i)))
       enddo
      enddo
c
      do i=1,npm+1
       do j=1,npm+1
        muuBL(j,i)=muBL(j,i)
        ro1BL(j,i)=rhoBL(j,i)
        ro2BL(j,i)=rhoBL(j,i)
        muuBR(j,i)=muBR(j,i)
        ro1BR(j,i)=rhoBR(j,i)
        ro2BR(j,i)=rhoBR(j,i)
        muuTL(j,i)=muTL(j,i)
        ro1TL(j,i)=rhoTL(j,i)
        ro2TL(j,i)=rhoTL(j,i)
        muuTR(j,i)=muTR(j,i)
        ro1TR(j,i)=rhoTR(j,i)
        ro2TR(j,i)=rhoTR(j,i)
       enddo
      enddo
c     
      do i=1,npm
       do j=2,npm+1
        muuBL(j,i)=4./(1/muBL(j,i)+ 1/muBL(j,i+1)+
     &             1/muBL(j-1,i+1)+ 1/muBL(j-1,i))
        ro1BL(j,i)=2./((1/rhoBL(j,i)+1/rhoBL(j,i+1)))
        ro2BL(j,i)=2./((1/rhoBL(j,i)+1/rhoBL(j-1,i)))
c
        muuBR(j,i)=4./(1/muBR(j,i)+ 1/muBR(j,i+1)+
     &             1/muBR(j-1,i+1)+ 1/muBR(j-1,i))
        ro1BR(j,i)=2./((1/rhoBR(j,i)+1/rhoBR(j,i+1)))
        ro2BR(j,i)=2./((1/rhoBR(j,i)+1/rhoBR(j-1,i)))
c        
        muuTL(j,i)=4./(1/muTL(j,i)+ 1/muTL(j,i+1)+
     &             1/muTL(j-1,i+1)+ 1/muTL(j-1,i))
        ro1TL(j,i)=2./((1/rhoTL(j,i)+1/rhoTL(j,i+1)))
        ro2TL(j,i)=2./((1/rhoTL(j,i)+1/rhoTL(j-1,i)))
c
        muuTR(j,i)=4./(1/muTR(j,i)+ 1/muTR(j,i+1)+
     &             1/muTR(j-1,i+1)+ 1/muTR(j-1,i))
        ro1TR(j,i)=2./((1/rhoTR(j,i)+1/rhoTR(j,i+1)))
        ro2TR(j,i)=2./((1/rhoTR(j,i)+1/rhoTR(j-1,i)))
       enddo
      enddo
      return
      end
ccccccccccccccccccccccccccccccccccc
       subroutine Initdistance
       include 'arrays.f'
c
c LEFT MATCHING
       do i=1,npm
        do j=1,ny-2
         xx=float(npm+1-i)
         vp=sqrt((lamL(j,i)+2*muL(j,i))/rhoL(j,i))
         ah=float(npm)*dx
         ah=1/ah
         xx=xx/npm
         ah=4.5*vp*ah
         DiXvL(j,i)=ah*xx**6
         xx=float(npm+1-i)-.5
         xx=xx/npm 
         DiXoL(j,i)=ah*xx**6
        enddo
       enddo
c RIGHT MATCHING
       do i=2,npm+1
        do j=3,ny-2
         xx=float(i-2)
         vp=sqrt((lamR(j,i)+2*muR(j,i))/rhoR(j,i))
         ah=float(npm)*dx
         ah=1/ah
         xx=xx/npm
         ah=4.5*vp*ah
         DiXvR(j,i)=ah*xx**6
         xx=float(i-2)+.5
         xx=xx/npm
         DiXoR(j,i)=ah*xx**6
        enddo
       enddo
c BOTTOM MATCHING
       do i=3,nx-2
        do j=2,npm+1
         xx=float(j-1)
         vp=sqrt((lamB(j,i)+2*muB(j,i))/rhoB(j,i))
         ah=float(npm)*dx
         ah=1/ah
         xx=xx/npm
         ah=4.5*vp*ah
         DiYoB(j,i)=ah*xx**6
         xx=float(j-1)-0.5
         xx=xx/npm
         DiYvB(j,i)=ah*xx**6
        enddo
       enddo
c TOP MATCHING
       do i=3,nx-2
        do j=1,npm
         xx=float(npm-j)
         vp=sqrt((lamT(j,i)+2*muT(j,i))/rhoT(j,i))
         ah=float(npm)*dx
         ah=1/ah
         xx=xx/npm
         ah=4.5*vp*ah
         DiYoT(j,i)=ah*xx**6
         xx=float(npm-j)+0.5
         xx=xx/npm
         DiYvT(j,i)=ah*xx**6
        enddo
       enddo
c BOTTOM LEFT MATCHING
       do i=1,npm 
        do j=2,npm+1
         xx=float(j-1)
         vp=sqrt((lamBL(j,i)+2*muBL(j,i))/rhoBL(j,i))
         ah=float(npm)*dx
         ah=1/ah    
         xx=xx/npm
         ah=4.5*vp*ah
         DiYoBL(j,i)=ah*xx**6
         xx=float(j-1)-0.5
         xx=xx/npm
         DiYvBL(j,i)=ah*xx**6
         xx=float(npm+1-i)
         xx=xx/npm
         DiXvBL(j,i)=ah*xx**6
         xx=float(npm+1-i)-.5
         xx=xx/npm 
         DiXoBL(j,i)=ah*xx**6
        enddo
       enddo
c BOTTOM RIGHT MATCHING
       do i=2,npm+1
        do j=2,npm+1
         xx=float(j-1)
         vp=sqrt((lamBR(j,i)+2*muBR(j,i))/rhoBR(j,i))
         ah=float(npm)*dx
         ah=1/ah   
         xx=xx/npm
         ah=4.5*vp*ah
         DiYoBR(j,i)=ah*xx**6
         xx=float(j-1)-0.5
         xx=xx/npm
         DiYvBR(j,i)=ah*xx**6
         xx=float(i-2)
         xx=xx/npm
         DiXvBR(j,i)=ah*xx**6
         xx=float(i-2)+.5
         xx=xx/npm
         DiXoBR(j,i)=ah*xx**6
        enddo
       enddo
c TOP LEFT MATCHING
       do i=1,npm 
        do j=1,npm
         xx=float(npm-j)
         vp=sqrt((lamTL(j,i)+2*muTL(j,i))/rhoTL(j,i))
         ah=float(npm)*dx
         ah=1/ah    
         xx=xx/npm
         ah=4.5*vp*ah
         DiYoTL(j,i)=ah*xx**6
         xx=float(npm-j)+0.5
         xx=xx/npm
         DiYvTL(j,i)=ah*xx**6
         xx=float(npm+1-i)
         xx=xx/npm
         DiXvTL(j,i)=ah*xx**6
         xx=float(npm+1-i)-.5
         xx=xx/npm 
         DiXoTL(j,i)=ah*xx**6
        enddo
       enddo
c TOP RIGHT MATCHING
       do i=2,npm+1
        do j=1,npm
         xx=float(npm-j)
         vp=sqrt((lamTR(j,i)+2*muTR(j,i))/rhoTR(j,i))
         ah=float(npm)*dx
         ah=1/ah   
         xx=xx/npm
         ah=4.5*vp*ah
         DiYoTR(j,i)=ah*xx**6
         xx=float(npm-j)+0.5
         xx=xx/npm
         DiYvTR(j,i)=ah*xx**6
         xx=float(i-2)
         xx=xx/npm
         DiXvTR(j,i)=ah*xx**6
         xx=float(i-2)+.5
         xx=xx/npm
         DiXoTR(j,i)=ah*xx**6
        enddo
       enddo 
       return
       end
