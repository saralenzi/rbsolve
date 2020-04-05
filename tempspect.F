C Horizontal temperature and coherence spectrum
#include "config.h"

      program main

#include "param.h"

      complex*16 v(nwax+1,2*nwaz+1,0:Nylmem+1)
      complex*16 t(nwax+1,2*nwaz+1,0:Nylmem+1)
      real*8 ss(Nx),ssx(Nx),ssy(Nx)
      real*8 umax,umin,dk
      real*8 vmax,vmin
      real xmin,xmax
      integer istart,iend,istep,it,jmax,jmin,jj,iargc
      character*12 file0
      character*45 file

      character cc*1,str*30

      real*8 ptan,qtan
      real*8 DD,Lx,Lz,Re,invRe,dt,dx,dz
      real*8 kxmin,kzmin,Ra,Kscalar,Pr,sum

      integer*4 nsave,ttot,nn(Nx),nnx(Nx),nny(Nx),nmin,nmax
      real*8 alfa(3),beta(3),gamma(3),zeta(3)

      common/ALLPARA/ptan,qtan,DD,Re,Lx,Lz,invRe,dt,dx,dz,
     +   kxmin,kzmin,nsave,ttot,alfa,beta,gamma,zeta,Ra,Kscalar,Pr

      real*8 kxm(nwax+1),kzm(2*nwaz+1),pwr
      complex*16 mm(nwax+1),nnn(2*nwaz+1)
      real*8 K2(nwax+1,2*nwaz+1),dum,rt,rs,rr
	integer ik



#ifndef NOMPI 
      print *,"Please compile with the NOMPI option!"
      stop
#endif      

      iam=1

      if(iargc().lt.3) then
         print *, 'Usage: spectrumh cc y istart [iend istep]'
	 print *,'Computes the horizontal spectrum of variable yy'
	 print *, 'y<=0 means all levels'
         stop
      endif
      call slice()


      call getarg(1,cc)

      call getarg(2,str)
      read(str,*) jj

      call getarg(3,str)
      read(str,*) istart

      if(iargc().gt.3) then
       call getarg(4,str)
       read(str,*) iend
       call getarg(5,str)
       read(str,*) istep
      else
        iend=istart
        istep=1
      endif

      call read_para(0)
      call Crea_wave(kxm,kzm,K2,mm,nnn,kxmin,kzmin)

c  use x scale as granularity	
	dk=kxm(2)
	print *,dk


	sum=0.

	    do i=1,4
              do k=1,4
   	        ik=sqrt(K2(i,k))/dk+1.5
	        print *,i-1,k-1,ik
	   
	       enddo
	    enddo

      do i=1,Nx
	ss(i)=0
      enddo
      do it=istart,iend,istep

      if(cc.eq.'t') then
             call readfield(t,'t',it)
         
	if(jj.le.0) then
c	Sum over shells
	 do j=-jj,Ny+1+jj
           do i=1,nwax+1
              do k=1,2*nwaz+1
   	 	 ik=sqrt(K2(i,k))/dk+1.5
                 pwr=dreal(dconjg(t(i,k,j))*t(i,k,j))
                 ss(ik)=ss(ik)+pwr
                 nn(ik)=nn(ik)+1
	      enddo
           enddo
	sum=sum+1
	  enddo
	else
	    do i=1,nwax+1
              do k=1,2*nwaz+1
                 ik=sqrt(K2(i,k))/dk+1.5
	         pwr=dreal(dconjg(t(i,k,jj))*t(i,k,jj))
                 ss(ik)=ss(ik)+pwr
                 nn(ik)=nn(ik)+1
              enddo
           enddo
	sum=sum+1
	endif
c if cc ne e MIXED: CO-SPECTRUM
        elseif(cc.eq.'m') then
        call readfield(t,'t',it)
        call readfield(v,'v',it)

        if(jj.le.0) then
c       Sum over shells
         do j=-jj,Ny+1+jj
           do i=1,nwax+1
              do k=1,2*nwaz+1
                 ik=sqrt(K2(i,k))/dk+1.5
                 pwr=dreal(dconjg(v(i,k,j))*t(i,k,j))
                 ss(ik)=ss(ik)+pwr
                 nn(ik)=nn(ik)+1
              enddo
           enddo
        sum=sum+1
          enddo
        else
            do i=1,nwax+1
              do k=1,2*nwaz+1
                 ik=sqrt(K2(i,k))/dk+1.5
                 pwr=dreal(dconjg(v(i,k,jj))*t(i,k,jj))
                 ss(ik)=ss(ik)+pwr
                 nn(ik)=nn(ik)+1
              enddo
           enddo
        sum=sum+1
        endif

	endif

      enddo

      call makenamebob(file0,cc,istart,'.spe')
      open(unit=20,file=file0,form='formatted')
      do i=2,Nx/2+1
	sum=sum+ss(i)
      enddo
      sum=sum*Nx*Nx
      do i=1,Nx/2+1
           write(20,'(i5,2e16.8)') (i-1),ss(i)/nn(i)*(i-1),sum
      enddo
      close(20)

      end



