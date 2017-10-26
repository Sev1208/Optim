      program bmo
      implicit double precision (a-h,o-z)
      parameter(ndimm=10000)
      double precision v(ndimm),v0(ndimm)
      double precision xmin(ndimm),xmax(ndimm)
      double precision fpx(ndimm),x1(ndimm),fpx0(ndimm)
      double precision hgc(ndimm)
      double precision vinit(ndimm),temp(ndimm)
      double precision xsave(ndimm),xsave0(ndimm)
      double precision cstr,cstropt
      common/comcstr/cstr(100),cstropt(100),fseul,fseulopt,ncstr
      common/comsave/costsaveming,gnormsave,xoptg(ndimm),costsaveming0
      common/nbeval/nbeval,nbevalp
c
      ncstr=0
c
      nbeval=0
      nbevalp=0
c
ccccccccccccccccccccccccccccccccccccccccccccccccc
c
      open(1,file='DATA_BMO.txt',status='unknown')
      read(1,*) nbrestart
      read(1,*) nbext1
      read(1,*) nbbvp
      read(1,*) nbgrad
      read(1,*) ifd,epsfd000
      read(1,*) rho000
      read(1,*) epsij
      close(1)
c       nbrestart=10
c       nbext1=4
c       nbbvp=4
c       nbgrad=4
c       ifd=1
c       epsfd000=0.9
c       rho000=1.
c       epsij=1.e-4
c
      call init(ndim,vinit,xmin,xmax)
      if(ndim.gt.ndimm) stop 'ndim'
c
       call func(ndim,vinit,finit)
       f=finit
       f0=finit
c
          do ii=1,ndim
           xoptg(ii)=vinit(ii)
          enddo
          if(finit.lt.epsij) then
           do ii=1,ncstr
             cstropt(ii)=cstr(ii)
           enddo
           fseulopt=fseul
           goto 4001
          endif
c
          epsij=epsij*finit
c
          print *,' F_init = ',finit
        if(ncstr.gt.0) print *,'CSTR = ',(cstr(ii),ii=1,ncstr)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
        open(2,file='hist.J',status='unknown')
        open(3,file='hist.JC',status='unknown')
c
        write(3,999) fseul,(cstr(ii),ii=1,ncstr)
        write(2,*) finit,1.,xoptg(1),xoptg(2),xoptg(3),xoptg(4)
        close(2)
c
      costsaveming = finit
      costsaveming0 = finit
      do 4000 irestart=1,nbrestart
       do ii=1,ndim
        xsave(ii)=xoptg(ii)
       enddo
        rho00=max(rho000/2**(irestart-1),1.e-5)
        epsfd00=max(epsfd000/2**(irestart-1),1.e-4)

      do 3000 iter1=1,nbext1
       do ii=1,ndim
        v(ii)=xsave(ii)
       enddo
        rho0=max(rho00/2**(iter1-1),1.e-5)
        epsfd0=epsfd00
       do 2000 iterbvp=1,nbbvp
c
       do ii=1,ndim
        x1(ii)=v(ii)
       enddo
c
        rho=max(rho0/2**(iterbvp-1),1.e-5)
        epsfd=max(epsfd0/2**(iterbvp-1),1.e-4)
c
        call gradopt(ifd,epsfd,nbgrad,ndim,x1,fpx,
     1  xmin,xmax,temp,rho,f,gnorm,epsij,fpx0,hgc,irestart,iterbvp)
c
         if(costsaveming.lt.epsij) goto 2001
c
         call tir(ndim,iterbvp,f,f0,v,v0,temp,xmin,xmax,fpx,itir)
c
        if(itir.eq.0) goto 2001
        do ii=1,ndim
         v0(ii)=v(ii)
         v(ii)=temp(ii)
        enddo
         f0=f
c
 2000   continue
 2001   continue
c
        if(costsaveming.lt.epsij) goto 3001
c
        costsave=f
c
        call tir(ndim,iter1,costsave,costsave0,
     1            xsave,xsave0,temp,xmin,xmax,fpx,itir)
c
        if(itir.eq.0) goto 3001
          costsave0=costsave
          do ii=1,ndim
          xsave0(ii)=xsave(ii)
          xsave(ii)=temp(ii)
          enddo
c
 3000   continue
 3001    continue
        xtest=abs(costsaveming-costsaveming0)/
     1     abs(costsaveming+costsaveming0)
      print *,'ITER = ',irestart,' FM = ',costsaveming,' RES = ',xtest
        if(costsaveming.lt.epsij) goto 4001
       if(irestart.gt.1.and.xtest.lt.1.e-4) goto 4001
4000    continue
4001    continue
c
       call result(ndim,xoptg,vinit)
c
        close(2)
        close(3)
c
        if(ncstr.gt.0) then
        print *,'-------------------------------------------'
         print *,'F seul = ',fseulopt
        print *,'-------------------------------------------'
         print *,'CSTR = ',(cstropt(ii),ii=1,ncstr)
        endif
        print *,'-------------------------------------------'
        if(ndim.lt.20) print *,' x = ',(xoptg(ii),ii=1,ndim)
        print *,'-------------------------------------------'
 999     format(100(e19.7,1x))
         print *,'NB FUNCTIONAL J EVALUATIONS '
         print *,'nbeval (sequential, for J, rho_opt)= ',nbeval
         print *,'nbevalp (parallel, for Grad J) = ',nbevalp
         print *,'SEQUENTIAL COMPLEXITY (# eval J) = '
     1            ,int(nbeval+nbevalp/(1.*ndim))
        print *,'-------------------------------------------'

c
c        call echant(ndim,xoptg,xmin,xmax)
c
       end
c
         subroutine tir(ndim,iter,f,f0,v,v0,temp,xmin,xmax,fpx,itir)
         implicit double precision (a-h,o-z)  
         double precision v(*),v0(*),xmin(*),xmax(*),fpx(*),temp(*)
         common/nbeval/nbeval,nbevalp
c
            itir=1
         if(iter.ge.2) then
c
           if(abs(f-f0)/abs(f+f0).gt.1.e-5) then
            do ii=1,ndim
             fpx(ii)=-f/(f-f0)*(v(ii)-v0(ii))
            temp(ii)=v(ii)+fpx(ii)
            temp(ii)=max(min(temp(ii),xmax(ii)),xmin(ii))
            enddo

           else
            itir=0
            print *,'tir degenere'
            return
           endif
c
         else
c
            do ii=1,ndim
          call xrandme(iter+nbeval+nbevalp,xrdran)
             temp(ii)=xmin(ii)+xrdran*(xmax(ii)-xmin(ii))
             temp(ii)=max(min(temp(ii),xmax(ii)),xmin(ii))
            enddo

         endif
c
         end
c
      subroutine xrandme(ii,xrd)
      implicit double precision (a-h,o-z)
c
       xilim=dble(2**30)
       xrd2=dabs(irand(ii)/xilim)
       xrd3=dabs(irand(irand(ii))/xilim)
       xrd4=dabs(irand(irand(irand(ii)))/xilim)
       xrd5=dabs(irand(irand(irand(irand(ii))))/xilim)
       xrd6=dabs(irand(irand(irand(irand(irand(ii)))))/xilim)
       xrd=(xrd2+xrd3+xrd4+xrd5+xrd6)/10.
      end
c
      subroutine gradopt(ifd,epsfd,nbgrad,ndim,x1,fpx,xmin,xmax,
     1       temp,rho,f,gnorm,epsij,fpx0,hgc,irestart,iterbvp)
      parameter(ndimm=10000)
      implicit double precision (a-h,o-z)
      double precision fpx(*),xmin(*),xmax(*),temp(*),x1(*)
      double precision fpx0(*),hgc(*)
      double precision cstr,cstropt
       common/nbeval/nbeval,nbevalp
      common/comcstr/cstr(100),cstropt(100),fseul,fseulopt,ncstr
      common/comsave/costsaveming,gnormsave,xoptg(ndimm),costsaveming0
c
      if(irestart.gt.1.or.iterbvp.gt.1) call func(ndim,x1,f)
      epscv=1.e-3/2**(irestart-1)
c
      f00=f
      igc=1
      do igr=1,nbgrad
       f0=f
c
       xnorm=0.
      do ii=1,ndim
         fpx0(ii)=fpx(ii)
         xnorm=xnorm+fpx0(ii)**2
      enddo
c
       nbeval=-nbeval
      call funcp(ifd,epsfd,ndim,x1,xmin,xmax,f,fpx)
       nbeval=-nbeval
c
      gamgc=0.
      if(igc.eq.1.and.igr.ge.2.and.xnorm.gt.1.e-10) then
      do ii=1,ndim
         gamgc=gamgc+(fpx(ii)-fpx0(ii))*fpx(ii)/xnorm
      enddo
      endif
c
      do ii=1,ndim
         hgc(ii)=fpx(ii)+gamgc*hgc(ii)
      enddo
c
       call ropt_dicho(ndim,x1,xmin,xmax,temp,rho,hgc,f)
c
       do ii=1,ndim
       x0=x1(ii)
       x1(ii)=x1(ii)-rho*hgc(ii)
       x1(ii)=min(x1(ii),xmax(ii))
       x1(ii)=max(x1(ii),xmin(ii))
       enddo
c
       call func(ndim,x1,f)
c
       gnorm=0.
       do ii=1,ndim
       gnorm=gnorm+fpx(ii)**2
       enddo
       gnorm=sqrt(gnorm)
       if(igr.eq.1)   gnorm0=gnorm   
       if(gnorm0.lt.1.e-6) return
       gnorm=gnorm/gnorm0
c
        open(2,file='hist.J',access='append')
        write(2,*) f,gnorm*gnorm0,x1(1),x1(2),x1(3),x1(4)
        close(2)
        write(3,999) fseul,(cstr(ii),ii=1,ncstr)
c
        if(f.lt.costsaveming) then
          costsaveming0=costsaveming    
          costsaveming=f
          gnormsave=gnorm
          do ii=1,ncstr
             cstropt(ii)=cstr(ii)
          enddo
          fseulopt=fseul
          do ii=1,ndim
           xoptg(ii)=x1(ii)
          enddo
         endif
c
!        if(f.lt.epsij.or.abs(f/f00).lt.epscv
!     1        .or.abs(f-f0)/f00.lt.1.e-10) goto 888

c
999      format(9(e19.7,1x))
!       if(gnorm.lt.1.e-10.or.gnorm*gnorm0.lt.1.e-10) goto 888
c
       enddo

888     return       
        end
c
      subroutine ropt_dicho(n,x,xmin,xmax,temp,ro,g,cout)
      implicit double precision (a-h,o-z)
      double precision romin(3),fmin(3)
      double precision x(*),g(*),temp(*),xmin(*),xmax(*)
c
      numi=0
      numimax=5
 240  continue
      romin(1)=ro*0.5
      romin(2)=ro
      romin(3)=ro*2.
      l=0
  300 l=l+1
      numi=numi+1
      call fun(n,x,xmin,xmax,temp,g,romin(l),fmin(l))
      
      if (l.ne.2) goto 352
      if (fmin(1).lt.fmin(2)) goto 380
      goto 355
  352 continue
      if (l.ne.1) goto 355
      if (fmin(1).lt.cout) goto 355
      ro=ro*0.5
c******  test d'arret
      if(abs(ro).lt.1.e-5.or.numi.gt.numimax) go to 500
c******
      goto 240
  355 continue
      if (l.lt.3) goto 300
  360 if (fmin(1).lt.fmin(2)) goto 380
  370 if (fmin(2)-fmin(3)) 450,450,420
  380 romin(3)=romin(2)
      fmin(3)=fmin(2)
      romin(2)=romin(1)
      fmin(2)=fmin(1)
      romin(1)=romin(1)*0.5
      numi=numi+1
      call fun(n,x,xmin,xmax,temp,g,romin(1),fmin(1))
      goto 360
 420  continue
  425 romin(1)=romin(2)
      fmin(1)=fmin(2)
      romin(2)=romin(3)
      fmin(2)=fmin(3)
      romin(3)=romin(3)*2.
      numi=numi+1
      call fun(n,x,xmin,xmax,temp,g,romin(3),fmin(3))
      goto 370
  450 ro=romin(2)
       if(2*abs(fmin(2)-fmin(3))/(fmin(2)+fmin(3)).lt.1.e-4
     1       .or.numi.gt.numimax) go to 500
  451 continue
c****** calcul de ro interpole
      sn=0.
      sd=0.
      do 454 i=1,3,1
      s=0.
      pr=1.
      do 452 j=1,3
      if (i.eq.j) goto 452
      s=s+romin(j)
      pr=pr*(romin(i)-romin(j))
  452 continue
      sn=sn+fmin(i)*s/pr
  454 sd=sd+fmin(i)/pr
      ro=sn/sd/2.
  460 continue
  468 continue
 500  continue
      call fun(n,x,xmin,xmax,temp,g,ro,fm)
      cout=fm
      if(fm.gt.fmin(2)) then
              ro=romin(2)
              cout=fmin(2)
      endif
c
      end
c
      subroutine fun(n,x,xmin,xmax,temp,g,ro,f)
      implicit double precision (a-h,o-z)
      double precision x(*),g(*),temp(*),xmin(*),xmax(*)
      do ii=1,n
        temp(ii)=x(ii)-ro*g(ii)
        temp(ii)=max(min(temp(ii),xmax(ii)),xmin(ii))
      enddo
      call func(n,temp,f)
      end
c
      subroutine funcp(ifd,epsfd,ndim,x,xmin,xmax,f,fpx)
      implicit double precision (a-h,o-z)
      double precision fpx(*),x(*),xmin(*),xmax(*)
      double precision fpx0(1000),fpx1(1000)
c
       if(ifd.eq.0) then
       call funcpex(ndim,x,fpx)
       elseif(ifd.eq.1) then
c
c boucle à paralleliser
c
        do ii=1,ndim       
         fpx(ii)=0
        enddo
        do ii=1,ndim       
        x00=x(ii)
        dx0=xmax(ii)-xmin(ii)
        dx=max(1.e-4*dx0,dx0*epsfd)
c        epsifd=max(min(abs(x00)*0.5,dx,dx0*0.2),dx0*1.e-4,5.e-5)		
		epsifd=epsfd
c        if(x00+epsifd.le.xmax(ii)) then
        x(ii)=x00+epsifd
        call func(ndim,x,fp)
c		else
        x(ii)=x00-2*epsifd
        call func(ndim,x,fm)
c        epsifd=-epsifd
c        endif
        fpx(ii)=(fp-fm)/(2*epsifd)
        x(ii)=x00
        enddo
       endif
c
      open(9,file='df.data',status='unknown')
         do ii=1,ndim
        write(9,*) fpx(ii)
        enddo
      close(9)
c
      end
c
      subroutine echant(ndim,xinit,xmin,xmax)
      implicit double precision (a-h,o-z)
      double precision xmin(*),xmax(*),xx(1000),xinit(*)
c
      nx=20
      i1=1
      i2=2
c
      do ii=1,ndim
         xx(ii)=xinit(ii)
      enddo
c
       x0=xmin(i1)
       y0=xmin(i2)
       x1=xmax(i1)
       y1=xmax(i2)
       dx=(x1-x0)
       dy=(y1-y0)
       xh=dx/(nx-1)
       yh=dy/(nx-1)

c
       ifile=7
       do ix=1,nx
          do iy=1,nx
           xx(i1)=x0+xh*(ix-1)
           xx(i2)=y0+yh*(iy-1)
           call func(ndim,xx,f)
         write(ifile,*) xx(i1),xx(i2),f
c         write(ifile,*) (xx(i1)-x0)/dx,
c     1                  (xx(i2)-y0)/dy,log(abs(f))
           enddo
         write(ifile,'()')
       enddo
c
      end
