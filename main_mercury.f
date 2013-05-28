
      implicit none
      include 'mercury.inc'
c
      integer j,algor,nbod,nbig,opt(8),stat(NMAX),lmem(NMESS)
      integer opflag,ngflag,ndump,nfun
      real*8 m(NMAX),xh(3,NMAX),vh(3,NMAX),s(3,NMAX),rho(NMAX)
      real*8 rceh(NMAX),epoch(NMAX),ngf(4,NMAX),rmax,rcen,jcen(3)
      real*8 cefac,time,tstart,tstop,dtout,h0,tol,en(3),am(3)
      character*8 id(NMAX)
      character*80 outfile(3), dumpfile(4), mem(NMESS)
      external mdt_mvs, mdt_bs1, mdt_bs2, mdt_ra15, mdt_hy, mdt_cb
      external mdt_wb,mco_dh2h,mco_h2dh,mco_wb2h,mco_h2wb,mco_cb2h
      external mco_h2cb,mco_b2h,mco_h2b,mco_h2mvs,mco_mvs2h,mco_iden
c
      data opt/0,1,1,2,0,1,0,0/
c
c------------------------------------------------------------------------------
c
c Get initial conditions and integration parameters
      call mio_in (time,tstart,tstop,dtout,algor,h0,tol,rmax,rcen,jcen,
     %  en,am,cefac,ndump,nfun,nbod,nbig,m,xh,vh,s,rho,rceh,stat,id,
     %  epoch,ngf,opt,opflag,ngflag,outfile,dumpfile,lmem,mem)
c
c If this is a new integration, integrate all the objects to a common epoch.
      if (opflag.eq.-2) then
  20    open (23,file=outfile(3),status='old',access='append',err=20)
        write (23,'(/,a)') mem(55)(1:lmem(55))
        write (*,'(a)') mem(55)(1:lmem(55))
        call mxx_sync (time,tstart,h0,tol,jcen,nbod,nbig,m,xh,vh,s,rho,
     %    rceh,stat,id,epoch,ngf,opt,ngflag)
        write (23,'(/,a,/)') mem(56)(1:lmem(56))
        write (*,'(a)') mem(56)(1:lmem(56))
        opflag = -1
        close (23)
      end if
c
c Main integration
      if (algor.eq.1) call mal_hcon (time,tstart,tstop,dtout,algor,h0,
     %  tol,jcen,rcen,rmax,en,am,cefac,ndump,nfun,nbod,nbig,m,xh,vh,s,
     %  rho,rceh,stat,id,ngf,opt,opflag,ngflag,outfile,dumpfile,mem,
     %  lmem,mdt_mvs,mco_h2mvs,mco_mvs2h)
c
      if (algor.eq.9) call mal_hcon (time,tstart,tstop,dtout,algor,h0,
     %  tol,jcen,rcen,rmax,en,am,cefac,ndump,nfun,nbod,nbig,m,xh,vh,s,
     %  rho,rceh,stat,id,ngf,opt,opflag,ngflag,outfile,dumpfile,mem,
     %  lmem,mdt_mvs,mco_iden,mco_iden)
c
      if (algor.eq.2) call mal_hvar (time,tstart,tstop,dtout,algor,h0,
     %  tol,jcen,rcen,rmax,en,am,cefac,ndump,nfun,nbod,nbig,m,xh,vh,s,
     %  rho,rceh,stat,id,ngf,opt,opflag,ngflag,outfile,dumpfile,mem,
     %  lmem,mdt_bs1)
c
      if (algor.eq.3) call mal_hvar (time,tstart,tstop,dtout,algor,h0,
     %  tol,jcen,rcen,rmax,en,am,cefac,ndump,nfun,nbod,nbig,m,xh,vh,s,
     %  rho,rceh,stat,id,ngf,opt,opflag,ngflag,outfile,dumpfile,mem,
     %  lmem,mdt_bs2)
c
      if (algor.eq.4) call mal_hvar (time,tstart,tstop,dtout,algor,h0,
     %  tol,jcen,rcen,rmax,en,am,cefac,ndump,nfun,nbod,nbig,m,xh,vh,s,
     %  rho,rceh,stat,id,ngf,opt,opflag,ngflag,outfile,dumpfile,mem,
     %  lmem,mdt_ra15)
c
      if (algor.eq.10) call mal_hcon (time,tstart,tstop,dtout,algor,h0,
     %  tol,jcen,rcen,rmax,en,am,cefac,ndump,nfun,nbod,nbig,m,xh,vh,s,
     %  rho,rceh,stat,id,ngf,opt,opflag,ngflag,outfile,dumpfile,mem,
     %  lmem,mdt_hy,mco_h2dh,mco_dh2h)
c
      if (algor.eq.11) call mal_hcon (time,tstart,tstop,dtout,algor,h0,
     %  tol,jcen,rcen,rmax,en,am,cefac,ndump,nfun,nbod,nbig,m,xh,vh,s,
     %  rho,rceh,stat,id,ngf,opt,opflag,ngflag,outfile,dumpfile,mem,
     %  lmem,mdt_cb,mco_h2cb,mco_cb2h)
c
      if (algor.eq.12) call mal_hcon (time,tstart,tstop,dtout,algor,h0,
     %  tol,jcen,rcen,rmax,en,am,cefac,ndump,nfun,nbod,nbig,m,xh,vh,s,
     %  rho,rceh,stat,id,ngf,opt,opflag,ngflag,outfile,dumpfile,mem,
     %  lmem,mdt_wb,mco_h2wb,mco_wb2h)
c
c Do a final data dump
      do j = 2, nbod
        epoch(j) = time
      end do
      call mio_dump (time,tstart,tstop,dtout,algor,h0,tol,jcen,rcen,
     %  rmax,en,am,cefac,ndump,nfun,nbod,nbig,m,xh,vh,s,rho,rceh,stat,
     %  id,ngf,epoch,opt,opflag,dumpfile,mem,lmem)
c
c Calculate and record the overall change in energy and ang. momentum
  50  open  (23, file=outfile(3), status='old', access='append',
     %  err=50)
      write (23,'(/,a)') mem(57)(1:lmem(57))
      call mxx_en (jcen,nbod,nbig,m,xh,vh,s,en(2),am(2))
c
      write (23,231) mem(58)(1:lmem(58)),
     %  abs((en(2) + en(3) - en(1)) / en(1))
      write (23,232) mem(59)(1:lmem(59)),
     %  abs((am(2) + am(3) - am(1)) / am(1))
      write (23,231) mem(60)(1:lmem(60)), abs(en(3) / en(1))
      write (23,232) mem(61)(1:lmem(61)), abs(am(3) / am(1))
      close (23)
      write (*,'(a)') mem(57)(1:lmem(57))
c
c------------------------------------------------------------------------------
c
 231  format (/,a,1p1e12.5)
 232  format (a,1p1e12.5)
      stop
      end
