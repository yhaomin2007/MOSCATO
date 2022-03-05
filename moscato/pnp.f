c subroutines for pnp solver
c-----------------------------------------------------------------------
      subroutine cal_potential_terms
      include 'SIZE'
      include 'TOTAL'
      include "moscato/ECM"
      include "moscato/CASE"

      real phie_max,phie_min
      integer ips

      ntot = lx1*ly1*lz1*nelv

      if(istep.eq.0) call rzero(phie,ntot)
      call rzero(kappae,ntot)
      call rzero(rhs_phie,ntot)
     
      do ips = 1,nion
      phie_max = glmax(t(1,1,1,1,ips+1),ntot)
      phie_min = glmin(t(1,1,1,1,ips+1),ntot)
      if(nid.eq.0) write(6,*) 'c max/min : ', phie_max,'/',phie_min

c      call lap_usr1(t(1,1,1,1,ips+1),DiffCoeff(ips),lapc(1,1,1,1,ips)) 
      call lap_usr2(t(1,1,1,1,ips+1),DiffCoeff(1,1,1,1,ips),
     & lapc(1,1,1,1,ips)) 

      do i = 1,ntot 

        kappae(i,1,1,1) = kappae(i,1,1,1)
     &+F_const*F_const*charge(ips)*charge(ips)
     &*mobility(i,1,1,1,ips)*t(i,1,1,1,ips+1)

        rhs_phie(i,1,1,1) = rhs_phie(i,1,1,1) 
     & + F_const*charge(ips)*lapc(i,1,1,1,ips)
      enddo

      enddo

      phie_max = glmax(kappae,ntot)
      phie_min = glmin(kappae,ntot)
      if(nid.eq.0) then
       write(6,*) 'kappae max/min : ', phie_max,'/',phie_min
      endif
      phie_max = glmax(rhs_phie,ntot)
      phie_min = glmin(rhs_phie,ntot)
      if(nid.eq.0) write(6,*)'rhs_phie max/min:',phie_max,'/',phie_min

      return
      end
c-----------------------------------------------------------------------
      subroutine cal_pnp_migration()
c calcualte pnp migraion term for each ion
c mgr = div (zuFc) grad Phie

      include 'SIZE'
      include 'TOTAL'
      include "moscato/ECM"
      include "moscato/CASE"

      integer ips

      ntot = lx1*ly1*lz1*nelv

      do ips = 1,nion
        call lap_usr2(phie,t(1,1,1,1,ips+1),migration(1,1,1,1,ips))

        do i = 1,ntot
        migration(i,1,1,1,ips) =  migration(i,1,1,1,ips)
     & *charge(ips)*mobility(i,1,1,1,ips)*F_const
        enddo
		
      phie_max = glmax(migration(1,1,1,1,ips),ntot)
      phie_min = glmin(migration(1,1,1,1,ips),ntot)
      if(nid.eq.0) then
       write(6,*) 'migration max/min : ', phie_max,'/',phie_min
      endif
		
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine lap_usr1(aa,bb,cc)
c return cc = div(bb,grad(aa)), bb is constant

      include 'SIZE'
      include 'TOTAL'

      common /scrns/ w1(lx1,ly1,lz1,lelt)
     $              ,w2(lx1,ly1,lz1,lelt)
     $              ,w3(lx1,ly1,lz1,lelt)
     $              ,tx(lx1,ly1,lz1,lelt)
     $              ,ty(lx1,ly1,lz1,lelt)
     $              ,tz(lx1,ly1,lz1,lelt)

      real aa(1),cc(1)
      real bb

      ntot = lx1*ly1*lz1*nelv
      ifld_save = ifield
      ifield = 1
 
      call opgrad  (tx,ty,tz,aa)
      call opdssum (tx,ty,tz)
      call opcolv  (tx,ty,tz,binvm1)
	  
      do i = 1, ntot
       w1(i,1,1,1) = bb
      enddo

      call opcolv  (tx,ty,tz,w1)
      call opdiv   (w2,tx,ty,tz)

      call dssum   (w2,lx1,ly1,lz1)
      call col2    (w2,binvm1,ntot)

      call copy(cc,w2,ntot)

      ifield = ifield_save

      return
      end
c-----------------------------------------------------------------------
      subroutine lap_usr2(aa,bb,cc)
c return c = div(b,grad(a)), b is variable

      include 'SIZE'
      include 'TOTAL'
	  
       common /scrns/ w1(lx1,ly1,lz1,lelt)
     $              ,w2(lx1,ly1,lz1,lelt)
     $              ,w3(lx1,ly1,lz1,lelt)
     $              ,tx(lx1,ly1,lz1,lelt)
     $              ,ty(lx1,ly1,lz1,lelt)
     $              ,tz(lx1,ly1,lz1,lelt)
      real aa(1),bb(1),cc(1)
	 
      ntot = lx1*ly1*lz1*nelv
      ifld_save = ifield
      ifield = 1

      call opgrad  (tx,ty,tz,aa)
      call opdssum (tx,ty,tz)
      call opcolv  (tx,ty,tz,binvm1)

      call opcolv  (tx,ty,tz,bb)
      call opdiv   (w2,tx,ty,tz)

      call dssum   (w2,lx1,ly1,lz1)
      call col2    (w2,binvm1,ntot)

      call copy(cc,w2,ntot)

      ifield = ifield_save

      return
      end
c-----------------------------------------------------------------------
      subroutine cal_inertial()
c re-calculate inertial ion concentration to ensure electroneutrality
c
      include 'SIZE'
      include 'TOTAL'
      include "moscato/ECM"

      real en(lt)
      integer ips,i
      ntot = lx1*ly1*lz1*nelv

      call rzero(en,ntot)

      do i = 1,ntot  
      do ips = 1,nion
        if (t(i,1,1,1,ips+1).lt.0.0) t(i,1,1,1,ips+1) = 0.0
      enddo
      enddo
	  
      do i = 1,ntot  
      do ips = 1,nion
         en(i) = en(i) + charge(ips)*t(i,1,1,1,ips+1)
      enddo
      enddo

      en_max = glmax(en,ntot)
      en_min = glmin(en,ntot)
	  
      if(nid.eq.0) write(6,*) 'en max/min : ', en_max,'/',en_min
	  
      call rzero(en,ntot)

      do i = 1,ntot  
      do ips = 1,nion
        if (ips.ne.inertial) en(i)=en(i)+charge(ips)*t(i,1,1,1,ips+1)
      enddo
        t(i,1,1,1,inertial+1) = -en(i)/charge(inertial)
      enddo
 
      return
      end
c-----------------------------------------------------------------------
      subroutine bc_for_potential_equation()
c calculate bc for potential equation 
c flux of inertial ion is zero. 
c flux = -D grad c - z*m*F*c grad phie = 0 at boundaries.
c
      include 'SIZE'
      include 'TOTAL'

      include "moscato/ECM"

      integer ips

      integer i,j,k,e,f,idf,ie,ifc
      integer i0,i1,j0,j1,k0,k1

      real gradT(lx1,ly1,lz1,lelt,3)
      real nGradT(lx1,ly1,lz1,lelt)
      real srnl(3),alpha,beta
  
      ntot = lx1*ly1*lz1*nelv
  
      call rzero(gradT,3*ntot)
      call rzero(nGradT,ntot)
      call rzero(nGradPhie,ntot)
  
c calcualte gradient inertial ion
        call opgrad(gradT(1,1,1,1,1),gradT(1,1,1,1,2),
     &              gradT(1,1,1,1,3),t(1,1,1,1,inertial+1))
        call opdssum(gradT(1,1,1,1,1),gradT(1,1,1,1,2),
     &               gradT(1,1,1,1,3))
        call opcolv(gradT(1,1,1,1,1),gradT(1,1,1,1,2),
     &              gradT(1,1,1,1,3),binvm1)

c calcualte normal gradient t at all wall
      do ie = 1,nelt
       do ifc = 1,2*ndim
        if(cbc(ifc,ie,1).eq.'W  ') then
		call facind(i0,i1,j0,j1,k0,k1,nx1,ny1,nz1,ifc)
        do k=k0,k1
        do j=j0,j1
        do i=i0,i1
         call getSnormal(srnl,i,j,k,ifc,ie)

         ! srnl is outward of the element,
         ! revert the difnition here. nGradT intward the element

         nGradT(i,j,k,ie) =  - gradT(i,j,k,ie,1)*srnl(1) 
     & - gradT(i,j,k,ie,2)*srnl(2) 
         if(ndim.eq.3)  nGradT(i,j,k,ie) =  
     & nGradT(i,j,k,ie)- gradT(i,j,k,ie,3)*srnl(3)

        enddo
        enddo
        enddo
        endif
       enddo
      enddo

      nGradT_max = glmax(nGradT,ntot)
      nGradT_min = glmin(nGradT,ntot)
      if(nid.eq.0) then
      write(6,*) 'inertial ion nGradT max/min : ', 
     & nGradT_max,'/',nGradT_min
      endif


cc translate to normal gradiant of Phie
      do ie = 1,nelt
       do ifc = 1,2*ndim
        if(cbc(ifc,ie,1).eq.'W  ') then
		call facind(i0,i1,j0,j1,k0,k1,nx1,ny1,nz1,ifc)
        do k=k0,k1
        do j=j0,j1
        do i=i0,i1

       ! determine surface normal gradiant of phie based on 
       ! inertial ion flux zero on wall 
!       nGradPhie(i,j,k,ie) = DiffCoeff(i,j,k,ie,inertial)
!     & *nGradT(i,j,k,ie) 
!     & / (charge(inertial)*mobility(i,j,k,ie,inertial)
!     & *F_const*t(i,j,k,ie,inertial+1))

      nGradPhie(i,j,k,ie) = currentDensity(i,j,k,ie)/kappae(i,j,k,ie)
	 
!       nGradPhie(i,j,k,ie) = -DiffCoeff(i,j,k,ie,inertial)
!     & *nGradT(i,j,k,ie) 
!     & / (charge(inertial)*mobility(i,j,k,ie,inertial)
!     & *F_const*t(i,j,k,ie,inertial+1))
	 
	
! sum(z_i*grad(c_i)) = 
! 
!        alpha = 0.0	
!        beta =  0.0
!         do iion = 1,nion
!
!         alpha = alpha 
!     & + charge(iion)*charge(iion)*mobility(i,j,k,ie,iion)*F_const
!     & *t(i,j,k,ie,iion+1)/DiffCoeff(i,j,k,ie,iion)
!
!         beta = beta 
!     & + charge(iion)*ionflux(i,j,k,ie,iion)
!     & /DiffCoeff(i,j,k,ie,iion)
!
!         enddo		 
!
!        nGradPhie(i,j,k,ie)  = beta/alpha

        enddo
        enddo
        enddo
        endif
       enddo
      enddo
 
      nGradPhie_max = glmax(nGradPhie,ntot)
      nGradPhie_min = glmin(nGradPhie,ntot)
      if(nid.eq.0) then
      write(6,*) 'nGradPhie max/min : ', 
     & nGradPhie_max,'/',nGradPhie_min
      endif
 
      return
      end
c----------------------------------------------------------------------
      subroutine dump_ecm_data
      include 'SIZE'
      include 'TOTAL'
      include "moscato/ECM"

      iastep = param(68)
      if  (iastep.eq.0) iastep=param(15)   ! same as iostep

      if ( (mod(istep,iastep).eq.0.and.istep.gt.1) .or.lastep.eq.1) then
        ! for 2d case, vz position is not dumped.
        ! 
      call outpost2(phie,kappae,
     & rhs_phie,rhs_phie,migration,nion,'ecm')

c      call outpost(t(1,1,1,1,2),t(1,1,1,1,3),
c     & rhs_phie,t(1,1,1,1,4),t(1,1,1,1,5),'ccc')
c      call outpost(lapc(1,1,1,1,1),lapc(1,1,1,1,2),
c     & rhs_phie,lapc(1,1,1,1,3),lapc(1,1,1,1,4),'lpc')
c      call outpost(migration(1,1,1,1,1),migration(1,1,1,1,2),
c     & rhs_phie,migration(1,1,1,1,3),migration(1,1,1,1,4),'mgr')

      endif

      return
      end
c---------------------------------------------------------------------
      subroutine pnp_poisson_test
c     
      include 'SIZE'
      include 'TOTAL'
      include 'ORTHOT' ! hsolve proj
      include "moscato/ECM"

      real h1(lt)
     $   , h2(lt)
     $   , t1(lt)
     $   , t2(lt)
     $   , phi_exact(lt)

      n = lx1*ly1*lz1*nelt

      isd = 1
      imsh = 1
      maxit = 200
      idpss(ifldpot-1)=0 ! active solver
      ifield = ifldpot

c      call outpost(tmask(1,1,1,1,ifldpot-1),tmult(1,1,1,1,ifldpot-1)
c     $            ,t1,bpmask,pmask,'cbc')

      call vprops
      call sethlm(h1,h2,0)

      call copy(h1,kappae,n)

c      if(nid.eq.0) write(6,*)'h1 minmax',glmax(h1,n),glmin(h1,n)
c      if(nid.eq.0) write(6,*)'h2 minmax',glmax(h2,n),glmin(h2,n)

      call bcneusc (t2,-1) ! robin
      call Xaddcol3(h2,t2,h1,n) ! h2=h2+t2*h1

c      call rzero(phie,n) ! initial condition
c      call copy(phie,t(1,1,1,1,ifldpot-1),n)  ! from useric (exact solution)
      call bcdirsc(phie) ! dirichlet
      call axhelm (t1,phie,h1,h2,imsh,isd)
	  
c      phie_max = glmax(bq(1,1,1,1,ifldpot-1),n)
c      phie_min = glmin(bq(1,1,1,1,ifldpot-1),n)
c      if(nid.eq.0) then
c       write(6,*) 'bq(1..ifldpot-1) max/min : ', phie_max,'/',phie_min
c      endif
c	  
c      phie_max = glmax(rhs_phie,n)
c      phie_min = glmin(rhs_phie,n)
c      if(nid.eq.0) then
c       write(6,*) 'rhs_phie max/min : ', phie_max,'/',phie_min
c      endif
c
c      call setqvol(bq(1,1,1,1,ifldpot-1)) ! rhs
c      call col3(t2,bq(1,1,1,1,ifldpot-1),bm1,n)
	
      call setqvol(rhs_phie) ! rhs
      call col3(t2,rhs_phie,bm1,n)	
	
      call sub2(t2,t1,n)

      call bcneusc (t1,1) ! neumann
      call col2(t1,h1,n)
      call add2(t2,t1,n)

!     Solve helmholtz: [h1 A + h2 B] u_0 = rhs - h1 A u_b
      if (nio.eq.0) write(6,*) 'SOLVING MY POISSON:'
c      ifield = 1
c      call dssum(t2,lx1,ly1,lz1)
c      call col2 (t2,tmask(1,1,1,1,ifldpot-1),n)
c      ifield = ifldmhd
c      call hmh_gmres(t2,h1,h2,tmult(1,1,1,1,ifldpot-1),maxit)
c      call add2(phi,t2,n)

      ifield = ifldpot
      call hsolve   ('POTE',T1,T2,H1,H2
     $              ,tmask(1,1,1,1,ifield-1)
     $              ,tmult(1,1,1,1,ifield-1)
     $              ,imsh,tolht(ifield),maxit,isd
     $              ,approxt(1,0,ifield-1),napproxt(1,ifield-1)
     $              ,bintm1)
      call add2(phie,t1,n) ! u = u_0 + u_b

      ! explicit filter phie
      ncut = param(101)+1
      wght = param(103)
      call my_filter(phie,wght,ncut)

      call copy(t(1,1,1,1,ifldpot-1),phie,n)

      phie_max = glmax(phie,n)
      phie_min = glmin(phie,n)
      if(nid.eq.0) then
       write(6,*) 'phie max/min : ', phie_max,'/',phie_min
      endif

      idpss(ifldpot-1)=-1 ! disalbe in timestep

      return
      end
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      subroutine my_filter(u,wght,ncut)
      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer lxv,ncut,ifldt
      parameter(lxv=lx1-1)

      real u(lx1*ly1*lz1*lelt),wght
      real intv(lx1,lx1)
      save intv

c     working arrays
      common /ctmp0/ intw,intt
     $             , wk1,wk2
     $             , zgmv,wgtv,zgmp,wgtp,tmax(100),omax(103)
      real intw(lx1,lx1)
      real intt(lx1,lx1)
      real wk1  (lx1,lx1,lx1,lelt)
      real wk2  (lx1,lx1,lx1)
      real zgmv(lx1),wgtv(lx1),zgmp(lx1),wgtp(lx1)
      real tmax,omax

      integer icalld
      save    icalld
      data    icalld /0/

      if (icalld.eq.0) then
         icalld = 1
         call build_new_filter(intv,zgm1,lx1,ncut,wght,nio)
      endif

      ifldt  = ifield
      ifield = 2

      call filterq(u,intv,
     $                lx1,lz1,wk1,wk2,intt,if3d,tmax(100))

      ifield = ifldt   ! RESTORE ifield
      return
      end
c-----------------------------------------------------------------------