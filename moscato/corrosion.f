c corrosion model subroutines
c-------------------------------------------------------------------
      subroutine corrosion
c
c calcualte ion flux for oxi and red based on electrode kinetics of corrosion
c
      include 'SIZE'
      include 'TOTAL'
      include "moscato/ECM"
      include "moscato/ALLOY"
      logical ifverbose
      real E1,E2,E3,C1,C2,dcde
      real deltaE
	  
      ntot = lx1*ly1*lz1*nelv

      deltaE = 0.01

c use newton iteration to find E_alloy that will result total current=0

      niter = 20
      do iter = 1,niter

      E1 = E_alloy
      E2 = E_alloy + deltaE
      call total_current(E1,C1)
	  
      if (abs(C1).lt.1e-8) then
        E_alloy = E1
        if(nid.eq.0) then
         write(6,*) 'curent converged'
         write(6,*) 'E_alloy:',E_alloy
         write(6,*) 'current:',C1
        endif
      exit
      endif
	  
      call total_current(E2,C2)

      dcde = (C2-C1)/deltaE

      E3 = E1 - C1/dcde    ! find next input

      deltaE = abs(E1-E3)/100.0
      E_alloy = E3
	  
      ifverbose = .FALSE.
      if (ifverbose) then
	  if(nid.eq.0) then
         write(6,*) 'E1:',E1,' E2:',E2,' E3:', E3
         write(6,*) 'C1:',C1,' C2:',C2
      endif
      endif

      enddo

      ! after ion flux calcualted, determine mass loss rate
      ! calculate mass loss is meaningless, we should calculate mass loss rate
      call local_mass_loss_rate()  

      call metalDiffusion()
	  
      call dump_salt_surface_corrosion_data()

      if (ifverbose) then
      flux_max = glmax(ionflux(1,1,1,1,ired),ntot)
      flux_min = glmin(ionflux(1,1,1,1,ired),ntot)
      if(nid.eq.0) then
       write(6,*) 'Red flux max/min : ', flux_max,'/',flux_min
      endif

      flux_max = glmax(ionflux(1,1,1,1,ioxi),ntot)
      flux_min = glmin(ionflux(1,1,1,1,ioxi),ntot)
      if(nid.eq.0) then
       write(6,*) 'Oxi flux max/min : ', flux_max,'/',flux_min
      endif
      endif

      return
      end
c-------------------------------------------------------------------
      subroutine map_alloy_salt_interface
c
c  map interaface concentration from salt to alloy
c 
c  map flux from alloy to salt
c
      include 'SIZE'
      include 'TOTAL'
      include "moscato/ECM"
      include "moscato/ALLOY"
      
      logical ifverbose 
      integer e,f,i,j,k,i0,i1,j0,j1,k0,k1
      integer ix,iy
	  
      ntot = lx1*ly1*lz1*nelv
      ntot1 = lx1*ly1*nfe*lelt

      do e=1,nelv
      do f=1,2*ldim
        if (cbc(f,e,1).eq.'W  ') then

        if (ldim.eq.3) then ! 3d

        call facind(i0,i1,j0,j1,k0,k1,nx1,ny1,nz1,f)
        do k=k0,k1
        do j=j0,j1
        do i=i0,i1

        if (k0.eq.k1) then
        ix = i
        iy = j
        elseif (j0.eq.j1) then
        ix = i
        iy = k
        elseif (i0.eq.i1) then
        ix = j
        iy = k
        endif

        c_oxi_salt(ix,iy,f,e) = t(i,j,k,e,ioxi+1)
        c_red_salt(ix,iy,f,e) = t(i,j,k,e,ired+1)      ! map salt concentration to alloy

        pxyz_alloy(1,ix,iy,f,e) = xm1(i,j,k,e)
        pxyz_alloy(2,ix,iy,f,e) = ym1(i,j,k,e)
        pxyz_alloy(3,ix,iy,f,e) = zm1(i,j,k,e)

        T_salt(ix,iy,f,e) = t(i,j,k,e,1)

        ionflux(i,j,k,e,ired) = flux_red(ix,iy,f,e)    ! map alloy flux to salt
        ionflux(i,j,k,e,ioxi) = flux_oxi(ix,iy,f,e)

		COAOSS(i,j,k,e) = c_oxi_alloy(1,ix,iy,f,e)   ! map alloy surface concentration to salt domain
        CRAOSS(i,j,k,e) = c_red_alloy(1,ix,iy,f,e)
        MCROSS(i,j,k,e) = MCR(ix,iy,f,e)

        enddo
        enddo
        enddo
		
        else  ! 2d

        call facind(i0,i1,j0,j1,k0,k1,nx1,ny1,nz1,f)
        k = 1
        do j=j0,j1
        do i=i0,i1

        if (j0.eq.j1) then
        ix = i
        elseif (i0.eq.i1) then
        ix = j
        endif
        iy = 1

        c_oxi_salt(ix,iy,f,e) = t(i,j,k,e,ioxi+1)
        c_red_salt(ix,iy,f,e) = t(i,j,k,e,ired+1)      ! map salt concentration to alloy

        pxyz_alloy(1,ix,iy,f,e) = xm1(i,j,k,e)
        pxyz_alloy(2,ix,iy,f,e) = ym1(i,j,k,e)
        pxyz_alloy(3,ix,iy,f,e) = zm1(i,j,k,e)

        T_salt(ix,iy,f,e) = t(i,j,k,e,1)

        ionflux(i,j,k,e,ired) = flux_red(ix,iy,f,e)    ! map alloy flux to salt
        ionflux(i,j,k,e,ioxi) = flux_oxi(ix,iy,f,e)

		COAOSS(i,j,k,e) = c_oxi_alloy(1,ix,iy,f,e)   ! map alloy surface concentration to salt domain
        CRAOSS(i,j,k,e) = c_red_alloy(1,ix,iy,f,e)
        MCROSS(i,j,k,e) = MCR(ix,iy,f,e)

        enddo
        enddo

        endif
		
       endif
      enddo
      enddo

     
      ifverbose = .FALSE.
      if (ifverbose) then

      c_max = glmax(c_Fe_salt,ntot1)
      c_min = glmin(c_Fe_salt,ntot1)
      if(nid.eq.0) then
       write(6,*) 'interface: c_Fe_salt max/min : ', c_max,'/',c_min
      endif
	  
      c_max = glmax(c_Cr_salt,ntot1)
      c_min = glmin(c_Cr_salt,ntot1)
      if(nid.eq.0) then
       write(6,*) 'interface: c_Cr_salt max/min : ', c_max,'/',c_min
      endif
	  
      c_max = glmax(pxyz_alloy,ntot1*3)
      c_min = glmin(pxyz_alloy,ntot1*3)
      if(nid.eq.0) then
       write(6,*) 'interface: pxyz_alloy max/min : ', c_max,'/',c_min
      endif
	  
      
      f_max = glmax(flux_Fe,ntot1)
      f_min = glmin(flux_Fe,ntot1)
      if(nid.eq.0) then
       write(6,*) 'interface: flux_Fe max/min : ', f_max,'/',f_min
      endif
	  
      f_max = glmax(flux_Cr,ntot1)
      f_min = glmin(flux_Cr,ntot1)
      if(nid.eq.0) then
       write(6,*) 'interface: flux_Cr max/min : ', f_max,'/',f_min
      endif


      c_max = glmax(t(1,1,1,1,2),ntot)
      c_min = glmin(t(1,1,1,1,2),ntot)
      if(nid.eq.0) then
       write(6,*) 'c_Cr max/min : ', c_max,'/',c_min
      endif
	  
      c_max = glmax(t(1,1,1,1,3),ntot)
      c_min = glmin(t(1,1,1,1,3),ntot)
      if(nid.eq.0) then
       write(6,*) 'c_Fe max/min : ', c_max,'/',c_min
      endif

      endif

      return
      end
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      subroutine total_current(EA,totCurrent)
c calcualte total current based on electrode potential
c  
c
      include 'SIZE'
      include 'TOTAL'
      include "moscato/ECM"
      include "moscato/ALLOY"
      include "moscato/CASE"

      integer e,f,ix,iy
      real EA,totCurrent,phi_eq_oxi,phi_eq_red
      real i0
      real totFluxR,totFluxO,sint2,sarea2
      real totMassChange,mass_flux

      logical ifbv

      ifbv = .True. ! if using Butler-Volmer kinetics 
                     ! otherwise use linear kinetics
      ! coeffcient for Butler-Volmer kinetics
       
      ! obtained fromm Weizel's thesis
      alpha_oxi = 0.32
      alpha_red = 0.09
      k0_oxi = 2.17e-10
      k0_red = 5.77e-9

      ! for linear kinetics
      !i0 = 1.0e-3 !  dummy i0
      i0 = 1.0e-2 !  dummy i0

c	  
      call map_alloy_salt_interface
	  
c determine flux of each redux and oxide	
      do e=1,nelv
      do f=1,2*ldim
        if (cbc(f,e,1).eq.'W  ') then

        if (ldim.eq.3) then
        ny_corrosion = lx1
        else
        ny_corrosion = 1
        endif
  
        do ix = 1,lx1
        do iy = 1,ny_corrosion

        ! calculate local equilibrium potential
        ! 
        if(ifCorK.eq.0) then
        phi_eq_oxi = E_eq_oxi_0 
     & + (R_const*(T_salt(ix,iy,f,e)+273.15)/(F_const*nct))
     & *log(c_oxi_salt(ix,iy,f,e)/c_oxi_alloy(1,ix,iy,f,e))
        phi_eq_red = E_eq_red_0 
     & + (R_const*(T_salt(ix,iy,f,e)+273.15)/(F_const*nct))
     & *log(c_red_salt(ix,iy,f,e)/c_red_alloy(1,ix,iy,f,e))
        else 
        phi_eq_oxi = E_eq_oxi_0 
     & + (R_const*(T_salt(ix,iy,f,e))/(F_const*nct))
     & *log(c_oxi_salt(ix,iy,f,e)/c_oxi_alloy(1,ix,iy,f,e))
        phi_eq_red = E_eq_red_0 
     & + (R_const*(T_salt(ix,iy,f,e))/(F_const*nct))
     & *log(c_red_salt(ix,iy,f,e)/c_red_alloy(1,ix,iy,f,e))
        endif

          if (ifbv) then
          i0_oxi = nct*F_const*k0_oxi
     & *c_oxi_salt(ix,iy,f,e)**(1.0-alpha_oxi)
     & *c_oxi_alloy(1,ix,iy,f,e)**alpha_oxi

          i0_red = nct*F_const*k0_red
     & *c_red_salt(ix,iy,f,e)**(1.0-alpha_red)
     & *c_red_alloy(1,ix,iy,f,e)**alpha_red

        if(ifCorK.eq.0) then
          op_oxi = EA - phi_eq_oxi
          flux_oxi(ix,iy,f,e)  = (i0_oxi/(F_const*nct))
     & *(exp((1.0-alpha_oxi)*nct*F_const*op_oxi
     & /(R_const*(T_salt(ix,iy,f,e)+273.15)))
     & -exp((-alpha_oxi)*nct*F_const*op_oxi
     & /(R_const*(T_salt(ix,iy,f,e)+273.15))))

          op_red = EA - phi_eq_red
          flux_red(ix,iy,f,e)  = (i0_red/(F_const*nct))
     & *(exp((1.0-alpha_red)*nct*F_const*op_red
     & /(R_const*(T_salt(ix,iy,f,e)+273.15)))
     & -exp((-alpha_red)*nct*F_const*op_red
     & /(R_const*(T_salt(ix,iy,f,e)+273.15))))

        else
          op_oxi = EA - phi_eq_oxi
          flux_oxi(ix,iy,f,e)  = (i0_oxi/(F_const*nct))
     & *(exp((1.0-alpha_oxi)*nct*F_const*op_oxi
     & /(R_const*(T_salt(ix,iy,f,e))))
     & -exp((-alpha_oxi)*nct*F_const*op_oxi
     & /(R_const*(T_salt(ix,iy,f,e)))))

          op_red = EA - phi_eq_red
          flux_red(ix,iy,f,e)  = (i0_red/(F_const*nct))
     & *(exp((1.0-alpha_red)*nct*F_const*op_red
     & /(R_const*(T_salt(ix,iy,f,e))))
     & -exp((-alpha_red)*nct*F_const*op_red
     & /(R_const*(T_salt(ix,iy,f,e)))))
		
        endif

          else
          flux_oxi(ix,iy,f,e)  = (i0/(F_const*nct))*(EA-phi_eq_oxi)
          flux_red(ix,iy,f,e)  = (i0/(F_const*nct))*(EA-phi_eq_red)
          endif

        enddo
        enddo

        endif		

      enddo
      enddo


cc 
      call map_alloy_salt_interface

c 
      totFluxR = 0.0
      totFluxO = 0.0
      totMassChange = 0.0
	  
      do e=1,nelv
      do f=1,2*ldim
        if (cbc(f,e,1).eq.'W  ') then
        call surface_int(sint2,sarea2,ionflux(1,1,1,1,ired),e,f)
        totFluxR = totFluxR + sint2

        call surface_int(sint2,sarea2,ionflux(1,1,1,1,ioxi),e,f)
        totFluxO = totFluxO +  sint2

        endif
      enddo
      enddo

      totFluxR = glsum(totFluxR,1)
      totFluxO = glsum(totFluxO,1)
      totCurrent = 
     & (totFluxR*charge(ired) + totFluxO*charge(ioxi))*F_const

      if (nid.eq.0) write(6,*) 'totCurrent: ',totCurrent
      if (nid.eq.0) write(6,*) 'totFluxR: ',totFluxR
      if (nid.eq.0) write(6,*) 'totFluxO: ',totFluxO

      return
      end
c----------------------------------------------------------------------------
      subroutine local_mass_loss_rate
c calculate local and global mass loss
c
      include 'SIZE'
      include 'TOTAL'
      include "moscato/ECM"
      include "moscato/ALLOY"

      integer e,f,ix,iy
      real EA,totCurrent,phi_eq_oxi,phi_eq_red
      real i0
      real totFluxR,totFluxO,sint2,sarea2
      real totMassChange,mass_flux
	  
c determine flux of each redux and oxide	
      do e=1,nelv
      do f=1,2*ldim
        if (cbc(f,e,1).eq.'W  ') then

        if (ldim.eq.3) then
        ny_corrosion = lx1
        else
        ny_corrosion = 1
        endif
  
        do ix = 1,lx1
        do iy = 1,ny_corrosion

c convert to flux to mass change rate (MCR)
c MCR does not have Na for output
        mass_flux = 0.0
        mass_flux = mass_flux - flux_oxi(ix,iy,f,e)*MM_oxi
        mass_flux = mass_flux - flux_red(ix,iy,f,e)*MM_red

        MCR(ix,iy,f,e) = mass_flux   ! mass change rate, g/m2-s

        enddo
        enddo

        endif		

      enddo
      enddo


cc 
      call map_alloy_salt_interface
cc
      totMassChangeRate = 0.0

      do e=1,nelv
      do f=1,2*ldim
        if (cbc(f,e,1).eq.'W  ') then

        call surface_int(sint2,sarea2,MCROSS,e,f)
        totMassChangeRate = totMassChangeRate + sint2

        endif
      enddo
      enddo

      totMassChangeRate = glsum(totMassChangeRate,1)

      if (nid.eq.0) write(6,*) 'totMassChangeRate: ',totMassChangeRate

      return
      end
c----------------------------------------------------------------------------
      subroutine metalDiffusion
c metalDiffusion march in time
c  for both oxi and red
c
      include 'SIZE'
      include 'TOTAL'
      include "moscato/ECM"
      include "moscato/ALLOY"
      integer idd,e,f,ix,iy,ny_corrosion
      real lap,c1,c2,c3

      real flux_point

      real max_oxi_alloy,min_oxi_alloy
      real max_red_alloy,min_red_alloy

      min_oxi_alloy = 1e6
      max_oxi_alloy = -1e6
	  
      min_red_alloy = 1e6
      max_red_alloy = -1e6

      do e=1,nelv
      do f=1,2*ldim
        if (cbc(f,e,1).eq.'W  ') then
  
        if (ldim.eq.3) then
        ny_corrosion = lx1
        else
        ny_corrosion = 1
        endif
  
        do ix = 1,lx1
        do iy = 1,ny_corrosion
  
c solve diffusion equation ...
c for oxi

        flux_point = flux_oxi(ix,iy,f,e) 
        ! assuming ionflux >0, then oxi is leaving alloy

        idd = 1
        c_oxi_alloy(idd,ix,iy,f,e) = c_oxi_alloy(idd+1,ix,iy,f,e) - 
     & flux_point*dx_alloy/D_oxi_alloy

      min_oxi_alloy = min(min_oxi_alloy,c_oxi_alloy(idd,ix,iy,f,e))
      max_oxi_alloy = max(max_oxi_alloy,c_oxi_alloy(idd,ix,iy,f,e))

 
        do idd = 2,ndd-1
         c1 = c_oxi_alloy(idd-1,ix,iy,f,e)
         c2 = c_oxi_alloy(idd,ix,iy,f,e)
         c3 = c_oxi_alloy(idd+1,ix,iy,f,e)
         call oneDimensionlap(c1,c2,c3,dx_alloy,lap)
   	     c_oxi_alloy(idd,ix,iy,f,e) = 
     & c_oxi_alloy(idd,ix,iy,f,e) + dt*D_oxi_alloy*lap ! explicit time march

      min_oxi_alloy = min(min_oxi_alloy,c_oxi_alloy(idd,ix,iy,f,e))
      max_oxi_alloy = max(max_oxi_alloy,c_oxi_alloy(idd,ix,iy,f,e))
        enddo

        idd = ndd
        c_oxi_alloy(idd,ix,iy,f,e) = c0_oxi_alloy 
        ! for alloy bulk end
      min_oxi_alloy = min(min_oxi_alloy,c_oxi_alloy(idd,ix,iy,f,e))
      max_oxi_alloy = max(max_oxi_alloy,c_oxi_alloy(idd,ix,iy,f,e))

c for red 

       flux_point = flux_red(ix,iy,f,e) 
        ! assuming ionflux <0, then red is entering alloy

        idd = 1
        c_red_alloy(idd,ix,iy,f,e) = c_red_alloy(idd+1,ix,iy,f,e) - 
     & flux_point*dx_alloy/D_red_alloy

      min_red_alloy = min(min_red_alloy,c_red_alloy(idd,ix,iy,f,e))
      max_red_alloy = max(max_red_alloy,c_red_alloy(idd,ix,iy,f,e))

 
        do idd = 2,ndd-1
         c1 = c_red_alloy(idd-1,ix,iy,f,e)
         c2 = c_red_alloy(idd,ix,iy,f,e)
         c3 = c_red_alloy(idd+1,ix,iy,f,e)
         call oneDimensionlap(c1,c2,c3,dx_alloy,lap)
   	     c_red_alloy(idd,ix,iy,f,e) = 
     & c_red_alloy(idd,ix,iy,f,e) + dt*D_red_alloy*lap ! explicit time march

      min_red_alloy = min(min_red_alloy,c_red_alloy(idd,ix,iy,f,e))
      max_red_alloy = max(max_red_alloy,c_red_alloy(idd,ix,iy,f,e))
        enddo

        idd = ndd
        c_red_alloy(idd,ix,iy,f,e) = c0_red_alloy 
        ! for alloy bulk end
      min_red_alloy = min(min_red_alloy,c_red_alloy(idd,ix,iy,f,e))
      max_red_alloy = max(max_red_alloy,c_red_alloy(idd,ix,iy,f,e))

        enddo
        enddo
   
        endif
      enddo
      enddo
	  
      min_oxi_alloy = glmin(min_oxi_alloy,1)
      max_oxi_alloy = glmax(max_oxi_alloy,1)

      if (nid.eq.0) then
      write(6,*) 'min_oxi_alloy: ', min_oxi_alloy
      write(6,*) 'max_oxi_alloy: ', max_oxi_alloy
      endif
	  
      min_red_alloy = glmin(min_red_alloy,1)
      max_red_alloy = glmax(max_red_alloy,1)

      if (nid.eq.0) then
      write(6,*) 'min_red_alloy: ', min_red_alloy
      write(6,*) 'max_red_alloy: ', max_red_alloy
      endif
	  

      ! dump_alloy_profile 
      ! 	  
      if ((mod(istep,iostep).eq.0).and.(istep.gt.0)) then
      i_a_dump = i_a_dump + 1
        if (nid.eq.0) then
        write(6,*) "dumping: alloy_file_",i_a_dump
        endif
      call dump_alloy_profile(i_a_dump)
      endif

      return
      end
c----------------------------------------------------------------------------
      subroutine dump_alloy_profile(crfi)
c dump alloy profile
c for both oxi and red
c for plotting and restart
c 
c modify based on gen_rea_bc. 
c but gen_rea_bc is based on element face 
c here we need gll point based
c
c c_oxi_alloy(ndd,lx1,ly1,nfe,lelt)
c c_red_alloy(ndd,lx1,ly1,nfe,lelt)

      include 'SIZE'
      include 'TOTAL'
      include "moscato/ECM"
      include "moscato/ALLOY"

      integer e,eb,eg,i,ii,ix,ix0,iy,mid,kb
      parameter (lblock=500)
c	(ndd,lx1,ly1,nfe,lelt)
      common /scrns_alloy1/ oxia(ndd,lx1,ly1,nfe,lblock),
     & reda(ndd,lx1,ly1,nfe,lblock),
     & pxyza(3,lx1,ly1,nfe,lblock),
     & mca(lx1,ly1,nfe,lblock),
     & wka(ndd*lx1*ly1*nfe*lblock),
     & wka2(3*lx1*ly1*nfe*lblock),
     & wka3(lx1*ly1*nfe*lblock),
     & ibc(6,lblock),wk(5*6*lblock)

      character*1 s4(4)
      character*3 s3
      integer     i4
      equivalence(i4,s4)
      equivalence(s3,s4)
	  
      character*80 alloy_file_name
      character*80 dump_frmt

      integer crfi,ny_corrosion
      integer nface,nlg,nemax

      real xp,yp,zp

      nface = 2*ldim
      nlg = nelg(1)

      call blank(cr_alloy_file_name,80)      
      call blank(dump_frmt,80)
  
c  110 format('(a,i0)')

      if (nid.eq.0) then
	   write(alloy_file_name,'(a,i0)') "alloy_file_",crfi
	  
       write(6,*) "dumping ",alloy_file_name

       open(unit=510,file=alloy_file_name,status='unknown')
	   write(510,*) 'time ', time               ! write time tag
       write(510,*) 'eg f ix iy x y z mc oxi .. red ..'  ! header of file
       write(dump_frmt,'(a,i0,a,i0,a)') 
     & "(i9,3i3,4g16.8,",ndd,"g16.8,",ndd,"g16.8)"
      endif
	  


c gather information per block into nid=0 processor
c and only write at nid=0 processor
c 
c 
      do eb=1,nlg,lblock
         nemax = min(eb+lblock-1,nlg)

         call izero(ibc, 6*lblock)
         call rzero(oxia,ndd*lx1*ly1*nfe*lblock)
         call rzero(reda,ndd*lx1*ly1*nfe*lblock)
         call rzero(pxyza,3*lx1*ly1*nfe*lblock)
         call rzero(mca,lx1*ly1*nfe*lblock)
		 
         kb = 0
         do eg=eb,nemax    ! loop global element index for this block
            mid = gllnid(eg)   ! eg -> processors id
            e   = gllel (eg)   !  eg -> e, logal element index
            kb  = kb+1         ! local element account in this eb 
            if (mid.eq.nid) then ! if eg belong to this nid, then store info here.
               do i=1,nface
                  i4 = 0
                  call chcopy(s4,cbc(i,e,1),3)
                  ibc(i,kb) = i4
				  
	              ! store c_Cr_alloy to cra 
                  ny_corrosion = 1
                  if (ldim.eq.3) ny_corrosion = lx1

                  do ix = 1,lx1
	              do iy = 1,ny_corrosion
                  call copy(oxia(1,ix,iy,i,kb),
     & c_oxi_alloy(1,ix,iy,i,e),ndd)
	 
	              call copy(reda(1,ix,iy,i,kb),
     & c_red_alloy(1,ix,iy,i,e),ndd)

                  call copy(pxyza(1,ix,iy,i,kb),
     & pxyz_alloy(1,ix,iy,i,e),3)
	 
                  call copy(mca(ix,iy,i,kb),
     & MCR(ix,iy,i,e),1)
	 
                  enddo
                  enddo
				  
               enddo
            endif
         enddo

         call igop(ibc,wk,'+  ', 6*lblock)                ! Sum across all processors
         call gop(oxia,wka,'+  ',ndd*lx1*ly1*nfe*lblock)  ! Sum across all processors
         call gop(reda,wka,'+  ',ndd*lx1*ly1*nfe*lblock)  ! Sum across all processors
         call gop(pxyza,wka2,'+  ',3*lx1*ly1*nfe*lblock)  ! Sum across all processors
         call gop(mca,wka3,'+  ',lx1*ly1*nfe*lblock)      ! Sum across all processors

         ! write at nid=0
         if (nid.eq.0) then
            kb = 0
            do eg=eb,nemax
               kb  = kb+1

               do i=1,nface
                  i4 = ibc(i,kb)   ! equivalenced to s4 and s3
                  !write(6,*) s3
                  if (s3.eq.'W  ') then  ! only dump for wall
                    !write(6,*) 'writing face data'
                    !write(510,*) 'writing face data',lx1,ly1

                    ny_corrosion = 1
                    if (ldim.eq.3) ny_corrosion = lx1

                    do ix0 =1,lx1*ny_corrosion
                    ix = mod(ix0,lx1)
                    if (ix.eq.0) ix = lx1
                    iy = (ix0-ix)/lx1 + 1
					
                    !xp = pxyza(1,ix,iy,i,kb)
                    !yp = pxyza(2,ix,iy,i,kb)
                    !zp = pxyza(3,ix,iy,i,kb)
		
                    write(510,dump_frmt)
     &  eg,i,ix,iy,(pxyza(ii,ix,iy,i,kb),ii=1,3),mca(ix,iy,i,kb),
     & (oxia(ii,ix,iy,i,kb),ii=1,ndd),(reda(ii,ix,iy,i,kb),ii=1,ndd)
                    enddo
	
                  endif    


               enddo
            enddo
         endif
      enddo
	  
      if (nid.eq.0) close(510)
	  
      return
      end
c----------------------------------------------------------------------------
      subroutine read_alloy_profile(crfi)
c read cr alloy profile for restart
c this reader is done in parallel
c
      include 'SIZE'
      include 'TOTAL'
      include "moscato/ECM"
      include "moscato/ALLOY"
	  
      integer e,eb,eg,i,ii,ix,ix0,iy,mid,kb
      common /scrns_alloy2/ oxia2(ndd),reda2(ndd),mca2
	  
      character*80 alloy_file_name
      character*80 dump_frmt
      character*32 dummyline
	  
      integer crfi,ny_corrosion
      integer iostatus
	  
      real xp,yp,zp

      nface = 2*ldim
      nlg = nelg(1)

      if (nid.eq.0) then 
       write(6,*) 'CALL: read_alloy_profile for restart'
      endif
	 
c read is parallel
c so open file for all nid

      ierr  = 0
	  write(alloy_file_name,'(a,i0)') "alloy_file_",crfi
	  
      if (nid.eq.0) then
        open(unit=511,file=alloy_file_name,status='old',iostat=ierr)      
	    if (ierr .gt. 0) call exitti('Cannot open cr_alloy_file!$',1) 
      endif
	
      if (nid.ne.0) open(unit=511,file=alloy_file_name,status='old')      
      
	  ! skipper header line 
      read(511,*) dummyline
      read(511,*) dummyline

      if (nid.eq.0) write(6,*) dummyline

      write(dump_frmt,'(a,i0,a,i0,a)') 
     & "(i9,3i3,4g16.8,",ndd,"g16.8,",ndd,"g16.8)"

      iostatus = 0
      ! read line by line, until file ends

      do

	  read(511,dump_frmt,iostat=iostatus) 
     & ieg,ifc,ix,iy,xp,yp,zp,mca2,
     & (oxia2(ii),ii=1,ndd),(reda2(ii),ii=1,ndd)
   
      if (iostatus.ne.0) exit ! file ends  , exit loop

      if (gllnid(ieg).eq.nid) then ! if this elemenet belong to this nid 
         iel=gllel(ieg)
         call copy(c_oxi_alloy(1,ix,iy,ifc,iel),oxia2, ndd)
         call copy(c_red_alloy(1,ix,iy,ifc,iel),reda2, ndd)
         call copy(MCR(ix,iy,ifc,iel),mca2,1)
      endif

      enddo 

      close(511)

	  
      return
      end
c----------------------------------------------------------------------------
      subroutine oneDimensionlap(c1,c2,c3,dx,lap)
      real c1,c2,c3,dx,lap	  

	  lap = (c1 -2.0*c2 + c3)/(dx*dx)

      return
      end
c----------------------------------------------------------------------------
      subroutine dump_salt_surface_corrosion_data
      include 'SIZE'
      include 'TOTAL'
      include "moscato/ECM"
      include "moscato/ALLOY"

      iastep = param(68)
      if  (iastep.eq.0) iastep=param(15)   ! same as iostep

      if ( (mod(istep,iastep).eq.0.and.istep.gt.1) .or.lastep.eq.1) then
        ! for 2d case, vz position is not dumped.

      call outpost(ionflux(1,1,1,1,ired),ionflux(1,1,1,1,ioxi),
     & COAOSS,CRAOSS,MCROSS,'cro')

      endif

      return
      end
c----------------------------------------------------------------------