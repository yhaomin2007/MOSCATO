c for salazar paper, experimental
c
c------------------------------------------------------------------
      subroutine butler_volmer()
c special subroutine for slazzar paper.
c use butler-volmer to determin flux for ferri and ferro in anode and cathode
c
c      ioxi
c      ired 
c     anode_faceID
c     cathode_faceID

      include 'SIZE'
      include 'TOTAL'
      include "moscato/ECM"
      logical ifverbose
  
      real E1,E2,E3,deltaE,dcde,C1,C2
      integer iter, niter

      data deltaE /0.001/
  
      ntot = lx1*ly1*lz1*nelv

c use newton iteration to find E_anode that will result total current=0

      niter = 20
      do iter = 1,niter

      E1 = E_anode
      E2 = E_anode + deltaE
      call total_current_salazar(E1,E_ca,C1)
  
      if (abs(C1).lt.1e-8) then
        E_anode = E1
        if(nid.eq.0) then
         write(6,*) 'curent converged'
         write(6,*) 'E_anode:',E_anode
         write(6,*) 'current:',C1
        endif
      exit
      endif

      call total_current_salazar(E2,E_ca,C2)

      dcde = (C2-C1)/deltaE

      E3 = E1 - C1/dcde    ! find next input

      deltaE = abs(E1-E3)/100.0
      E_anode = E3
  
      ifverbose = .TRUE. !.FALSE.
      if (ifverbose) then
      if(nid.eq.0) then
         write(6,*) 'E1:',E1,' E2:',E2,' E3:', E3
         write(6,*) 'C1:',C1,' C2:',C2
      endif
      endif

      enddo

      E_cathode = E_anode + E_ca

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
	  
      flux_max = glmax(currentDensity,ntot)
      flux_min = glmin(currentDensity,ntot)
      if(nid.eq.0) then
       write(6,*) 'currentDensity max/min : ', flux_max,'/',flux_min
      endif
	  
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine total_current_salazar(Ea,Eca,totCurrent)
c calcualte total current based on electrode potential
c  
c
      include 'SIZE'
      include 'TOTAL'
      include "moscato/ECM"
      include "moscato/ALLOY"

      real lcro(lx1,ly1,lz1,lelt)

      real Ea,Eca,totCurrent
      real E_eq,OP,E_eq_a,E_eq_c
      real cr,co,tl
      real sint2,sarea2
      real totCurretAnode,totCurretCathode
      real anodeArea,cathodeArea
      real anodePhie,cathodePhie
      real anodeLCRO,cathodeLCRO
      real coa,cra,coc,crc
      

      integer i,j,k,e,f,idf
      integer i0,i1,j0,j1,k0,k1

! determine flux of each redux and oxide

      totCurrent = 0.0
      totCurretAnode = 0.0
      totCurretCathode = 0.0

      anodeArea = 0.0
      cathodeArea = 0.0
 
      anodePhie = 0.0
      cathodePhie = 0.0

!------------------------------------------------------
! get avg co,cr on anode and cathode
! coa,cra,coc,crc
       coa = 0.0
       cra = 0.0
       coc = 0.0 
       crc = 0.0

! obgain log(Cr/Co) 
       ntot = lx1*ly1*lz1*nelv

      do i = 1,ntot
       lcro(i,1,1,1) = log(t(i,1,1,1,ired+1)/t(i,1,1,1,ioxi+1)) 
      enddo


      do e=1,nelv
      do f=1,2*ldim
        idf = bc(5,f,e,1)
        if (idf.eq.anode_faceID) then ! for anode

        call surface_int(sint2,sarea2,t(1,1,1,1,ioxi+1),e,f)
        coa = coa + sint2
        anodeArea = anodeArea + sarea2

        call surface_int(sint2,sarea2,t(1,1,1,1,ired+1),e,f)
        cra = cra + sint2

        call surface_int(sint2,sarea2,phie,e,f)
        anodePhie = anodePhie + sint2

        call surface_int(sint2,sarea2,lcro,e,f)
        anodeLCRO = anodeLCRO + sint2

        endif

        if (idf.eq.cathode_faceID) then ! for cathode

        call surface_int(sint2,sarea2,t(1,1,1,1,ioxi+1),e,f)
        coc = coc + sint2
        cathodeArea = cathodeArea + sarea2

        call surface_int(sint2,sarea2,t(1,1,1,1,ired+1),e,f)
        crc = crc + sint2

        call surface_int(sint2,sarea2,phie,e,f)
        cathodePhie = cathodePhie + sint2

        call surface_int(sint2,sarea2,lcro,e,f)
        cathodeLCRO = cathodeLCRO + sint2

        endif

      enddo
      enddo
      coa = glsum(coa,1)
      cra = glsum(cra,1)
      coc = glsum(coc,1)
      crc = glsum(crc,1)

      anodeArea = glsum(anodeArea,1)
      cathodeArea = glsum(cathodeArea,1)

      coa = coa/anodeArea
      cra = cra/anodeArea
      coc = coc/cathodeArea
      crc = crc/cathodeArea

      if (nid.eq.0) write(6,*) 'coa: ',coa
      if (nid.eq.0) write(6,*) 'cra: ',cra
      if (nid.eq.0) write(6,*) 'coc: ',coc
      if (nid.eq.0) write(6,*) 'crc: ',crc

      anodePhie = glsum(anodePhie,1)
      cathodePhie = glsum(cathodePhie,1)

      anodePhie = anodePhie/anodeArea
      cathodePhie = cathodePhie/cathodeArea 

      anodeLCRO = glsum(anodeLCRO,1)
      cathodeLCRO = glsum(cathodeLCRO,1)

      anodeLCRO = anodeLCRO/anodeArea
      cathodeLCRO = cathodeLCRO/cathodeArea 

      E_eq_a = (T_anode-sT0)*sSrx/(nct*F_const)
!     & - (R_const*T_anode/(nct*F_const))*log(cra/coa) 
 
      E_eq_c = (T_cathode-sT0)*sSrx/(nct*F_const)
!     & - (R_const*T_cathode/(nct*F_const))*log(crc/coc)

      E_eq_a1 = (T_anode-sT0)*sSrx/(nct*F_const)
     & - (R_const*T_anode/(nct*F_const))*log(cra/coa)
 
      E_eq_c1 = (T_cathode-sT0)*sSrx/(nct*F_const)
     & - (R_const*T_cathode/(nct*F_const))*log(crc/coc)
 
      E_eq_a2 = (R_const*T_anode/(nct*F_const))*log(coa/cra) 
      E_eq_c2 = (R_const*T_cathode/(nct*F_const))*log(coc/crc)


!      E_eq_a = (R_const*T_anode/(nct*F_const))*log(coa/cra) 	 
!      E_eq_c = (R_const*T_cathode/(nct*F_const))*log(coc/crc)

      OP_a = Ea - E_eq_a
      OP_c = Ea + Eca - E_eq_c

      OP_a1 = Ea - E_eq_a1
      OP_c1 = Ea + Eca - E_eq_c1

      OP_a2 = Ea - E_eq_a2
      OP_c2 = Ea + Eca - E_eq_c2

!------------------------------------------------------
      anodeArea = 0.0
      cathodeArea = 0.0

      do e=1,nelv
      do f=1,2*ldim
        idf = bc(5,f,e,1)
        if (idf.eq.anode_faceID) then ! for anode
          call facind(i0,i1,j0,j1,k0,k1,nx1,ny1,nz1,f)
          do k=k0,k1
          do j=j0,j1
          do i=i0,i1

          cr = t(i,j,k,e,ired+1) ! local red
          co = t(i,j,k,e,ioxi+1) ! local oxi
          tl = t(i,j,k,e,1)  ! local temperature
          !may change definition accordinng to salazar paper definition
!          E_eq = (tl-sT0)*sSrx/(nct*F_const)
!     & +(R_const*tl/(nct*F_const))*log(co/cr) 

!           E_eq = (R_const*tl/(nct*F_const))*log(co/cr) 
!           OP = Ea - E_eq

          OP = OP_a1
          !OP = OP + (R_const*tl/(nct*F_const))*log(cr/co) 
          ! the direction of current is consistent with oxi ion flux ??

          ionflux(i,j,k,e,ioxi) =
     & -sk0*exp((sExa/R_const)*(1.0/sT0-1.0/tl))*
     & (co*exp((-nct*F_const*(1.0-stheta)*OP)/(R_const*tl)) 
     & -cr*exp((nct*F_const*stheta*OP)/(R_const*tl)))

!          ionflux(i,j,k,e,ioxi) =
!     & -sk0*exp((sExa/R_const)*(1.0/sT0-1.0/tl))*
!     & (co*exp((-nct*F_const*(1.0-stheta)*OP)/(R_const*tl)
!     & - (1.0-stheta)*log(cr/co)) 
!     & -cr*exp((nct*F_const*stheta*OP)/(R_const*tl) 
!     & + stheta*log(cr/co)))


!           ionflux(i,j,k,e,ioxi) =
!     & -sk0*exp((sExa/R_const)*(1.0/sT0-1.0/tl))*
!     & (coa*exp((-nct*F_const*(1.0-stheta)*OP)/(R_const*tl)) 
!     & -cra*exp((nct*F_const*stheta*OP)/(R_const*tl)))

          ionflux(i,j,k,e,ired) = - ionflux(i,j,k,e,ioxi)

          currentDensity(i,j,k,e) = nct*F_const*ionflux(i,j,k,e,ioxi)
  
!          if (ionflux(i,j,k,e,ioxi).lt.0.0) then
!          write(6,*) 
!     & 'ERROR: oxidized ionflux on anode should be positive'
!          endif

          enddo
          enddo
          enddo


        call surface_int(sint2,sarea2,ionflux(1,1,1,1,ioxi),e,f)
        totCurretAnode = totCurretAnode + sint2*nct*F_const
        anodeArea = anodeArea + sarea2

        endif

        if (idf.eq.cathode_faceID) then ! for cathode

          call facind(i0,i1,j0,j1,k0,k1,nx1,ny1,nz1,f)
          do k=k0,k1
          do j=j0,j1
          do i=i0,i1

          cr = t(i,j,k,e,ired+1) ! local redox 
          co = t(i,j,k,e,ioxi+1) ! local oxi
          tl = t(i,j,k,e,1)  ! local temperature
          !may change definition accordinng to salazar paper definition
!          E_eq = (tl-sT0)*sSrx/(nct*F_const)
!     & + (R_const*tl/(nct*F_const))*log(co/cr)
!          OP = Ea + Eca - E_eq
!          E_eq = (R_const*tl/(nct*F_const))*log(co/cr)
!          OP = Ea + Eca - E_eq

         OP = OP_c1
         !OP = OP + (R_const*tl/(nct*F_const))*log(cr/co) 

          ionflux(i,j,k,e,ioxi) =
     & -sk0*exp((sExa/R_const)*(1.0/sT0-1.0/tl))*
     & (co*exp((-nct*F_const*(1.0-stheta)*OP)/(R_const*tl)) 
     & -cr*exp((nct*F_const*stheta*OP)/(R_const*tl)))
 
!         ionflux(i,j,k,e,ioxi) =
!    & -sk0*exp((sExa/R_const)*(1.0/sT0-1.0/tl))*
!    & (co*exp((-nct*F_const*(1.0-stheta)*OP)/(R_const*tl)
!    & - (1.0-stheta)*log(cr/co)) 
!    & -cr*exp((nct*F_const*stheta*OP)/(R_const*tl)
!    & + stheta*log(cr/co)))

!          ionflux(i,j,k,e,ioxi) =
!     & -sk0*exp((sExa/R_const)*(1.0/sT0-1.0/tl))*
!     & (coc*exp((-nct*F_const*(1.0-stheta)*OP)/(R_const*tl)) 
!     & -crc*exp((nct*F_const*stheta*OP)/(R_const*tl)))

          ionflux(i,j,k,e,ired) = - ionflux(i,j,k,e,ioxi) 
          currentDensity(i,j,k,e) = nct*F_const*ionflux(i,j,k,e,ioxi)

!          if (ionflux(i,j,k,e,ioxi).gt.0.0) then
!          write(6,*) 
!     & 'ERROR: oxidized ionflux on cathode should be negative'
!          endif


          enddo
          enddo
          enddo

        call surface_int(sint2,sarea2,ionflux(1,1,1,1,ioxi),e,f)
        totCurretCathode = totCurretCathode + sint2*nct*F_const
        cathodeArea = cathodeArea + sarea2

        endif

      enddo
      enddo

      totCurretAnode = glsum(totCurretAnode,1)
      totCurretCathode = glsum(totCurretCathode,1)
      totCurrent = (totCurretAnode + totCurretCathode)

      if (nid.eq.0) write(6,*) 'totCurrent: ',totCurrent
      if (nid.eq.0) write(6,*) 'totCurretAnode: ',totCurretAnode
      if (nid.eq.0) write(6,*) 'totCurretCathode: ',totCurretCathode

      anodeArea = glsum(anodeArea,1)
      cathodeArea = glsum(cathodeArea,1)
      if (nid.eq.0) write(6,*) 'anodeArea: ',anodeArea
      if (nid.eq.0) write(6,*) 'cathodeArea: ',cathodeArea

      if (nid.eq.0) write(6,*) 'avg anode current density: '
     & ,totCurretAnode/anodeArea
      if (nid.eq.0) write(6,*) 'avg cathode current density: '
     & ,totCurretCathode/cathodeArea

      if (nid.eq.0) write(6,*) 'avg anode phie: '
     & ,anodePhie
      if (nid.eq.0) write(6,*) 'avg cathode phie: '
     & ,cathodePhie

      if (nid.eq.0) write(6,*) 'solution ohmic drop: '
     & ,anodePhie - cathodePhie

    
c      V_cell = V_oc + OP_a2 - OP_c2
c      if (nid.eq.0) write(6,*) 'V_cell1: ',V_cell
c
c      V_cell = V_oc - OP_a2 - OP_c2
c      if (nid.eq.0) write(6,*) 'V_cell2: ',V_cell
c	  
c      V_cell = V_oc + OP_a2 + OP_c2 
c      if (nid.eq.0) write(6,*) 'V_cell3: ',V_cell
c	  
c      V_cell = V_oc - OP_a2 + OP_c2 
c      if (nid.eq.0) write(6,*) 'V_cell4: ',V_cell
c	  
      V_cell = V_oc - OP_a + OP_c
      if (nid.eq.0) write(6,*) 'V_cell5: ',V_cell

      V_cell = V_oc - OP_a1 + OP_c1
      if (nid.eq.0) write(6,*) 'V_cell6: ',V_cell

      V_cell = V_oc - OP_a + OP_c - (anodePhie - cathodePhie)
      if (nid.eq.0) write(6,*) 'V_cell7: ',V_cell

      V_cell = V_oc - OP_a1 + OP_c1 - (anodePhie - cathodePhie)
      if (nid.eq.0) write(6,*) 'V_cell8: ',V_cell

      if (nid.eq.0) write(6,*) 'E_eq_a: ',E_eq_a
      if (nid.eq.0) write(6,*) 'E_eq_c: ',E_eq_c
      if (nid.eq.0) write(6,*) 'E_eq_a1: ',E_eq_a1
      if (nid.eq.0) write(6,*) 'E_eq_c1: ',E_eq_c1
      if (nid.eq.0) write(6,*) 'E_eq_a2: ',E_eq_a2
      if (nid.eq.0) write(6,*) 'E_eq_c2: ',E_eq_c2
 
      if (nid.eq.0) write(6,*) 'OP_a: ',OP_a
      if (nid.eq.0) write(6,*) 'OP_c: ',OP_c
      if (nid.eq.0) write(6,*) 'OP_a1: ',OP_a1
      if (nid.eq.0) write(6,*) 'OP_c1: ',OP_c1
      if (nid.eq.0) write(6,*) 'OP_a2: ',OP_a2
      if (nid.eq.0) write(6,*) 'OP_c2: ',OP_c2
      if (nid.eq.0) write(6,*) 'V_oc: ',V_oc

      return
      end
c----------------------------------------------------------------------------
