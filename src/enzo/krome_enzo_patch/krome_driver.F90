!KROME_DRIVER

subroutine krome_driver(d, e, ge, u, v, w, &
      De, HM, HI, HeI, &
      H2I, OI, OHI, H2OI, &
      O2I, CHI, COI, CI, &
      CO2I, CH2I, CH3I, CH4I, &
      HII, HeII, H2II, OII, &
      OHII, H2OII, H3OII, O2II, &
      CII, COII, HeIII,  in, jn, kn, imethod, &
      idual, idim, &
      is, js, ks, ie, je, ke, &
      dt, aye, &
      utem, uxyz, uaye, urho, utim, &
      gamma, fh, dtoh)

  !     SOLVE MULTI-SPECIES RATE EQUATIONS AND RADIATIVE COOLING
  !
  !     2014, KROME DEVELOPERS to interface the package with ENZO
  !
  !     PURPOSE:
  !     Solve the multi-species rate and cool equations via KROME.
  !
  !     INPUTS:
  !     in,jn,kn - dimensions of 3D fields
  !
  !     d        - total density field
  !
  !     is,ie    - start and end indices of active region (zero based)
  !     idual    - dual energy formalism flag (0 = off, 1 = on)
  !     idim     - dimensionality (rank) of problem
  !     imethod  - Hydro method (0 = PPMDE, 2 = ZEUS-type)
  !
  !     fh       - Hydrogen mass fraction (typically 0.76)
  !     dtoh     - Deuterium to H mass ratio
  !     dt       - timestep to integrate over
  !     aye      - expansion factor (in code units)
  !
  !     utim     - time units (i.e. code units to CGS conversion factor)
  !     uaye     - expansion factor conversion factor (uaye = 1/(1+zinit))
  !     urho     - density units
  !     uxyz     - length units
  !     utem     - temperature(-like) units
  !
  !     OUTPUTS:
  !     update chemical abundances densities (HI, HII, etc) and energy
  !
  !     PARAMETERS:
  !     mh      - H mass in cgs units
  !
  !-----------------------------------------------------------------------
  !     USE KROME
  use krome_main
  use krome_user
  use krome_constants

  implicit none

  real*8,parameter::mh=p_mass !mass_h
  real*8::dom,factor,tgas,tgasold,krome_x(krome_nmols)
  real*8::dt,dt_hydro,idom,edot,krome_tiny
  real*8::d(in,jn,kn),e(in,jn,kn),ge(in,jn,kn)
  real*8::u(in,jn,kn),v(in,jn,kn),w(in,jn,kn)
  real*8::aye,utem,uxyz,uaye,urho,utim,gamma,fh,dtoh
  integer::in,jn,kn,imethod,idual,is,js,ks,ie,je,ke,idim
  integer::i,j,k

  real*8::De(in,jn,kn)
  real*8::HM(in,jn,kn)
  real*8::HI(in,jn,kn)
  real*8::HeI(in,jn,kn)
  real*8::H2I(in,jn,kn)
  real*8::OI(in,jn,kn)
  real*8::OHI(in,jn,kn)
  real*8::H2OI(in,jn,kn)
  real*8::O2I(in,jn,kn)
  real*8::CHI(in,jn,kn)
  real*8::COI(in,jn,kn)
  real*8::CI(in,jn,kn)
  real*8::CO2I(in,jn,kn)
  real*8::CH2I(in,jn,kn)
  real*8::CH3I(in,jn,kn)
  real*8::CH4I(in,jn,kn)
  real*8::HII(in,jn,kn)
  real*8::HeII(in,jn,kn)
  real*8::H2II(in,jn,kn)
  real*8::OII(in,jn,kn)
  real*8::OHII(in,jn,kn)
  real*8::H2OII(in,jn,kn)
  real*8::H3OII(in,jn,kn)
  real*8::O2II(in,jn,kn)
  real*8::CII(in,jn,kn)
  real*8::COII(in,jn,kn)
  real*8::HeIII(in,jn,kn)

  !******************************

  !set units
  dom = urho*(aye**3)/mh

  !scaling factor for comoving->proper
  factor = aye**(-3)

  !check minimal value and comoving->proper
  krome_tiny = 1d-30
  do k = ks+1, ke+1
    do j = js+1, je+1
      do i = is+1, ie+1
        d(i,j,k) = d(i,j,k) * factor
        !scale comoving->proper
        De(i,j,k) = De(i,j,k) * factor
        HM(i,j,k) = HM(i,j,k) * factor
        HI(i,j,k) = HI(i,j,k) * factor
        HeI(i,j,k) = HeI(i,j,k) * factor
        H2I(i,j,k) = H2I(i,j,k) * factor
        OI(i,j,k) = OI(i,j,k) * factor
        OHI(i,j,k) = OHI(i,j,k) * factor
        H2OI(i,j,k) = H2OI(i,j,k) * factor
        O2I(i,j,k) = O2I(i,j,k) * factor
        CHI(i,j,k) = CHI(i,j,k) * factor
        COI(i,j,k) = COI(i,j,k) * factor
        CI(i,j,k) = CI(i,j,k) * factor
        CO2I(i,j,k) = CO2I(i,j,k) * factor
        CH2I(i,j,k) = CH2I(i,j,k) * factor
        CH3I(i,j,k) = CH3I(i,j,k) * factor
        CH4I(i,j,k) = CH4I(i,j,k) * factor
        HII(i,j,k) = HII(i,j,k) * factor
        HeII(i,j,k) = HeII(i,j,k) * factor
        H2II(i,j,k) = H2II(i,j,k) * factor
        OII(i,j,k) = OII(i,j,k) * factor
        OHII(i,j,k) = OHII(i,j,k) * factor
        H2OII(i,j,k) = H2OII(i,j,k) * factor
        H3OII(i,j,k) = H3OII(i,j,k) * factor
        O2II(i,j,k) = O2II(i,j,k) * factor
        CII(i,j,k) = CII(i,j,k) * factor
        COII(i,j,k) = COII(i,j,k) * factor
        HeIII(i,j,k) = HeIII(i,j,k) * factor

        !mimimal value check
        De(i,j,k) = max(De(i,j,k), krome_tiny)
        HM(i,j,k) = max(HM(i,j,k), krome_tiny)
        HI(i,j,k) = max(HI(i,j,k), krome_tiny)
        HeI(i,j,k) = max(HeI(i,j,k), krome_tiny)
        H2I(i,j,k) = max(H2I(i,j,k), krome_tiny)
        OI(i,j,k) = max(OI(i,j,k), krome_tiny)
        OHI(i,j,k) = max(OHI(i,j,k), krome_tiny)
        H2OI(i,j,k) = max(H2OI(i,j,k), krome_tiny)
        O2I(i,j,k) = max(O2I(i,j,k), krome_tiny)
        CHI(i,j,k) = max(CHI(i,j,k), krome_tiny)
        COI(i,j,k) = max(COI(i,j,k), krome_tiny)
        CI(i,j,k) = max(CI(i,j,k), krome_tiny)
        CO2I(i,j,k) = max(CO2I(i,j,k), krome_tiny)
        CH2I(i,j,k) = max(CH2I(i,j,k), krome_tiny)
        CH3I(i,j,k) = max(CH3I(i,j,k), krome_tiny)
        CH4I(i,j,k) = max(CH4I(i,j,k), krome_tiny)
        HII(i,j,k) = max(HII(i,j,k), krome_tiny)
        HeII(i,j,k) = max(HeII(i,j,k), krome_tiny)
        H2II(i,j,k) = max(H2II(i,j,k), krome_tiny)
        OII(i,j,k) = max(OII(i,j,k), krome_tiny)
        OHII(i,j,k) = max(OHII(i,j,k), krome_tiny)
        H2OII(i,j,k) = max(H2OII(i,j,k), krome_tiny)
        H3OII(i,j,k) = max(H3OII(i,j,k), krome_tiny)
        O2II(i,j,k) = max(O2II(i,j,k), krome_tiny)
        CII(i,j,k) = max(CII(i,j,k), krome_tiny)
        COII(i,j,k) = max(COII(i,j,k), krome_tiny)
        HeIII(i,j,k) = max(HeIII(i,j,k), krome_tiny)

      end do
    end do
  end do

  !loop over zones
  do k = ks+1, ke+1
    do j = js+1, je+1
      do i = is+1, ie+1

        !rhogas = #KROME_sum to be removed

        !convert to number densities
        krome_x(krome_idx_E) = De(i,j,k) * dom
        krome_x(krome_idx_Hk) = HM(i,j,k) * dom
        krome_x(krome_idx_H) = HI(i,j,k) * dom
        krome_x(krome_idx_HE) = HeI(i,j,k) * dom * 0.25d0
        krome_x(krome_idx_H2) = H2I(i,j,k) * dom * 0.5d0
        krome_x(krome_idx_O) = OI(i,j,k) * dom * 0.0625d0
        krome_x(krome_idx_OH) = OHI(i,j,k) * dom * 0.0588235294118d0
        krome_x(krome_idx_H2O) = H2OI(i,j,k) * dom * 0.0555555555556d0
        krome_x(krome_idx_O2) = O2I(i,j,k) * dom * 0.03125d0
        krome_x(krome_idx_CH) = CHI(i,j,k) * dom * 0.0769230769231d0
        krome_x(krome_idx_CO) = COI(i,j,k) * dom * 0.0357142857143d0
        krome_x(krome_idx_C) = CI(i,j,k) * dom * 0.0833333333333d0
        krome_x(krome_idx_CO2) = CO2I(i,j,k) * dom * 0.0227272727273d0
        krome_x(krome_idx_CH2) = CH2I(i,j,k) * dom * 0.0714285714286d0
        krome_x(krome_idx_CH3) = CH3I(i,j,k) * dom * 0.0666666666667d0
        krome_x(krome_idx_CH4) = CH4I(i,j,k) * dom * 0.0625d0
        krome_x(krome_idx_Hj) = HII(i,j,k) * dom
        krome_x(krome_idx_HEj) = HeII(i,j,k) * dom * 0.25d0
        krome_x(krome_idx_H2j) = H2II(i,j,k) * dom * 0.5d0
        krome_x(krome_idx_Oj) = OII(i,j,k) * dom * 0.0625d0
        krome_x(krome_idx_OHj) = OHII(i,j,k) * dom * 0.0588235294118d0
        krome_x(krome_idx_H2Oj) = H2OII(i,j,k) * dom * 0.0555555555556d0
        krome_x(krome_idx_H3Oj) = H3OII(i,j,k) * dom * 0.0526315789474d0
        krome_x(krome_idx_O2j) = O2II(i,j,k) * dom * 0.03125d0
        krome_x(krome_idx_Cj) = CII(i,j,k) * dom * 0.0833333333333d0
        krome_x(krome_idx_COj) = COII(i,j,k) * dom * 0.0357142857143d0
        krome_x(krome_idx_HEjj) = HeIII(i,j,k) * dom * 0.25d0

        call evaluate_tgas(d(i,j,k), e(i,j,k), ge(i,j,k),&
            u(i,j,k), v(i,j,k), w(i,j,k),&
            krome_x(:),imethod,idual,idim,tgas,&
            utem)

        !store old tgas
        tgasold = tgas

        dt_hydro = utim*dt !dt*time_conversion

        !call KROME solver
        call krome(krome_x(:),tgas,dt_hydro)

        idom = 1.d0/dom
        !convert back to code units
        De(i,j,k) = krome_x(krome_idx_E) * idom
        HM(i,j,k) = krome_x(krome_idx_Hk) * idom
        HI(i,j,k) = krome_x(krome_idx_H) * idom
        HeI(i,j,k) = krome_x(krome_idx_HE) * idom * 4d0
        H2I(i,j,k) = krome_x(krome_idx_H2) * idom * 2d0
        OI(i,j,k) = krome_x(krome_idx_O) * idom * 16d0
        OHI(i,j,k) = krome_x(krome_idx_OH) * idom * 17d0
        H2OI(i,j,k) = krome_x(krome_idx_H2O) * idom * 18d0
        O2I(i,j,k) = krome_x(krome_idx_O2) * idom * 32d0
        CHI(i,j,k) = krome_x(krome_idx_CH) * idom * 13d0
        COI(i,j,k) = krome_x(krome_idx_CO) * idom * 28d0
        CI(i,j,k) = krome_x(krome_idx_C) * idom * 12d0
        CO2I(i,j,k) = krome_x(krome_idx_CO2) * idom * 44d0
        CH2I(i,j,k) = krome_x(krome_idx_CH2) * idom * 14d0
        CH3I(i,j,k) = krome_x(krome_idx_CH3) * idom * 15d0
        CH4I(i,j,k) = krome_x(krome_idx_CH4) * idom * 16d0
        HII(i,j,k) = krome_x(krome_idx_Hj) * idom
        HeII(i,j,k) = krome_x(krome_idx_HEj) * idom * 4d0
        H2II(i,j,k) = krome_x(krome_idx_H2j) * idom * 2d0
        OII(i,j,k) = krome_x(krome_idx_Oj) * idom * 16d0
        OHII(i,j,k) = krome_x(krome_idx_OHj) * idom * 17d0
        H2OII(i,j,k) = krome_x(krome_idx_H2Oj) * idom * 18d0
        H3OII(i,j,k) = krome_x(krome_idx_H3Oj) * idom * 19d0
        O2II(i,j,k) = krome_x(krome_idx_O2j) * idom * 32d0
        CII(i,j,k) = krome_x(krome_idx_Cj) * idom * 12d0
        COII(i,j,k) = krome_x(krome_idx_COj) * idom * 28d0
        HeIII(i,j,k) = krome_x(krome_idx_HEjj) * idom * 4d0

        !evaluate energy from temperature difference
        edot = (tgas - tgasold) * d(i,j,k) &
            / ((gamma - 1.d0) * utem * dt)

        !update internal energy
        e(i,j,k)  = e(i,j,k) + edot / d(i,j,k) * dt
        !when using dual
        if (idual .eq. 1) ge(i,j,k) = ge(i,j,k)+ edot / d(i,j,k) * dt
      end do
    end do
  end do

  !scale comoving<-proper
  factor = aye**3
  do k = ks+1, ke+1
    do j = js+1, je+1
      do i = is+1, ie+1
        d(i,j,k) = d(i,j,k) * factor
        !scale comoving->proper
        De(i,j,k) = De(i,j,k) * factor
        HM(i,j,k) = HM(i,j,k) * factor
        HI(i,j,k) = HI(i,j,k) * factor
        HeI(i,j,k) = HeI(i,j,k) * factor
        H2I(i,j,k) = H2I(i,j,k) * factor
        OI(i,j,k) = OI(i,j,k) * factor
        OHI(i,j,k) = OHI(i,j,k) * factor
        H2OI(i,j,k) = H2OI(i,j,k) * factor
        O2I(i,j,k) = O2I(i,j,k) * factor
        CHI(i,j,k) = CHI(i,j,k) * factor
        COI(i,j,k) = COI(i,j,k) * factor
        CI(i,j,k) = CI(i,j,k) * factor
        CO2I(i,j,k) = CO2I(i,j,k) * factor
        CH2I(i,j,k) = CH2I(i,j,k) * factor
        CH3I(i,j,k) = CH3I(i,j,k) * factor
        CH4I(i,j,k) = CH4I(i,j,k) * factor
        HII(i,j,k) = HII(i,j,k) * factor
        HeII(i,j,k) = HeII(i,j,k) * factor
        H2II(i,j,k) = H2II(i,j,k) * factor
        OII(i,j,k) = OII(i,j,k) * factor
        OHII(i,j,k) = OHII(i,j,k) * factor
        H2OII(i,j,k) = H2OII(i,j,k) * factor
        H3OII(i,j,k) = H3OII(i,j,k) * factor
        O2II(i,j,k) = O2II(i,j,k) * factor
        CII(i,j,k) = CII(i,j,k) * factor
        COII(i,j,k) = COII(i,j,k) * factor
        HeIII(i,j,k) = HeIII(i,j,k) * factor

        !mimimal value check
        De(i,j,k) = max(De(i,j,k), krome_tiny)
        HM(i,j,k) = max(HM(i,j,k), krome_tiny)
        HI(i,j,k) = max(HI(i,j,k), krome_tiny)
        HeI(i,j,k) = max(HeI(i,j,k), krome_tiny)
        H2I(i,j,k) = max(H2I(i,j,k), krome_tiny)
        OI(i,j,k) = max(OI(i,j,k), krome_tiny)
        OHI(i,j,k) = max(OHI(i,j,k), krome_tiny)
        H2OI(i,j,k) = max(H2OI(i,j,k), krome_tiny)
        O2I(i,j,k) = max(O2I(i,j,k), krome_tiny)
        CHI(i,j,k) = max(CHI(i,j,k), krome_tiny)
        COI(i,j,k) = max(COI(i,j,k), krome_tiny)
        CI(i,j,k) = max(CI(i,j,k), krome_tiny)
        CO2I(i,j,k) = max(CO2I(i,j,k), krome_tiny)
        CH2I(i,j,k) = max(CH2I(i,j,k), krome_tiny)
        CH3I(i,j,k) = max(CH3I(i,j,k), krome_tiny)
        CH4I(i,j,k) = max(CH4I(i,j,k), krome_tiny)
        HII(i,j,k) = max(HII(i,j,k), krome_tiny)
        HeII(i,j,k) = max(HeII(i,j,k), krome_tiny)
        H2II(i,j,k) = max(H2II(i,j,k), krome_tiny)
        OII(i,j,k) = max(OII(i,j,k), krome_tiny)
        OHII(i,j,k) = max(OHII(i,j,k), krome_tiny)
        H2OII(i,j,k) = max(H2OII(i,j,k), krome_tiny)
        H3OII(i,j,k) = max(H3OII(i,j,k), krome_tiny)
        O2II(i,j,k) = max(O2II(i,j,k), krome_tiny)
        CII(i,j,k) = max(CII(i,j,k), krome_tiny)
        COII(i,j,k) = max(COII(i,j,k), krome_tiny)
        HeIII(i,j,k) = max(HeIII(i,j,k), krome_tiny)

      end do
    end do
  end do

end subroutine krome_driver
