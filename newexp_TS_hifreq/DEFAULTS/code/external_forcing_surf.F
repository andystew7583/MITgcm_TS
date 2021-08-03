C $Header: /u/gcmpack/MITgcm/model/src/external_forcing_surf.F,v 1.57 2013/07/04 23:05:09 jmc Exp $
C $Name: checkpoint64s $

#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"
#undef CHECK_OVERLAP_FORCING

CBOP
C     !ROUTINE: EXTERNAL_FORCING_SURF
C     !INTERFACE:
      SUBROUTINE EXTERNAL_FORCING_SURF(
     I             iMin, iMax, jMin, jMax,
     I             myTime, myIter, myThid )
C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE EXTERNAL_FORCING_SURF
C     | o Determines forcing terms based on external fields
C     |   relaxation terms etc.
C     *==========================================================*
C     \ev

C     !USES:
      IMPLICIT NONE
C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "FFIELDS.h"
#include "DYNVARS.h"
#include "GRID.h"
#include "SURFACE.h"
#ifdef ALLOW_SHELFICE
# include "SHELFICE.h"
#endif /* ALLOW_SHELFICE */
#ifdef ALLOW_AUTODIFF_TAMC
# include "tamc.h"
# include "tamc_keys.h"
#endif

C     !INPUT/OUTPUT PARAMETERS:
C     === Routine arguments ===
C     iMin,iMax, jMin,jMax :: Range of points for calculation
C     myTime :: Current time in simulation
C     myIter :: Current iteration number in simulation
C     myThid :: Thread no. that called this routine.
      INTEGER iMin, iMax
      INTEGER jMin, jMax
      _RL myTime
      INTEGER myIter
      INTEGER myThid

C     !LOCAL VARIABLES:
C     === Local variables ===
C     bi,bj  :: tile indices
C     i,j    :: loop indices
C     ks     :: index of surface interface layer
      INTEGER bi,bj
      INTEGER i,j
      INTEGER ks

C---------------------------------
C--   BEGIN STEWART CODE
C---------------------------------

C--   Constants 
      _RL L0_als      !-- Latent heat of melting at 0 deg C (J/kg)      
      _RL ci_als      !-- Specific heat capacity of ice (J/kg/C)
      _RL cw_als      !-- Specific heat capacity of water (J/kg/C)
      _RL rhow_als    !-- Reference density of water (kg/m^3)
      _RL rhoi_als    !-- Reference density of ice (kg/m^3)
      _RL kapi_als    !-- Molecular thermal diffusivity of ice (m^2/s)
      _RL gamT_als    !-- Thermal turbulent exchange velocity (m/s)
      _RL delh_als    !-- Ice thickness (m)
      _RL frS_als     !-- Ice salinity retention (dimensionless)
      _RL mu_als      !-- Linear freezing temperature coefficient (deg C/psu)
      _RL Ti_als      !-- Ice temperature (deg C)
      _RL Si_als      !-- Ice salinity (psu)
      _RL lam_als     !-- Vertical temperature diffusivity in ice (J/deg C/m/s)

      _RL Lshelf_als  !-- Width of ice shelf (m)
C--   TODO Lshelf should be an external parameter

C--   Reciprocals
      _RL recip_cw_als
      _RL recip_delh_als 
      _RL recip_rhow_als

C--   Variables
      LOGICAL isMelting_als    !-- Determines whether ice is melting or freezing
      _RL To_als               !-- Ocean temperature (deg C)
      _RL So_als               !-- Ocean salinity (psu)
      _RL Sb_als               !-- Boundary layer salinity (psu)
      _RL Tb_als               !-- Boundary layer temperature (deg C)
      _RL Tib_als              !-- Temperature of source (ice or ocean)
      _RL Sib_als              !-- Salinity of source (ice or ocean)
      _RL Eo_als               !-- Ocean internal energy
      _RL Eb_als               !-- Boundary layer internal energy
      _RL Eib_als              !-- Internal energy of source (ice or ocean)
      _RL Fdiff_als            !-- Vertical diffusive heat transport through sea ice
      _RL Fturb_als            !-- Turbulent heat transport into ocean-ice boundary layer
      _RL Fm_als               !-- Downward mass flux due to melting/freezing
      _RL FH_als               !-- Downward heat flux due to melting/freezing
      _RL FS_als               !-- Downward salt flux due to melting/freezing
      _RL FHeff_als            !-- Effective downward heat flux
      _RL FSeff_als            !-- Effective downward salt flux

C---------------------------------
C--   END STEWART CODE
C---------------------------------


#ifdef ALLOW_BALANCE_FLUXES
      _RS tmpVar(1)
#endif
#ifdef CHECK_OVERLAP_FORCING
      _RS fixVal
#endif
CEOP

      IF ( usingPCoords ) THEN
       ks        = Nr
      ELSE
       ks        = 1
      ENDIF

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

C--   Apply adjustment (balancing forcing) and exchanges
C     to oceanic surface forcing

#ifdef ALLOW_BALANCE_FLUXES
C     balance fluxes
      tmpVar(1) = oneRS
      IF ( balanceEmPmR .AND. (.NOT.useSeaice .OR. useThSIce) ) THEN
        CALL REMOVE_MEAN_RS( 1, EmPmR, maskInC, maskInC, rA, tmpVar,
     &                       'EmPmR', myTime, myThid )
      ENDIF
      IF ( balanceQnet  .AND. (.NOT.useSeaice .OR. useThSIce) ) THEN
        CALL REMOVE_MEAN_RS( 1, Qnet,  maskInC, maskInC, rA, tmpVar,
     &                       'Qnet ', myTime, myThid )
      ENDIF
#endif /* ALLOW_BALANCE_FLUXES */

C-    Apply exchanges (if needed)

#ifdef CHECK_OVERLAP_FORCING
C     Put large value in overlap of forcing array to check if exch is needed
c     IF ( .NOT. useKPP ) THEN
       fixVal = 1.
       CALL RESET_HALO_RS ( EmPmR, fixVal, 1, myThid )
       fixVal = 400.
       CALL RESET_HALO_RS ( Qnet, fixVal, 1, myThid )
       fixVal = -200.
       CALL RESET_HALO_RS ( Qsw, fixVal, 1, myThid )
       fixVal = 40.
       CALL RESET_HALO_RS ( saltFlux, fixVal, 1, myThid )
c     ENDIF
#endif /* CHECK_OVERLAP_FORCING */

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

#ifdef EXACT_CONSERV
C NB: synchronous time step: PmEpR lag 1 time step behind EmPmR
C     to stay consitent with volume change (=d/dt etaH).
# ifdef ALLOW_AUTODIFF_TAMC
CADJ STORE PmEpR = comlev1, key = ikey_dynamics,  kind = isbyte
# endif
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        IF ( staggerTimeStep ) THEN
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
           PmEpR(i,j,bi,bj) = -EmPmR(i,j,bi,bj)
          ENDDO
         ENDDO
        ENDIF
       ENDDO
      ENDDO
#endif /* EXACT_CONSERV */

C--   set surfaceForcingT,S to zero.
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
           surfaceForcingT(i,j,bi,bj) = 0. _d 0
           surfaceForcingS(i,j,bi,bj) = 0. _d 0
         ENDDO
        ENDDO
       ENDDO
      ENDDO

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

C--   Start with surface restoring term :
      IF ( doThetaClimRelax .OR. doSaltClimRelax ) THEN
       CALL FORCING_SURF_RELAX(
     I              iMin, iMax, jMin, jMax,
     I              myTime, myIter, myThid )
      ENDIF

#ifdef ALLOW_PTRACERS
C--   passive tracer surface forcing:
#ifdef ALLOW_AUTODIFF_TAMC
CADJ STORE surfaceForcingS = comlev1, key = ikey_dynamics,
CADJ &    kind = isbyte
#endif
      IF ( usePTRACERS ) THEN
       DO bj=myByLo(myThid),myByHi(myThid)
        DO bi=myBxLo(myThid),myBxHi(myThid)
         CALL PTRACERS_FORCING_SURF(
     I        surfaceForcingS(1-OLx,1-OLy,bi,bj),
     I        bi, bj, iMin, iMax, jMin, jMax,
     I        myTime,myIter,myThid )
        ENDDO
       ENDDO
      ENDIF
#endif /* ALLOW_PTRACERS */

C- Notes: setting of PmEpR and pTracers surface forcing could have been
C         moved below, inside a unique bi,bj block. However this results
C         in tricky dependencies for TAF (and recomputations).
C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)

C--   Surface Fluxes :
        DO j = jMin, jMax
         DO i = iMin, iMax

C     Zonal wind stress fu:
          surfaceForcingU(i,j,bi,bj) = fu(i,j,bi,bj)*mass2rUnit
C     Meridional wind stress fv:
          surfaceForcingV(i,j,bi,bj) = fv(i,j,bi,bj)*mass2rUnit
C     Net heat flux Qnet:
          surfaceForcingT(i,j,bi,bj) = surfaceForcingT(i,j,bi,bj)
     &       - ( Qnet(i,j,bi,bj)
#ifdef SHORTWAVE_HEATING
     &          -Qsw(i,j,bi,bj)
#endif
     &         ) *recip_Cp*mass2rUnit
C     Net Salt Flux :
          surfaceForcingS(i,j,bi,bj) = surfaceForcingS(i,j,bi,bj)
     &      -saltFlux(i,j,bi,bj)*mass2rUnit

         ENDDO
        ENDDO




C---------------------------------
C--   BEGIN STEWART CODE
C---------------------------------

C--   Calculate heat/ salt fluxes using the 2-equation model of 
C--   Schmidt et al. (2004) or Holland & Jenkins (1999)

C--   Set constants 
      L0_als = 3.34 _d 5                  !-- Latent heat of melting at 0 deg C (J/kg)      
      ci_als = 2.06 _d 3                  !-- Specific heat capacity of ice (J/kg/C)
      cw_als = HeatCapacity_Cp            !-- Specific heat capacity of water (J/kg/C)
      rhow_als = rhoConst                 !-- Reference density of water (kg/m^3)
      rhoi_als = 9.2 _d 2                 !-- Reference density of ice (kg/m^3)
      kapi_als = 1.14 _d -6               !-- Molecular thermal diffusivity of ice (m^2/s)
      gamT_als = 6. _d -5                 !-- Thermal turbulent exchange velocity (m/s)
      delh_als = 1. _d 0                  !-- Ice thickness (m)
      frS_als = 1.4 _d -1                 !-- Ice salinity retention (dimensionless)
      mu_als = 5.4 _d -2                  !-- Linear freezing temperature coefficient (deg C/psu)
      Ti_als = - 5. _d 0                  !-- Ice temperature (deg C)
      Si_als = 5. _d 0                    !-- Ice salinity (psu)
      lam_als = rhoi_als*ci_als*kapi_als  !-- Vertical temperature diffusivity in ice (J/deg C/m/s)
      Lshelf_als = 5. _d 4                !-- Width of ice shelf (m)      

C--   Reciprocals
      recip_cw_als = recip_Cp
      recip_delh_als = 1. _d 0 / delh_als
      recip_rhow_als = recip_rhoConst

      DO j = jMin, jMax

        IF (yC(i,j,bi,bj) > Lshelf_als) THEN

          DO i = iMin, iMax

C--         Ocean and boundary layer temperature/salinity
            To_als = theta(i,j,ks,bi,bj)     !-- Ocean temperature (deg C)
            So_als = salt(i,j,ks,bi,bj)      !-- Ocean salinity (psu)
            Sb_als = So_als                  !-- Boundary layer salinity (psu)
            Tb_als = -mu_als*Sb_als          !-- Boundary layer temperature (deg C)

C--         Compute heat fluxes due to diffusion within ice and
C--         turbulent exchance with the ocean-ice boundary layer, 
C--         in order to determine whether melting or freezing is 
C--         taking place
            Fdiff_als = - lam_als*(Ti_als-Tb_als)*recip_delh_als
            Fturb_als = rhow_als*cw_als*gamT_als*(Tb_als-To_als)       
            isMelting_als = Fdiff_als+Fturb_als <= 0

C--         Determine source temperature and salinity depending on 
C--         whether melting or freezing is taking place
            IF (isMelting_als) THEN
              Tib_als = Ti_als
              Sib_als = Si_als
            ELSE         
              Tib_als = Tb_als
              Sib_als = frS_als*Sb_als
            ENDIF

C--         Internal energies
            Eo_als = cw_als*To_als 
            Eb_als = cw_als*Tb_als
            Eib_als = - L0_als + ci_als*Tib_als

C--         Downward mass, heat and salt fluxes
            Fm_als = - (Fdiff_als+Fturb_als) / (Eb_als-Eib_als)
            FH_als = - Fdiff_als + Fm_als*Eib_als
            FS_als = Fm_als*Sib_als

C--         "Virtual" heat and salt fluxes into the ocean
            FHeff_als = FH_als  - Fm_als*Eo_als
            FSeff_als = FS_als - Fm_als*So_als

C--         Add fluxes to surface forcing
            surfaceForcingT(i,j,bi,bj) = surfaceForcingT(i,j,bi,bj)
     &         + FHeff_als*recip_cw_als*recip_rhow_als
            surfaceForcingS(i,j,bi,bj) = surfaceForcingS(i,j,bi,bj)
     &         + FSeff_als*recip_rhow_als

          ENDDO

        ELSE

          To_als = theta(i,j,ks,bi,bj)     !-- Ocean temperature (deg C)
          So_als = salt(i,j,ks,bi,bj)      !-- Ocean salinity (psu)
          Tb_als = -mu_als*So_als          !-- Boundary layer (freezing) temperature (deg C)

C--       Restore to freezing temperature, using turbulent exchange velocity
          surfaceForcingT(i,j,bi,bj) = surfaceForcingT(i,j,bi,bj)
     &         + gamT_als*(Tb_als-To_als)          

        ENDIF

      ENDDO

C---------------------------------
C--   END STEWART CODE
C---------------------------------




#ifdef ALLOW_SALT_PLUME
C saltPlume is the amount of salt rejected by ice while freezing;
C it is here subtracted from surfaceForcingS and will be redistributed
C to multiple vertical levels later on as per Duffy et al. (GRL 1999)
        IF ( useSALT_PLUME ) THEN
         CALL SALT_PLUME_FORCING_SURF(
     I        bi, bj, iMin, iMax, jMin, jMax,
     I        myTime,myIter,myThid )
        ENDIF
#endif /* ALLOW_SALT_PLUME */

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C--   Fresh-water flux

C-    Apply mask on Fresh-Water flux (if useRealFreshWaterFlux)
C     <== removed: maskInC is applied directly in S/R SOLVE_FOR_PRESSURE

#ifdef EXACT_CONSERV
      IF ( (nonlinFreeSurf.GT.0 .OR. usingPCoords)
     &     .AND. useRealFreshWaterFlux ) THEN

C--   NonLin_FrSurf and RealFreshWaterFlux : PmEpR effectively changes
C     the water column height ; temp., salt, (tracer) flux associated
C     with this input/output of water is added here to the surface tendency.

       IF (temp_EvPrRn.NE.UNSET_RL) THEN
        DO j = jMin, jMax
         DO i = iMin, iMax
          surfaceForcingT(i,j,bi,bj) = surfaceForcingT(i,j,bi,bj)
     &      + PmEpR(i,j,bi,bj)
     &          *( temp_EvPrRn - theta(i,j,ks,bi,bj) )
     &          *mass2rUnit
         ENDDO
        ENDDO
       ENDIF

       IF (salt_EvPrRn.NE.UNSET_RL) THEN
        DO j = jMin, jMax
         DO i = iMin, iMax
          surfaceForcingS(i,j,bi,bj) = surfaceForcingS(i,j,bi,bj)
     &      + PmEpR(i,j,bi,bj)
     &          *( salt_EvPrRn - salt(i,j,ks,bi,bj) )
     &          *mass2rUnit
         ENDDO
        ENDDO
       ENDIF

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      ELSE
#else /* EXACT_CONSERV */
      IF (.TRUE.) THEN
#endif /* EXACT_CONSERV */

C--   EmPmR does not really affect the water column height (for tracer budget)
C     and is converted to a salt tendency.

       IF (convertFW2Salt .EQ. -1.) THEN
C-    use local surface tracer field to calculate forcing term:

        IF (temp_EvPrRn.NE.UNSET_RL) THEN
C     account for Rain/Evap heat content (temp_EvPrRn) using local SST
         DO j = jMin, jMax
          DO i = iMin, iMax
           surfaceForcingT(i,j,bi,bj) = surfaceForcingT(i,j,bi,bj)
     &       + EmPmR(i,j,bi,bj)
     &           *( theta(i,j,ks,bi,bj) - temp_EvPrRn )
     &           *mass2rUnit
          ENDDO
         ENDDO
        ENDIF
        IF (salt_EvPrRn.NE.UNSET_RL) THEN
C     converts EmPmR to salinity tendency using surface local salinity
         DO j = jMin, jMax
          DO i = iMin, iMax
           surfaceForcingS(i,j,bi,bj) = surfaceForcingS(i,j,bi,bj)
     &       + EmPmR(i,j,bi,bj)
     &           *( salt(i,j,ks,bi,bj) - salt_EvPrRn )
     &           *mass2rUnit
          ENDDO
         ENDDO
        ENDIF

       ELSE
C-    use uniform tracer value to calculate forcing term:

        IF (temp_EvPrRn.NE.UNSET_RL) THEN
C     account for Rain/Evap heat content (temp_EvPrRn) assuming uniform SST (=tRef)
         DO j = jMin, jMax
          DO i = iMin, iMax
           surfaceForcingT(i,j,bi,bj) = surfaceForcingT(i,j,bi,bj)
     &       + EmPmR(i,j,bi,bj)
     &           *( tRef(ks) - temp_EvPrRn )
     &           *mass2rUnit
          ENDDO
         ENDDO
        ENDIF
        IF (salt_EvPrRn.NE.UNSET_RL) THEN
C     converts EmPmR to virtual salt flux using uniform salinity (default=35)
         DO j = jMin, jMax
          DO i = iMin, iMax
           surfaceForcingS(i,j,bi,bj) = surfaceForcingS(i,j,bi,bj)
     &       + EmPmR(i,j,bi,bj)
     &           *( convertFW2Salt - salt_EvPrRn )
     &           *mass2rUnit
          ENDDO
         ENDDO
        ENDIF

C-    end local-surface-tracer / uniform-value distinction
       ENDIF

      ENDIF

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

#ifdef ATMOSPHERIC_LOADING
C-- Atmospheric surface Pressure loading : added to phi0surf when using Z-coord;
C   Not yet implemented for Ocean in P: would need to be applied to the other end
C   of the column, as a vertical velocity (omega); (meaningless for Atmos in P).
C- Note:
C   Using P-coord., a hack (now directly applied from S/R INI_FORCING)
C   is sometime used to read phi0surf from a file (pLoadFile) instead
C   of computing it from bathymetry & density ref. profile.

        IF ( usingZCoords ) THEN
C   The true atmospheric P-loading is not yet implemented for P-coord
C   (requires time varying dP(Nr) like dP(k-bottom) with NonLin FS).
         IF ( useRealFreshWaterFlux ) THEN
          DO j = jMin, jMax
           DO i = iMin, iMax
            phi0surf(i,j,bi,bj) = ( pLoad(i,j,bi,bj)
     &                          +sIceLoad(i,j,bi,bj)*gravity
     &                            )*recip_rhoConst
           ENDDO
          ENDDO
         ELSE
          DO j = jMin, jMax
           DO i = iMin, iMax
            phi0surf(i,j,bi,bj) = pLoad(i,j,bi,bj)*recip_rhoConst
           ENDDO
          ENDDO
         ENDIF
c       ELSEIF ( usingPCoords ) THEN
C-- This is a hack used to read phi0surf from a file (pLoadFile)
C   instead of computing it from bathymetry & density ref. profile.
C   ==> now done only once, in S/R INI_FORCING
c         DO j = jMin, jMax
c          DO i = iMin, iMax
c           phi0surf(i,j,bi,bj) = pLoad(i,j,bi,bj)
c          ENDDO
c         ENDDO
        ENDIF
#endif /* ATMOSPHERIC_LOADING */

#ifdef ALLOW_SHELFICE
        IF ( usingZCoords ) THEN
         IF ( useSHELFICE) THEN
          DO j = jMin, jMax
           DO i = iMin, iMax
            phi0surf(i,j,bi,bj) = phi0surf(i,j,bi,bj)
     &         + shelficeLoadAnomaly(i,j,bi,bj)*recip_rhoConst
           ENDDO
          ENDDO
         ENDIF
        ENDIF
#endif /* ALLOW_SHELFICE */

C--   end bi,bj loops.
       ENDDO
      ENDDO

      RETURN
      END
