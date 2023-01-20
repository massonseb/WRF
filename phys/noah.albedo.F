  SUBROUTINE ALCALC (ALB,SNOALB,EMBRD,ALBEDO,EMISSI,   &
       &             DT,SNOWNG,SNOTIME1)

! ----------------------------------------------------------------------
! CALCULATE ALBEDO INCLUDING SNOW EFFECT (0 -> 1)
!   ALB     SNOWFREE ALBEDO
!   SNOALB  MAXIMUM (DEEP) SNOW ALBEDO
!   ALBEDO  SURFACE ALBEDO INCLUDING SNOW EFFECT
!   TSNOW   SNOW SURFACE TEMPERATURE (K)
! ----------------------------------------------------------------------
    IMPLICIT NONE

! ----------------------------------------------------------------------
! SNOALB IS ARGUMENT REPRESENTING MAXIMUM ALBEDO OVER DEEP SNOW,
! AS PASSED INTO SFLX, AND ADAPTED FROM THE SATELLITE-BASED MAXIMUM
! SNOW ALBEDO FIELDS PROVIDED BY D. ROBINSON AND G. KUKLA
! (1985, JCAM, VOL 24, 402-411)
! ----------------------------------------------------------------------
    REAL,    INTENT(IN)    :: ALB, SNOALB, EMBRD
    REAL,    INTENT(IN)    :: DT
    LOGICAL, INTENT(IN)    :: SNOWNG
    REAL,    INTENT(INOUT) :: SNOTIME1
    REAL,    INTENT(OUT)   :: ALBEDO, EMISSI
    REAL                   :: SNOALB2
    REAL                   :: TM,SNOALB1
    REAL,    PARAMETER     :: SNACCA=0.94,SNACCB=0.58,SNTHWA=0.82,SNTHWB=0.46
    ALBEDO = ALB + (SNOALB-ALB)
    EMISSI = EMBRD + (EMISSI_S - EMBRD)
    SNOALB1 = SNOALB+LVCOEF_DATA*(0.85-SNOALB)
    SNOALB2=SNOALB1
    IF (SNOWNG) THEN
       SNOTIME1 = 0.
    ELSE
       SNOTIME1=SNOTIME1+DT
       SNOALB2=SNOALB1*(SNACCA**((SNOTIME1/86400.0)**SNACCB))
    ENDIF

    SNOALB2 = MAX ( SNOALB2, ALB )
    ALBEDO = ALB + (SNOALB2-ALB)
    IF (ALBEDO .GT. SNOALB2) ALBEDO=SNOALB2
! ----------------------------------------------------------------------
  END SUBROUTINE ALCALC

