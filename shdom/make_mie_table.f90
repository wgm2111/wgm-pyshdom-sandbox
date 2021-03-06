PROGRAM MAKE_MIE_TABLE
! 
! Does Mie computations to create a scattering table as a function of
! effective radius for gamma, modified gamma, or lognormal size distributions 
! of spherical particles.  The particles may be water or ice (in which case
! the program provides the index of refraction depending on wavelength) or
! "aerosols" (in which case the index of refraction is user specified).
! For water or ice particles the scattering properties may be averaged 
! over the desired spectral range with Planck function weighting.  
! The phase functions in the output scattering table are represented 
! with Legendre series.  For polarized output, the six unique elements of
! the phase matrix are represented with Wigner d-function expansions 
! (Doicu et al., 2013, JQSRT, http://dx.doi.org/10.1016/j.jqsrt.2012.12.009).
!
!  compile: pgf90 -fast -o make_mie_table  make_mie_table.f90 
!                 indexwatice.f  miewig.f
!    
!    Frank Evans    University of Colorado       May 2003
!  Modified for polarization                     March 2012
!  Modified for modified gamma distribution      September 2012
!  Modified for Wigner d-function coefficients   July 2012
  IMPLICIT NONE
  INTEGER :: NRETAB
  LOGICAL :: POLTAB, LOGRE
  REAL    :: WAVELEN1, WAVELEN2, DELTAWAVE, PARDENS
  REAL    :: SRETAB, ERETAB, ALPHA, GAMMA, MAXRADIUS
  COMPLEX :: RINDEX
  CHARACTER(LEN=1) :: PARTYPE, AVGFLAG, DISTFLAG
  CHARACTER(LEN=80) :: MIETABFILE
  INTEGER :: NSIZE, MAXRANK, I, J, L, NL
  REAL    :: PI, WAVELENCEN, XMAX, REFF, SCATTER
  INTEGER, ALLOCATABLE :: NRANK1(:), NRANK(:)
  REAL, ALLOCATABLE :: RADII(:), ND(:)
  REAL, ALLOCATABLE :: QEXT(:), QSCA(:)
  REAL, ALLOCATABLE :: EXTINCT1(:), SCATTER1(:), WIGCOEF1(:,:,:)
  REAL, ALLOCATABLE :: EXTINCT(:), SSALB(:), WIGCOEF(:,:,:)


   ! Get the input parameters
  CALL USER_INPUT (POLTAB, WAVELEN1, WAVELEN2, PARTYPE, RINDEX, PARDENS, &
                   AVGFLAG, DELTAWAVE, DISTFLAG, ALPHA, GAMMA, &
                   NRETAB, SRETAB, ERETAB, LOGRE, MAXRADIUS,  MIETABFILE)
  IF (PARTYPE == 'W') THEN
    PARDENS = 1.0
  ELSE IF (PARTYPE == 'I') THEN
    PARDENS = 0.916
  ENDIF
  
   ! Calculate the maximum size parameter and the max number of Legendre terms
  CALL GET_CENTER_WAVELEN (WAVELEN1, WAVELEN2, WAVELENCEN)
  PI = ACOS(-1.0)
  IF (AVGFLAG == 'A') THEN
    XMAX = 2*PI*MAXRADIUS/WAVELEN1
  ELSE
    XMAX = 2*PI*MAXRADIUS/WAVELENCEN
  ENDIF  
  MAXRANK = NINT(2*(XMAX + 4.0*XMAX**0.3334 + 2))

   ! Get the average index of refraction for water or ice
  IF (PARTYPE /= 'A') THEN
    CALL GET_REFRACT_INDEX (PARTYPE, WAVELEN1, WAVELEN2, RINDEX)
  ENDIF

   ! Figure the number of radii there will be
  CALL GET_NSIZE (SRETAB, MAXRADIUS, WAVELENCEN, NSIZE)

   ! Allocate all the arrays here
  ALLOCATE (RADII(NSIZE), ND(NSIZE), NRANK1(NSIZE))
  ALLOCATE (EXTINCT1(NSIZE), SCATTER1(NSIZE), WIGCOEF1(6,0:MAXRANK,NSIZE))
  ALLOCATE (EXTINCT(NRETAB), SSALB(NRETAB))
  ALLOCATE (NRANK(NRETAB), WIGCOEF(6,0:MAXRANK,NRETAB))


   ! Make up the discrete particle radii to use
  CALL GET_SIZES (SRETAB, MAXRADIUS, WAVELENCEN, NSIZE, RADII)

   ! Do the Mie computations for each radius, which may involve several
   !   Mie calculation over the wavelength integration
  CALL COMPUTE_MIE_ALL_SIZES (AVGFLAG, WAVELEN1, WAVELEN2, DELTAWAVE, PARTYPE, &
                              WAVELENCEN, RINDEX, NSIZE, RADII, MAXRANK, &
                              EXTINCT1, SCATTER1, NRANK1, WIGCOEF1)


  ! Loop over the number of output tabulated effective radii
  DO I = 1, NRETAB
    ! Set tabulated effective radius
    IF (NRETAB <= 1) THEN
      REFF = SRETAB
    ELSE IF (LOGRE) THEN
      REFF = EXP((LOG(ERETAB)-LOG(SRETAB))*FLOAT(I-1)/(NRETAB-1) +LOG(SRETAB))
    ELSE
      REFF = (ERETAB-SRETAB)*FLOAT(I-1)/(NRETAB-1) + SRETAB
    ENDIF

    ! Calculate the discrete size number concentrations (ND), which vary
    ! according to a truncated gamma, modified gamma, or lognormal 
    ! distribution that gives the desired effective radius (REFF) and LWC (1 g/m^3).
    CALL MAKE_SIZE_DIST (DISTFLAG, PARDENS, NSIZE, RADII, REFF, ALPHA, GAMMA, &
                         ND)

    ! Sum the scattering properties over the discrete size distribution
    EXTINCT(I) = 0.0
    SCATTER = 0.0
    WIGCOEF(:,:,I) = 0.0
    NL = 1
    DO J = 1, NSIZE
      EXTINCT(I) = EXTINCT(I) + ND(J)*EXTINCT1(J)
      SCATTER = SCATTER + ND(J)*SCATTER1(J)
      NL = MAX(NL,NRANK1(J))
      WIGCOEF(:,0:NL,I) = WIGCOEF(:,0:NL,I) + ND(J)*WIGCOEF1(:,0:NL,J)
    ENDDO
    DO L = 0, NL
      WIGCOEF(:,L,I) = WIGCOEF(:,L,I)/SCATTER
      IF (WIGCOEF(1,L,I) .GT. 0.5E-5) NRANK(I) = L
    ENDDO
    IF (ABS(WIGCOEF(1,0,I)-1.0) > 0.0001) THEN
      PRINT *,'Phase function not normalized for Reff=',REFF,WIGCOEF(1,0,I)
      STOP
    ENDIF 
    IF (EXTINCT(I) > 0.0) THEN
      SSALB(I) = SCATTER/EXTINCT(I)
    ENDIF
    EXTINCT(I) = 0.001*EXTINCT(I)

  ENDDO  ! end of effective radius loop


  CALL WRITE_MIE_TABLE (MIETABFILE, POLTAB, WAVELEN1, WAVELEN2, DELTAWAVE, &
                        PARTYPE, PARDENS, RINDEX, &
                        DISTFLAG, ALPHA, GAMMA, NRETAB, SRETAB, ERETAB, LOGRE,&
                        EXTINCT, SSALB, MAXRANK, NRANK, WIGCOEF)

  DEALLOCATE (EXTINCT, SSALB, NRANK, WIGCOEF)
  DEALLOCATE (ND, EXTINCT1, SCATTER1, NRANK1, WIGCOEF1)
END




SUBROUTINE USER_INPUT (POLTAB, WAVELEN1,WAVELEN2, PARTYPE, RINDEX, PARDENS, &
                       AVGFLAG, DELTAWAVE, DISTFLAG, ALPHA, GAMMA, &
                       NRETAB, SRETAB, ERETAB, LOGRE, MAXRADIUS,  MIETABFILE)
 ! Reads the input parameters from the standard input
  IMPLICIT NONE
  INTEGER, INTENT(OUT) :: NRETAB
  LOGICAL, INTENT(OUT) :: POLTAB, LOGRE
  REAL,    INTENT(OUT) :: WAVELEN1, WAVELEN2, DELTAWAVE, PARDENS
  REAL,    INTENT(OUT) :: SRETAB, ERETAB, ALPHA, GAMMA, MAXRADIUS
  COMPLEX, INTENT(OUT) :: RINDEX
  CHARACTER(LEN=1), INTENT(OUT) :: PARTYPE, AVGFLAG, DISTFLAG
  CHARACTER(LEN=*), INTENT(OUT) :: MIETABFILE

  WRITE(*,*) 'Making Mie scattering tables for spherical particles'

  WRITE (*,*) 'Make polarized Mie table (T or F)'
  READ (*,*) POLTAB
    WRITE (*,'(L)') POLTAB

  WRITE (*,*) 'Wavelength range (micron)'
  READ (*,*) WAVELEN1, WAVELEN2
    WRITE (*,'(2(1X,F9.3))') WAVELEN1, WAVELEN2
  IF (WAVELEN1 > WAVELEN2) STOP 'USER_INPUT: wavelength1 must be <= wavelength2'

  WRITE(*,*) 'Water, Ice, or Aerosol spherical particles (W,I,A)'
  READ(*,'(A1)') PARTYPE
    WRITE(*,'(1X,A1)') PARTYPE

  IF (PARTYPE == 'W' .OR. PARTYPE == 'I') THEN
    WRITE(*,*) 'Average Mie properties over wavelength or use Planck weighted'
    WRITE(*,*) '  index of refraction at center wavelength (A or C)'
    READ(*,'(A1)') AVGFLAG
    WRITE (*,'(1X,A1)') AVGFLAG   
    IF (AVGFLAG == 'A') THEN
      WRITE(*,*) 'Wavelength interval for averaging (micron)'
      READ(*,*) DELTAWAVE
      WRITE(*,'(1X,F6.3)') DELTAWAVE
    ELSE
      AVGFLAG = 'C'
      DELTAWAVE = -1.0
    ENDIF
  ELSE
    PARTYPE = 'A'
    WRITE (*,*) 'Aerosol complex index of refraction (negative imaginary part)'
    READ (*,*) RINDEX 
    WRITE (*,*) RINDEX
    WRITE (*,*) 'Aerosol particle bulk density (g/cm^3)'
    READ (*,*) PARDENS
    WRITE (*,'(1X,F5.3)') PARDENS
    AVGFLAG = 'C'
  ENDIF

  WRITE (*,*) 'Particle size distribution type: G = Gamma, M = modified gamma, L = Lognormal'
  READ (*,*) DISTFLAG
  WRITE (*,*) DISTFLAG
  IF (DISTFLAG == 'L') THEN
    WRITE (*,*) 'Log normal size distribution log standard deviation'
    READ (*,*) ALPHA
    WRITE (*,'(1X,F6.3)') ALPHA
  ELSE IF (DISTFLAG == 'G') THEN
    WRITE(*,*) 'Gamma size distribution shape parameter (alpha)'
    READ (*,*) ALPHA
    WRITE (*,'(1X,F6.3)') ALPHA
  ELSE IF (DISTFLAG == 'M') THEN
    WRITE(*,*) 'Modified gamma size distribution shape parameters (alpha & gamma)'
    READ (*,*) ALPHA, GAMMA
    WRITE (*,'(2(1X,F6.3))') ALPHA, GAMMA
  ELSE
    WRITE (*,*) 'Unrecognized size distribution type'
    STOP
  ENDIF

  WRITE(*,*) 'Number, starting, and ending tabulated effective radius (micron)'
  READ(*,*) NRETAB, SRETAB, ERETAB
    WRITE (*,'(1X,I3,2(1X,F7.2))') NRETAB, SRETAB, ERETAB

  WRITE(*,*) 'Log-spaced effective radius (T or F) (F for evenly spaced)'
  READ(*,*) LOGRE
    WRITE (*,'(L)') LOGRE

  WRITE(*,*) 'Maximum particle radius in size distribution (micron)'
  READ(*,*) MAXRADIUS
    WRITE (*,'(1X,F7.2)') MAXRADIUS

  WRITE (*,*) 'Output Mie scattering table name'
  READ (*,*) MIETABFILE
    WRITE(*,'(1X,A70)') MIETABFILE
END SUBROUTINE USER_INPUT



SUBROUTINE GET_NSIZE (SRETAB, MAXRADIUS, WAVELEN, NSIZE)
 ! Calculates the number of radii for which the Mie computation will be run.
 ! The formula and spacing in size parameter can be changed to trade
 ! off size distribution integration accuracy vs. computer time.
  IMPLICIT NONE
  REAL,    INTENT(IN)  :: SRETAB, MAXRADIUS, WAVELEN
  INTEGER, INTENT(OUT) :: NSIZE
  REAL    :: TWOPI, RADMIN, RAD, X, DELX, DELRAD, DELRADMAX

  TWOPI = 2.0*ACOS(-1.0)
  RADMIN = 0.02*SRETAB
  RAD = RADMIN
  DELRADMAX = 0.002*MAXRADIUS
  NSIZE = 1
  DO WHILE (RAD < MAXRADIUS)
    X = TWOPI*RAD/WAVELEN
    DELX = MAX(0.01,0.03*X**0.5)    ! coarser spacing at large size parameters
!    DELX = 0.1                     ! One alternative method
    DELRAD = MIN(DELX*WAVELEN/TWOPI,DELRADMAX)
    RAD = RAD + DELRAD
    NSIZE = NSIZE + 1
  ENDDO
END SUBROUTINE GET_NSIZE


SUBROUTINE GET_SIZES (SRETAB, MAXRADIUS, WAVELEN, NSIZE, RADII)
 ! Calculates the radii for which the Mie computation will be run and
 ! from which all the size distributions will be computed.
 ! The formula and spacing in size parameter can be changed to trade
 ! off size distribution integration accuracy vs. computer time.
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NSIZE
  REAL,    INTENT(IN) :: SRETAB, MAXRADIUS, WAVELEN
  REAL,    INTENT(OUT) :: RADII(NSIZE)
  INTEGER :: N
  REAL    :: TWOPI, RADMIN, RAD, X, DELX, DELRAD, DELRADMAX

  TWOPI = 2.0*ACOS(-1.0)
  RAD = 0.02*SRETAB
  DELRADMAX = 0.002*MAXRADIUS
  RADII(1) = RAD
  DO N = 2, NSIZE
    X = TWOPI*RAD/WAVELEN
    DELX = MAX(0.01,0.03*X**0.5)    ! coarser spacing at large size parameters
!    DELX = 0.1                     ! One alternative method
    DELRAD = MIN(DELX*WAVELEN/TWOPI,DELRADMAX)
    RAD = RAD + DELRAD
    RADII(N) = RAD
  ENDDO
END SUBROUTINE GET_SIZES



SUBROUTINE GET_CENTER_WAVELEN (WAVELEN1, WAVELEN2, WAVELENCEN)
!  Returns the Planck weighted center wavelength averaged over the 
! wavelength interval (WAVELEN1 < WAVELEN2 [microns]).  A solar
! blackbody temperature of 5800 K is used for the Planck weighting
! if the average wavelength is less than 3 microns, no Planck weighting
! is done for an average wavelength between 3 and 5 microns, and a 
! blackbody temperature of 270 K is done for an average wavelength 
! greater than 5 microns.
  IMPLICIT NONE
  REAL, INTENT(IN)  :: WAVELEN1, WAVELEN2
  REAL, INTENT(OUT) :: WAVELENCEN
  REAL :: WAVECEN, DELWAVE, WAVE, BBTEMP, PLANCK, SUMP, SUMW

  IF (WAVELEN1 == WAVELEN2) THEN
    WAVELENCEN = WAVELEN1  
  ELSE
    WAVECEN = 0.5*(WAVELEN1+WAVELEN2)
    IF (WAVECEN < 3.0) THEN
      BBTEMP = 5800.0
    ELSE IF (WAVECEN > 5.0) THEN
      BBTEMP = 270.0
    ELSE
      BBTEMP = -1.0
      PLANCK = 1.0
    ENDIF 
    DELWAVE = MIN(WAVECEN/100.,0.1*ABS(WAVELEN2-WAVELEN1))
    SUMP = 0.0
    SUMW = 0.0
    WAVE = WAVELEN1   
    DO WHILE (WAVE <= WAVELEN2)
      IF (BBTEMP > 0) PLANCK = (1.19E8/WAVE**5)/(EXP(1.439E4/(WAVE*BBTEMP))-1)
      SUMP = SUMP + PLANCK
      SUMW = SUMW + PLANCK*WAVE
      WAVE = WAVE + DELWAVE
    ENDDO
    WAVELENCEN = 0.001*NINT(1000*SUMW/SUMP)
  ENDIF
END SUBROUTINE GET_CENTER_WAVELEN




SUBROUTINE GET_REFRACT_INDEX (PARTYPE, WAVELEN1, WAVELEN2, RINDEX)
 ! Returns the index of refraction for water or ice averaged over
 ! the wavelength interval (WAVELEN1 < WAVELEN2 [microns]).   The
 ! averaging is done at 0.05 micron intervals and is weighted by
 ! a Planck function at 5800 K temperature for central wavelengths
 ! less than 3 microns, a flat function between 3 and 5 microns, and
 ! 270 K Planck function beyond 5 microns.  
 ! The index of refraction is at -30 C for ice and +10 C for water
 ! (the temperature dependence is important in the microwave).
  IMPLICIT NONE
  CHARACTER(LEN=1), INTENT(IN) :: PARTYPE
  REAL, INTENT(IN) :: WAVELEN1, WAVELEN2
  COMPLEX, INTENT(OUT) :: RINDEX
  REAL :: WAVECEN, WAVECUT, DELWAVE, WAVE, BBTEMP, PLANCK
  REAL :: MRE, MIM, SUMP, SUMMR, SUMMI, A

  WAVECEN = 0.5*(WAVELEN1+WAVELEN2)
  IF (WAVECEN < 3.0) THEN
    BBTEMP = 5800.0
  ELSE IF (WAVECEN > 5.0) THEN
    BBTEMP = 270.0
  ELSE
    BBTEMP = -1.0
    PLANCK = 1.0
  ENDIF 
  DELWAVE = MIN(WAVECEN/100.,0.1*ABS(WAVELEN2-WAVELEN1))
  DELWAVE = MAX(DELWAVE,WAVECEN*1.0E-5)
  SUMP = 0.0
  SUMMR = 0.0
  SUMMI = 0.0
  WAVE = WAVELEN1
  DO WHILE (WAVE <= WAVELEN2)
    IF (BBTEMP > 0) PLANCK = (1.19E8/WAVE**5)/(EXP(1.439E4/(WAVE*BBTEMP))-1)
    SUMP = SUMP + PLANCK
    IF (PARTYPE == 'I') THEN
      CALL REFICE (0, WAVE, 243.0, MRE, MIM, A, A)
    ELSE
      CALL REFWAT (0, WAVE, 283.0, MRE, MIM, A, A)
    ENDIF
    SUMMR = SUMMR + PLANCK*MRE
    SUMMI = SUMMI + PLANCK*MIM
    WAVE = WAVE + DELWAVE
  ENDDO
  MRE = SUMMR/SUMP
  MIM = SUMMI/SUMP
  RINDEX = CMPLX(MRE,-MIM)
END SUBROUTINE GET_REFRACT_INDEX


 

SUBROUTINE COMPUTE_MIE_ALL_SIZES (AVGFLAG, WAVELEN1, WAVELEN2, DELTAWAVE, &
                                PARTYPE, WAVELENCEN, RINDEX, NSIZE, RADII, &
                                MAXRANK, EXTINCT1, SCATTER1, NRANK1, WIGCOEF1)
 ! Does a Mie computation for each particle radius in RADII and returns the
 ! optical properties in arrays EXTINCT1, SCATTER1, NRANK1, and WIGCOEF1.
 ! For AVGFLAG='C' the computation is done at a single wavelength (WAVELENCEN),
 ! using the input index of refraction (RINDEX).  For AVGFLAG='A' an
 ! integration of the Mie properties over wavelength is performed for
 ! each radius.  For each wavelength, with spacing DELTAWAVE, the water 
 ! or ice (depending on PARTYPE) index of refraction is obtained and
 ! used in the Mie computation for that wavelength, and the Mie optical
 ! properties are averaged with Planck function weighting (blackbody
 ! temperature depends on wavelength).  The Wigner d coefficients are
 ! returned with the product of the phase function times the scattering
 ! coefficient.
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NSIZE, MAXRANK
  REAL,    INTENT(IN) :: WAVELEN1, WAVELEN2, DELTAWAVE, WAVELENCEN
  REAL,    INTENT(IN) :: RADII(NSIZE)
  COMPLEX, INTENT(IN) :: RINDEX
  CHARACTER(LEN=1), INTENT(IN) :: AVGFLAG, PARTYPE
  INTEGER, INTENT(OUT) :: NRANK1(NSIZE)
  REAL,    INTENT(OUT) :: EXTINCT1(NSIZE), SCATTER1(NSIZE)
  REAL,    INTENT(OUT) :: WIGCOEF1(6,0:MAXRANK,NSIZE)
  INTEGER :: I, NL
  REAL    :: WAVECEN, WAVE, BBTEMP, PLANCK, SUMP, A
  REAL    :: MRE, MIM, EXT, SCAT, COEF(6,0:MAXRANK)
  COMPLEX :: REFIND

  IF (AVGFLAG == 'C') THEN
     ! For using one central wavelength: just call Mie routine for each radius
    DO I = 1, NSIZE
      CALL MIE_ONE (WAVELENCEN, RINDEX, RADII(I), MAXRANK, &
                    EXTINCT1(I), SCATTER1(I), NRANK1(I), WIGCOEF1(1,0,I) )
    ENDDO

  ELSE
     ! For averaging over wavelength range:
    WAVECEN = 0.5*(WAVELEN1+WAVELEN2)
    IF (WAVECEN < 3.0) THEN
      BBTEMP = 5800.0
    ELSE IF (WAVECEN > 5.0) THEN
      BBTEMP = 270.0
    ELSE
      BBTEMP = -1.0
      PLANCK = 1.0
    ENDIF 
    EXTINCT1(:) = 0.0
    SCATTER1(:) = 0.0
    NRANK1(:) = 1
    WIGCOEF1(:,:,:) = 0.0
    SUMP = 0.0
    WAVE = WAVELEN1
    DO WHILE (WAVE <= WAVELEN2)   ! Loop over the wavelengths
      IF (BBTEMP > 0) PLANCK = (1.19E8/WAVE**5)/(EXP(1.439E4/(WAVE*BBTEMP))-1)
      SUMP = SUMP + PLANCK
      IF (PARTYPE == 'I') THEN   ! Get the index of refraction of water or ice
        CALL REFICE (0, WAVE, 243.0, MRE, MIM, A, A)
      ELSE
        CALL REFWAT (0, WAVE, 283.0, MRE, MIM, A, A)
      ENDIF
      REFIND = CMPLX(MRE,-MIM)
      DO I = 1, NSIZE
        CALL MIE_ONE (WAVE, REFIND, RADII(I), MAXRANK, EXT, SCAT, NL, COEF)
        EXTINCT1(I) = EXTINCT1(I) + PLANCK*EXT
        SCATTER1(I) = SCATTER1(I) + PLANCK*SCAT
        NRANK1(I) = MAX(NRANK1(I),NL)
        WIGCOEF1(:,0:NL,I) = WIGCOEF1(:,0:NL,I) + PLANCK*COEF(:,0:NL)
      ENDDO
      WAVE = WAVE + DELTAWAVE
    ENDDO
    EXTINCT1(:) = EXTINCT1(:)/SUMP
    SCATTER1(:) = SCATTER1(:)/SUMP
    WIGCOEF1(:,:,:) = WIGCOEF1(:,:,:)/SUMP
  ENDIF
END SUBROUTINE COMPUTE_MIE_ALL_SIZES





SUBROUTINE MAKE_SIZE_DIST (DISTFLAG, PARDENS, NSIZE, RADII, REFF, ALPHA, GAMMA,&
                           ND)
 ! Calculates the number concentrations (ND in cm^-3) for the NSIZE
 ! discrete particle radii (micron) of a gamma or lognormal size distribution
 ! with an effective radius of REFF (micron), gamma shape parameter or 
 ! lognormal standard deviation of ALPHA, and mass content of 1 g/m^3.
 ! Also handles modified gamma distributions specified by ALPHA and GAMMA.
  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: NSIZE
  REAL,    INTENT(IN)  :: RADII(NSIZE), REFF, ALPHA, GAMMA, PARDENS
  REAL,    INTENT(OUT) :: ND(NSIZE)
  CHARACTER(LEN=1), INTENT(IN) :: DISTFLAG
  REAL, PARAMETER :: TOL=0.001  ! fractional tolerance in achieving Reff
  INTEGER :: I
  REAL    :: TRUERE, F, REHI, RELO, REMID

   ! See if the true effective radius is already close enough
  CALL DO_SIZE_DIST (PARDENS, DISTFLAG, ALPHA, GAMMA, REFF, &
                     NSIZE, RADII, ND, TRUERE)
  IF (ABS(TRUERE-REFF) < TOL*REFF) RETURN
  F = REFF/TRUERE

  IF (TRUERE < REFF) THEN
    ! Find Reff that gives true Reff above desired value
    RELO = REFF
    REHI = F*REFF
    I = 0
    TRUERE = REFF/F
    DO WHILE (TRUERE <= REFF .AND. I < 8)
      REHI = F*REHI
      I = I + 1
      CALL DO_SIZE_DIST (PARDENS,DISTFLAG, ALPHA, GAMMA, REHI, &
                         NSIZE, RADII, ND, TRUERE)
    ENDDO
    IF (TRUERE <= REFF) THEN
      PRINT *, 'MAKE_SIZE_DIST: effective radius cannot be achieved',REFF,TRUERE
      STOP
    ENDIF
  ELSE
    ! Find Reff that gives true Reff below desired value
    REHI = REFF
    RELO = F*REFF
    I = 0
    TRUERE = REFF/F
    DO WHILE (TRUERE >= REFF .AND. I < 8)
      RELO = F*RELO
      I = I + 1
      CALL DO_SIZE_DIST (PARDENS,DISTFLAG, ALPHA, GAMMA, RELO, &
                         NSIZE, RADII, ND, TRUERE)
    ENDDO
    IF (TRUERE >= REFF) THEN
      PRINT *, 'MAKE_SIZE_DIST: effective radius cannot be achieved',REFF,TRUERE
      STOP
    ENDIF
  ENDIF
  ! Do bisection to get correct effective radius
  DO WHILE (ABS(TRUERE-REFF) > TOL*REFF)
    REMID = 0.5*(RELO+REHI)
    CALL DO_SIZE_DIST (PARDENS, DISTFLAG, ALPHA, GAMMA, REMID, &
                       NSIZE, RADII, ND, TRUERE)
    IF (TRUERE < REFF) THEN
      RELO = REMID
    ELSE
      REHI = REMID
    ENDIF
  ENDDO  
END SUBROUTINE MAKE_SIZE_DIST



SUBROUTINE DO_SIZE_DIST (PARDENS, DISTFLAG, ALPHA, GAMMA, RE, NSIZE, RADII, &
                         ND, TRUERE)
 ! For the input effective radius (RE) [um], returns the number concentrations 
 ! ND [cm^-3] and the calculated effective radius TRUERE [um] for a 
 ! gamma, modified gamma, or lognormal size distribution with mass content 
 ! of 1 g/m^3.
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NSIZE
  REAL,    INTENT(IN) :: PARDENS, ALPHA, GAMMA, RE, RADII(NSIZE)
  CHARACTER(LEN=1), INTENT(IN) :: DISTFLAG
  REAL,    INTENT(OUT) :: ND(NSIZE), TRUERE
  INTEGER :: J
  REAL    :: DENS, PI, A, B, LWC, R, DELR, SUM2, SUM3, GAMMLN

  PI = ACOS(-1.0)
  IF (DISTFLAG == 'G') THEN
    B = (ALPHA+3)/RE
    A = 1.E6/( (4*PI/3.)*PARDENS *B**(-ALPHA-4) *EXP(GAMMLN(ALPHA+4.)) )
  ELSE IF (DISTFLAG == 'M') THEN
    B = (EXP(GAMMLN((ALPHA+4.)/GAMMA)-GAMMLN((ALPHA+3.)/GAMMA)) /RE)**GAMMA
    A = 1.E6*GAMMA *B**((ALPHA+4)/GAMMA) &
        /( (4*PI/3.)*PARDENS *EXP(GAMMLN((ALPHA+4.)/GAMMA)) )
  ELSE IF (DISTFLAG == 'L') THEN
    B = RE*EXP(-2.5*ALPHA**2)
    A = 1.E6/( (4*PI/3.)*PARDENS *SQRT(2*PI)*ALPHA * B**3 *EXP(4.5*ALPHA**2) )
  ENDIF
  LWC = 0.0
  SUM2 = 0.0
  SUM3 = 0.0
  DO J = 1, NSIZE
    R = RADII(J)
    DELR = SQRT(RADII(J)*RADII(MIN(NSIZE,J+1))) &
         - SQRT(RADII(J)*RADII(MAX(1,J-1)))
    IF (DISTFLAG == 'G') THEN
      ND(J) = A* R**ALPHA *EXP(-B*R) *DELR
    ELSE IF (DISTFLAG == 'M') THEN
      ND(J) = A* R**ALPHA *EXP(-B*R**GAMMA) *DELR
    ELSE IF (DISTFLAG == 'L') THEN
      ND(J) = (A/R)*EXP(-0.5*(LOG(R/B))**2/ALPHA**2) *DELR
    ENDIF
    LWC = LWC + 1.0E-6*PARDENS*ND(J)*(4*PI/3)*R**3
    SUM2 = SUM2 + ND(J)*R**2
    SUM3 = SUM3 + ND(J)*R**3
  ENDDO
  ND(:) = (1.0/LWC)*ND(:)
  TRUERE = SUM3/SUM2
END SUBROUTINE DO_SIZE_DIST




SUBROUTINE WRITE_MIE_TABLE (MIETABFILE, POLTAB, WAVELEN1,WAVELEN2, DELTAWAVE, &
                            PARTYPE, PARDENS, RINDEX, &
                            DISTFLAG, ALPHA, GAMMA, NRETAB,SRETAB,ERETAB,LOGRE,&
                            EXTINCT, SSALB, MAXRANK, NRANK, WIGCOEF)
 ! Writes the table of Mie scattering properties as a function of 
 ! effective radius.  There are two types of tables output: 
 ! unpolarized with Legendre (Wigner d_l_00) series of a single phase 
 ! matrix element (P11) or polarized with six Wigner d-function series 
 ! of the Mie phase matrix elements.
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NRETAB, MAXRANK, NRANK(NRETAB)
  LOGICAL, INTENT(IN) :: POLTAB, LOGRE
  COMPLEX, INTENT(IN) :: RINDEX
  REAL,    INTENT(IN) :: WAVELEN1, WAVELEN2, DELTAWAVE
  REAL,    INTENT(IN) :: PARDENS, SRETAB, ERETAB, ALPHA, GAMMA
  REAL,    INTENT(IN) :: EXTINCT(NRETAB), SSALB(NRETAB)
  REAL,    INTENT(IN) :: WIGCOEF(6,0:MAXRANK,NRETAB)
  CHARACTER(LEN=1), INTENT(IN) :: PARTYPE, DISTFLAG
  CHARACTER(LEN=*), INTENT(IN) :: MIETABFILE
  INTEGER :: I, J, JP, K, L, NL, NPHASEPOL
  REAL    REFF

  OPEN (UNIT=3, FILE=MIETABFILE, STATUS='REPLACE')
  IF (POLTAB) THEN
    NPHASEPOL = 6
    WRITE (3,'(A)') '! Polarized Mie scattering table vs. effective radius (LWC=1 g/m^3)'
  ELSE
    NPHASEPOL = 1
    WRITE (3,'(A)') '! Mie scattering table vs. effective radius (LWC=1 g/m^3)'
  ENDIF
  IF (DELTAWAVE < 0.0) THEN
    WRITE (3,'(2(1X,F8.3),A)') WAVELEN1, WAVELEN2, '  wavelength range (micron)'
  ELSE
    WRITE (3,'(3(1X,F8.3),A)') WAVELEN1, WAVELEN2, DELTAWAVE, '  wavelength range and averaging step (micron)'
  ENDIF
  WRITE (3,'(1X,F5.3,2X,A1,A)') PARDENS, PARTYPE, '   particle density (g/cm^3) and type (Water, Ice, Aerosol)'
  WRITE (3,'(2(1X,E13.6),A)') RINDEX, '  particle index of refraction'
  IF (DISTFLAG == 'L') THEN
    WRITE (3,'(F7.5,A)') ALPHA, '  lognormal log standard deviation'
  ELSE IF (DISTFLAG == 'G') THEN
    WRITE (3,'(F7.5,A)') ALPHA, ' gamma size distribution shape parameter'
  ELSE IF (DISTFLAG == 'M') THEN
    WRITE (3,'(2(1X,F7.5),A)') ALPHA, GAMMA, ' modified gamma size distribution shape parameters'
  ENDIF
  WRITE (3,'(1X,I3,2(1X,F8.3),A)') NRETAB, SRETAB, ERETAB, &
        '  number, starting, ending effective radius'

  DO I = 1, NRETAB
    IF (NRETAB <= 1) THEN
      REFF = SRETAB
    ELSE IF (LOGRE) THEN
      REFF = EXP((LOG(ERETAB)-LOG(SRETAB))*FLOAT(I-1)/(NRETAB-1) +LOG(SRETAB))
    ELSE
      REFF = (ERETAB-SRETAB)*FLOAT(I-1)/(NRETAB-1) + SRETAB
    ENDIF
    WRITE (3,'(1X,F8.4,1X,E12.5,1X,F8.6,1X,I6,A)') &
        REFF, EXTINCT(I), SSALB(I), NRANK(I), '  Reff  Ext  Alb  Nrank'
    DO J = 1, NPHASEPOL
      IF (POLTAB) THEN
        WRITE (3,'(I1,1X,201(1X,F10.5))') J, (WIGCOEF(J,L,I), L=0,MIN(NRANK(I),200))
      ELSE
        WRITE (3,'(2X,201(1X,F10.5))') (WIGCOEF(J,L,I), L=0,MIN(NRANK(I),200))
      ENDIF
      DO K = 200, NRANK(I)-1, 200
        WRITE (3,'(2X,200(1X,F10.5))') (WIGCOEF(J,K+L,I),L=1,MIN(200,NRANK(I)-K))
      ENDDO
    ENDDO
  ENDDO
  CLOSE (3)
END SUBROUTINE WRITE_MIE_TABLE

