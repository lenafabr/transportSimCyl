MODULE PARTICLEUTIL
! utilities for defining parameters and states of groups of particles

  TYPE PARTICLEGROUP
  ! object describing a group of particles and microtubules

    INTEGER :: NPART ! number of particles in the group

    INTEGER :: NMT ! number of microtubules

    ! particle position
    ! first index is x,y,z, second is particle number
    DOUBLE PRECISION, POINTER :: PARTPOS(:,:)

    ! microtubule tip position
    ! first index is x,y,z; second is MT number
    DOUBLE PRECISION, POINTER :: MTPOS(:,:)

    ! microtubule nucleation point
    DOUBLE PRECISION, POINTER :: MTNUC(:,:)

    ! which track is each particle on
    INTEGER, POINTER :: MTINDS(:)

    ! microtubule state
    ! 1: growing, 2: shrinking, 3: stable
    INTEGER, POINTER :: MTSTATE(:)

    ! current direction
    ! +1 = processive to the right
    ! -1 = processive to the left
    ! 0 = diffusive
    INTEGER, POINTER :: CURDIR(:)

    ! absorption times
    DOUBLE PRECISION, POINTER :: ABSTIMES(:)

    ! absorbing position
    DOUBLE PRECISION, POINTER :: ABSPOS(:,:)

    ! is the particle active?
    LOGICAL, POINTER :: ISACTIVE(:)

    ! arrays have been allocated
    LOGICAL :: ARRAYSET = .FALSE.
      
  END TYPE PARTICLEGROUP

  TYPE PARAMLIST
  ! object to store relevant parameters

    ! number of particles
    INTEGER :: NPART

    ! number of microtubules
    INTEGER :: NMT

    ! walk velocity
    DOUBLE PRECISION :: VEL

    ! diffusion coefficient
    DOUBLE PRECISION :: DIFFCONST

    ! domain length
    DOUBLE PRECISION :: DOMLEN

    ! domain radius
    DOUBLE PRECISION :: DOMRAD

    ! binding rate
    DOUBLE PRECISION :: KBIND

    ! direction reversal rate
    DOUBLE PRECISION :: KREV

    ! capture radius
    DOUBLE PRECISION :: CRAD

    ! jump rate
    DOUBLE PRECISION :: KJUMP

    ! time step
    DOUBLE PRECISION :: DELT

    ! parameters for microtubule dynamics
    DOUBLE PRECISION :: V_GROWTH,V_SHRINK
    DOUBLE PRECISION :: K_CAT,K_RESCUE,K_PAUSE

    ! starting state for particles
    ! 1 = start on MT
    ! 2 = start at domain end
    ! 3 = start uniformly within the domain (not on microtubules)
    INTEGER :: STARTSTATE

    ! are particles only captured at + ends?
    LOGICAL :: TIPCAPTURE

    ! which state are particles captured in
    ! 1 = any state
    ! 2 = growing and paused
    ! 3 = only paused
    INTEGER :: CAPTUREINSTATE

    ! fpt to capture or to reach cell body?
    LOGICAL :: STOPONCAPTURE

    ! do particles fall off on reaching MT ends?
    LOGICAL :: STAYONMT

    ! do particles switch directions while jumping on another microtubule?
    LOGICAL :: SWITCHONJUMP 

    ! do MTs start at radial center of domain
    LOGICAL :: STARTATCENTER

    ! do particles get absorbed in the soma (instead of reflecting back)
    ! relevant for tip capture
    LOGICAL :: ABSBOUND

    ! arrays have been allocated
    LOGICAL :: ARRAYSET = .FALSE.
  END TYPE PARAMLIST

CONTAINS

  SUBROUTINE OUTPUTSNAPSHOT(GROUPLIST,FILENAME,INFO,NTRIALS,APPEND)
  ! Output a snapshot of current particle configuration
  ! ------------------
  ! input parameters:
  ! GROUPLIST: object describing multiple particle groups
  ! FILENAME: output file
  ! INFO: additional informational float (eg: time)
  ! APPEND: append to file or rewrite?
  ! -------------
  ! output format:
  ! 
    IMPLICIT NONE

    TYPE(PARTICLEGROUP), TARGET :: GROUPLIST(NTRIALS)
    TYPE(PARTICLEGROUP), POINTER :: PGROUP
    INTEGER, INTENT(IN) :: NTRIALS
    CHARACTER(LEN=*), INTENT(IN) :: FILENAME
    DOUBLE PRECISION, INTENT(IN) :: INFO(:)
    LOGICAL, INTENT(IN) :: APPEND

    INTEGER :: PC,MC,NINFO,TRC
    CHARACTER(LEN=2) :: X1
    CHARACTER(LEN=19) :: FMT1

    IF (APPEND) THEN
      OPEN(UNIT=99,FILE=FILENAME,STATUS='UNKNOWN',POSITION='APPEND')
    ELSE
        OPEN(UNIT=99,FILE=FILENAME,STATUS='UNKNOWN')
    END IF

    NINFO = SIZE(INFO)
    WRITE (X1,'(I2)') NINFO
    FMT1 = '('//TRIM(X1)//'ES15.6E2)'
    WRITE(99,FMT1) INFO

    DO TRC = 1,NTRIALS
      PGROUP=>GROUPLIST(TRC)

      DO MC = 1,PGROUP%NMT
        WRITE(99,'(2I8,6ES18.4E2)') MC, PGROUP%MTSTATE(MC), PGROUP%MTPOS(:,MC), PGROUP%MTNUC(:,MC)
      END DO

      DO PC = 1,PGROUP%NPART
        WRITE(99,'(2I8,3ES18.4E2)') PC, PGROUP%MTINDS(PC), PGROUP%PARTPOS(:,PC)
      END DO
    END DO

    CLOSE(99)

  END SUBROUTINE OUTPUTSNAPSHOT

  SUBROUTINE OUTPUTFPT(GROUPLIST,FILENAME,NPART,NTRIALS)
  ! Output first passage times
  ! ------------------
  ! input parameters:
  ! GROUPLIST: object describing multiple particle groups
  ! FILENAME: output file
  ! -------------
    IMPLICIT NONE

    TYPE(PARTICLEGROUP), TARGET :: GROUPLIST(NTRIALS)
    TYPE(PARTICLEGROUP), POINTER :: PGROUP
    INTEGER, INTENT(IN) :: NTRIALS,NPART
    CHARACTER(LEN=*), INTENT(IN) :: FILENAME

    INTEGER :: TRC
    CHARACTER(LEN=8) :: X1
    CHARACTER(LEN=19) :: FMT1

    OPEN(UNIT=99,FILE=FILENAME,STATUS='UNKNOWN')

    WRITE (X1,'(I8)') NPART
    FMT1 = '('//TRIM(X1)//'ES15.6E2)'

    WRITE(99,'(2I8)') NPART,NTRIALS
    DO TRC = 1,NTRIALS
      PGROUP=>GROUPLIST(TRC)

      WRITE(99,FMT1) PGROUP%ABSTIMES

    END DO

    CLOSE(99)
  
  END SUBROUTINE OUTPUTFPT

  SUBROUTINE OUTPUTBINDPOS(GROUPLIST,FILENAME,NPART,NTRIALS)
    ! Output first passage times
    ! ------------------
    ! input parameters:
    ! GROUPLIST: object describing multiple particle groups
    ! FILENAME: output file
    ! -------------

      IMPLICIT NONE
  
      TYPE(PARTICLEGROUP), TARGET :: GROUPLIST(NTRIALS)
      TYPE(PARTICLEGROUP), POINTER :: PGROUP
      INTEGER, INTENT(IN) :: NTRIALS,NPART
      CHARACTER(LEN=*), INTENT(IN) :: FILENAME
  
      INTEGER :: TRC,PC
      CHARACTER(LEN=8) :: X1
      CHARACTER(LEN=19) :: FMT1
  
      OPEN(UNIT=99,FILE=FILENAME,STATUS='UNKNOWN')
  
      WRITE (X1,'(I8)') NPART
      FMT1 = '('//TRIM(X1)//'ES15.6E2)'
  
      WRITE(99,'(2I8)') NPART,NTRIALS
      DO TRC = 1,NTRIALS
        PGROUP=>GROUPLIST(TRC)
        
        DO PC = 1,3
          WRITE(99,FMT1) PGROUP%ABSPOS(PC,:)
        END DO
  
      END DO
  
      CLOSE(99)
    
    END SUBROUTINE OUTPUTBINDPOS

  SUBROUTINE READMTLENS(MTLENS,MTLENFILE,NVALS,NMT)
  ! Read in microtubule lengths
    USE INPUTPARAMS, ONLY : READLINE, READA, READF, READI, READO

    IMPLICIT NONE 

    CHARACTER(LEN=*),INTENT(IN) :: MTLENFILE
    INTEGER,INTENT(IN) :: NVALS,NMT
    DOUBLE PRECISION,INTENT(OUT) :: MTLENS(NMT,NVALS)
    LOGICAL :: LDUM,FILEEND
    ! file id number
    INTEGER, PARAMETER :: MF = 51
    INTEGER :: NC,MC,NITEMS
    
    ! attempt to open file
    PRINT*, 'Reading MT lengths from: ', MTLENFILE
    INQUIRE(FILE=MTLENFILE,EXIST=LDUM)
    IF (.NOT.LDUM) THEN
       PRINT*, 'ERROR in READMTLEN: file ', TRIM(ADJUSTL(MTLENFILE)), ' does not exist.'
       STOP 1
    ENDIF
    OPEN(UNIT=MF, FILE=MTLENFILE, STATUS='OLD')

    DO NC = 1,NVALS
      CALL READLINE(MF,FILEEND,NITEMS)
      IF (FILEEND.and.NITEMS.eq.0) EXIT
      DO MC = 1,MIN(NMT,NITEMS)
        CALL READF(MTLENS(MC,NC))
      END DO
    END DO

    CLOSE(MF)

  END SUBROUTINE READMTLENS

  SUBROUTINE READMTNUC(MTNUC,MTNUCFILE,NVALS,NMT)
  ! Read in microtubule nucleation points
    USE INPUTPARAMS, ONLY : READLINE, READA, READF, READI, READO

    IMPLICIT NONE 

    CHARACTER(LEN=*),INTENT(IN) :: MTNUCFILE
    INTEGER,INTENT(IN) :: NVALS,NMT
    DOUBLE PRECISION,INTENT(OUT) :: MTNUC(NMT,NVALS)
    LOGICAL :: LDUM,FILEEND
    ! file id number
    INTEGER, PARAMETER :: MF = 51
    INTEGER :: NC,MC,NITEMS
    
    ! attempt to open file
    PRINT*, 'Reading MT nucleation points from: ', MTNUCFILE
    INQUIRE(FILE=MTNUCFILE,EXIST=LDUM)
    IF (.NOT.LDUM) THEN
       PRINT*, 'ERROR in READMTNUC: file ', TRIM(ADJUSTL(MTNUCFILE)), ' does not exist.'
       STOP 1
    ENDIF
    OPEN(UNIT=MF, FILE=MTNUCFILE, STATUS='OLD')

    DO NC = 1,NVALS
      CALL READLINE(MF,FILEEND,NITEMS)
      IF (FILEEND.and.NITEMS.eq.0) EXIT
      DO MC = 1,MIN(NMT,NITEMS)
        CALL READF(MTNUC(MC,NC))
      END DO
    END DO

    CLOSE(MF)

  END SUBROUTINE READMTNUC
 

  SUBROUTINE SETUPPARTICLEGROUP(PGROUP,NPART,NMT,PARAMP,MTLENS,RANDNUC,MTNUC)
  ! Initialize arrays for a group of particles and microtubules
    USE MT19937, ONLY: GRND
    
    IMPLICIT NONE

    DOUBLE PRECISION, PARAMETER :: PI = 3.14149265
    TYPE(PARTICLEGROUP), POINTER :: PGROUP
    TYPE(PARAMLIST), POINTER :: PARAMP
    INTEGER, INTENT(IN) :: NPART,NMT
    DOUBLE PRECISION,INTENT(IN) :: MTLENS(NMT),MTNUC(NMT)
    LOGICAL, INTENT(IN) :: RANDNUC
    INTEGER :: MC,PC
    DOUBLE PRECISION :: TH,R0,Z,Z0
    DOUBLE PRECISION :: P_GROW,P_PAUSE,P_SET,DENOM
    
    Z0 = 0D0

    PGROUP%NMT = NMT
    ALLOCATE(PGROUP%MTPOS(3,NMT))
    ALLOCATE(PGROUP%MTNUC(3,NMT))
    ALLOCATE(PGROUP%MTSTATE(NMT))

    DENOM = (PARAMP%K_CAT+PARAMP%K_PAUSE)/PARAMP%V_GROWTH+&
                PARAMP%K_CAT/PARAMP%V_SHRINK
    P_GROW = PARAMP%K_CAT/PARAMP%V_GROWTH/DENOM
    P_PAUSE = PARAMP%K_PAUSE/PARAMP%V_GROWTH/DENOM
    
    ! choose starting position for MTs
    IF(NMT.GT.0) THEN
      DO MC = 1,NMT
        IF(PARAMP%STARTATCENTER) THEN
          TH = 0D0
          R0 = 0D0
        ELSE
          TH = 2*PI*GRND()
          R0 = PARAMP%DOMRAD*SQRT(GRND())
        END IF
        
        IF(MTLENS(MC).GT.0D0) THEN
          Z = MTLENS(MC)
        ELSE
          Z = PARAMP%DOMLEN*GRND()
        END IF

        IF(RANDNUC) THEN
          IF(MTNUC(MC).GE.0D0) THEN
            Z0 = MTNUC(MC)
          ELSE 
            Z0 = GRND()*(PARAMP%DOMLEN-Z)
          END IF
        END IF

        PGROUP%MTPOS(:,MC) = (/R0*COS(TH),R0*SIN(TH),Z+Z0/)
        PGROUP%MTNUC(:,MC) = (/R0*COS(TH),R0*SIN(TH),Z0/)
        
        
        IF(PGROUP%MTPOS(3,MC).GE.PARAMP%DOMLEN) THEN
          PGROUP%MTSTATE(MC) = 3 ! mt starts in paused state
        ELSE
          P_SET = GRND()
          IF(P_SET.LE.P_GROW) THEN
            PGROUP%MTSTATE(MC) = 1 ! MT starts in growing state
          ELSEIF(P_SET.GT.P_GROW .AND. P_SET.LE.(P_GROW+P_PAUSE)) THEN
            PGROUP%MTSTATE(MC) = 3 ! MT starts in paused state
          ELSE
            PGROUP%MTSTATE(MC) = 2 ! MT starts in shrinking state
          END IF
        END IF

      END DO

    END IF

    PGROUP%NPART = NPART
    ALLOCATE(PGROUP%PARTPOS(3,NPART))
    ALLOCATE(PGROUP%MTINDS(NPART))
    ALLOCATE(PGROUP%CURDIR(NPART))
    SELECT CASE(PARAMP%STARTSTATE)
      CASE(1)
        DO PC = 1,NPART
          ! distribute particles over MTs
          PGROUP%MTINDS(PC) = 1+FLOOR(NMT*GRND())

          !set initial position for particles
          PGROUP%PARTPOS(1:2,PC) = PGROUP%MTPOS(1:2,PGROUP%MTINDS(PC))
          PGROUP%PARTPOS(3,PC) = PGROUP%MTNUC(3,PGROUP%MTINDS(PC))

          !set initial direction of motion
          PGROUP%CURDIR(PC) = 1
        END DO
      CASE(2)
        DO PC = 1,NPART
          ! particles start diffusive
          PGROUP%MTINDS(PC) = 0

          !set initial position for particles
          TH = 2*PI*GRND()
          R0 = PARAMP%DOMRAD*SQRT(GRND())
          PGROUP%PARTPOS(:,PC) = (/R0*COS(TH),R0*SIN(TH),PARAMP%DOMLEN/)

          !set initial direction of motion
          PGROUP%CURDIR(PC) = -1
        END DO
      CASE(3)
        DO PC = 1,NPART
          ! particles start diffusive
          PGROUP%MTINDS(PC) = 0

          !set initial position for particles
          TH = 2*PI*GRND()
          R0 = PARAMP%DOMRAD*SQRT(GRND())
          Z = PARAMP%DOMLEN*GRND()
          PGROUP%PARTPOS(:,PC) = (/R0*COS(TH),R0*SIN(TH),Z/)

          !set initial direction of motion
          PGROUP%CURDIR(PC) = -1
        END DO
    END SELECT
    ALLOCATE(PGROUP%ABSTIMES(NPART))
    PGROUP%ABSTIMES = -1D0
    ALLOCATE(PGROUP%ISACTIVE(NPART))
    PGROUP%ISACTIVE = .TRUE.
    ALLOCATE(PGROUP%ABSPOS(3,NPART))
    PGROUP%ABSPOS = -1D0

    PGROUP%ARRAYSET= .TRUE.
  END SUBROUTINE SETUPPARTICLEGROUP

  SUBROUTINE SETUPPARAMS(PARAMP)
  ! set up and initialize object to track particle and MT parameters
    USE KEYS, ONLY: ABSBOUND,CAPTUREINSTATE,CRAD,DELT,DIFFCONST,DOMLEN,DOMRAD,&
                    KBIND,KREV,KJUMP,K_CAT,K_PAUSE,K_RESCUE,NPART,NMT,&
                    STOPONCAPTURE,STAYONMT,SWITCHONJUMP,&
                    STARTATCENTER,TIPCAPTURE,V_GROWTH,V_SHRINK,VEL
    IMPLICIT NONE
    TYPE(PARAMLIST),POINTER :: PARAMP

    PARAMP%NPART = NPART

    PARAMP%NMT = NMT

    PARAMP%VEL = VEL

    PARAMP%DIFFCONST = DIFFCONST

    PARAMP%DOMLEN = DOMLEN

    PARAMP%DOMRAD = DOMRAD

    PARAMP%KBIND = KBIND

    PARAMP%K_CAT = K_CAT

    PARAMP%K_PAUSE = K_PAUSE

    PARAMP%KREV = KREV

    PARAMP%K_RESCUE = K_RESCUE

    PARAMP%CRAD = CRAD

    PARAMP%KJUMP = KJUMP

    PARAMP%TIPCAPTURE = TIPCAPTURE

    PARAMP%CAPTUREINSTATE = CAPTUREINSTATE

    PARAMP%DELT = DELT

    PARAMP%STARTSTATE = 0

    PARAMP%STOPONCAPTURE = STOPONCAPTURE

    PARAMP%STAYONMT = STAYONMT

    PARAMP%SWITCHONJUMP = SWITCHONJUMP

    PARAMP%STARTATCENTER = STARTATCENTER

    PARAMP%V_GROWTH = V_GROWTH

    PARAMP%V_SHRINK = V_SHRINK

    PARAMP%ABSBOUND = ABSBOUND

    PARAMP%ARRAYSET = .TRUE.

  END SUBROUTINE SETUPPARAMS


  SUBROUTINE CLEANUPPARTICLEGROUP(PGROUP)
  ! deallocate arrays
    IMPLICIT NONE
    TYPE(PARTICLEGROUP), POINTER :: PGROUP

    DEALLOCATE(PGROUP%MTPOS,PGROUP%MTSTATE)
    DEALLOCATE(PGROUP%PARTPOS,PGROUP%MTINDS,PGROUP%CURDIR)
    DEALLOCATE(PGROUP%ABSTIMES,PGROUP%ABSPOS,PGROUP%ISACTIVE)
    PGROUP%NMT = 0
    PGROUP%NPART = 0

    PGROUP%ARRAYSET = .FALSE.
  END SUBROUTINE CLEANUPPARTICLEGROUP

  SUBROUTINE CLEANUPPARAMS(PARAMP)
  ! deallocate arrays
    IMPLICIT NONE
    TYPE(PARAMLIST), POINTER:: PARAMP

    PARAMP%VEL = 0D0

    PARAMP%DIFFCONST = 0D0

    PARAMP%DOMLEN = 0D0

    PARAMP%DOMRAD = 0D0

    PARAMP%KBIND = 0D0

    PARAMP%K_CAT = 0D0

    PARAMP%K_PAUSE = 0D0

    PARAMP%K_RESCUE = 0D0

    PARAMP%KREV = 0D0

    PARAMP%CRAD = 0D0

    PARAMP%KJUMP = 0D0

    PARAMP%DELT = 0

    PARAMP%TIPCAPTURE = .FALSE.

    PARAMP%CAPTUREINSTATE = 1

    PARAMP%STARTATCENTER = .FALSE.

    PARAMP%ARRAYSET = .FALSE.

  END SUBROUTINE CLEANUPPARAMS

END MODULE PARTICLEUTIL
