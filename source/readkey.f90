SUBROUTINE READKEY
! this subroutine reads in keywords from a parameter file
! it sets the various global variables defined in KEYS module
! name of the parameter file is param.* where * is a keyword argument
! if no keyword argument is supplied, the default is just a file called param
! The EXTRAPARAMFILES keyword will allow extra parameter files to be 
! read in as well

  USE KEYS
  USE INPUTPARAMS, ONLY : READLINE, READA, READF, READI, READO
  USE GENUTIL

  IMPLICIT NONE

  ! ---- stuff for inputing the parameter file in free format --------
  CHARACTER*100 :: ARG ! command line argument
  INTEGER :: NUMARG ! number of command line arguments
  INTEGER :: NITEMS ! number of items on the line in the parameter file
  INTEGER :: PF ! input file unit
  LOGICAL :: FILEEND=.FALSE. ! done reading file?
  CHARACTER*100 :: WORD ! keyword
  ! -------------- for reading multiple parameter files --------  
  INTEGER, PARAMETER :: MAXNFILES = 10
  CHARACTER*100 :: PARAMFILES(MAXNFILES)
  INTEGER :: NPARAMFILES, NPARAMREAD
  ! ------ for initializing random number generator
  INTEGER :: TIMEVAL(8), SEED
  ! ---------------- temporary variables ---------------
  INTEGER :: I
  CHARACTER*100 :: DUMSTR
  CHARACTER*2 :: X1
  LOGICAL :: LDUM
  DOUBLE PRECISION :: TMP

  ! ------------------------
  ! set variable defaults
  ! ------------------------

  ! general program control
  ACTION = 'NONE'
  RNGSEED = 0
  VERBOSE = .FALSE.

  ! simulation info
  NTRIALS = 1 ! number of trials

  ! particle properties
  NPART = 1 ! number of particles
  VEL = 1D0 ! velocity of particles
  KREV = 1D0 ! reversal rate (per second)
  KJUMP = 1D0 ! jump rate (per second)
  KBIND = 1D0 ! bind rate (per second)
  CRAD = 0.25D0 ! capture radius (microns)
  DIFFCONST = 0.015 ! diffusion coefficient
  STARTSTATE  = 1 ! starting state of the particle, 1: start on MT, 2: start at end,
                  ! 3: start uniformly in the domain
  STOPONCAPTURE = .TRUE. ! stop recording FPT when particle is captured by MT
  STAYONMT = .TRUE. ! do particles stay on MT after reaching the tip?
  SWITCHONJUMP = .FALSE. ! do particles switch directions after jumping to another MT?
  
  ! domain info
  DOMRAD = 1D0 ! domain radius (micron)
  DOMLEN = 20D0 ! domain length (micron)
  ABSBOUND = .FALSE. ! absorbing boundary at cell body?

  ! microtubule properties
  NMT = 5             ! number of microtubules
  MTDYN = .FALSE.     ! implement MT dynamics
  RESAMPLE = .FALSE.  ! resample MTs instead of continuous dynamics
  V_GROWTH = 1D0      ! growth velocity
  V_SHRINK = 1D0      ! shrink velocity
  K_CAT = 1D0         ! catastrophe rate
  K_RESCUE  = 1D0     ! rescue rate
  K_PAUSE = 1D0       ! pause rate
  TIPCAPTURE = .FALSE.! particles can be captured anywhere
  CAPTUREINSTATE = 1  ! particles are captured in any MT state
  READLENS = .FALSE.  ! do not read MT lengths from file
  RANDNUC = .FALSE.   ! microtubules nucleate at 0
  READNUC = .FALSE.   ! do not read MT nucleation points from file
  TOTMTLEN = 0D0      ! total available MT length
  STARTATCENTER = .FALSE. ! start MTs at radial center of domain
  
  ! input/output  
  SNAPFILE = '*.snap.out' ! snapshot file
  DUMPSNAPSHOTS = .FALSE. ! periodically dump chain snapshots
  APPENDSNAPSHOTS = .FALSE. ! append snapshots to file rather than replacing
  SNAPEVERY = 0 ! how often to dump snapshots
  USEPARAMFILENAME = .TRUE. ! use param file name to replace substrings while reading files
  FPTFILE = '*.fpt.out'
  BINDPOSFILE = '*.bindpos.out'
  MTLENFILE = '*.txt'
  MTNUCFILE = '*_nuc.txt'
  
  ! dynamics
  DELT = 1D-4 ! time step
  NSTEPS = INT(1D6) ! number of brownian steps to run
  PRINTEVERY = 0 ! how often to print output
  DOBROWN = .TRUE.
    
  ! -------------------------
  ! Read in all parameter files, starting with the ones specified on command line
  ! --------------------------

  PF = 55 ! i/o unit number to be used for parameter files

  ! get input parameter files from command line
  NPARAMFILES = 0
  NUMARG = COMMAND_ARGUMENT_COUNT()  
  IF (NUMARG==0) THEN
    NPARAMFILES = 1
    PARAMFILES(1) = 'param'
    ARG = ''
  ELSE
    DO I = 1,NUMARG
      CALL GETARG(I, ARG)
      NPARAMFILES = NPARAMFILES + 1
      WRITE(DUMSTR,'(A)') 'param.' //TRIM(ADJUSTL(ARG))
      PARAMFILES(NPARAMFILES) = DUMSTR
    END DO
    ! reset arg to its original value
    IF (NUMARG.GT.1) CALL GETARG(1,ARG)
  ENDIF

  NPARAMREAD = 0 ! keep track of how many files have been read
  DO WHILE (NPARAMREAD.LT.NPARAMFILES)
    NPARAMREAD = NPARAMREAD + 1

    PRINT*, 'Reading parameter file: ', PARAMFILES(NPARAMREAD)
    INQUIRE(FILE=PARAMFILES(NPARAMREAD),EXIST=LDUM)
    IF (.NOT.LDUM) THEN
      PRINT*, 'ERROR in READKEY: Parameter file ', TRIM(ADJUSTL(PARAMFILES(NPARAMREAD))), ' does not exist.'
      STOP 1
    ENDIF
    OPEN(UNIT=PF, FILE=PARAMFILES(NPARAMREAD), STATUS='OLD')

    ! read in the keywords one line at a time
    DO 
      CALL READLINE(PF,FILEEND,NITEMS)
      IF (FILEEND.and.nitems.eq.0) EXIT

      ! skip empty lines
      IF (NITEMS.EQ.0) CYCLE

      ! Read in the keyword for this line
      CALL READA(WORD,CASESET=1)

      ! Skip any empty lines or any comment lines
      IF (WORD(1:1).EQ.'#') CYCLE

      SELECT CASE(WORD) ! pick which keyword
        CASE('ACTION')
          CALL READA(ACTION, CASESET=1)
        CASE('ABSBOUND')
          CALL READO(ABSBOUND)
        CASE('BINDPOS')
          CALL READA(BINDPOSFILE)
        CASE('CAPTUREINSTATE')
          CALL READI(CAPTUREINSTATE)
        CASE('CRAD')
          CALL READF(CRAD)
        CASE('DIFFCONST')
          CALL READF(DIFFCONST)
        CASE('DOBROWN')
          CALL READO(DOBROWN)
        CASE('DOMRAD')
          CALL READF(DOMRAD)
        CASE('DOMLEN')
          CALL READF(DOMLEN)
        CASE('DELT')
          CALL READF(DELT)
        CASE('FPT')
          CALL READA(FPTFILE)
        CASE('KBIND')
          CALL READF(KBIND)
        CASE('K_CAT')
          CALL READF(K_CAT)
        CASE('KJUMP')
          CALL READF(KJUMP)
        CASE('K_PAUSE')
          CALL READF(K_PAUSE)
        CASE('K_RESCUE')
          CALL READF(K_RESCUE)
        CASE('KREV')
          CALL READF(KREV)
        CASE('MTDYN')
          CALL READO(MTDYN)
        CASE('NMT')
          CALL READF(TMP)
          NMT = INT(TMP)
        CASE('NPART')
          CALL READI(NPART)
        CASE('NSTEPS')
          CALL READI(NSTEPS)
        CASE('NTRIALS')
          CALL READI(NTRIALS)
        CASE('PRINTEVERY')
          CALL READI(PRINTEVERY)
          IF(PRINTEVERY.LT.0D0) PRINTEVERY = NSTEPS/10 
        CASE('MTLENFILE') 
          CALL READA(MTLENFILE)
        CASE('MTNUCFILE')
          CALL READA(MTNUCFILE)
        CASE('RANDNUC')
          CALL READO(RANDNUC)
        CASE('READNUC')
          CALL READO(READNUC)
        CASE('READLENS')
          CALL READO(READLENS)
        CASE('RESAMPLE')
          CALL READO(RESAMPLE)
        CASE('RNGSEED')
          CALL READI(RNGSEED)
        CASE('SNAPFILE')
          CALL READA(SNAPFILE)  
        CASE('SNAPSHOTS')
          IF(NITEMS.GT.1) CALL READI(SNAPEVERY)
          IF(SNAPEVERY.NE.0) DUMPSNAPSHOTS = .TRUE.
          IF(SNAPEVERY.LT.0) SNAPEVERY = NSTEPS
          IF(NITEMS.GT.2) CALL READA(SNAPFILE)
          IF(NITEMS.GT.3) CALL READO(APPENDSNAPSHOTS)
        CASE('STARTATCENTER')
          CALL READO(STARTATCENTER)
        CASE('STARTSTATE')
          CALL READI(STARTSTATE)
        CASE('STAYONMT')
          CALL READO(STAYONMT)
        CASE('STOPONCAPTURE')
          CALL READO(STOPONCAPTURE)
        CASE('SWITCHONJUMP')
          CALL READO(SWITCHONJUMP)
        CASE('TIPCAPTURE')
          CALL READO(TIPCAPTURE)
        CASE('TOTMTLEN')
          CALL READF(TOTMTLEN)
        CASE('USEPARAMFILENAME')
          CALL READO(USEPARAMFILENAME)
        CASE('V_GROWTH')
          CALL READF(V_GROWTH)
        CASE('V_SHRINK')
          CALL READF(V_SHRINK)
        CASE('VEL')
          CALL READF(VEL)
        CASE('VERBOSE')
          CALL READO(VERBOSE)
        CASE DEFAULT
          print*, 'ERROR: unidentified keyword ', TRIM(WORD), " Will ignore."
      END SELECT
    ENDDO
    CLOSE(PF)
  ENDDO

  ! -----------------
  ! check validity of some values, raise errors or adjust as necessary
  ! -----------------  
  IF(RANDNUC.AND..NOT.STOPONCAPTURE.AND.STAYONMT) THEN
    PRINT*, 'EXIT CONDITION MAY NOT BE REACHED, CHECK PARAM FILE'
  END IF

  ! ----------- fix file names -----------
  CALL REPLACESUBSTR(SNAPFILE,'*',TRIM(ADJUSTL(ARG)))
  CALL REPLACESUBSTR(FPTFILE,'*',TRIM(ADJUSTL(ARG)))
  CALL REPLACESUBSTR(BINDPOSFILE,'*',TRIM(ADJUSTL(ARG)))
  IF(USEPARAMFILENAME) THEN
    CALL REPLACESUBSTR(MTLENFILE,'*',TRIM(ADJUSTL(ARG)))
    CALL REPLACESUBSTR(MTNUCFILE,'*',TRIM(ADJUSTL(ARG)))
  ELSE
    WRITE(X1,'(I2)') NMT
    CALL REPLACESUBSTR(MTLENFILE,'*',TRIM(ADJUSTL(X1)))
    CALL REPLACESUBSTR(MTNUCFILE,'*',TRIM(ADJUSTL(X1)))
  END IF
  ! ---------------------------

  ! Initiate random number generator 
  IF (RNGSEED.EQ.0) THEN
    ! use the current time of day in milliseconds
    CALL DATE_AND_TIME(VALUES=TIMEVAL)
    SEED = TIMEVAL(5)*3600*1000 + TIMEVAL(6)*60*1000 + TIMEVAL(7)*1000 + TIMEVAL(8)
  ELSEIF (RNGSEED.EQ.-1) THEN
    ! use the last 5 characters in the command-line argument
    SEED = STRING2NUM(TRIM(ADJUSTL(ARG)))    
  ELSEIF (RNGSEED.EQ.-2) THEN
    ! use the last 4 characters in the command-line argument 
    ! and additionally the millisecond time 
    CALL DATE_AND_TIME(VALUES=TIMEVAL)
    SEED = STRING2NUM(TRIM(ADJUSTL(ARG)),TIMEVAL(8))
  ELSE
    ! use this seed directly
    SEED = RNGSEED
  ENDIF

  print*, 'Initiating Mersenne twister random number generator with seed:', SEED
  CALL SGRND(SEED)

  print*, '------------Parameter values : -------------------'
  print*, 'ACTION: ', TRIM(ADJUSTL(ACTION))
  IF (DUMPSNAPSHOTS) THEN
    PRINT*, 'Dumping snapshot every', SNAPEVERY,'steps. In file:', TRIM(ADJUSTL(SNAPFILE))
  ENDIF
  
  print*, 'number of particles:', NPART
  print*, 'number of microtubules:', NMT
      
  IF(DOMRAD.GT.0D0) PRINT('(A15,F12.4)'), 'domain radius: ',DOMRAD
  
  IF(DOMLEN.GT.0D0) PRINT('(A15,F12.4)'), 'domain length: ',DOMLEN

  IF(MTDYN) PRINT('(A24)'), 'implementing MT dynamics'
  
  PRINT*, '----------------------------------------------------'
	
END SUBROUTINE READKEY
