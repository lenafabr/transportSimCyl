MODULE KEYS
  ! keyword parameters that are globally used in many different places in the code
  IMPLICIT NONE

  ! ----------------------------------------------------------------------------
  ! general program control
  ! ----------------------------------------------------------------------------
  CHARACTER*100 :: ACTION
  INTEGER :: RNGSEED
  LOGICAL :: VERBOSE
  
  ! ----------------------------------------------------------------------------
  ! simulation info
  ! ----------------------------------------------------------------------------
  INTEGER:: NTRIALS

  ! ----------------------------------------------------------------------------
  ! particle properties
  ! ----------------------------------------------------------------------------
  INTEGER :: NPART
  DOUBLE PRECISION :: VEL
  DOUBLE PRECISION :: KREV
  DOUBLE PRECISION :: KJUMP
  DOUBLE PRECISION :: CRAD
  DOUBLE PRECISION :: DIFFCONST
  DOUBLE PRECISION :: KBIND
  INTEGER :: STARTSTATE ! 1 = start on MT, 2 = start at domain end, 3 = start uniformly
  LOGICAL :: STOPONCAPTURE ! stop recording FPT when captured
  LOGICAL :: STAYONMT ! do particles stay on MT when reaching the endpoint?
  LOGICAL :: SWITCHONJUMP ! do particles randomize direction on jumping?

  ! ----------------------------------------------------------------------------
  ! microtubule properties
  ! ----------------------------------------------------------------------------
  INTEGER :: NMT
  LOGICAL :: MTDYN
  LOGICAL :: RESAMPLE
  DOUBLE PRECISION :: V_GROWTH,V_SHRINK
  DOUBLE PRECISION :: K_CAT,K_RESCUE,K_PAUSE
  LOGICAL :: TIPCAPTURE ! are particles only captured at MT +ends?
  
  ! what state does capture occur in
  ! 1 = any state
  ! 2 = growing and paused
  ! 3 = only paused
  INTEGER :: CAPTUREINSTATE 

  LOGICAL :: READLENS ! read microtubule lengths from file
  LOGICAL :: RANDNUC ! do microtubules nucleate randomly?
  LOGICAL :: READNUC ! read microtubule nucleation points from file
  DOUBLE PRECISION :: TOTMTLEN ! total available MT length
  LOGICAL :: STARTATCENTER ! start MTs at radial center

  ! ----------------------------------------------------------------------------
  ! domain info
  ! ----------------------------------------------------------------------------
  DOUBLE PRECISION :: DOMRAD
  DOUBLE PRECISION :: DOMLEN
  LOGICAL :: ABSBOUND
  
  ! ----------------------------------------------------------------------------
  ! input/output
  ! ----------------------------------------------------------------------------
  CHARACTER*100 :: SNAPFILE,FPTFILE,BINDPOSFILE,MTLENFILE,MTNUCFILE
  LOGICAL :: USEPARAMFILENAME
  LOGICAL :: DUMPSNAPSHOTS, APPENDSNAPSHOTS
  INTEGER :: SNAPEVERY
  
  ! ----------------------------------------------------------------------------
  ! dynamics
  ! ----------------------------------------------------------------------------
  LOGICAL :: DOBROWN
  DOUBLE PRECISION :: DELT
  INTEGER :: NSTEPS, PRINTEVERY
	
END MODULE KEYS