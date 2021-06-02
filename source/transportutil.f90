MODULE TRANSPORTUTIL
! subroutines to simulate transport on dynamic microtubules in a cylinder
! VERSION BY: SAURABH MOGRE 
! DATE: APR 2020

USE MT19937, ONLY: GRND,RNORM,GRND
  
IMPLICIT NONE

CONTAINS

  SUBROUTINE TRANSPORTSIM(GROUPLIST,PARAMP,NTRIALS)

	!------------------------------------------------------
  !simulates particles undergoing multi-modal transport in a 3d cylinder          
  !particles move on microtubules distributed uniformly within the cross section
  !length of microtubules is uniformly distributed                    
	!NPART - NUMBER OF PARTICLES
	!NMT - NUMBER OF MICROTUBULES
	!NSTEP - NUMBER OF STEPS
	!DELT - TIMESTEP
  !V - PARTICLE VELOCITY ON A MICROTUBULE
  !KR - RATE OF REVERSAL OF A MOVING PARTICLE
	!KJ - RATE OF JUMPING TO ANOTHER MICROTUBULE WITHIN REACH
	!CRAD - RADIUS WITHIN WHICH A MICROTUBULE CAN CAPTURE PARTICLES
  !DOMRAD - RADIUS OF THE DOMAIN
  !DOMLEN - LENGTH OF THE DOMAIN
  !PRINTEVERY - HOW OFTEN STATUS IS DISPLAYED WHILE RUNNING SIMS
  !PARTFILE - FILE TO WHICH PARTICLE POSITIONS ARE WRITTEN
  !SNAPSHOTEVERY - FREQUENCY AT WHICH OUTPUTS ARE WRITTEN TO FILE
  !------------------------------------------------------
    USE KEYS, ONLY : NSTEPS,DELT,MTDYN,RESAMPLE,DOBROWN,PRINTEVERY,SNAPFILE,&
                      SNAPEVERY,FPTFILE,BINDPOSFILE,MTLENFILE
    USE PARTICLEUTIL, ONLY: PARTICLEGROUP,OUTPUTSNAPSHOT,PARAMLIST,OUTPUTFPT,&
                            OUTPUTBINDPOS,READMTLENS
    USE MT19937, ONLY: RNORM,GRND
    USE GENUTIL,ONLY: RANDARRAY

    INTEGER, INTENT(IN) :: NTRIALS
    TYPE(PARTICLEGROUP),TARGET :: GROUPLIST(NTRIALS)
    TYPE(PARTICLEGROUP),POINTER :: PGROUP
    TYPE(PARAMLIST),POINTER :: PARAMP
    DOUBLE PRECISION,PARAMETER :: PI = 3.14159265
  	
    DOUBLE PRECISION,ALLOCATABLE :: INFO(:)
    INTEGER :: TC,PC,MC,NJUMP,TRC
    
    INTEGER :: NEWIND, DIRVEC(2)
    INTEGER :: REACHEDMT(NTRIALS),REACHEDSOMA(NTRIALS),TOTREACHEDMT,TOTREACHEDSOMA
    INTEGER :: ACTIVEPART(NTRIALS),TOTACTIVEPART
    INTEGER,ALLOCATABLE :: MTJUMPS(:)
    DOUBLE PRECISION :: VMT(3), NEWPOS(3), TMPPOS(3)
    LOGICAL :: APPENDTOFILE,JUMPFLAG,CANBIND
    DOUBLE PRECISION :: TIPDIST,NUCDIST,POS,PBIND
    DOUBLE PRECISION :: NX1,NX2,NX3,R,NXNORM
    DOUBLE PRECISION :: P_PAUSE,P_RESCUE,P_CAT
    DOUBLE PRECISION, ALLOCATABLE :: TMPARRAY(:)
    DOUBLE PRECISION, ALLOCATABLE :: MTLENLIST(:,:)
    
    IF(PRINTEVERY.NE.0) THEN
      PRINT'(A5,I8,A4,I8)', 'Step ',0,' of ',NSTEPS
    END IF

    ! initial snapshot
    APPENDTOFILE = .FALSE.
    IF(SNAPEVERY.GT.0) THEN
      ALLOCATE(INFO(9))
      INFO = (/0D0,DFLOAT(NSTEPS/SNAPEVERY),DELT*SNAPEVERY,DFLOAT(NTRIALS),&
              DFLOAT(PARAMP%NPART),DFLOAT(PARAMP%NMT),PARAMP%DOMRAD,&
              PARAMP%DOMLEN,PARAMP%CRAD/)
      CALL OUTPUTSNAPSHOT(GROUPLIST,SNAPFILE,INFO,NTRIALS,APPENDTOFILE)
    END IF

    ! initialize capture arrays
    REACHEDMT = 0
    REACHEDSOMA = 0
    TOTREACHEDMT = 0
    TOTREACHEDSOMA = 0
    ACTIVEPART = PARAMP%NPART
    TOTACTIVEPART = PARAMP%NPART*NTRIALS

    ! directions available
    DIRVEC = (/-1,1/)

    ! array of microtubule velocities
    VMT = (/PARAMP%V_GROWTH,-PARAMP%V_SHRINK,0D0/)

    ALLOCATE(TMPARRAY(GROUPLIST(1)%NMT))

    ! get pool of configurations to sample from
    IF(MTDYN.AND.RESAMPLE) THEN
      ALLOCATE(MTLENLIST(GROUPLIST(1)%NMT,INT(1D6)))
      CALL READMTLENS(MTLENLIST,MTLENFILE,INT(1D6),GROUPLIST(1)%NMT)
    END IF

    ! run simulation steps
    DO TC = 1,NSTEPS
      IF((PRINTEVERY.NE.0).AND.MOD(TC,PRINTEVERY).EQ.0) THEN
      ! print information on particle status
        PRINT'(A5,I8,A4,I8)', 'Step ',TC,' of ',NSTEPS
        IF(PARAMP%STOPONCAPTURE) THEN
          PRINT'(I8,A4,I8,A41)', TOTREACHEDMT,' of ',NTRIALS*PARAMP%NPART,&
                                                  ' particles have attached to a microtubule'
        ELSE
          PRINT'(I8,A4,I8,A32)', TOTREACHEDSOMA,' of ',NTRIALS*PARAMP%NPART,&
                                                  ' particles have reached the soma'
        END IF
        PRINT'(I8,A27)',TOTACTIVEPART, ' active particles remaining'
      END IF

      ! stop if all particles have reached the target
      IF(TOTACTIVEPART.EQ.0) THEN
        PRINT*, 'NO ACTIVE PARTICLE REMAINING'
        EXIT
      END IF
      
      ! pre-calculate probabilities
            
      ! probability to pause while growing
      P_PAUSE = 1-EXP(-PARAMP%K_PAUSE*DELT)
      
      ! probability to rescue while shrinking
      P_RESCUE = 1-EXP(-PARAMP%K_RESCUE*DELT)

      ! probability to catastrophe while paused
      P_CAT = 1-EXP(-(PARAMP%K_CAT)*DELT)

      ! go through all trials
      DO TRC = 1,NTRIALS
        PGROUP=>GROUPLIST(TRC)    

        IF(ACTIVEPART(TRC).EQ.0) CYCLE

        ALLOCATE(MTJUMPS(PGROUP%NMT))

        ! microtubule dynamics 
        IF(MTDYN) THEN
          IF(RESAMPLE) THEN
          ! resample MT distribution on a given time scale
            IF(GRND().LE.P_CAT) THEN
              PGROUP%MTPOS(3,:) = MTLENLIST(:,FLOOR(GRND()*SIZE(MTLENLIST,2))+1)
            END IF
          ELSE
          ! explicitly run dynamics
            DO MC = 1,PGROUP%NMT
              SELECT CASE (PGROUP%MTSTATE(MC))
                CASE(1) ! growing MT
                  ! check if change of state occurs within next timestep
                  IF(GRND().LE.P_CAT) THEN
                    PGROUP%MTSTATE(MC) = 2
                  END IF

                CASE(2) !shrinking MT
                  ! check if change of state occurs within next timestep
                  IF(GRND().LE.P_RESCUE) THEN
                    PGROUP%MTSTATE(MC) = 1
                  END IF
                CASE(3) !paused MT
                  ! check if change of state occurs within next timestep
                  IF(GRND().LE.P_CAT) THEN
                    PGROUP%MTSTATE(MC) = 2
                  END IF
                CASE DEFAULT
              END SELECT

              POS = PGROUP%MTPOS(3,MC)+VMT(PGROUP%MTSTATE(MC))*DELT

              IF(POS.GT.PARAMP%DOMLEN) THEN
                PGROUP%MTPOS(3,MC) = PARAMP%DOMLEN
                PGROUP%MTSTATE(MC) = 3
              ELSEIF(POS.LT.0D0) THEN
                PGROUP%MTPOS(3,MC) = 0D0
                PGROUP%MTSTATE(MC) = 1
              ELSE
                PGROUP%MTPOS(3,MC) = POS
              END IF
              
            END DO
          END IF
        END IF

        ! particle motion
        DO PC = 1,PGROUP%NPART
          IF(.NOT.PGROUP%ISACTIVE(PC)) CYCLE
          ! check if particle is on a microtubule
          IF(PGROUP%MTINDS(PC).NE.0) THEN
            
            ! particle is on a microtubule, count number of possible jumps
            NJUMP = 0
            IF(PGROUP%NMT.GT.1.AND.PARAMP%KJUMP.GT.0D0) THEN
              MTJUMPS = 0
              ! check which jumps are possible
              DO MC = 1,PGROUP%NMT
                IF(MC.EQ.PGROUP%MTINDS(PC)) CYCLE
                IF((((PGROUP%PARTPOS(1,PC)-PGROUP%MTPOS(1,MC))**2+&
                  (PGROUP%PARTPOS(2,PC)-PGROUP%MTPOS(2,MC))**2).LE.PARAMP%CRAD**2)&
                  .AND.(PGROUP%PARTPOS(3,PC).LE.PGROUP%MTPOS(3,MC))&
                  .AND.(PGROUP%PARTPOS(3,PC).GE.PGROUP%MTNUC(3,MC))) THEN
                  MTJUMPS(MC) = 1
                  NJUMP = NJUMP+1
                END IF
              END DO
            END IF

            ! check if reversal or jump event occurs
            IF(PARAMP%KREV.GT.0D0.AND.PARAMP%KJUMP.GT.0D0) THEN
              IF(GRND().LE.(1-EXP(-(PARAMP%KREV+NJUMP*PARAMP%KJUMP)*DELT))) THEN
                ! determine whether event is reversal or jump
                IF(GRND().LE.(PARAMP%KREV/(PARAMP%KREV+NJUMP*PARAMP%KJUMP))) THEN
                  !reverse direction
                  PGROUP%CURDIR(PC) = -PGROUP%CURDIR(PC)
                ELSE
                  !choose microtubule to jump to
                  NEWIND = MAXLOC(MTJUMPS*RANDARRAY(PGROUP%NMT),1)
                  !update xy position and mtinds                  
                  PGROUP%MTINDS(PC) = NEWIND
                  !choose direction
                  IF(PARAMP%SWITCHONJUMP) PGROUP%CURDIR(PC) = DIRVEC(FLOOR(1+2*GRND()))
                END IF
              END IF
            END IF

            ! attempt to move particle
            POS = PGROUP%PARTPOS(3,PC)+PGROUP%CURDIR(PC)*PARAMP%VEL*DELT
            TMPPOS = (/PGROUP%PARTPOS(1,PC),PGROUP%PARTPOS(2,PC),POS/)
            ! check if particle falls off MT  
            IF(POS.LE.PGROUP%MTPOS(3,PGROUP%MTINDS(PC)).AND.&
                                            POS.GE.PGROUP%MTNUC(3,PGROUP%MTINDS(PC))) THEN
              PGROUP%PARTPOS(3,PC) = POS
            ELSE
              TIPDIST = SQRT(SUM((TMPPOS-PGROUP%MTPOS(:,PGROUP%MTINDS(PC)))**2))
              NUCDIST = SQRT(SUM((TMPPOS-PGROUP%MTNUC(:,PGROUP%MTINDS(PC)))**2))
              IF(TIPDIST.LE.PARAMP%CRAD.OR.NUCDIST.LE.PARAMP%CRAD) THEN
                PGROUP%PARTPOS(3,PC) = POS
              ELSE
                IF(PARAMP%STAYONMT) THEN
                  ! equilibrate position around MT tip
                  NX1 = RNORM()
                  NX2 = RNORM()
                  NX3 = RNORM()
                  R = PARAMP%CRAD*GRND()**(1./3.)
                  NXNORM = SQRT(NX1**2+NX2**2+NX3**2)
                  PGROUP%PARTPOS(1,PC) = R*NX1/NXNORM
                  PGROUP%PARTPOS(2,PC) = R*NX2/NXNORM
                  PGROUP%PARTPOS(3,PC) = R*NX3/NXNORM
                ELSE
                  PGROUP%MTINDS(PC) = 0
                END IF

              END IF

            END IF
            
          END IF
          
          ! if particle is not on a microtubule, execute Brownian motion
          IF(DOBROWN) THEN
            IF(PGROUP%MTINDS(PC).EQ.0) THEN
              CANBIND = .TRUE.
              ! Brownian step
              CALL BDPROPAGATE(PGROUP%PARTPOS(:,PC),PARAMP,DELT,NEWPOS,CANBIND)
              PGROUP%PARTPOS(:,PC) = NEWPOS

              NJUMP = 0
              MTJUMPS = 0
              
              ! check if the particle can bind to a microtubule
              IF(CANBIND) THEN
                ! look for microtubules to bind
                DO MC = 1,PGROUP%NMT

                  ! does capture only occur in specific states?
                  SELECT CASE (PARAMP%CAPTUREINSTATE)
                    CASE(1)
                      CONTINUE
                    CASE(2)
                      IF(PGROUP%MTSTATE(MC).EQ.2) CYCLE !do not bind to shrinking MTs
                    CASE(3)
                      IF(PGROUP%MTSTATE(MC).NE.3) CYCLE !only bind to paused MTs
                    CASE DEFAULT
                      CONTINUE
                  END SELECT
                  
                  IF((((PGROUP%PARTPOS(1,PC)-PGROUP%MTPOS(1,MC))**2+&
                    (PGROUP%PARTPOS(2,PC)-PGROUP%MTPOS(2,MC))**2).LE.PARAMP%CRAD**2)) THEN
                    TIPDIST = SQRT(SUM((PGROUP%PARTPOS(:,PC)-PGROUP%MTPOS(:,MC))**2))
                    NUCDIST = SQRT(SUM((PGROUP%PARTPOS(:,PC)-PGROUP%MTNUC(:,MC))**2))
                    IF(PARAMP%TIPCAPTURE) THEN
                      ! only capture at plus ends
                      JUMPFLAG = (TIPDIST.LE.PARAMP%CRAD)
                    ELSE
                      JUMPFLAG = .TRUE.
                      IF(PGROUP%PARTPOS(3,PC).GT.PGROUP%MTPOS(3,MC)) THEN
                        JUMPFLAG = (TIPDIST.LE.PARAMP%CRAD)
                      ELSEIF(PGROUP%PARTPOS(3,PC).LT.PGROUP%MTNUC(3,MC)) THEN
                        JUMPFLAG = (NUCDIST.LE.PARAMP%CRAD)
                      END IF
                    END IF
                    IF(JUMPFLAG) THEN
                      MTJUMPS(MC) = 1
                      NJUMP = NJUMP + 1
                    END IF
                  END IF
                END DO

                ! check if binding event occurs
                IF(NJUMP.GT.0) THEN
                  IF(NJUMP*PARAMP%KBIND*DELT.LT.8) THEN
                    PBIND = 1-EXP(-NJUMP*PARAMP%KBIND*DELT)
                  ELSE
                    PBIND = 1
                  END IF
                  IF(GRND().LE.PBIND) THEN
                    ! randomly choose MT to bind
                    TMPARRAY = RANDARRAY(PGROUP%NMT)
                    NEWIND = MAXLOC(MTJUMPS*TMPARRAY,1)
                    PGROUP%MTINDS(PC) = NEWIND
  
                    ! switch movement direction if specified
                    IF(PARAMP%SWITCHONJUMP) THEN
                      PGROUP%CURDIR(PC) = DIRVEC(FLOOR(1+2*GRND()))
                    ELSE
                      PGROUP%CURDIR(PC) = -1 ! particle moves towards the cell body
                    END IF

                    ! record capture at MT tip
                    IF(PARAMP%STOPONCAPTURE) THEN
                      REACHEDMT(TRC) = REACHEDMT(TRC)+1
                      TOTREACHEDMT = TOTREACHEDMT+1
                      PGROUP%ABSTIMES(PC) = TC*DELT
                      PGROUP%ABSPOS(:,PC) = PGROUP%PARTPOS(:,PC)
                      PGROUP%ISACTIVE(PC) = .FALSE.
                      ACTIVEPART(TRC) = ACTIVEPART(TRC)-1
                      TOTACTIVEPART = TOTACTIVEPART-1
                    END IF

                  END IF
                END IF
              END IF  

            END IF
          END IF
          
          ! check if particle has reached cell body
          IF(PGROUP%PARTPOS(3,PC).LE.0D0) THEN
            IF(PARAMP%ABSBOUND.OR.(.NOT.PARAMP%STOPONCAPTURE)) THEN
              ! record capture at soma
              REACHEDSOMA(TRC) = REACHEDSOMA(TRC)+1
              TOTREACHEDSOMA = TOTREACHEDSOMA+1
              PGROUP%ABSTIMES(PC) = TC*DELT
              PGROUP%ABSPOS(:,PC) = PGROUP%PARTPOS(:,PC)
              PGROUP%ISACTIVE(PC) = .FALSE.
              ACTIVEPART(TRC) = ACTIVEPART(TRC)-1
              TOTACTIVEPART = TOTACTIVEPART-1
            END IF
          END IF

        END DO
        DEALLOCATE(MTJUMPS)
      END DO

      ! output snapshot to file
      IF(SNAPEVERY.GT.0) THEN
        IF(MOD(TC,SNAPEVERY).EQ.0) THEN
          INFO(1) = TC
          CALL OUTPUTSNAPSHOT(GROUPLIST,SNAPFILE,INFO,NTRIALS,.TRUE.)
        END IF
      END IF

    END DO

    ! output fpt distribution
    CALL OUTPUTFPT(GROUPLIST,FPTFILE,PARAMP%NPART,NTRIALS)

    ! output binding position
    CALL OUTPUTBINDPOS(GROUPLIST,BINDPOSFILE,PARAMP%NPART,NTRIALS)

  END SUBROUTINE TRANSPORTSIM

  SUBROUTINE BDPROPAGATE(PREVPOS,PARAMP,DELT,NEWPOS,CANBIND)
  ! do a brownian dynamics step to propagate the particle over a time interval
  ! R0, Z0 are initial radial and axial positions
  ! DDT = D*dt (diffusion coefficient times time interval)
  ! RNEW, ZNEQ are final radial and axial positions
    USE PARTICLEUTIL,ONLY: PARAMLIST
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: PREVPOS(3)
    TYPE(PARAMLIST), POINTER, INTENT(IN) :: PARAMP
    DOUBLE PRECISION, INTENT(IN) :: DELT
    DOUBLE PRECISION, INTENT(OUT) :: NEWPOS(3)
    LOGICAL, INTENT(OUT) :: CANBIND
    DOUBLE PRECISION :: STEPXYZ(3), COEF
    INTEGER :: I
    
    CANBIND = .TRUE.

    COEF = SQRT(2*PARAMP%DIFFCONST*DELT);
 
    DO I = 1,3
       STEPXYZ(I) = COEF*RNORM()
    END DO
    
    NEWPOS = PREVPOS+STEPXYZ
    
    IF((NEWPOS(3).GT.PARAMP%DOMLEN).OR.(NEWPOS(3).LT.0D0)&
                                .OR.(NEWPOS(1)**2+NEWPOS(2)**2).GT.PARAMP%DOMRAD**2) THEN
      IF(NEWPOS(3).LT.0D0.AND.PARAMP%ABSBOUND) THEN
        CANBIND = .FALSE.
        RETURN
      ELSE
        CALL REFLECTPART(PREVPOS,PARAMP%DOMRAD,PARAMP%DOMLEN,NEWPOS)                            
      END IF
    END IF
    
  END SUBROUTINE BDPROPAGATE

  SUBROUTINE REFLECTPART(POS0,DOMRAD,DOMLEN,POS1)
  ! subroutine to reflect particles off the boundary
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: POS0(3),DOMRAD,DOMLEN
    DOUBLE PRECISION, INTENT(OUT) :: POS1(3)
    DOUBLE PRECISION :: X1,Y1,X2,Y2,X0,Y0,A,B,C,C1,M,DISC
      
    DO WHILE ((POS1(3).GT.DOMLEN).OR.(POS1(3).LT.0D0))
      IF(POS1(3).GT.DOMLEN) POS1(3) = 2*DOMLEN-POS1(3)
      IF(POS1(3).LT.0D0) POS1(3) = -POS1(3)
    END DO
    
        
    DO WHILE (SQRT((POS1(1)**2+POS1(2)**2)).GT.DOMRAD)
      
      !find co-ordinates of the intersection point
      X1 = POS0(1)
      Y1 = POS0(2)
      X2 = POS1(1)
      Y2 = POS1(2)
      IF(ABS(X2-X1).LE.1D-6) THEN
        ! do something
        X0 = X1
        Y0 = SQRT(DOMRAD**2-X1**2)
      ELSE
        M = (Y2-Y1)/(X2-X1)
        C1 = Y1-M*X1
          
        A = 1+M**2
        B = 2*M*C1
        C = C1**2-DOMRAD**2
        
        IF(B**2-4*A*C .LT.0) THEN
          PRINT*, 'TESTXA:', B**2-4*A*C 
          PRINT*, 'TESTXA:',POS0,POS1
          PRINT*, 'TESTXA:',SQRT((POS0(1)**2+POS0(2)**2)),SQRT((POS1(1)**2+POS1(2)**2))
          STOP 1
        END IF
        DISC = SQRT(B**2-4*A*C)
        
          
        IF(((-B-DISC)/(2*A).GT.MIN(X1,X2)).AND.((-B-DISC)/(2*A).LT.MAX(X1,X2))) THEN
          X0 = (-B-DISC)/(2*A)
        ELSE
          X0 = (-B+DISC)/(2*A)
        END IF
          
        Y0 = M*X0+C1
      END IF

      !get reflection off the tangent at intersection point
      POS1(1) = X2-2*X0*(X0*X2+Y0*Y2-DOMRAD**2)/DOMRAD**2
      POS1(2) = Y2-2*Y0*(X0*X2+Y0*Y2-DOMRAD**2)/DOMRAD**2
        
    END DO
  END SUBROUTINE REFLECTPART
 
END MODULE TRANSPORTUTIL