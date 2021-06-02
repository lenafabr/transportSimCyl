MODULE GENUTIL
  ! generally useful utilities
  USE MT19937 ! mersenne random number generator

  IMPLICIT NONE
  DOUBLE PRECISION, PARAMETER :: PI = 3.141592653589793D0
  
CONTAINS  

  INTEGER FUNCTION STRING2NUM(STRINGIN,APPENDNUM)
    ! convert the string to a unique number based on ascii characters
    ! the characters SPACE, {,},(,),[,],",`,<,> and all nonprintable characters are ignored
    ! at most the last five characters (ignoring the unacceptable characters above) at the end of the string are used
    ! any leading "!" do not affect the final number (these map to 0)
    ! if APPENDNUM is specificied, only use the last 4 characters of the string as well as the additional number modulo 84

    IMPLICIT NONE
    CHARACTER(LEN=*) :: STRINGIN
    CHARACTER(LEN=5) :: STRING
    INTEGER, OPTIONAL :: APPENDNUM
    INTEGER :: DIGARRAY(5)
    INTEGER :: ALLOWED(84)
    INTEGER :: NVAL, I, D, COUNT
    CHARACTER*84 :: ALLOWEDSTR

    ! set the allowed characters
    ALLOWED(1:6) = (/33,35,36,37,38,39/)
    ALLOWED(7:24) = (/(I,I=42,59)/)
    ALLOWED(25:27) = (/61,63,64/)
    ALLOWED(28:53) = (/(I, I=65,90)/)
    ALLOWED(54:56) = (/92,94,95/)
    ALLOWED(57:82) = (/(I, I=97,122)/)
    ALLOWED(83:84) = (/124,126/)

    NVAL = LEN(STRINGIN)
    IF (PRESENT(APPENDNUM)) THEN
       STRING(1:4) = STRINGIN(NVAL-3:NVAL)
       STRING(5:5) = ACHAR(ALLOWED(MOD(APPENDNUM,84)+1))
    ELSE
       STRING = STRINGIN(NVAL-4:NVAL)
    ENDIF
    NVAL =  5


    DO I = 1,84
       ALLOWEDSTR(I:I) = ACHAR(ALLOWED(I))
    ENDDO

    DIGARRAY = 0
    COUNT = 0
    DO I = 0,NVAL-1
       D = INDEX(ALLOWEDSTR,STRING(NVAL-I:NVAL-I),.FALSE.)
       IF (D.EQ.0) THEN
          print*, 'Ignoring character:', D
          CYCLE
       ENDIF

       DIGARRAY(5-COUNT) = D-1
       COUNT = COUNT + 1
       IF (COUNT.GE.5) EXIT
    ENDDO

    STRING2NUM = BASE2DEC(DIGARRAY,5,84)
  END FUNCTION STRING2NUM

  INTEGER FUNCTION BASE2DEC(DIGARRAY,NVAL,BASE)
  ! given a number in some integer base (specified as a list of digits)
  ! convert that number to a decimal integer
  ! N is the size of the list
  ! if resulting number is too large, wrap around to negative numbers
  ! starting from the right, only use as many of the digits as 
  ! will fit into the resulting integer between -HUGE and HUGE  
  ! if any digit is greater than base-1, print error and stop

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NVAL, BASE
  INTEGER, DIMENSION(NVAL) :: DIGARRAY
  INTEGER :: MAXDIG, I, D

  MAXDIG = INT(LOG(2*DBLE(HUGE(BASE))+2)/LOG(DBLE(BASE)))

  BASE2DEC = 0
  DO I = 0, MIN(NVAL-1,MAXDIG-1)
     D = DIGARRAY(NVAL-I)
     IF (D.EQ.0) CYCLE
     IF (D.GT.BASE-1) THEN
        PRINT*, 'ERROR in BASE2DEC: digit is bigger than base.', I, D, BASE
        STOP 1
     ENDIF
     
     BASE2DEC = BASE2DEC + D*BASE**I
  ENDDO
  
  END FUNCTION BASE2DEC

  SUBROUTINE INTERPARRAY(ARRAY,NA,COL,VAL,IND,INTERP)
    ! for an 2D array with dimensions NA
    ! use the values in column COL to interpolate for the value VAL
    ! return the index IND such that ARRAY(IND,COL)<VAL<ARRAY(IND+1,COL)
    ! and the interpolation of all other columns in INTERP
    ! If val is out of bounds returns IND=0 or IND=NL
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NA(2), COL
    DOUBLE PRECISION, INTENT(IN) :: ARRAY(NA(1),NA(2))
    DOUBLE PRECISION, INTENT(IN) :: VAL
    INTEGER, INTENT(OUT) :: IND
    DOUBLE PRECISION, INTENT(OUT) :: INTERP(NA(2))
    DOUBLE PRECISION :: FRAC

    CALL INTERP1(ARRAY(:,COL),NA(1),VAL,IND)

    IF (IND.EQ.0.OR.IND.GE.NA(1)) RETURN

    FRAC = (VAL-ARRAY(IND,COL))/(ARRAY(IND+1,COL)-ARRAY(IND,COL))
    INTERP = (1-FRAC)*ARRAY(IND,:)+FRAC*ARRAY(IND+1,:)
       
  END SUBROUTINE INTERPARRAY

  SUBROUTINE INTERP1(LIST,NL,VAL,IND)
    ! for a monotonically increasing, double precision list
    ! find the index where LIST(IND) <VAL<LIST(IND+1)
    INTEGER, INTENT(IN) :: NL
    DOUBLE PRECISION, INTENT(IN) :: LIST(NL)
    DOUBLE PRECISION, INTENT(IN) :: VAL
    INTEGER, INTENT(OUT) :: IND
    DOUBLE PRECISION :: MINL, MAXL,PREVMAXL
    INTEGER :: MINI, MAXI,PREVMAXI
    LOGICAL :: VERBOSE = .FALSE.

    MINI = 1; MAXI = NL; MINL = LIST(1); MAXL = LIST(NL)
    IF (VAL.LT.MINL) THEN
       IND = 0; RETURN
    ELSEIF (VAL.EQ.MINL) THEN
       IND = 1; RETURN       
    ELSEIF (VAL.GT.MAXL) THEN
       IND = NL; RETURN
    ELSEIF (VAL.EQ.MAXL) THEN
       IND = NL-1; RETURN
    ENDIF

    DO WHILE (MAXI-MINI.GT.1.OR.MAXL.LT.VAL)
       IF (MAXL.GT.VAL) THEN
          PREVMAXL = MAXL; PREVMAXI = MAXI
          MAXI = MINI + (MAXI-MINI)/2
          MAXL = LIST(MAXI)
       ELSE
          MINI = MAXI; MAXI = PREVMAXI
          MINL = MAXL; MAXL = PREVMAXL
       ENDIF
       IF (VERBOSE) PRINT*, 'MINI, MAXI, MINL, MAXL', MINI, MAXI, MINL, MAXL,VAL
       if (maxi.eq.mini) then
          print*, 'something weird in interp1:', list(1), list(nl), val          
          stop 1
       endif
    ENDDO

    IF (.NOT.(MAXI.EQ.MINI+1.AND.LIST(MINI).LE.VAL.AND.LIST(MAXI).GE.VAL)) THEN
       PRINT*, 'SOMETHING IS WEIRD IN INTERP1', val, mini, maxi, list(mini), list(maxi)
       STOP 1
    ENDIF

    IND = MINI
  END SUBROUTINE INTERP1

  FUNCTION RANDARRAY(NVAL)
  ! returns an array of size nval filled with uniform random numbers between 0 and 1
    IMPLICIT NONE
    INTEGER :: NVAL
    DOUBLE PRECISION :: RANDARRAY(NVAL)
    INTEGER :: NC

    RANDARRAY = 0D0
    DO NC = 1,NVAL
      RANDARRAY(NC) = GRND()
    END DO
  END FUNCTION RANDARRAY

  SUBROUTINE REPLACESUBSTR(INSTRING,C,REPL)
    ! replace the last instance of the substring C in INSTRING with REPL
    IMPLICIT NONE
    CHARACTER*100 :: INSTRING
    CHARACTER(LEN=*) :: C
    CHARACTER(LEN=*) :: REPL
    INTEGER :: LENC, IND

    INSTRING = ADJUSTL(INSTRING)

    LENC = LEN_TRIM(C)

    IND = INDEX(INSTRING,C,.TRUE.)
    IF (IND.GT.0) THEN! if * was found in the string

       INSTRING = INSTRING(1:IND-1) // TRIM(ADJUSTL(REPL)) // INSTRING(IND+LENC:100)   
    END IF
  END SUBROUTINE REPLACESUBSTR

  SUBROUTINE NORMALIZE(X)
    ! normalize a 3 dimensional vector

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(INOUT) :: X(3)
    DOUBLE PRECISION :: DX

    DX = SQRT(DOT_PRODUCT(X,X))
    X(:) = X(:)/DX

    RETURN
  END SUBROUTINE NORMALIZE

  DOUBLE PRECISION FUNCTION NORM(X)
    ! norm of 3D vector X

    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: X(3)

    NORM = sqrt(DOT_PRODUCT(X,X))

  END FUNCTION NORM

  SUBROUTINE CROSS_PRODUCT(A, B, C)
    ! take the cross product of 3D vectors A and B; return result in C

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN) :: A(3), B(3)
    DOUBLE PRECISION, INTENT(OUT) :: C(3)

    C(1) = A(2)*B(3) - A(3)*B(2)
    C(2) = A(3)*B(1)-A(1)*B(3)
    C(3) = A(1)*B(2) - A(2)*B(1)

    RETURN
  END SUBROUTINE CROSS_PRODUCT

  SUBROUTINE RANDOMAXIS(REFAX,CTRANGE,RANDAX)
    ! generate a random axis, within a certain range in cos(theta) 
    ! relative to the reference axis    
    DOUBLE PRECISION, INTENT(IN) :: REFAX(3), CTRANGE
    DOUBLE PRECISION, INTENT(OUT) :: RANDAX(3)
    DOUBLE PRECISION :: THETA, pHI, E1(3), E2(3), E3(3), X, Y, Z

    THETA = 1.0D0
    DO WHILE (THETA.EQ.1.0D0)
       ! get random number; MT19937 algorithm uses closed interval [0,1], 
       ! so ignore when R is exactly 1
       THETA = GRND() !get a random number
    ENDDO
    !CALL RANDOM_NUMBER(THETA)
    THETA = acos(1D0 - THETA*MAX(CTRANGE,2D0))

    PHI = 1.0D0
    DO WHILE (PHI.EQ.1.0D0)
       PHI = GRND() !get a random number
    ENDDO
    !CALL RANDOM_NUMBER(PHI)
    PHI = PHI*2*PI

    ! axis system relative to which angles are defined
    E3 = REFAX
    IF (E3(2) == 0 .AND. E3(3) == 0) THEN
       E2 = (/0D0,1D0,0D0/)
    ELSE
       CALL CROSS_PRODUCT(E3, (/1D0,0D0,0D0/),E2)
    END IF
    CALL CROSS_PRODUCT(E2, E3, E1)

    CALL NORMALIZE(E1); CALL NORMALIZE(E2); CALL NORMALIZE(E3)

    ! generate the axis around which to rotate
    X = sin(THETA)*cos(PHI)
    Y = sin(THETA)*sin(PHI)
    Z = cos(THETA)

    RANDAX = X*E1 + Y*E2 + Z*E3

  END SUBROUTINE RANDOMAXIS

  SUBROUTINE GETPERP(V1,V2)
    ! get some unit vector perpendicular to V1 and store it in V2
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: V1(3)
    DOUBLE PRECISION, INTENT(OUT) :: V2(3)

    IF (V1(2).EQ.0.AND.V1(3).EQ.0) THEN
       V2 = (/0D0,1D0,0D0/)
    ELSE
       CALL CROSS_PRODUCT(V1,(/0D0,1D0,0D0/),V2)
       CALL NORMALIZE(V2)
    ENDIF
  END SUBROUTINE GETPERP
  
  
  SUBROUTINE GETORTHOPROJ(V1,V2,V0,MINDIST,MVAL,VOUT)
  !get the orthogonal projection of V0 on the line joining v1 and v2
  !the orthogonal distance is sored in mindist
  !fraction of distance from v1 is stored in mval
  !vout gives the location of the projected point
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: V1(3),V2(3),V0(3)
    DOUBLE PRECISION, INTENT(OUT) :: MINDIST,MVAL,VOUT(3)
    DOUBLE PRECISION :: UVEC(3),UMAG
    UVEC = V2-V1
    UMAG = NORM2(UVEC)
    MVAL = DOT_PRODUCT(V0-V1,UVEC)/UMAG**2
    VOUT = V1+MVAL*UVEC
    MINDIST = NORM2(V0-VOUT)
  END SUBROUTINE GETORTHOPROJ

  SUBROUTINE GETCLOSESTPT(V1,V2,V0,MINDIST,MVAL,VOUT)
  !returns the point on the segment joining v1 and v2 that is closest to v0
  !if MVAL is between 0 and 1, the closest point is on the segment
  !if MVAL = 0, the closest point is V1
  !if MVAL = 1, the closet point is V2
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: V1(3),V2(3),V0(3)
    DOUBLE PRECISION, INTENT(OUT) :: MINDIST,MVAL,VOUT(3)

    CALL GETORTHOPROJ(V1,V2,V0,MINDIST,MVAL,VOUT)
    IF(MVAL.LT.0D0) THEN
      MINDIST = NORM2(V1-V0)
      MVAL = 0D0
      VOUT = V1
    ELSE IF(MVAL.GT.1D0)THEN
      MINDIST = NORM2(V2-V0)
      MVAL = 1D0
      VOUT = V2
    END IF

    IF(MVAL.GT.1D0.OR.MVAL.LT.0D0) THEN
      PRINT*, 'CHECK INPUT VALUES TO GETCLOSESTPT'
      STOP 1
    END IF
  END SUBROUTINE GETCLOSESTPT

END MODULE GENUTIL
