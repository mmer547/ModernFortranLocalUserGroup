C**********************************************************************
C     UNSTEADY FLOW AROUND CIRCULAR CYLINDER
C           PSI-OMEGA METHOD
C**********************************************************************
C
      program main
      PARAMETER(MX=51,MY=51)
      DIMENSION PSI(MX,MY),OMG(MX,MY),TMP(MX,MY)
C
C***  READ AND CALCULATE PARAMETERS
  99  WRITE(*,*) 'INPUT NUMBER OF MESH--AZIMUSAL & RADIAL(<51) (40,40)'
       READ(*,*) NA,NB
         NX = NA + 1
         NY = NB + 1
      WRITE(*,*) 'INPUT REYNOLDS NUMBER RE (40)'
       READ(*,*) RE
      WRITE(*,*) 'INPUT TIME & SPACE(RADIAL) INCREMENT DT&DY (0.01,0.1)'
       READ(*,*) DT,DY
      WRITE(*,*) 'INPUT NUMBER OF TIME STEP (500)'
       READ(*,*) NMAX
      WRITE(*,*) 'INPUT MAX. NUMBERS OF ITERATION FOR POISSON EQ. (40)'
       READ(*,*) KK
      WRITE(*,*) 'INPUT ACCELARATION PARAMETER (1.0)'
       READ(*,*) CONST1
      WRITE(*,*) 'INPUT MAXNUM ERROR (0.01)'
       READ(*,*) EPS
C
         PAI = ATAN(1.)*4.
         DX = PAI/FLOAT(NX-1)
         DXI = 1./DX
         DYI = 1./DY
         REI =  1./RE
         DX2 = DXI*DXI
         DY2 = DYI*DYI
         FCT = 1./(2.*DX2+2.*DY2)
C
C***  INITIAL CONDITION FOR PSI AND OMEGA
      DO 10 J = 1,NY
      DO 10 I = 1,NX
        PSI(I,J) = EXP((J-1)*DY)*SIN(DX*(I-1))
        OMG(I,J) = 0.0
   10 CONTINUE
C
C***  MAIN LOOP
C
      DO 100 N = 1,NMAX
         FFF = (N-1)/30.
         IF(FFF.GE.1) FFF=1.
C
C***  BOUNDARY CONDITON (STEP1)
C***  ON THE CYLINDER
        DO 20 I = 1,NX
          OMG(I,1) = -2.*PSI(I,2)*DYI*DYI*FFF
          PSI(I,1) = 0.
   20   CONTINUE
C*** ON THE FAR BOUNDARY
        DO 30 I = 1,NX
          PSI(I,NY) = EXP((NY-1)*DY)*SIN(DX*(I-1))
          OMG(I,NY) = 0.
   30   CONTINUE
C*** ALONG THE SYMMETRY LINE
        DO 40 J = 1,NY
          PSI(1,J) = 0.
          OMG(1,J) = 0.
          PSI(MX,J)=0.
          OMG(MX,J)=0.
   40   CONTINUE
C
C*** SOLVE POISSON EQUATION FOR PSI (STEP2)
        FCT = 1./(2.*DX2+2.*DY2)
        DO 50 K = 1,KK
            ERR=0.
          DO 60  J = 2,NY-1
          DO 60  I = 2,NX-1
               RHS = ((PSI(I+1,J)+PSI(I-1,J))*DX2
     1               +(PSI(I,J+1)+PSI(I,J-1))*DY2
     2               +OMG(I,J)*EXP(2.*(J-1)*DY))*FCT
               ERR = ERR+(RHS-PSI(I,J))**2
          PSI(I,J) = PSI(I,J)*(1.-CONST1)+RHS*CONST1
   60     CONTINUE
        IF(ERR.LT .0.00001) GO TO 65
   50   CONTINUE
   65   IF(MOD(N,5).EQ.0)
     1  WRITE(*,*) 'ITERATION NO. =',K,'   ERROR(L2) =',ERR
C
C***  CALCURATE NEW OMEGA (STEP3)
        DO 70 J = 2,NY-1
        DO 70 I = 2,NX-1
C
          TMP(I,J) = OMG(I,J)
C
          RHS = ((OMG(I+1,J)-2.*OMG(I,J)+OMG(I-1,J))*DX2
     1         +(OMG(I,J+1)-2.*OMG(I,J)+OMG(I,J-1))*DY2)*REI
     2         +((PSI(I+1,J)-PSI(I-1,J))*(OMG(I,J+1)-OMG(I,J-1))
     3         -(PSI(I,J+1)-PSI(I,J-1))*(OMG(I+1,J)-OMG(I-1,J)))
     4         *DXI*DYI/4.
          OMG(I,J) = OMG(I,J)+DT*RHS*EXP(-2.*(J-1)*DY)
   70   CONTINUE
C
        ERR1 = 0.
        DO 80 J = 2,NY-1
        DO 80 I = 2,NX-1
          BB = ABS(OMG(I,J)-TMP(I,J))
          IF(BB.GE.ERR1) ERR1 = BB
   80   CONTINUE
C
        IF(MOD(N,5).EQ.0)
     1  WRITE(*,*) N,' *** ERROR(OMG)=' ,ERR1, '  ***'
        IF(N.GT.10.AND.ERR1.LE.EPS) GO TO 90
C
  100 CONTINUE
C***  END OF MAIN LOOP
C
      WRITE(*,*) 'NOT CONVERGE!  DO YOU WANT CONTINUE? (YES=1)'
      READ(*,*) II
      IF(II.EQ.1) GO TO 99
   90 CALL OUT2(PSI,MX,MY,NX,NY,DY)
C     
      WRITE(*,*) 'Save data? Yes=1, No=0'
        READ(*,*) ISAVE
        IF(ISAVE.EQ.1) THEN
          OPEN(8,FIle='Result.txt')
          DO 95 J = 1,NY
          DO 95 I = 1,NX
            WRITE(8,*) I,J,PSI(I,J),OMG(I,J)
   95     CONTINUE
        END IF
      STOP
      end program main
C
      SUBROUTINE OUT2(A,MX,MY,NX,NY,DY)
      DIMENSION A(MX,MY),INDEX(39,15)
C
      PAI=4.*ATAN(1.)
      DX=PAI/FLOAT(NX-1)
C
      AMIN=A(1,1)
      DO 10 J=1,NY
      DO 10 I=1,NX
        IF(A(I,J).LT.AMIN) AMIN=A(I,J)
   10 CONTINUE
      DO 20 J=1,NY
      DO 20 I=1,NX
        A(I,J)=A(I,J)-AMIN
   20 CONTINUE
      AMAX=A(1,1)
      DO 30 J=1,NY
      DO 30 I=1,NX
        IF(A(I,J).GT.AMAX) AMAX=A(I,J)
   30 CONTINUE
C
      DO 40 J=1,15
      DO 40 I=1,39
        IND=0
        IF(I.NE.25) RT=FLOAT(J-1)/ABS(FLOAT(I-25))
        TET=PAI/2.
          IF(I.LE.24) TET=PAI-ATAN(RT)
          IF(I.GE.26) TET=ATAN(RT)
        RR=SQRT(FLOAT((I-25)**2+(J-1)**2))/3.5
          IF(RR.NE.0.) JJ=ALOG(RR)/DY+1
        II=TET/DX+1.5
        IF((II.GE.1.AND.II.LE.NX).AND.(JJ.GE.1.AND.JJ.LE.NY)) THEN
          AA=A(II,JJ)*100./AMAX
          IND=AA+2
            IF(AA.LT.0.) IND=8
        END IF
        INDEX(I,J)=MOD(IND,10)*11
   40 CONTINUE
      DO 50 J=15,1,-1
        WRITE(*,600) (INDEX(I,J),I=39,1,-1)
   50 CONTINUE
      DO 60 J=2,15
        WRITE(*,600) (INDEX(I,J),I=39,1,-1)
   60 CONTINUE
  600   FORMAT(1H,39I2)
      RETURN
      END     


        