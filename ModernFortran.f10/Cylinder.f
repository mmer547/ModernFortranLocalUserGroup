C**********************************************************************
C     UNSTEADY FLOW AROUND CIRCULAR CYLINDER
C           PSI-OMEGA METHOD
C**********************************************************************
C
      program main
      parameter(MX=51,MY=51)
      dimension PSI(MX,MY),OMG(MX,MY),TMP(MX,MY)
C
C***  READ AND CALCULATE PARAMETERS
  99  print *, 'INPUT NUMBER OF MESH--AZIMUSAL & RADIAL(<51) (40,40)'
       read *, NA,NB
         NX = NA + 1
         NY = NB + 1
      print *, 'INPUT REYNOLDS NUMBER RE (40)'
       read *, RE
      print *, 'INPUT TIME & SPACE(RADIAL) INCREMENT DT&DY (0.01,0.1)'
       read *, DT,DY
      print *, 'INPUT NUMBER OF TIME STEP (500)'
       read *, NMAX
      print *, 'INPUT MAX. NUMBERS OF ITERATION FOR POISSON EQ. (40)'
       read *, KK
      print *, 'INPUT ACCELARATION PARAMETER (1.0)'
       read *, CONST1
      print *, 'INPUT MAXNUM ERROR (0.01)'
       read *, EPS
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
      do 10 J = 1,NY
      do 10 I = 1,NX
        PSI(I,J) = EXP((J-1)*DY)*SIN(DX*(I-1))
        OMG(I,J) = 0.0
   10 continue
C
C***  MAIN LOOP
C
      do 100 N = 1,NMAX
         FFF = (N-1)/30.
         if(FFF.GE.1) FFF=1.
C
C***  BOUNDARY CONDITON (STEP1)
C***  ON THE CYLINDER
        do 20 I = 1,NX
          OMG(I,1) = -2.*PSI(I,2)*DYI*DYI*FFF
          PSI(I,1) = 0.
   20   continue
C*** ON THE FAR BOUNDARY
        do 30 I = 1,NX
          PSI(I,NY) = EXP((NY-1)*DY)*SIN(DX*(I-1))
          OMG(I,NY) = 0.
   30   continue
C*** ALONG THE SYMMETRY LINE
        do 40 J = 1,NY
          PSI(1,J) = 0.
          OMG(1,J) = 0.
          PSI(MX,J)=0.
          OMG(MX,J)=0.
   40   continue
C
C*** SOLVE POISSON EQUATION FOR PSI (STEP2)
        FCT = 1./(2.*DX2+2.*DY2)
        do 50 K = 1,KK
            ERR=0.
          do 60  J = 2,NY-1
          do 60  I = 2,NX-1
               RHS = ((PSI(I+1,J)+PSI(I-1,J))*DX2
     1               +(PSI(I,J+1)+PSI(I,J-1))*DY2
     2               +OMG(I,J)*EXP(2.*(J-1)*DY))*FCT
               ERR = ERR+(RHS-PSI(I,J))**2
          PSI(I,J) = PSI(I,J)*(1.-CONST1)+RHS*CONST1
   60     continue
        if(ERR.LT .0.00001) go to 65
   50   continue
   65   if(MOD(N,5).EQ.0)
     1  print *, 'ITERATION NO. =',K,'   ERROR(L2) =',ERR
C
C***  CALCURATE NEW OMEGA (STEP3)
        do 70 J = 2,NY-1
        do 70 I = 2,NX-1
C
          TMP(I,J) = OMG(I,J)
C
          RHS = ((OMG(I+1,J)-2.*OMG(I,J)+OMG(I-1,J))*DX2
     1         +(OMG(I,J+1)-2.*OMG(I,J)+OMG(I,J-1))*DY2)*REI
     2         +((PSI(I+1,J)-PSI(I-1,J))*(OMG(I,J+1)-OMG(I,J-1))
     3         -(PSI(I,J+1)-PSI(I,J-1))*(OMG(I+1,J)-OMG(I-1,J)))
     4         *DXI*DYI/4.
          OMG(I,J) = OMG(I,J)+DT*RHS*EXP(-2.*(J-1)*DY)
   70   continue
C
        ERR1 = 0.
        do 80 J = 2,NY-1
        do 80 I = 2,NX-1
          BB = ABS(OMG(I,J)-TMP(I,J))
          if(BB.GE.ERR1) ERR1 = BB
   80   continue
C
        if(MOD(N,5).EQ.0)
     1  print *, N,' *** ERROR(OMG)=' ,ERR1, '  ***'
        if(N.GT.10.AND.ERR1.LE.EPS) go to 90
C
  100 continue
C***  END OF MAIN LOOP
C
      print *, 'NOT CONVERGE!  DO YOU WANT CONTINUE? (YES=1)'
      read *, II
      if(II.EQ.1) go to 99
   90 call OUT2(PSI,MX,MY,NX,NY,DY)
C     
      print *, 'Save data? Yes=1, No=0'
        read *, ISAVE
        if(ISAVE.EQ.1) THEN
          open(8,FIle='Result.txt')
          do 95 J = 1,NY
          do 95 I = 1,NX
            write(8,*) I,J,PSI(I,J),OMG(I,J)
   95     continue
        end if
      stop
      end program main
C
      SUBROUTINE OUT2(A,MX,MY,NX,NY,DY)
      dimension A(MX,MY),INDEX(39,15)
C
      PAI=4.*ATAN(1.)
      DX=PAI/FLOAT(NX-1)
C
      AMIN=A(1,1)
      do 10 J=1,NY
      do 10 I=1,NX
        if(A(I,J).LT.AMIN) AMIN=A(I,J)
   10 continue
      do 20 J=1,NY
      do 20 I=1,NX
        A(I,J)=A(I,J)-AMIN
   20 continue
      AMAX=A(1,1)
      do 30 J=1,NY
      do 30 I=1,NX
        if(A(I,J).GT.AMAX) AMAX=A(I,J)
   30 continue
C
      do 40 J=1,15
      do 40 I=1,39
        IND=0
        if(I.NE.25) RT=FLOAT(J-1)/ABS(FLOAT(I-25))
        TET=PAI/2.
          if(I.LE.24) TET=PAI-ATAN(RT)
          if(I.GE.26) TET=ATAN(RT)
        RR=SQRT(FLOAT((I-25)**2+(J-1)**2))/3.5
          if(RR.NE.0.) JJ=ALOG(RR)/DY+1
        II=TET/DX+1.5
        if((II.GE.1.AND.II.LE.NX).AND.(JJ.GE.1.AND.JJ.LE.NY)) THEN
          AA=A(II,JJ)*100./AMAX
          IND=AA+2
            if(AA.LT.0.) IND=8
        end if
        INDEX(I,J)=MOD(IND,10)*11
   40 continue
      do 50 J=15,1,-1
        write(*,600) (INDEX(I,J),I=39,1,-1)
   50 continue
      do 60 J=2,15
        write(*,600) (INDEX(I,J),I=39,1,-1)
   60 continue
  600   format(1H,39I2)
      return
      END     


        