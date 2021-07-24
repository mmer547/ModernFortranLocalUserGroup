C**********************************************************************
C     UNSTEADY FLOW AROUND CIRCULAR CYLINDER
C           PSI-OMEGA METHOD
C**********************************************************************
C
      program main
        implicit none
        integer, parameter :: MX=51, MY=51
        real PSI(MX, MY), OMG(MX, MY), TMP(MX, MY)
        integer NA, NB, NX, NY
        integer RE
        real DT, DY
        integer NMAX, K, KK
        real CONST1, EPS
        real PAI
        real DX, DXI, DYI, REI, DX2, DY2, FCT
        integer J, I, N, II, ISAVE
        real FFF
        real ERR, ERR1
        real RHS
        real BB
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
      do J = 1,NY
        do I = 1,NX
          PSI(I,J) = EXP((J-1)*DY)*SIN(DX*(I-1))
          OMG(I,J) = 0.0
        end do
      end do
C
C***  MAIN LOOP
C
      do N = 1,NMAX
         FFF = (N-1)/30.
         if(FFF.GE.1) FFF=1.
C
C***  BOUNDARY CONDITON (STEP1)
C***  ON THE CYLINDER
        do I = 1,NX
          OMG(I,1) = -2.*PSI(I,2)*DYI*DYI*FFF
          PSI(I,1) = 0.
        end do
C*** ON THE FAR BOUNDARY
        do I = 1,NX
          PSI(I,NY) = EXP((NY-1)*DY)*SIN(DX*(I-1))
          OMG(I,NY) = 0.
        end do
C*** ALONG THE SYMMETRY LINE
        do J = 1,NY
          PSI(1,J) = 0.
          OMG(1,J) = 0.
          PSI(MX,J)=0.
          OMG(MX,J)=0.
        end do
C
C*** SOLVE POISSON EQUATION FOR PSI (STEP2)
        FCT = 1./(2.*DX2+2.*DY2)
        do K = 1,KK
            ERR=0.
          do J = 2,NY-1
            do I = 2,NX-1
              RHS = ((PSI(I+1,J)+PSI(I-1,J))*DX2
     1             +(PSI(I,J+1)+PSI(I,J-1))*DY2
     2             +OMG(I,J)*EXP(2.*(J-1)*DY))*FCT
              ERR = ERR+(RHS-PSI(I,J))**2
              PSI(I,J) = PSI(I,J)*(1.-CONST1)+RHS*CONST1
            end do
          end do
        if(ERR.LT .0.00001) exit
        end do
        if(MOD(N,5).EQ.0)
     1  print *, 'ITERATION NO. =',K,'   ERROR(L2) =',ERR
C
C***  CALCURATE NEW OMEGA (STEP3)
        do J = 2,NY-1
          do I = 2,NX-1
C
            TMP(I,J) = OMG(I,J)
C
            RHS = ((OMG(I+1,J)-2.*OMG(I,J)+OMG(I-1,J))*DX2
     1           +(OMG(I,J+1)-2.*OMG(I,J)+OMG(I,J-1))*DY2)*REI
     2           +((PSI(I+1,J)-PSI(I-1,J))*(OMG(I,J+1)-OMG(I,J-1))
     3           -(PSI(I,J+1)-PSI(I,J-1))*(OMG(I+1,J)-OMG(I-1,J)))
     4           *DXI*DYI/4.
            OMG(I,J) = OMG(I,J)+DT*RHS*EXP(-2.*(J-1)*DY)
          end do
        end do
C
        ERR1 = 0.
        do J = 2,NY-1
          do I = 2,NX-1
            BB = ABS(OMG(I,J)-TMP(I,J))
            if(BB.GE.ERR1) ERR1 = BB
          end do
        end do
C
        if(MOD(N,5).EQ.0)
     1  print *, N,' *** ERROR(OMG)=' ,ERR1, '  ***'
        if(N.GT.10.AND.ERR1.LE.EPS) go to 90
C
      end do
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
          do J = 1,NY
            do I = 1,NX
              write(8,*) I,J,PSI(I,J),OMG(I,J)
            end do
          end do
        end if
      stop
      end program main
C
      SUBROUTINE OUT2(A,MX,MY,NX,NY,DY)
      implicit none
      real A(MX,MY)
      integer INDEX(39,15)
      real PAI
      integer MX, MY, NX, NY
      real DX, DY
      real AMIN, AMAX
      integer J, I, IND, II, JJ
      real RT, TET, RR, AA
C
      PAI=4.*ATAN(1.)
      DX=PAI/FLOAT(NX-1)
C
      AMIN=A(1,1)
      do J=1,NY
        do I=1,NX
          if(A(I,J).LT.AMIN) AMIN=A(I,J)
        end do
      end do
      do J=1,NY
        do I=1,NX
          A(I,J)=A(I,J)-AMIN
        end do
      end do
      AMAX=A(1,1)
      do J=1,NY
        do I=1,NX
          if(A(I,J).GT.AMAX) AMAX=A(I,J)
        end do
      end do
C
      do J=1,15
        do I=1,39
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
        end do
      end do
      do J=15,1,-1
        write(*,600) (INDEX(I,J),I=39,1,-1)
      end do
      do J=2,15
        write(*,600) (INDEX(I,J),I=39,1,-1)
      end do
  600   format(1H,39I2)
      return
      END     


        