      PROGRAM ROS
C YPARXEI KAI GROUND EFFECT.
C YPOLOGIZETAI H DIND
      IMPLICIT REAL*4 (A-H,O-Z)
      COMMON/AVV1/ANX(4000),ANY(4000),ANZ(4000),GX(4000),GY(4000),
     1GZ(4000)
      COMMON/AVV2/ATX(4000),ATY(4000),ATZ(4000),ALX(4000),ALY(4000),
     1ALZ(4000)
      COMMON/GAVV/X(4000),Y(4000),Z(4000)
      COMMON/GW/XW(60,4000),YW(60,4000),ZW(60,4000)
      COMMON/MP/AL(4000)
      COMMON/ANPORT/DS(4000)
      COMMON/GVV/C(4000)
      COMMON/VW/CW(60,4000),ICW(4000,2)
      COMMON/GRID/IC(4000,4),MARK(4000),MARKW(4000)
      COMMON/PPG/IPROP(4000)
      COMMON/WAICOR/NC(4000)
      COMMON/WIN/IPW(4000)
      COMMON/WINTER/IPROPW(4000)
      COMMON/NEIB/KN(4000,4000)
      COMMON/VG1/A(4000,4000)
      COMMON/VG2/B(4000)
      
      
      RHO=1.225
      ITER=0
      ALIFT1=0.


      OPEN(UNIT=21,FILE='WN.TST',STATUS='UNKNOWN')
      OPEN(UNIT=10,FILE='WN.RES',STATUS='UNKNOWN')
      OPEN(UNIT=20,FILE='WN.WAK',STATUS='UNKNOWN')
      OPEN(UNIT=30,FILE='WN.PAN',STATUS='UNKNOWN')
      OPEN(UNIT=40,FILE='WN.DAT',STATUS='UNKNOWN')
      OPEN(UNIT=7,FILE='DELTA.PAN',STATUS='UNKNOWN')
      OPEN(UNIT=8,FILE='DELTA.DAT',STATUS='UNKNOWN')
      OPEN(UNIT=9,FILE='DELTA.INP',STATUS='UNKNOWN')
      open(unit=41,file='tt',status='unknown')
      OPEN(UNIT=46,FILE='SVWN',STATUS='UNKNOWN')
      OPEN(UNIT=42,FILE='WVWN',STATUS='UNKNOWN')
      OPEN(UNIT=43,FILE='NWN',STATUS='UNKNOWN')
      OPEN(UNIT=45,FILE='IWWN',STATUS='UNKNOWN')
      OPEN(UNIT=47,FILE='LINEWN.INP',STATUS='UNKNOWN') 
      OPEN(UNIT=48,FILE='LIWN',STATUS='UNKNOWN')
C 
      VAIP=0.
C
      CALL GEOM(NPAN,NGRID,ALF,BET,GAM,VINIT,EPS,DT,NSYM,NGRND,HFL)
      write (*,'(a)') "--------------------------"
      write (*,'(a)') "Printing ROS Configuration"
      write (*,'(a)') "Angles  (rads)"
      write (*,'(a)') "Lengths (meters)"
      write (*,'(a)') "Max Iterations : 60"
      write (*,'(a)') "--------------------------"
      write (*,'(a)') "Angles"
      write (*,'(a,f10.3)') "   Alpha : ", ALF
      write (*,'(a,f10.3)') "   Beta  : ", BET
      write (*,'(a,f10.3)') "   Gamma : ", GAM
      write (*,'(a)') "--------------------------"
      write (*,'(a,f10.3)') "Vinit        : ", VINIT
      write (*,'(a,f10.8)') "Epsilon      : ", EPS
      write (*,'(a,f10.3)') "DT           : ", DT
      write (*,'(a)') "--------------------------"
      write (*,'(a)') "Indices"
      write (*,'(a,I1)') "   Symmetry      : ", NSYM
      write (*,'(a,I1)') "   Ground Effect : ", NGRND
      write (*,'(a)') "--------------------------"
      write (*,'(a)') " "
      write (*,'(a)') "STARTING SOLVER"


      CALL ANALGEO(NPAN)
      IF(ITER.EQ.0) GO TO 94
   93 CALL WAKE(ITER,VINIT,DT,NPAN,NPW,NGW,NGRID)
      IF(ITER.EQ.1) CALL NEIBORG(NPAN,NGW)
      IF(ITER.EQ.1) GO TO 96
      CALL WAKINT(NPW,ITER)
      CALL WAKCOR(NPAN,NGW,ITER)
   96 CALL WAKREL(ITER,NPAN,NGW,NPW,DT,NSYM,NGRND)
      CALL WAKINT(NPW,ITER)
C      CALL IWCOR(ITER,NGW,NPAN)
      IF(ITER.EQ.1) GO TO 94
      CALL WAKCOR(NPAN,NGW,ITER)
   94 CALL VORCALC(VINIT,NPAN,ITER,NPW,NSYM,VAIP,NGRND)
      IF(ITER.EQ.0) GO TO 103
      ALIFT1 = ALIFT
	   CALL AIRLOAD1(NPAN,VINIT,RHO,ALIFT,DRAG,SIDE,NSYM,ITER,
	1NPW,NGRND)
      write (*,'(a)') "--------------------------"
      write (*,'(a,i2)') "Iteration : ", ITER
      write (*,'(a,f15.3)') "Lift   : ", ALIFT
      write (*,'(a,f15.3)') "Drag   : ", DRAG 
      write (*,'(a,f12.8)') "Error  : ", ABS((ALIFT1-ALIFT)/ALIFT)
      write (*,'(a)') "--------------------------"
      WRITE (48,*) ITER,ALIFT
      IF(ABS((ALIFT1-ALIFT)/ALIFT).LE.EPS.OR.ITER.EQ.60) THEN
      GO TO 95
      END IF
  103 ITER=ITER+1
      GO TO 93
 1002 FORMAT (5X,I8,F15.4)
   95 DO 99 I=1,ITER+1
      DO 98 J=1,NGW
      WRITE(20,1001) I,J,XW(I,J),YW(I,J),ZW(I,J)
   98 CONTINUE
   99 CONTINUE
 1001 FORMAT (5X,2I8,3F15.4)
      DO 102 I=1,NPAN
      WRITE(10,1002) I,AL(I)
  102 CONTINUE
      DO 104 I=1,NGRID
      WRITE(40,110) I,X(I),Y(I),Z(I),MARK(I)
  104 CONTINUE
  110 FORMAT(I10,3F15.5,I10)
      DO 105 I=1,NPAN
      WRITE(30,120) I,IC(I,1),IC(I,2),IC(I,3),IC(I,4)
  105 CONTINUE
  120 FORMAT(5I10)
      IF(VAIP.EQ.0.) GO TO 121
      CALL CPAIP(ITER,NPAN,NPW,NSYM,VINIT,NGRND) 
  121 CONTINUE
      DO 122 I=1,NPAN
      WRITE(46,*) C(I)
  122 CONTINUE
      DO 124 I=1,ITER
      DO 123 J=1,NPW
      WRITE(42,*) CW(I,J)
  123 CONTINUE
  124 CONTINUE
      DO 125 I=1,NPAN
      WRITE(43,*) ANX(I),ANY(I),ANZ(I),ALX(I),ALY(I),ALZ(I),
     1ATX(I),ATY(I),ATZ(I),GX(I),GY(I),GZ(I),DS(I)
  125 CONTINUE
      DO 126 I=1,NPW
      WRITE(45,*) ICW(I,1),ICW(I,2)
  126 CONTINUE
      WRITE(47,*) NPAN,NGRID,VINIT,ITER,NPW,NGW,NSYM,NGRND,
     1ALF,BET,GAM,HFL
      END


      SUBROUTINE GEOM(NPAN,NGRID,ALF,BET,GAM,VINIT,EPS,DT,NSYM,NGRND,
	1HFL)
      IMPLICIT REAL*4 (A-H,O-Z)
      COMMON/GAVV/X(4000),Y(4000),Z(4000)
      COMMON/GRID/IC(4000,4),MARK(4000),MARKW(4000)
      COMMON/PPG/IPROP(4000)
      PI=3.14159
      READ(9,100) NPAN,NGRID,ALF,BET,GAM,VINIT,EPS,AL,NSYM,INCH,
     1NGRND,HFL
      ALF=ALF*PI/180.
      BET=BET*PI/180.
      GAM=GAM*PI/180.
C      IF(INCH.EQ.1) AL=AL/1000.
C      DT=0.25*AL/VINIT
C      DT=0.03
      DO 10 I=1,NGRID
      READ(8,110) J,X(I),Y(I),Z(I),MARK(I)
      IF(INCH.EQ.0) GO TO 5
      X(I)=X(I)/1000.
      Y(I)=Y(I)/1000.
      Z(I)=Z(I)/1000.
    5 XX=X(I)*COS(ALF)*COS(BET)+Y(I)*SIN(BET)+Z(I)*SIN(ALF)*COS(BET)
      YY=-X(I)*(COS(ALF)*SIN(BET)*COS(GAM)-SIN(ALF)*SIN(GAM))+
     1Y(I)*COS(BET)*COS(GAM)-Z(I)*(SIN(ALF)*SIN(BET)*COS(GAM)-
     2COS(ALF)*SIN(GAM))
      ZZ=X(I)*(COS(ALF)*SIN(BET)*SIN(GAM)-SIN(ALF)*COS(GAM))-
     1Y(I)*COS(BET)*SIN(GAM)+Z(I)*(SIN(ALF)*SIN(BET)*SIN(GAM)+
     2COS(ALF)*COS(GAM))
      X(I)=XX
      Y(I)=YY
      Z(I)=ZZ
      IF(NGRND.EQ.1) Z(I)=Z(I)+HFL
   10 CONTINUE
      AL = ABS(MINVAL(X) - MAXVAL(X))
      print *, "###################################"
      print *, "TROUBLESHOOTING PRINT :"
      print *, "MIN X :", MINVAL(X)
      print *, "MAX X :", MAXVAL(X)
      print *, "AL ABS(MIN(X) - MAX(X)) :", AL
      print *, "###################################"
      DT=0.25*AL/VINIT

      DO 20 I=1,NPAN
      READ(7,120) J,IC(I,1),IC(I,2),IC(I,3),IC(I,4),IPROP(I)
   20 CONTINUE 
  100 FORMAT(I10,/,I10,/,F10.3,/,F10.3,/,F10.3,/,F10.3,/,F10.3,/,F10.3
     1,/,I10,/,I10,/,I10,/,F10.3)
  110 FORMAT(I10,3F15.5,I10)
  120 FORMAT(6I10)
      CLOSE(7)
      CLOSE(8)
      CLOSE(9)
      RETURN
      END


      SUBROUTINE ANALGEO(NPAN)
      IMPLICIT REAL*4 (A-H,O-Z)
      COMMON/AVV1/ANX(4000),ANY(4000),ANZ(4000),GX(4000),GY(4000),
     1GZ(4000)
      COMMON/AVV2/ATX(4000),ATY(4000),ATZ(4000),ALX(4000),ALY(4000),
     1ALZ(4000)
      COMMON/GAVV/X(4000),Y(4000),Z(4000)
      COMMON/GRID/IC(4000,4),MARK(4000),MARKW(4000)
      COMMON/ANPORT/DS(4000)
      DO 10 I=1,NPAN
      IF(IC(I,3).NE.IC(I,4)) GO TO 5
      GX(I)=(X(IC(I,1))+X(IC(I,2))+X(IC(I,3)))/3.
      GY(I)=(Y(IC(I,1))+Y(IC(I,2))+Y(IC(I,3)))/3.
      GZ(I)=(Z(IC(I,1))+Z(IC(I,2))+Z(IC(I,3)))/3.
      GO TO 6
    5 GX(I)=(X(IC(I,1))+X(IC(I,2))+X(IC(I,3))+X(IC(I,4)))/4.
      GY(I)=(Y(IC(I,1))+Y(IC(I,2))+Y(IC(I,3))+Y(IC(I,4)))/4.
      GZ(I)=(Z(IC(I,1))+Z(IC(I,2))+Z(IC(I,3))+Z(IC(I,4)))/4.
    6 X13=X(IC(I,1))-X(IC(I,3))
      Y13=Y(IC(I,1))-Y(IC(I,3))
      Z13=Z(IC(I,1))-Z(IC(I,3))
      X42=X(IC(I,4))-X(IC(I,2))
      Y42=Y(IC(I,4))-Y(IC(I,2))
      Z42=Z(IC(I,4))-Z(IC(I,2))
      XX=Y13*Z42-Z13*Y42
      YY=Z13*X42-X13*Z42
      ZZ=X13*Y42-Y13*X42
      ANORM=SQRT(XX**2+YY**2+ZZ**2)
      ANX(I)=-XX/ANORM
      ANY(I)=-YY/ANORM
      ANZ(I)=-ZZ/ANORM
      X12=X(IC(I,2))-X(IC(I,1))
      Y12=Y(IC(I,2))-Y(IC(I,1))
      Z12=Z(IC(I,2))-Z(IC(I,1))
      X14=X(IC(I,4))-X(IC(I,1))
      Y14=Y(IC(I,4))-Y(IC(I,1))
      Z14=Z(IC(I,4))-Z(IC(I,1))
      XX=Y13*Z12-Z13*Y12
      YY=Z13*X12-X13*Z12
      ZZ=X13*Y12-Y13*X12
      DS1=0.5*SQRT(XX**2+YY**2+ZZ**2)
      XX=Y13*Z14-Z13*Y14
      YY=Z13*X14-X13*Z14
      ZZ=X13*Y14-Y13*X14
      DS2=0.5*SQRT(XX**2+YY**2+ZZ**2)
      DS(I)=DS1+DS2
      XM=(X(IC(I,1))+X(IC(I,2)))/2.
      YM=(Y(IC(I,1))+Y(IC(I,2)))/2.
      ZM=(Z(IC(I,1))+Z(IC(I,2)))/2.
      XX=XM-GX(I)
      YY=YM-GY(I)
      ZZ=ZM-GZ(I)
      ANORM=SQRT(XX**2+YY**2+ZZ**2)
      ALX(I)=XX/ANORM
      ALY(I)=YY/ANORM
      ALZ(I)=ZZ/ANORM
      XX=ALY(I)*ANZ(I)-ALZ(I)*ANY(I)
      YY=ALZ(I)*ANX(I)-ALX(I)*ANZ(I)
      ZZ=ALX(I)*ANY(I)-ALY(I)*ANX(I)
      ANORM=SQRT(XX**2+YY**2+ZZ**2)
      ATX(I)=XX/ANORM
      ATY(I)=YY/ANORM
      ATZ(I)=ZZ/ANORM
   10 CONTINUE
      RETURN
      END


      SUBROUTINE VORCALC(VINIT,NPAN,ITER,NPW,NSYM,VAIP,NGRND)
      IMPLICIT REAL*4 (A-H,O-Z)
      COMMON/AVV1/ANX(4000),ANY(4000),ANZ(4000),GX(4000),GY(4000),
     1GZ(4000)
      COMMON/GAVV/X(4000),Y(4000),Z(4000)
      COMMON/VOW/XX(4000),YY(4000),ZZ(4000)
      COMMON/VG1/A(4000,4000)
      COMMON/VG2/B(4000)
      COMMON/VVOR/XG,YG,ZG,X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4,
     1VVX,VVY,VVZ
      COMMON/VW/CW(60,4000),ICW(4000,2)
      COMMON/GW/XW(60,4000),YW(60,4000),ZW(60,4000)
      COMMON/GRID/IC(4000,4),MARK(4000),MARKW(4000)
      COMMON/GVV/XP(4000)
      COMMON/PPG/IPROP(4000)
      IF(ITER.NE.0) GO TO 25
      DO 22 I=1,NPAN
      DO 21 J=1,NPAN
      A(I,J)=0.
   21 CONTINUE
   22 CONTINUE
      DO 20 I=1,NPAN
      XG=GX(I)
      YG=GY(I)
      ZG=GZ(I)
      NG=0
    6 DO 10 J=1,NPAN
      NA=0
      X1=X(IC(J,1))
      Y1=Y(IC(J,1))
      Z1=Z(IC(J,1))
      X2=X(IC(J,2))
      Y2=Y(IC(J,2))
      Z2=Z(IC(J,2))
      X3=X(IC(J,3))
      Y3=Y(IC(J,3))
      Z3=Z(IC(J,3))
      X4=X(IC(J,4))
      Y4=Y(IC(J,4))
      Z4=Z(IC(J,4))
      IF(NG.EQ.0) GO TO 11
      Z1=-Z1
      Z2=-Z2
      Z3=-Z3
      Z4=-Z4
   11 CALL VORTEX   
      A(I,J)=A(I,J)+(VVX*ANX(I)+VVY*ANY(I)+VVZ*ANZ(I))*((-1)**NA)*
     1(-1)**NG
      IF(NSYM.EQ.0) GO TO 10
      IF(NA.EQ.1) GO TO 10
      NA=1
      Y1=-Y1
      Y2=-Y2
      Y3=-Y3
      Y4=-Y4
      GO TO 11
   10 CONTINUE
      IF(NGRND.EQ.0) GO TO 20
      IF(NG.EQ.1) GO TO 20
      NG=1
      GO TO 6
   20 CONTINUE

   25 DO 30 II=1,NPAN
      XG=GX(II)
      YG=GY(II)
      ZG=GZ(II)
      VXV=0.
      VYV=0.
      VZV=0.
      KPV=1
      IF(IPROP(II).EQ.14) KPV=0
      IF(ITER.EQ.0) GO TO 26
      CALL VELWAK(ITER,NPW,NSYM,NGRND,XG,YG,ZG,VXV,
     1VYV,VZV)
   26 B(II)=-(VINIT*KPV+VXV)*ANX(II)-VYV*ANY(II)-VZV*
     1ANZ(II)+VAIP
   30 CONTINUE
      CALL SVDBK(NPAN,ITER)
      RETURN
      END
      
      SUBROUTINE SVDBK(N,ITER)
      IMPLICIT REAL*4 (A-H,O-Z)
      COMMON/VG1/U(4000,4000) ! A
      COMMON/SSV/V(4000,4000),W(4000) ! VT , S
      COMMON/VG2/B(4000)
      COMMON/GVV/X(4000)
      DIMENSION ipiv(N), work(N)
      
      if (iter==0) then
         print *, "SGETRF"
         call SGETRF(n,n,U(1:N,1:N),n,ipiv,info)
         print *, "SGETRF DONE"
         if (info.ne.0) stop 'Matrix is numerically singular!'
         ! SGETRI computes the inverse of a matrix using the LU factorization
         ! computed by SGETRF.
         print *, "SGETRI"
         call SGETRI(n,U(1:N,1:N),n,ipiv,work,n,info)
         if (info.ne.0) stop 'Matrix inversion failed!'
         print *, "SGETRI DONE"
      end if
      
      X = matmul(U,B)
      END


      subroutine solve_axb(n, a, b, x)
         implicit none
         integer, intent(in) :: n
         real(4), dimension(n, n), intent(in)    :: a
         real(4), dimension(n)   , intent(in)    :: b
         real(4), dimension(n)   , intent(inout) :: x
       
         integer :: info
         real(4), dimension(n) :: ipiv(n)
         ! Call DGESV to solve the system AX = B. X is saved in B
         x = b
         call sgesv(n,1,a,n,ipiv,x,n,info)
         if (info .ne. 0) then
           print *, 'Error in dgesv: ', info
           ! Handle the error appropriately
         end if
       end subroutine solve_axb



      SUBROUTINE VORTEX
      IMPLICIT REAL*4 (A-H,O-Z)
      COMMON/VVOR/GX,GY,GZ,X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4,
     1VX,VY,VZ
      DIMENSION XX(4),YY(4),ZZ(4)
      PI=3.14159
      ZERO= 0.000000001
      XX(1)=X1
      YY(1)=Y1
      ZZ(1)=Z1
      XX(2)=X2
      YY(2)=Y2
      ZZ(2)=Z2
      XX(3)=X3
      YY(3)=Y3
      ZZ(3)=Z3
      XX(4)=X4
      YY(4)=Y4
      ZZ(4)=Z4
      VX=0.
      VY=0.
      VZ=0.
      DO 80 KK=1,4
      KIK=KK+1
      IF (KK.EQ.4) KIK=1
      ABX=XX(KIK)-XX(KK)
      ABY=YY(KIK)-YY(KK)
      ABZ=ZZ(KIK)-ZZ(KK)
      AB=SQRT(ABX**2+ABY**2+ABZ**2)
      IF(AB.LE.ZERO) GO TO 173
      APX=GX-XX(KK)
      APY=GY-YY(KK)
      APZ=GZ-ZZ(KK)
      AP=SQRT(APX**2+APY**2+APZ**2)
      IF(AP.LE.ZERO) GO TO 173
      COSTH1=(ABX*APX+ABY*APY+ABZ*APZ)/(AB*AP)
      BPX=GX-XX(KIK)
      BPY=GY-YY(KIK) 
      BPZ=GZ-ZZ(KIK)
      BP=SQRT(BPX**2+BPY**2+BPZ**2)
      IF(BP.LE.ZERO) GO TO 173
      COSTH2=-(ABX*BPX+ABY*BPY+ABZ*BPZ)/(AB*BP)
      V1=APY*BPZ-APZ*BPY
      V2=APZ*BPX-APX*BPZ
      V3=APX*BPY-APY*BPX
      H=SQRT(V1**2+V2**2+V3**2)/AB
C      IF(ABS(COSTH1+COSTH2).GT.ZERO.AND.H.GT.ZERO) GO TO 73
      IF(H.GT.ZERO) GO TO 73
  173 VPX=0.
      VPY=0.
      VPZ=0.
      GO TO 72
   73 VP=(COSTH1+COSTH2)/(4*PI*H)
      ABPX=ABY*APZ-ABZ*APY
      ABPY=ABZ*APX-ABX*APZ
      ABPZ=ABX*APY-ABY*APX
      ABP=SQRT(ABPX**2+ABPY**2+ABPZ**2)
      VPX=VP*ABPX/ABP
      VPY=VP*ABPY/ABP
      VPZ=VP*ABPZ/ABP
   72 VX=VX+VPX
      VY=VY+VPY
      VZ=VZ+VPZ
   80 CONTINUE
      RETURN
      END



         
      SUBROUTINE SVDCMP(N)
      IMPLICIT REAL*4 (A-H,O-Z)
      COMMON/VG1/A(4000,4000)
      COMMON/SSV/V(4000,4000),W(4000)
      DIMENSION RV1(4000)
      M=N
      G=0.
      SCALE=0.
      ANORM=0.
      DO 25 I=1,N    
      L=I+1
      RV1(I)=SCALE*G
      G=0.
      S=0.
      SCALE=0.
      IF(I.GT.M) GO TO 161
      DO 11 K=I,M
      SCALE=SCALE+ABS(A(K,I))
   11 CONTINUE
      IF(SCALE.EQ.0.) GO TO 161
      DO 12 K=I,M
      A(K,I)=A(K,I)/SCALE
      S=S+A(K,I)*A(K,I)
   12 CONTINUE
      F=A(I,I)
      G=-SIGN(SQRT(S),F)
      H=F*G-S
      A(I,I)=F-G
      DO 15 J=L,N
      S=0.
      DO 13 K=I,M
      S=S+A(K,I)*A(K,J)
   13 CONTINUE
      F=S/H
      DO 14 K=I,M
      A(K,J)=A(K,J)+F*A(K,I)
   14 CONTINUE
   15 CONTINUE
      DO 16 K=I,M
      A(K,I)=SCALE*A(K,I)
   16 CONTINUE
  161 CONTINUE
      W(I)=SCALE*G
      G=0.
      S=0.
      SCALE=0.
      IF((I.LE.M).AND.(I.NE.N)) GO TO 162
      GO TO 124
  162 DO 17 K=L,N
      SCALE=SCALE+ABS(A(I,K))
   17 CONTINUE
      IF(SCALE.EQ.0.) GO TO 124
      DO 18 K=L,N
      A(I,K)=A(I,K)/SCALE
      S=S+A(I,K)*A(I,K)
   18 CONTINUE
      F=A(I,L)
      G=-SIGN(SQRT(S),F)
      H=F*G-S
      A(I,L)=F-G
      DO 19 K=L,N
      RV1(K)=A(I,K)/H
   19 CONTINUE
      DO 23 J=L,M
      S=0.
      DO 21 K=L,N
      S=S+A(J,K)*A(I,K)
   21 CONTINUE
      DO 22 K=L,N
      A(J,K)=A(J,K)+S*RV1(K)
   22 CONTINUE
   23 CONTINUE
      DO 24 K=L,N
      A(I,K)=SCALE*A(I,K)
   24 CONTINUE
  124 CONTINUE
      ANORM=MAX(ANORM,(ABS(W(I))+ABS(RV1(I))))
   25 CONTINUE
      DO 32 I=N,1,-1  
      IF(I.GE.N) GO TO 131
      IF(G.EQ.0.) GO TO 129
      DO 26 J=L,N
      V(J,I)=(A(I,J)/A(I,L))/G
   26 CONTINUE
      DO 29 J=L,N
      S=0.
      DO 27 K=L,N
      S=S+A(I,K)*V(K,J)
   27 CONTINUE
      DO 28 K=L,N
      V(K,J)=V(K,J)+S*V(K,I)
   28 CONTINUE
   29 CONTINUE
  129 CONTINUE
      DO 31 J=L,N
      V(I,J)=0.
      V(J,I)=0.
   31 CONTINUE
  131 CONTINUE
      V(I,I)=1.
      G=RV1(I)
      L=I
   32 CONTINUE
      DO 39 I=MIN(M,N),1,-1      
      L=I+1
      G=W(I)
      DO 33 J=L,N
      A(I,J)=0.
   33 CONTINUE
      IF(G.EQ.0.) GO TO 137
      G=1./G
      DO 36 J=L,N
      S=0.
      DO 34 K=L,M
      S=S+A(K,I)*A(K,J)
   34 CONTINUE
      F=(S/A(I,I))*G
      DO 35 K=I,M
      A(K,J)=A(K,J)+F*A(K,I)
   35 CONTINUE
   36 CONTINUE
      DO 37 J=I,M
      A(J,I)=A(J,I)*G
   37 CONTINUE
      GO TO 138
  137 DO 38 J=I,M
      A(J,I)=0.
   38 CONTINUE
  138 A(I,I)=A(I,I)+1.
   39 CONTINUE
      DO 49 K=N,1,-1
      DO 48 ITS=1,30
      DO 41 L=K,1,-1
      NM=L-1
      IF((ABS(RV1(L))+ANORM).EQ.ANORM) GO TO 2
      IF((ABS(W(NM))+ANORM).EQ.ANORM) GO TO 1
   41 CONTINUE
    1 C=0.
      S=1.
      DO 43 I=L,K
      F=S*RV1(I)
      RV1(I)=C*RV1(I)
      IF((ABS(F)+ANORM).EQ.ANORM) GO TO 2
      G=W(I)
      H=PYTHAG(F,G)
      W(I)=H
      H=1./H
      C=G*H
      S=-F*H
      DO 42 J=1,M
      Y=A(J,NM)
      Z=A(J,I)
      A(J,NM)=Y*C+Z*S
      A(J,I)=-Y*S+Z*C
   42 CONTINUE
   43 CONTINUE
    2 Z=W(K)
      IF(L.NE.K) GO TO 244
      IF(Z.GE.0.) GO TO 144
      W(K)=-Z
      DO 44 J=1,N
      V(J,K)=-V(J,K)
   44 CONTINUE
  144 CONTINUE
      GO TO 3
  244 CONTINUE
      IF(ITS.EQ.30.) STOP 'NO CONVERGENCE IN SVD'
      X=W(L)
      NM=K-1
      Y=W(NM)
      G=RV1(NM)
      H=RV1(K)
      F=((Y-Z)*(Y+Z)+(G-H)*(G+H))/(2.*H*Y) 
      abc=1.
      G=PYTHAG(F,abc)
      F=((X-Z)*(X+Z)+H*((Y/(F+SIGN(G,F)))-H))/X
      C=1.
      S=1.
      DO 47 J=L,NM
      I=J+1
      G=RV1(I)
      Y=W(I)
      H=S*G
      G=C*G
      Z=PYTHAG(F,H)
      RV1(J)=Z
      C=F/Z
      S=H/Z
      F=X*C+G*S
      G=-X*S+G*C
      H=Y*S
      Y=Y*C
      DO 45 JJ=1,N
      X=V(JJ,J)
      Z=V(JJ,I)
      V(JJ,J)=X*C+Z*S
      V(JJ,I)=-X*S+Z*C
   45 CONTINUE
      Z=PYTHAG(F,H)
      W(J)=Z
      IF(Z.EQ.0.) GO TO 145
      Z=1./Z
      C=F*Z
      S=H*Z
  145 CONTINUE
      F=C*G+S*Y
      X=-S*G+C*Y
      DO 46 JJ=1,M
      Y=A(JJ,J)
      Z=A(JJ,I)
      A(JJ,J)=Y*C+Z*S
      A(JJ,I)=-Y*S+Z*C
   46 CONTINUE
   47 CONTINUE
      RV1(L)=0.
      RV1(K)=F
      W(K)=X
   48 CONTINUE
    3 CONTINUE
   49 CONTINUE
      RETURN
      END


      FUNCTION PYTHAG(A,B)
      IMPLICIT REAL*4 (A-H,O-Z)
      ABSA=ABS(A)
      ABSB=ABS(B)
      IF(ABSA.LE.ABSB) GO TO 10
      PYTHAG=ABSA*SQRT(1.+(ABSB/ABSA)**2)
      GO TO 30
   10 IF(ABSB.NE.0.) GO TO 20
      PYTHAG=0.
      GO TO 30
   20 PYTHAG=ABSB*SQRT(1.+(ABSA/ABSB)**2)
   30 CONTINUE
      RETURN
      END


      SUBROUTINE SVDBK_old(N,ITER)
      IMPLICIT REAL*4 (A-H,O-Z)
      COMMON/VG1/U(4000,4000)
      COMMON/SSV/V(4000,4000),W(4000)
      COMMON/VG2/B(4000)
      COMMON/GVV/X(4000)
      DIMENSION TMP(4000)
      M=N
      IF(ITER.NE.0) GO TO 6
      CALL SVDCMP(N)
      WMAX=0.
      DO 2 J=1,N
      IF(W(J).GT.WMAX) WMAX=W(J)
    2 CONTINUE
      WMIN=WMAX*1.D-12
      DO 3 J=1,N
      IF(W(J).LT.WMIN) W(J)=0.
    3 CONTINUE

    
    6 DO 12 J=1,N
      S=0.
      IF(W(J).EQ.0.) GO TO 111
      DO 11 I=1,M
      S=S+U(I,J)*B(I)
   11 CONTINUE
      S=S/W(J)
  111 CONTINUE
      TMP(J)=S
   12 CONTINUE
      DO 14 J=1,N
      S=0.
      DO 13 JJ=1,N
      S=S+V(J,JJ)*TMP(JJ)
   13 CONTINUE
      X(J)=S
   14 CONTINUE
      RETURN
      END




      SUBROUTINE WAKE(ITER,VINIT,DT,NPAN,NPW,NGW,NGRID)
      IMPLICIT REAL*4 (A-H,O-Z)
      COMMON/GW/XW(60,4000),YW(60,4000),ZW(60,4000)
      COMMON/GAVV/X(4000),Y(4000),Z(4000)
      COMMON/VW/CW(60,4000),ICW(4000,2)
      COMMON/GVV/C(4000)
      COMMON/GRID/IC(4000,4),MARK(4000),MARKW(4000)
      DIMENSION ICHECK(4000)
      COMMON/AVV1/ANX(4000),ANY(4000),ANZ(4000),GX(4000),GY(4000),
     1GZ(4000)
      COMMON/AVV2/ATX(4000),ATY(4000),ATZ(4000),ALX(4000),ALY(4000),
     1ALZ(4000)
      COMMON/ANPORT/DS(4000)
      COMMON/CORECT/ X1X(4000),Y1Y(4000),Z1Z(4000)
      COMMON/WAICOR/NC(4000)
      COMMON/WIN/IPW(4000)
      COMMON/WINTER/IPROPW(4000)
      COMMON/PPG/IPROP(4000)
      DO 5 I=1,NGRID
      MARKW(I)=0
    5 CONTINUE
      DO 7 I=1,NPAN
      ICHECK(I)=0
    7 CONTINUE
      IT1=ITER+1
      II=0
      NPW=0
      DO 20 J=1,NPAN
      IF (ICHECK(J).EQ.1) GO TO 20
      MK=4
      IF(IC(J,3).EQ.IC(J,4)) MK=3
      DO 10 K=1,MK
      K1=K
      K2=K+1
      IF (K.NE.MK) GO TO 11
      K1=1
      K2=MK
   11 IF((MARK(IC(J,K1)).NE.0).AND.(MARK(IC(J,K2)).NE.0)) GO TO 21
      GOTO 10
   21 IF((MARKW(IC(J,K1)).NE.0).AND.(MARKW(IC(J,K2)).NE.0)) GOTO 10
      GOTO 31
   10 CONTINUE
      GO TO 20
   31 DO 17 L=J+1,NPAN
      LPAN=0
      KK1=0
      KK2=0
      MKK=4
      IF(IC(L,3).EQ.IC(L,4)) MKK=3
      DO 16 M=1,MKK
      IF (IC(L,M).EQ.IC(J,K1)) KK1=M
      IF (IC(L,M).EQ.IC(J,K2)) KK2=M
   16 CONTINUE
      IF ((KK1*KK2).EQ.0) GOTO 17
      LPAN=L
      ICHECK(LPAN)=1
      GOTO 18
   17 CONTINUE
   18 CONTINUE
      IF(MARKW(IC(J,K1)).NE.0) GO TO 15
      II=II+1
      IPW(II)=1
      IF(IPROP(J).EQ.2) IPW(II)=2
      MARKW(IC(J,K1))=II
      X1X(II)=X(IC(J,K1))
      Y1Y(II)=Y(IC(J,K1))
      Z1Z(II)=Z(IC(J,K1))
      NC(II)=IC(J,K1)
      XW(IT1,II)=X(IC(J,K1))
      YW(IT1,II)=Y(IC(J,K1))
      ZW(IT1,II)=Z(IC(J,K1))
   15 IF(MARKW(IC(J,K2)).NE.0) GO TO 19
      II=II+1
      IPW(II)=1
      IF(IPROP(J).EQ.2) IPW(II)=2
      MARKW(IC(J,K2))=II
      X1X(II)=X(IC(J,K2))
      Y1Y(II)=Y(IC(J,K2))
      Z1Z(II)=Z(IC(J,K2))
      NC(II)=IC(J,K2)
      XW(IT1,II)=X(IC(J,K2))
      YW(IT1,II)=Y(IC(J,K2))
      ZW(IT1,II)=Z(IC(J,K2))
   19 NPW=NPW+1
      CW(ITER,NPW)=C(J)
      IF ((K2-K1).EQ.(MK-1)) CW(ITER,NPW)=-CW(ITER,NPW)
      WRITE(21,*) 'J,C(J)',J,C(J)
      WRITE(21,*) IC(J,1),IC(J,2),IC(J,3),IC(J,4)
      WRITE(21,*) 'K1,K2,NPW,CW',K1,K2,NPW,CW(ITER,NPW)
      IF (LPAN.EQ.0) GOTO 90
      IF ((KK2.LT.KK1).AND.((KK1-KK2).EQ.1)) CW(ITER,NPW)=CW(ITER,NPW)
     *-C(LPAN)
      IF ((KK2.LT.KK1).AND.((KK1-KK2).EQ.(MKK-1))) CW(ITER,NPW)=
     *CW(ITER,NPW)+C(LPAN)
      IF ((KK2.GT.KK1).AND.((KK2-KK1).EQ.1)) CW(ITER,NPW)=CW(ITER,NPW)
     *+C(LPAN)
      IF ((KK2.GT.KK1).AND.((KK2-KK1).EQ.(MKK-1))) CW(ITER,NPW)=
     *CW(ITER,NPW)-C(LPAN)
      WRITE(21,*) 'LPAN,C(LPAN)',LPAN,C(LPAN)
      WRITE(21,*) IC(LPAN,1),IC(LPAN,2),IC(LPAN,3),IC(LPAN,4)
      WRITE(21,*) 'KK1,KK2,CW',KK1,KK2,CW(ITER,NPW)
   90 WRITE(21,*)
      ICW(NPW,1)=MARKW(IC(J,K1))
      ICW(NPW,2)=MARKW(IC(J,K2))
      IPROPW(NPW)=0
      IF(IPROP(J).EQ.2) IPROPW(NPW)=1
   20 CONTINUE
      NGW=II
      DO 30 I=1,NPAN
      IF (ICHECK(I).EQ.1) GO TO 30
      MK=4
      IF(IC(I,3).EQ.IC(I,4)) MK=3
      DO 40 K=1,MK
      K1=K
      K2=K+1
      IF (K.NE.MK) GO TO 35
      K1=1
      K2=MK
   35 IF (MARKW(IC(I,K1))*MARKW(IC(I,K2)).NE.0) GO TO 38
      GOTO 40
   38 DO 37 L=1,NPW
      IF ((ICW(L,1).EQ.MARKW(IC(I,K1))).AND.(ICW(L,2).EQ.MARKW
     1(IC(I,K2)))) GO TO 40
   37 CONTINUE
      NPW=NPW+1
      DO 117 LL=I+1,NPAN
      LPAN=0
      KK1=0
      KK2=0
      MKK=4
      IF(IC(LL,3).EQ.IC(LL,4)) MKK=3
      DO 116 M=1,MKK
      IF (IC(LL,M).EQ.IC(I,K1)) KK1=M
      IF (IC(LL,M).EQ.IC(I,K2)) KK2=M
  116 CONTINUE
      IF (KK1*KK2.EQ.0) GOTO 117
      LPAN=LL
      ICHECK(LPAN)=1
      GOTO 118
  117 CONTINUE
  118 CONTINUE
      CW(ITER,NPW)=C(I)
      IF ((K2-K1).EQ.(MK-1)) CW(ITER,NPW)=-CW(ITER,NPW)
      WRITE(21,*) 'J,C(J)',I,C(I)
      WRITE(21,*) IC(I,1),IC(I,2),IC(I,3),IC(I,4)
      WRITE(21,*) 'K1,K2,NPW,CW',K1,K2,NPW,CW(ITER,NPW)
      IF (LPAN.EQ.0) GOTO 190
      IF ((KK2.LT.KK1).AND.((KK1-KK2).EQ.1)) CW(ITER,NPW)=CW(ITER,NPW)
     *-C(LPAN)
      IF ((KK2.LT.KK1).AND.((KK1-KK2).EQ.(MKK-1))) CW(ITER,NPW)=
     *CW(ITER,NPW)+C(LPAN)
      IF ((KK2.GT.KK1).AND.((KK2-KK1).EQ.1)) CW(ITER,NPW)=CW(ITER,NPW)
     *+C(LPAN)
      IF ((KK2.GT.KK1).AND.((KK2-KK1).EQ.(MKK-1))) CW(ITER,NPW)=
     *CW(ITER,NPW)-C(LPAN)
      WRITE(21,*) 'LPAN,C(LPAN)',LPAN,C(LPAN)
      WRITE(21,*) IC(LPAN,1),IC(LPAN,2),IC(LPAN,3),IC(LPAN,4)
      WRITE(21,*) 'KK1,KK2,CW',KK1,KK2,CW(ITER,NPW)
  190 WRITE(21,*)
      ICW(NPW,1)=MARKW(IC(I,K1))
      ICW(NPW,2)=MARKW(IC(I,K2))
      IPROPW(NPW)=0
      IF(IPROP(I).EQ.2) IPROPW(NPW)=1
   40 CONTINUE
   30 CONTINUE
      DO 121 KL=1,NGW
      X1=X1X(KL)
      Y1=Y1Y(KL)
      Z1=Z1Z(KL)
      XW(ITER,KL)=X1+VINIT*DT
      YW(ITER,KL)=Y1
      ZW(ITER,KL)=Z1
  121 CONTINUE
C      CALL IWCOR(ITER,NGW,NPAN)
      IF(ITER.EQ.1) GO TO 124
      DO 123 I=1,ITER-1
      DO 122 J=1,NGW
      XW(I,J)=XW(I,J)+VINIT*DT
      YW(I,J)=YW(I,J)
      ZW(I,J)=ZW(I,J)
  122 CONTINUE
  123 CONTINUE
  124 CONTINUE
      RETURN
      END



 


      SUBROUTINE WAKCOR(NPAN,NGW,ITER)
      IMPLICIT REAL*4 (A-H,O-Z)
      COMMON/AVV1/ANX(4000),ANY(4000),ANZ(4000),GX(4000),GY(4000),
     1GZ(4000)
      COMMON/AVV2/ATX(4000),ATY(4000),ATZ(4000),ALX(4000),ALY(4000),
     1ALZ(4000)
      COMMON/GAVV/X(4000),Y(4000),Z(4000)
      COMMON/GRID/IC(4000,4),MARK(4000),MARKW(4000)
      COMMON/GW/XW(60,4000),YW(60,4000),ZW(60,4000)
      COMMON/ANPORT/DS(4000)
      COMMON/CORECT/ X1X(4000),Y1Y(4000),Z1Z(4000)
      COMMON/WAICOR/NC(4000)
      NP1=0
      DO 130 J=1,NGW
      DO 129 I=1,ITER
      II=ITER+1-I
      XXW=XW(II,J)
      YYW=YW(II,J)
      ZZW=ZW(II,J)
      XXW1=XW(II+1,J)
      YYW1=YW(II+1,J)
      ZZW1=ZW(II+1,J)
      NCC=NC(J)
      CALL PENETR(NPAN,XXW,YYW,ZZW,XXW1,YYW1,ZZW1,NP1,NCC)
      XW(II,J)=XXW
      YW(II,J)=YYW
      ZW(II,J)=ZZW
  129 CONTINUE
  130 CONTINUE
      RETURN
      END



 
      SUBROUTINE CPAIP(ITER,NPAN,NPW,NSYM,VINIT,NGRND) 
      IMPLICIT REAL*4 (A-H,O-Z)
      COMMON/VVOR/XG,YG,ZG,X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4,VVX,
     1VVY,VVZ
      COMMON/AVV1/ANX(4000),ANY(4000),ANZ(4000),GX(4000),GY(4000),
     1GZ(4000)
      COMMON/ANPORT/DS(4000)
      COMMON/VW/CW(60,4000),ICW(4000,2)
      COMMON/GAVV/X(4000),Y(4000),Z(4000)
      COMMON/GW/XW(60,4000),YW(60,4000),ZW(60,4000)
      COMMON/GVV/C(4000)
      COMMON/GRID/IC(4000,4),MARK(4000),MARKW(4000)
      COMMON/PPG/IPROP(4000)
      DIMENSION XVV(4000),YVV(4000),ZVV(4000)
      DO 80 I=1,NPAN
      IF(IPROP(I).NE.2) GO TO 80
      XGP=GX(I)
      YGP=GY(I)
      ZGP=GZ(I)
      XG=XGP+1.3*ANX(I)
      YG=YGP+1.3*ANY(I)
      ZG=ZGP+1.3*ANZ(I)
      XVV(I)=0.
      YVV(I)=0.
      ZVV(I)=0.
      NG=0
   20 DO 30 IK=1,NPAN
      NA=0
      X1=X(IC(IK,1))
      Y1=Y(IC(IK,1))
      Z1=Z(IC(IK,1))
      X2=X(IC(IK,2))
      Y2=Y(IC(IK,2))
      Z2=Z(IC(IK,2))
      X3=X(IC(IK,3))
      Y3=Y(IC(IK,3))
      Z3=Z(IC(IK,3))
      X4=X(IC(IK,4))
      Y4=Y(IC(IK,4))
      Z4=Z(IC(IK,4))
      IF(NG.EQ.0) GO TO 31
      Z1=-Z1
      Z2=-Z2
      Z3=-Z3
      Z4=-Z4
   31 CALL VORTEX
      XVV(I)=XVV(I)+VVX*C(IK)*((-1)**NA)*(-1)**NG
      YVV(I)=YVV(I)+VVY*C(IK)*((-1)**NA)*(-1)**NG
      ZVV(I)=ZVV(I)+VVZ*C(IK)*((-1)**NA)*(-1)**NG
      IF(NSYM.EQ.0) GO TO 30
      IF(NA.EQ.1) GO TO 30
      NA=1
      Y1=-Y1
      Y2=-Y2
      Y3=-Y3
      Y4=-Y4
      GO TO 31
   30 CONTINUE
      IF(NGRND.EQ.0) GO TO 40
      IF(NG.EQ.1) GO TO 40
      NG=1
      GO TO 20
   40 NG=0
   41 DO 60 IK=1,ITER
      DO 50 JK=1,NPW
      NA=0
      X1=XW(IK,ICW(JK,1))
      Y1=YW(IK,ICW(JK,1))
      Z1=ZW(IK,ICW(JK,1))
      X2=XW(IK,ICW(JK,2))
      Y2=YW(IK,ICW(JK,2))
      Z2=ZW(IK,ICW(JK,2))
      X3=XW(IK+1,ICW(JK,2))
      Y3=YW(IK+1,ICW(JK,2))
      Z3=ZW(IK+1,ICW(JK,2))
      X4=XW(IK+1,ICW(JK,1))
      Y4=YW(IK+1,ICW(JK,1))
      Z4=ZW(IK+1,ICW(JK,1))
      IF(NG.EQ.0) GO TO 51
      Z1=-Z1
      Z2=-Z2
      Z3=-Z3
      Z4=-Z4
   51 CALL VORTEX
      XVV(I)=XVV(I)+VVX*CW(IK,JK)*((-1)**NA)*(-1)**NG
      YVV(I)=YVV(I)+VVY*CW(IK,JK)*((-1)**NA)*(-1)**NG
      ZVV(I)=ZVV(I)+VVZ*CW(IK,JK)*((-1)**NA)*(-1)**NG
      IF(NSYM.EQ.0) GO TO 50
      IF(NA.EQ.1) GO TO 50
      NA=1
      Y1=-Y1
      Y2=-Y2
      Y3=-Y3
      Y4=-Y4
      GO TO 51
   50 CONTINUE
   60 CONTINUE
      IF(NGRND.EQ.0) GO TO 70
      IF(NG.EQ.1) GO TO 70
      NG=1
      GO TO 41
   70 CONTINUE              
   80 CONTINUE  
      DO 90 I=1,NPAN
      IF(IPROP(I).NE.2) GO TO 90
      VRX=XVV(I)+VINIT
      VRY=YVV(I)
      VRZ=ZVV(I)
      VR=VRX**2+VRY**2+VRZ**2
      CAIP=1.-VR/(VINIT**2)
      WRITE(41,*) I,CAIP,VR
   90 CONTINUE
      RETURN
      END




      SUBROUTINE PENETR(NPAN,XW,YW,ZW,XW1,YW1,ZW1,NP1,NCC)
      IMPLICIT REAL*4 (A-H,O-Z)
      COMMON/AVV1/ANX(4000),ANY(4000),ANZ(4000),GX(4000),GY(4000),
     1GZ(4000)
      COMMON/AVV2/ATX(4000),ATY(4000),ATZ(4000),ALX(4000),ALY(4000),
     1ALZ(4000)
      COMMON/GAVV/X(4000),Y(4000),Z(4000)
      COMMON/GRID/IC(4000,4),MARK(4000),MARKW(4000)
      DIMENSION XX(4),YY(4),ZZ(4),PX(4),PY(4),PV(4)
  124 INT=0
      DO 128 K=1,NPAN
      INTR=0
      IF(K.EQ.NP1) GO TO 131
      KKI=4
      IF(IC(K,3).EQ.IC(K,4)) KKI=3
      DO 123 IK=1,KKI
      IF(IC(K,IK).EQ.NCC) GO TO 131
  123 CONTINUE
      INTR=1                         
      PGX=XW-GX(K)
      PGY=YW-GY(K)
      PGZ=ZW-GZ(K)
      WT=PGX*ATX(K)+PGY*ATY(K)+PGZ*ATZ(K)
      WL=PGX*ALX(K)+PGY*ALY(K)+PGZ*ALZ(K)
      WN=PGX*ANX(K)+PGY*ANY(K)+PGZ*ANZ(K)
      PGX1=XW1-GX(K)
      PGY1=YW1-GY(K)
      PGZ1=ZW1-GZ(K)
      WT1=PGX1*ATX(K)+PGY1*ATY(K)+PGZ1*ATZ(K)
      WL1=PGX1*ALX(K)+PGY1*ALY(K)+PGZ1*ALZ(K)
      WN1=PGX1*ANX(K)+PGY1*ANY(K)+PGZ1*ANZ(K)
c      IF(WN1.LT.0.) GO TO 127
c      IF((WN1.GT.0.).AND.(WN.GT.0.)) GO TO 127
      IF((WN1*WN).GT.0.OR.WN1.EQ.WN) GO TO 127
      XC=(WN1*WT-WN*WT1)/(WN1-WN)
      YC=(WN*WL1-WN1*WL)/(WN-WN1)
      DO 125 L=1,KKI
      XX(L)=X(IC(K,L))
      YY(L)=Y(IC(K,L))
      ZZ(L)=Z(IC(K,L))
      XP=XX(L)-GX(K)
      YP=YY(L)-GY(K)
      ZP=ZZ(L)-GZ(K)
      PX(L)=XP*ATX(K)+YP*ATY(K)+ZP*ATZ(K)
      PY(L)=XP*ALX(K)+YP*ALY(K)+ZP*ALZ(K)
  125 CONTINUE
      DO 126 L=1,KKI
      LL=L+1
      IF(L.EQ.KKI) LL=1
      XS=PX(LL)-PX(L)
      YS=PY(LL)-PY(L)
      XS1=PX(L)-XC
      YS1=PY(L)-YC
      PV(L)=XS*YS1-YS*XS1
  126 CONTINUE
      DO 132 L=1,KKI
      LL=L+1
      IF(L.EQ.KKI) LL=1
      PR=PV(L)*PV(LL)
      IF(PR.LT.0.) GO TO 127
  132 CONTINUE
      X1=XW-XW1
      Y1=YW-YW1
      Z1=ZW-ZW1
      DIST=-X1*ANX(K)-Y1*ANY(K)-Z1*ANZ(K)
      XW=XW+DIST*ANX(K)
      YW=YW+DIST*ANY(K)
      ZW=ZW+DIST*ANZ(K)
      GO TO 131
  127 INTR=0
  131 INT=INT+INTR
  128 CONTINUE        
      IF(INT.NE.0) GO TO 124
      RETURN
      END
      
      
      SUBROUTINE WAKINT(NPW,ITER)
      IMPLICIT REAL*4 (A-H,O-Z)
      COMMON/GW/XW(60,4000),YW(60,4000),ZW(60,4000)
      COMMON/VW/CW(60,4000),ICW(4000,2)
      COMMON/GRID/IC(4000,4),MARK(4000),MARKW(4000)
      COMMON/WINTER/IPROPW(4000)
      DIMENSION WX(60,4000),WY(60,4000),WZ(60,4000) 
      DO 9 I=1,ITER+1
      DO 5 J=1,NPW
      WX(I,ICW(J,1))=XW(I,ICW(J,1))
      WY(I,ICW(J,1))=YW(I,ICW(J,1))
      WZ(I,ICW(J,1))=ZW(I,ICW(J,1))
      WX(I,ICW(J,2))=XW(I,ICW(J,2))
      WY(I,ICW(J,2))=YW(I,ICW(J,2))
      WZ(I,ICW(J,2))=ZW(I,ICW(J,2))
    5 CONTINUE
    9 CONTINUE
      DO 30 II=1,ITER
      IK=ITER+1-II
      DO 20 I=1,NPW
      IF(IPROPW(I).NE.0) GO TO 20
      X1=WX(IK,ICW(I,1))   
      Y1=WY(IK,ICW(I,1)) 
      Z1=WZ(IK,ICW(I,1))   
      X2=WX(IK,ICW(I,2))   
      Y2=WY(IK,ICW(I,2))   
      Z2=WZ(IK,ICW(I,2))   
      X3=WX(IK+1,ICW(I,2))   
      Y3=WY(IK+1,ICW(I,2))   
      Z3=WZ(IK+1,ICW(I,2))   
      X4=WX(IK+1,ICW(I,1))   
      Y4=WY(IK+1,ICW(I,1))   
      Z4=WZ(IK+1,ICW(I,1)) 
      DO 12 IJ=1,ITER
      IJ1=ITER+1-IJ 
      DO 11 J=1,NPW
      IF(IPROPW(J).NE.1) GO TO 11
      XX1=WX(IJ1,ICW(J,1))   
      YY1=WY(IJ1,ICW(J,1))   
      ZZ1=WZ(IJ1,ICW(J,1))   
      XX2=WX(IJ1,ICW(J,2))   
      YY2=WY(IJ1,ICW(J,2))   
      ZZ2=WZ(IJ1,ICW(J,2))   
      XX3=WX(IJ1+1,ICW(J,2))   
      YY3=WY(IJ1+1,ICW(J,2))   
      ZZ3=WZ(IJ1+1,ICW(J,2))   
      XX4=WX(IJ1+1,ICW(J,1))   
      YY4=WY(IJ1+1,ICW(J,1))   
      ZZ4=WZ(IJ1+1,ICW(J,1)) 
      CALL CROSS(XX3,YY3,ZZ3,XX2,YY2,ZZ2,X1,Y1,Z1,X2,Y2,Z2,
     1X3,Y3,Z3,X4,Y4,Z4,ANEWXW2,ANEWYW2,ANEWZW2)
      CALL CROSS(XX4,YY4,ZZ4,XX1,YY1,ZZ1,X1,Y1,Z1,X2,Y2,Z2,
     1X3,Y3,Z3,X4,Y4,Z4,ANEWXW1,ANEWYW1,ANEWZW1)
      WX(IJ1,ICW(J,1))=ANEWXW1
      WY(IJ1,ICW(J,1))=ANEWYW1
      WZ(IJ1,ICW(J,1))=ANEWZW1
      WX(IJ1,ICW(J,2))=ANEWXW2
      WY(IJ1,ICW(J,2))=ANEWYW2
      WZ(IJ1,ICW(J,2))=ANEWZW2
   11 CONTINUE
   12 CONTINUE
   20 CONTINUE
   30 CONTINUE
      DO 50 JK=1,ITER
      DO 40 J=1,NPW
      IF(IPROPW(J).NE.1) GO TO 40
      XW(JK,ICW(J,1))=WX(JK,ICW(J,1))
      YW(JK,ICW(J,1))=WY(JK,ICW(J,1)) 
      ZW(JK,ICW(J,1))=WZ(JK,ICW(J,1))
      XW(JK,ICW(J,2))=WX(JK,ICW(J,2))
      YW(JK,ICW(J,2))=WY(JK,ICW(J,2)) 
      ZW(JK,ICW(J,2))=WZ(JK,ICW(J,2))
   40 CONTINUE
   50 CONTINUE
      RETURN
      END
      
      
      SUBROUTINE CROSS(XXA,YYA,ZZA,XXB,YYB,ZZB,X1,Y1,Z1,X2,Y2,Z2,
     1X3,Y3,Z3,X4,Y4,Z4,ANEWXW,ANEWYW,ANEWZW)
      IMPLICIT REAL*4 (A-H,O-Z)
      DIMENSION XI(4),YI(4),ZI(4),PR(4)
      ABX=XXB-XXA
      ABY=YYB-YYA
      ABZ=ZZB-ZZA
      ANORM=SQRT(ABX**2+ABY**2+ABZ**2)
      ANEWXW=XXB
      ANEWYW=YYB
      ANEWZW=ZZB
      XI(1)=X1
      YI(1)=Y1
      ZI(1)=Z1
      XI(2)=X2
      YI(2)=Y2
      ZI(2)=Z2
      XI(3)=X3
      YI(3)=Y3
      ZI(3)=Z3
      XI(4)=X4
      YI(4)=Y4
      ZI(4)=Z4
      DO 10 I=1,4
      II=I+1
      IF(I.EQ.4) II=1
      DAX1=XI(I)-XXA
      DAY1=YI(I)-YYA
      DAZ1=ZI(I)-ZZA
      DAX2=XI(II)-XXA
      DAY2=YI(II)-YYA
      DAZ2=ZI(II)-ZZA
      PRAX12=DAY1*DAZ2-DAZ1*DAY2
      PRAY12=DAZ1*DAX2-DAX1*DAZ2
      PRAZ12=DAX1*DAY2-DAY1*DAX2
      PR(I)=ABX*PRAX12+ABY*PRAY12+ABZ*PRAZ12
      IF(I.EQ.3.OR.I.EQ.4) PR(I)=-PR(I)
   10 CONTINUE
      DXA=X1-XXA
      DYA=Y1-YYA
      DZA=Z1-ZZA
      DXB=X1-XXB
      DYB=Y1-YYB
      DZB=Z1-ZZB
      PRA=ABX*DXA+ABY*DYA+ABZ*DZA
      PRB=ABX*DXB+ABY*DYB+ABZ*DZB
      IF((PRA*PRB).GT.0.) GO TO 20
      PR13=PR(1)*PR(3)
      PR24=PR(2)*PR(4)
      IF(PR13.GT.0.OR.PR24.GT.0.) GO TO 20
      DX1=X3-X1
      DY1=Y3-Y1
      DZ1=Z3-Z1
      DX2=X2-X4
      DY2=Y2-Y4
      DZ2=Z2-Z4
      PRX=DY1*DZ2-DZ1*DY2
      PRY=DZ1*DX2-DX1*DZ2
      PRZ=DX1*DY2-DY1*DX2
      WNORM=SQRT(PRX**2+PRY**2+PRZ**2)
      WNX=PRX/WNORM
      WNY=PRY/WNORM
      WNZ=PRZ/WNORM
      ABPR=WNX*ABX+WNY*ABY+WNZ*ABZ
      ANEWXW=XXB-WNX*ABPR
      ANEWYW=YYB-WNY*ABPR
      ANEWZW=ZZB-WNZ*ABPR
   20 CONTINUE
      RETURN
      END
      
     
    
      SUBROUTINE NEIBORG(NPAN,NGW)    
      IMPLICIT REAL*4 (A-H,O-Z)
      COMMON/GRID/IC(4000,4),MARK(4000),MARKW(4000)
      COMMON/NEIB/KN(4000,4000)
      COMMON/WAICOR/NC(4000)
      DO 6 I=1,NPAN
      KKI=4
      IF(IC(I,3).EQ.IC(I,4)) KKI=3
      DO 5 IK=1,KKI
      KN(I,IC(I,IK))=0
    5 CONTINUE
    6 CONTINUE
      DO 30 L=1,NGW
      DO 20 I=1,NPAN
      KKI=4
      IF(IC(I,3).EQ.IC(I,4)) KKI=3
      DO 10 IK=1,KKI
      IF(IC(I,IK).NE.NC(L)) GO TO 10
      KN(I,IC(I,IK))=NC(L)
   10 CONTINUE
   20 CONTINUE
   30 CONTINUE
      RETURN
      END


      SUBROUTINE IWCOR(ITER,NGW,NPAN)
      IMPLICIT REAL*4 (A-H,O-Z)
      COMMON/AVV1/ANX(4000),ANY(4000),ANZ(4000),GX(4000),GY(4000),
     1GZ(4000)
      COMMON/AVV2/ATX(4000),ATY(4000),ATZ(4000),ALX(4000),ALY(4000),
     1ALZ(4000)
      COMMON/GRID/IC(4000,4),MARK(4000),MARKW(4000)
      COMMON/GAVV/X(4000),Y(4000),Z(4000)
      COMMON/GW/XW(60,4000),YW(60,4000),ZW(60,4000)
      COMMON/WAICOR/NC(4000)
      COMMON/NEIB/KN(4000,4000)
      DIMENSION ICHECK(4000)
      DO 80 J=1,NGW
      DX=XW(ITER,J)-X(NC(J))
      DY=YW(ITER,J)-Y(NC(J))
      DZ=ZW(ITER,J)-Z(NC(J))
      DO 7 I=1,NPAN
      ICHECK(I)=0
      KKI=4
      IF(IC(I,3).EQ.IC(I,4)) KKI=3
      DO 5 IK=1,KKI
      IF(KN(I,IC(I,IK)).EQ.NC(J)) GO TO 6
    5 CONTINUE
      GO TO 7
    6 PR=DX*ANX(I)+DY*ANY(I)+DZ*ANZ(I)
      IF(PR.GT.0.) GO TO 80
      ICHECK(I)=1
    7 CONTINUE
      DO 70 I=1,NPAN
      IF(ICHECK(I).EQ.0) GO TO 70
      KKI=4
      IF(IC(I,3).EQ.IC(I,4)) KKI=3
      IF (IK == 2) THEN
         AX = X(IC(I, KKI)) - X(IC(I, IK))
         AY = Y(IC(I, KKI)) - Y(IC(I, IK))
         AZ = Z(IC(I, KKI)) - Z(IC(I, IK))
         BX = X(IC(I, 2)) - X(IC(I, IK))
         BY = Y(IC(I, 2)) - Y(IC(I, IK))
         BZ = Z(IC(I, 2)) - Z(IC(I, IK))
         GO TO 60
      ELSE IF (IK == 1) THEN
         AX = X(IC(I, 1)) - X(IC(I, IK))
         AY = Y(IC(I, 1)) - Y(IC(I, IK))
         AZ = Z(IC(I, 1)) - Z(IC(I, IK))
         BX = X(IC(I, 3)) - X(IC(I, IK))
         BY = Y(IC(I, 3)) - Y(IC(I, IK))
         BZ = Z(IC(I, 3)) - Z(IC(I, IK))
         GO TO 60
      ELSE IF (IK /= 3) THEN
         GO TO 40
      END IF      
      AX=X(IC(I,2))-X(IC(I,IK))
      AY=Y(IC(I,2))-Y(IC(I,IK))
      AZ=Z(IC(I,2))-Z(IC(I,IK))
      BX=X(IC(I,1))-X(IC(I,IK))
      BY=Y(IC(I,1))-Y(IC(I,IK))
      BZ=Z(IC(I,1))-Z(IC(I,IK))
      GO TO 60
   40 IF(IK.NE.3) GO TO 50
      AX=X(IC(I,2))-X(IC(I,IK))
      AY=Y(IC(I,2))-Y(IC(I,IK))
      AZ=Z(IC(I,2))-Z(IC(I,IK))
      BX=X(IC(I,4))-X(IC(I,IK))
      BY=Y(IC(I,4))-Y(IC(I,IK))
      BZ=Z(IC(I,4))-Z(IC(I,IK))
      GO TO 60
   50 AX=X(IC(I,3))-X(IC(I,IK))
      AY=Y(IC(I,3))-Y(IC(I,IK))
      AZ=Z(IC(I,3))-Z(IC(I,IK))
      BX=X(IC(I,1))-X(IC(I,IK))
      BY=Y(IC(I,1))-Y(IC(I,IK))
      BZ=Z(IC(I,1))-Z(IC(I,IK))
   60 PRADX=AY*DZ-AZ*DY
      PRADY=AZ*DX-AX*DZ
      PRADZ=AX*DY-AY*DX
      PRABX=AY*BZ-AZ*BY
      PRABY=AZ*BX-AX*BZ
      PRABZ=AX*BY-AY*BX
      PRODV=PRADX*PRABX+PRADY*PRABY+PRADZ*PRABZ
      IF(PRODV.LT.0.) GO TO 70
      IF(PR.EQ.0.) PR=-0.1
      XW(ITER,J)=XW(ITER,J)-1.1*PR*ANX(I)
      YW(ITER,J)=YW(ITER,J)-1.1*PR*ANY(I)
      ZW(ITER,J)=ZW(ITER,J)-1.1*PR*ANZ(I)
      NP1=I
      XX=XW(ITER,J)
      YY=YW(ITER,J)
      ZZ=ZW(ITER,J)
      XX1=X(NC(J))
      YY1=Y(NC(J))
      ZZ1=Z(NC(J))
      NCC=NC(J)
      CALL PENETR(NPAN,XX,YY,ZZ,XX1,YY1,ZZ1,NP1,NCC)
      XW(ITER,J)=XX
      YW(ITER,J)=YY
      ZW(ITER,J)=ZZ
   70 CONTINUE
   80 CONTINUE
      RETURN
      END


      
      SUBROUTINE VELWAK(ITER,NPW,NSYM,NGRND,GX,GY,GZ,XVV,
     1YVV,ZVV)
      IMPLICIT REAL*4 (A-H,O-Z)
      COMMON/VVOR/XG,YG,ZG,X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4,VVX,
     1VVY,VVZ
      COMMON/VW/CW(60,4000),ICW(4000,2)
      COMMON/GW/WX(60,4000),WY(60,4000),WZ(60,4000)
      NG=0
      XVV=0.
      YVV=0.
      ZVV=0.
      XG=GX
      YG=GY
      ZG=GZ
   32 DO 60 IK=1,ITER
      DO 50 JK=1,NPW
      NA=0
      X1=WX(IK,ICW(JK,1))
      Y1=WY(IK,ICW(JK,1))
      Z1=WZ(IK,ICW(JK,1))
      X2=WX(IK,ICW(JK,2))
      Y2=WY(IK,ICW(JK,2))
      Z2=WZ(IK,ICW(JK,2))
      X3=WX(IK+1,ICW(JK,2))
      Y3=WY(IK+1,ICW(JK,2))
      Z3=WZ(IK+1,ICW(JK,2))
      X4=WX(IK+1,ICW(JK,1))
      Y4=WY(IK+1,ICW(JK,1))
      Z4=WZ(IK+1,ICW(JK,1))
      IF(NG.EQ.0) GO TO 51
      Z1=-Z1
      Z2=-Z2
      Z3=-Z3
      Z4=-Z4
   51 CALL VORTEX
      XVV=XVV+VVX*CW(IK,JK)*((-1)**NA)*(-1)**NG
      YVV=YVV+VVY*CW(IK,JK)*((-1)**NA)*(-1)**NG
      ZVV=ZVV+VVZ*CW(IK,JK)*((-1)**NA)*(-1)**NG
      IF(NSYM.EQ.0) GO TO 50
      IF(NA.EQ.1) GO TO 50
      NA=1
      Y1=-Y1
      Y2=-Y2
      Y3=-Y3
      Y4=-Y4
      GO TO 51
   50 CONTINUE
   60 CONTINUE
      IF(NGRND.EQ.0) GO TO 70
      IF(NG.EQ.1) GO TO 70
      NG=1
      GO TO 32
   70 CONTINUE
      RETURN
      END


      
      SUBROUTINE VELPAN(NPAN,NSYM,NGRND,GX,GY,GZ,XVV,YVV,ZVV) 
      IMPLICIT REAL*4 (A-H,O-Z)
      COMMON/VVOR/XG,YG,ZG,X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4,VVX,
     1VVY,VVZ
      COMMON/GAVV/X(4000),Y(4000),Z(4000)
      COMMON/GVV/C(4000)
      COMMON/GRID/IC(4000,4),MARK(4000),MARKW(4000)
      XG=GX
      YG=GY
      ZG=GZ
      XVV=0.
      YVV=0.
      ZVV=0.
      NG=0
  501 DO 30 IK=1,NPAN
      NA=0
      X1=X(IC(IK,1))
      Y1=Y(IC(IK,1))
      Z1=Z(IC(IK,1))
      X2=X(IC(IK,2))
      Y2=Y(IC(IK,2))
      Z2=Z(IC(IK,2))
      X3=X(IC(IK,3))
      Y3=Y(IC(IK,3))
      Z3=Z(IC(IK,3))
      X4=X(IC(IK,4))
      Y4=Y(IC(IK,4))
      Z4=Z(IC(IK,4))
      IF(NG.EQ.0) GO TO 31
      Z1=-Z1
      Z2=-Z2
      Z3=-Z3
      Z4=-Z4
   31 CALL VORTEX
      XVV=XVV+VVX*C(IK)*((-1)**NA)*(-1)**NG
      YVV=YVV+VVY*C(IK)*((-1)**NA)*(-1)**NG
      ZVV=ZVV+VVZ*C(IK)*((-1)**NA)*(-1)**NG
      IF(NSYM.EQ.0) GO TO 30
      IF(NA.EQ.1) GO TO 30
      NA=1
      Y1=-Y1
      Y2=-Y2
      Y3=-Y3
      Y4=-Y4
      GO TO 31
   30 CONTINUE
      IF(NGRND.EQ.0) GO TO 311
      IF(NG.EQ.1) GO TO 311
      NG=1
      GO TO 501
  311 CONTINUE
      RETURN
      END

      

      SUBROUTINE WAKREL(ITER,NPAN,NGW,NPW,DT,NSYM,NGRND) 
      IMPLICIT REAL*4 (A-H,O-Z)
      COMMON/VVOR/XG,YG,ZG,X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4,VVX,
     1VVY,VVZ
      COMMON/AVV1/ANX(4000),ANY(4000),ANZ(4000),GX(4000),GY(4000),
     1GZ(4000)
      COMMON/ANPORT/DS(4000)
      COMMON/VW/CW(60,4000),ICW(4000,2)
      COMMON/GAVV/X(4000),Y(4000),Z(4000)
      COMMON/GW/WX(60,4000),WY(60,4000),WZ(60,4000)
      COMMON/GVV/C(4000)
      COMMON/GRID/IC(4000,4),MARK(4000),MARKW(4000)
      COMMON/WIN/IPW(4000)
      DIMENSION XVV(60,4000),YVV(60,4000),ZVV(60,4000),XW(60,4000),
     1YW(60,4000),ZW(60,4000)
      DO 401 I=1,ITER  
      DO 301 J=1,NGW
      XW(I,J)=WX(I,J)
      YW(I,J)=WY(I,J)
      ZW(I,J)=WZ(I,J)
  301 CONTINUE               
  401 CONTINUE
      DO 80 I=1,ITER
      DO 70 J=1,NGW
      XG=XW(I,J)
      YG=YW(I,J)
      ZG=ZW(I,J)
      CALL VELPAN(NPAN,NSYM,NGRND,XG,YG,ZG,XVP,YVP,ZVP) 
      CALL VELWAK(ITER,NPW,NSYM,NGRND,XG,YG,ZG,XVW,
     1YVW,ZVW)
      XVV(I,J)=XVP+XVW
      YVV(I,J)=YVP+YVW
      ZVV(I,J)=ZVP+ZVW
   70 CONTINUE              
   80 CONTINUE
      DO 100 I=1,ITER
      DO 90 J=1,NGW
      XW(I,J)=XW(I,J)+XVV(I,J)*DT
      YW(I,J)=YW(I,J)+YVV(I,J)*DT
      ZW(I,J)=ZW(I,J)+ZVV(I,J)*DT
   90 CONTINUE
  100 CONTINUE
      DO 106 I=1,ITER  
      DO 115 J=1,NGW
      WX(I,J)=XW(I,J)
      WY(I,J)=YW(I,J)
      WZ(I,J)=ZW(I,J)
  115 CONTINUE               
  106 CONTINUE
      DO 157 I=1,ITER
      DO 155 J=1,NGW
      XW(I,J)=0.
      YW(I,J)=0.
      ZW(I,J)=0.
  155 CONTINUE               
  157 CONTINUE
      RETURN
      END      
      
 
	SUBROUTINE AIRLOAD1(NPAN,VINIT,RHO,ALIFT,DRAG,SIDE,NSYM,ITER,
	1NPW,NGRND)
      IMPLICIT REAL*4 (A-H,O-Z)
      COMMON/MWKVV/VXV(4000),VYV(4000),VZV(4000)
      COMMON/AVV1/ANX(4000),ANY(4000),ANZ(4000),GX(4000),GY(4000),
     1GZ(4000)
      COMMON/GAVV/X(4000),Y(4000),Z(4000)
      COMMON/GVV/C(4000)
      COMMON/GRID/IC(4000,4),MARK(4000),MARKW(4000)
      COMMON/ANPORT/DS(4000)
      COMMON/MP/CP(4000)
      COMMON/GW/WX(60,4000),WY(60,4000),WZ(60,4000)
      COMMON/VW/CW(60,4000),ICW(4000,2)
      DIMENSION IP(2)
      COEF=0.5*RHO*VINIT**2 
      KN=0
      DRAG=0.
      SIDE=0.
      ALIFT=0.
    5 DO 80 I=1,NPAN
      KKI=4
      IF(IC(I,3).EQ.IC(I,4)) KKI=3
      DO 30 J=1,NPAN
      IF(J.EQ.I) GO TO 30
      KKJ=4
      IF(IC(J,3).EQ.IC(J,4)) KKJ=3 
      KP=0 
      DO 20 IK=1,KKI 
      DO 10 JK=1,KKJ
      IF(IC(I,IK).NE.IC(J,JK)) GO TO 10 
      KP=KP+1              
      IP(KP)=IK
      IF(KP.NE.2) GO TO 20  
      IF(MARK(IC(I,IP(1))).EQ.1.AND.MARK(IC(I,IP(2))).EQ.1) GO TO 30
      DX=X(IC(I,IP(2)))-X(IC(I,IP(1))) 
      DY=Y(IC(I,IP(2)))-Y(IC(I,IP(1))) 
      IF(KN.EQ.1) DY=-DY
      DZ=Z(IC(I,IP(2)))-Z(IC(I,IP(1)))
      IF(ABS(IP(2)-IP(1)).NE.(KKI-1)) GO TO 11
      IF(IP(2).GT.IP(1)) GO TO 12
      GO TO 13
   11 IF(IP(2).GT.IP(1)) GO TO 13 
   12 DY=-DY
      DZ=-DZ     
   13 ROC=RHO*(C(I)-C(J))*VINIT/2.
      RC=RHO*(C(I)-C(J))
      IF(KN.EQ.1) ROC=-ROC
      SIDE=SIDE-ROC*DZ
      ALIFT=ALIFT+ROC*DY  
	GGX=GX(I)                        
	GGY=GY(I)                        
	GGZ=GZ(I)                        
      CALL VELWAK(ITER,NPW,NSYM,NGRND,GGX,GGY,GGZ,XVV,YVV,ZVV)
      DRAG=DRAG-RC*(XVV*DY-YVV*DX)
      GO TO 30
   10 CONTINUE   
   20 CONTINUE
   30 CONTINUE
   80 CONTINUE   
      IF(NSYM.EQ.0) GO TO 90
      IF(KN.EQ.1) GO TO 90
      KN=1
      GO TO 5
   90 CONTINUE
      RETURN
      END
      
      
      
      
    
         

