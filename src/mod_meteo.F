!------------------------------------------------------------------------------------
!
!      FILE mod_meteo.F
!
!      This file is part of the FUNWAVE-TVD program under the Simplified BSD license
!
!-------------------------------------------------------------------------------------
! 
!    Copyright (c) 2017, FUNWAVE Development Team
!
!    (See http://www.udel.edu/kirby/programs/funwave/funwave.html
!     for Development Team membership)
!
!    All rights reserved.
!
!    FUNWAVE_TVD is free software: you can redistribute it and/or modify
!    it under the terms of the Simplified BSD License as released by
!    the Berkeley Software Distribution (BSD).
!
!    Redistribution and use in source and binary forms, with or without
!    modification, are permitted provided that the following conditions are met:
!
!    1. Redistributions of source code must retain the above copyright notice, this
!       list of conditions and the following disclaimer.
!    2. Redistributions in binary form must reproduce the above copyright notice,
!    this list of conditions and the following disclaimer in the documentation
!    and/or other materials provided with the distribution.
!
!    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
!    ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
!    WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
!    DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
!    ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
!    (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
!    LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
!    ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
!    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
!    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!  
!    The views and conclusions contained in the software and documentation are those
!    of the authors and should not be interpreted as representing official policies,
!    either expressed or implied, of the FreeBSD Project.
!  
!-------------------------------------------------------------------------------------
!
!  meteo is a module to model atmospherically forced meteo-tsunami waves
!
!  HISTORY :
!    1/25/2017  Jim Kirby  Copied initial code from mod_METEO.F
!
!-------------------------------------------------------------------------------------

# if defined (METEO)

MODULE METEO_MODULE
  USE PARAM
  USE GLOBAL,ONLY : Mloc,Nloc,Nghost,Ibeg,Iend,Jbeg,Jend,DX,DY, &
                    H
  USE INPUT_READ
#if defined (PARALLEL)
  USE GLOBAL,ONLY : myid,ier, npx,npy,PX,PY
  USE MPI
# endif
  IMPLICIT NONE
  SAVE

    INTEGER :: NumMETEO,Kves
    REAL(SP),DIMENSION(:),ALLOCATABLE :: XMETEO1,YMETEO1,XMETEO2,YMETEO2, &
                                       LengthMETEO,WidthMETEO, &
                                       AlphaMETEO,BetaMETEO,PMETEO, &
                                       TimeMETEO1,TimeMETEO2,ThetaMETEO
    REAL(SP),DIMENSION(:,:),ALLOCATABLE :: METEOPressureTotal,METEOPressureEach, &
                                       METEOPressureX,METEOPressureY
    REAL(SP),DIMENSION(:),ALLOCATABLE ::  Xco,Yco
    LOGICAL :: OUT_METEO = .TRUE.


!INTERFACE READ_FOUR_TYPE_VALUES
!  Module Procedure METEO_INITIAL
!  Module Procedure METEO_FORCING
!END INTERFACE

CONTAINS
  
! READ METEO

SUBROUTINE METEO_INITIAL
  USE GLOBAL,ONLY : itmp1,itmp2,itmp3,itmp4
                    
  USE INPUT_READ
  IMPLICIT NONE
  CHARACTER(LEN=80)::FILE_NAME=' '
  CHARACTER(LEN=80)::FILE_METEO=' '
  CHARACTER(LEN=80) :: METEO_FOLDER=' '
  CHARACTER(LEN=80)::TMP_NAME=' '
  INTEGER :: Ifile,ierr
  CHARACTER(LEN=80):: METEOName

! read METEO number and folder from input.txt
      FILE_NAME='input.txt'

! METEO folder
      CALL READ_STRING(METEO_FOLDER,FILE_NAME,'METEO_FOLDER',ierr)
# if defined (PARALLEL)
      if (myid.eq.0) WRITE(3,'(A15,A50)')'METEO_FOLDER:', METEO_FOLDER
# else
      WRITE(3,'(A15,A50)')'METEO_FOLDER:', METEO_FOLDER
# endif

      CALL READ_INTEGER(NumMETEO,FILE_NAME,'NumMETEO',ierr)
# if defined (PARALLEL)
      if (myid.eq.0) WRITE(3,'(A12,I3)') 'NumMETEO = ',NumMETEO
# else
      WRITE(3,'(A12,I3)') 'NumMETEO = ',NumMETEO
# endif

      CALL READ_LOGICAL(OUT_METEO,FILE_NAME,'OUT_METEO',ierr)

      ALLOCATE (XMETEO1(NumMETEO),YMETEO1(NumMETEO),  &
                XMETEO2(NumMETEO),YMETEO2(NumMETEO),  &
                TimeMETEO1(NumMETEO),TimeMETEO2(NumMETEO), &
                LengthMETEO(NumMETEO),WidthMETEO(NumMETEO), &
                AlphaMETEO(NumMETEO),BetaMETEO(NumMETEO),  &
                PMETEO(NumMETEO),ThetaMETEO(NumMETEO) )
      ALLOCATE (METEOPressureTotal(Mloc,Nloc), METEOPressureEach(Mloc,Nloc),&
                 Xco(Mloc),Yco(Nloc), METEOPressureX(Mloc,Nloc), &
                 METEOPressureY(Mloc,Nloc) )


! Xco, and Yco

# if defined (PARALLEL)
     Xco(Ibeg) = npx*(Mloc-2*Nghost)*DX
# else
     Xco(Ibeg) = ZERO
# endif
     DO I = Ibeg+1,Mloc
       Xco(I) = Xco(I-1)+DX
     ENDDO
     DO I = Ibeg-1,Ibeg-Nghost,-1
       Xco(I) = Xco(I+1)-DX
     ENDDO

# if defined (PARALLEL)
     Yco(Jbeg) = npy*(Nloc-2*Nghost)*DY
# else
     Yco(Jbeg) = ZERO
# endif
     DO J = Jbeg+1,Nloc
       Yco(J) = Yco(J-1)+DY
     ENDDO
     DO J = Jbeg-1,Jbeg-Nghost,-1
       Yco(J) = Yco(J+1)-DY
     ENDDO

  DO Kves = 1, NumMETEO

!  file name
    itmp1=mod(Kves/1000,10)
    itmp2=mod(Kves/100,10)
    itmp3=mod(Kves/10,10)
    itmp4=mod(Kves,10)
    write(FILE_METEO(1:1),'(I1)')itmp1
    write(FILE_METEO(2:2),'(I1)')itmp2
    write(FILE_METEO(3:3),'(I1)')itmp3
    write(FILE_METEO(4:4),'(I1)')itmp4

    TMP_NAME = TRIM(METEO_FOLDER)//'METEO_'//TRIM(FILE_METEO)

! check existing

 INQUIRE(FILE=TRIM(TMP_NAME),EXIST=FILE_EXIST)
  IF(.NOT.FILE_EXIST)THEN
# if defined (PARALLEL)
   IF(MYID==0)  &
   WRITE(*,*) TRIM(TMP_NAME), ' specified in ', TRIM(METEO_FOLDER), ' but CANNOT BE FOUND. STOP'
   CALL MPI_FINALIZE (ier)
   STOP
# else
    WRITE(*,*) TRIM(TMP_NAME), ' specified in ', TRIM(METEO_FOLDER), ' but CANNOT BE FOUND. STOP'
    STOP
# endif
  ENDIF

! open file
  Ifile=Kves+200
  OPEN(Ifile,FILE=TRIM(TMP_NAME))

! read file
         READ(Ifile,*)  ! title
         READ(Ifile,*)  METEOName  ! METEO name
         READ(Ifile,*)  ! length and width
         READ(Ifile,*)  LengthMETEO(Kves), WidthMETEO(Kves), &
                      AlphaMETEO(Kves),BetaMETEO(Kves),PMETEO(Kves)
         READ(Ifile,*)  ! t, x, y
         READ(Ifile,*)  TimeMETEO2(Kves),XMETEO2(Kves),YMETEO2(Kves)

         TimeMETEO1(Kves) = TimeMETEO2(Kves)
         XMETEO1(Kves) = XMETEO2(Kves)
         YMETEO1(Kves) = YMETEO2(Kves)

         AlphaMETEO(Kves) = Max(SMALL, AlphaMETEO(Kves))
         BetaMETEO(Kves) = Max(SMALL, BetaMETEO(Kves))
         AlphaMETEO(Kves) = Min(1.0_SP, AlphaMETEO(Kves))
         BetaMETEO(Kves) = Min(1.0_SP, BetaMETEO(Kves))

# if defined (PARALLEL)
   IF(MYID==0)THEN
   WRITE(3,*) 'METEO Name: ',  TRIM(METEOName)
   WRITE(3,*) 'METEO Length', LengthMETEO(Kves)
   WRITE(3,*) 'METEO Width', WidthMETEO(Kves)
   WRITE(3,*) 'METEO Alpha', AlphaMETEO(Kves)
   WRITE(3,*) 'METEO Beta', BetaMETEO(Kves)
   WRITE(3,*) 'METEO P', PMETEO(Kves)
   WRITE(3,*) 'Initial Time, X, Y', TimeMETEO2(Kves),XMETEO2(Kves),YMETEO2(Kves)
   ENDIF
# else
   WRITE(3,*) 'METEO Name: ',  TRIM(METEOName)
   WRITE(3,*) 'METEO Length', LengthMETEO(Kves)
   WRITE(3,*) 'METEO Width', WidthMETEO(Kves)
   WRITE(3,*) 'METEO Alpha', AlphaMETEO(Kves)
   WRITE(3,*) 'METEO Beta', BetaMETEO(Kves)
   WRITE(3,*) 'METEO P', PMETEO(Kves)
   WRITE(3,*) 'Initial Time, X, Y', TimeMETEO2(Kves),XMETEO2(Kves),YMETEO2(Kves)
# endif

  ENDDO  ! end Kves

End SUBROUTINE METEO_INITIAL

SUBROUTINE METEO_FORCING
  USE GLOBAL,ONLY : Mloc,Nloc,tmp1,tmp2,SMALL,TIME,ZERO
  USE INPUT_READ
  IMPLICIT NONE
  INTEGER :: Ifile,ierr,I,J
  REAL(SP) :: Xves,Yves,Lves,Wves

  METEOPressureTotal = ZERO

  DO Kves = 1,NumMETEO

    IF(TIME>TimeMETEO1(Kves).AND.TIME>TimeMETEO2(Kves)) THEN

         TimeMETEO1(Kves)=TimeMETEO2(Kves)
         XMETEO1(Kves) = XMETEO2(Kves)
         YMETEO1(Kves) = YMETEO2(Kves)

    Ifile = 200 + Kves

    READ(Ifile,*,END=120)  TimeMETEO2(Kves),XMETEO2(Kves),YMETEO2(Kves)

    ThetaMETEO(Kves) = ATAN2(YMETEO2(Kves)-YMETEO1(Kves),  &
                              XMETEO2(Kves)-XMETEO1(Kves))

# if defined (PARALLEL)
   IF(MYID==0)THEN
     WRITE(3,*)'Read METEO # ', Kves
     WRITE(3,*)'T,X,Y = ', TimeMETEO2(Kves),XMETEO2(Kves),YMETEO2(Kves)
   ENDIF
# else
     WRITE(3,*)'Read METEO # ', Kves
     WRITE(3,*)'T,X,Y = ', TimeMETEO2(Kves),XMETEO2(Kves),YMETEO2(Kves)
# endif

    ENDIF ! end time > timeMETEO2

! calculate force
    tmp2=ZERO
    tmp1=ZERO

    IF(TIME>TimeMETEO1(Kves))THEN
      IF(TimeMETEO1(Kves).EQ.TimeMETEO2(Kves))THEN
        ! no more data
        tmp2=ZERO
        tmp1=ZERO
      ELSE
      tmp2=(TimeMETEO2(Kves)-TIME) &
            /MAX(SMALL, ABS(TimeMETEO2(Kves)-TimeMETEO1(Kves)))
      tmp1=1.0_SP - tmp2;
      ENDIF  ! no more data?
    ENDIF ! time>time_1

    Xves = XMETEO2(Kves)*tmp1 +XMETEO1(Kves)*tmp2
    Yves = YMETEO2(Kves)*tmp1 +YMETEO1(Kves)*tmp2

! rectangular
    METEOPressureEach = ZERO
    DO J=1,Nloc
    DO I=1,Mloc
!     Lves=(Xco(I)-Xves)*COS(ThetaMETEO(Kves)) + (Yco(J)-Yves)*SIN(ThetaMETEO(Kves))
!     Wves=-(Xco(I)-Xves)*SIN(ThetaMETEO(Kves)) + (Yco(J)-Yves)*COS(ThetaMETEO(Kves))			    !     IF(ABS(Lves)<=0.5_SP*LengthMETEO(Kves).AND. &
!        ABS(Wves)<=0.5_SP*WidthMETEO(Kves)) THEN
         METEOPressureEach(I,J) = ((1/(2*pi*(sigma**2)))*exp(-((Xco(I)-Xves)**2+(Yco(J)-Yves)**2)/(2*(sigma**2))))**2                 
      ENDIF
    ENDDO
    ENDDO

120 CONTINUE  ! no more data for METEO Kves

    METEOPressureTotal = METEOPressureTotal+METEOPressureEach

  ENDDO  ! end Kves

! sourceX and sourceY

    DO J=Jbeg,Jend
    DO I=Ibeg,Iend

!   I modified the term to negative 11/22/2016 fyshi

       METEOPressureX(I,J) = -Grav*H(I,J)*  &
               (METEOPressureTotal(I+1,J)-METEOPressureTotal(I-1,J))/2.0_SP  &
               /DX
       METEOPressureY(I,J) = -Grav*H(I,J)*  &
               (METEOPressureTotal(I,J+1)-METEOPressureTotal(I,J-1))/2.0_SP  &
               /DY

    ENDDO
    ENDDO


END SUBROUTINE METEO_FORCING

END MODULE METEO_MODULE

# endif 
! end METEO
