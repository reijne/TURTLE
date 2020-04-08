      SUBROUTINE OUT(ICH,MODE,COORD) 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'param'
      COMMON /GEOKST/ NATOMS,LABELS(NUMATM),
     .NA(NUMATM),NB(NUMATM),NC(NUMATM)
      COMMON /ELEMTS/ ELEMNT(107)
      CHARACTER*2 ELEMNT              
      DIMENSION COORD(3,NUMATM)
      COMMON /ATINDEX/ IADEX(NUMATM)
     
      F=1.0d0/0.529177260d0
      if(ich.eq.1) then
      open (unit=1,file='coord')
      rewind 1
      if(mode.le.1.or.mode.gt.4) write(1,'(''$coord'')')
      if(mode.eq.2) write(1,'(''geometry'')')
      if(mode.eq.3) then
         write(ich,*)' setup'
         write(ich,*)
         write(ich,*)
      endif
      ENDIF

      
      DO I=1,NATOMS
         DO J=1,3
            C=COORD(J,I)
            IF(DABS(C).LT.1.0D-6) C=0.0D0
            COORD(J,I)=C
         ENDDO
      ENDDO

cmmg begin of additional zindo input
cmmg counting number of valence electrons (nvel),
cmmg basis sets (nbas?) etc. in order to complete
cmmg a zindo input

             IF(MODE.EQ.4) then
     		write(ich,'(a40)') '!number of atoms you fight with         '
		nbas=0
		nbass=0
		nbasp=0
		nbasd=0
		nbasf=0
		nbasg=0
		nbash=0
		nvel=0
		m=0
		do i=1,natoms
		   if((labels(i).le.2).and.(labels(i).gt.0)) then
		   nvel=nvel+labels(i)
		   m=m+1
		   nbass=nbass+1
		   else if((labels(i).ge.3).and.(labels(i).le.10)) then
      		   nvel=nvel+labels(i)-2
		   m=m+1
		   nbasp=nbasp+1
		   else if((labels(i).ge.11).and.(labels(i).le.18)) then
      		   nvel=nvel+labels(i)-10
		   m=m+1
		   nbasp=nbasp+1
		   else if((labels(i).ge.19).and.(labels(i).le.36)) then
      		   nvel=nvel+labels(i)-18
		   m=m+1
		   nbasd=nbasd+1
		   else if((labels(i).ge.37).and.(labels(i).le.54)) then
      		   nvel=nvel+labels(i)-36
		   m=m+1
		   nbasd=nbasd+1
		   else if((labels(i).ge.55).and.(labels(i).le.86)) then
      		   nvel=nvel+labels(i)-54
		   m=m+1
		   nbasf=nbasf+1
		   endif
                enddo
		write(ich,*) 'NAT',m
		write(ich,'(a40)') '!number of valence electrons             '
		write(ich,*) 'NEL',nvel
     		write(ich,'(a40)') '! basis set and ci info (latter ignored) '
		write(ich,'(a9,7i3)') 'DYNAL(1)=',nbaso,nbass,nbasp,nbasd,nbasf,
     .          nbasg,nbash
		write(ich,'(a40)') '!**** name of the outputfile ****        '
		write(ich,'(a40)') 'ONAME=zindo.out                          '
		write(ich,*) '  '
                write(ich,'('' $END'')')
		write(ich,*) '  '
                write(ich,'('' $DATAIN'')')
		write(ich,*) '  '
                endif 

cmmg end of additional zindo input


      IDUM=0
      DO 240 I=1,NATOMS
         IF(LABELS(I).LT.99) THEN
             IF(MODE.LE.1.or.mode.gt.4)
     .          WRITE(ICH,'(8X,3F14.6,5X,A2)') (F*COORD(J,I),J=1,3),
     .          ELEMNT(LABELS(I))
             IF(MODE.EQ.2) THEN
                if(iadex(i).eq.0) then
                write(ICH,'(1x,a2,'','',i1,'','',
     .          F16.6,'','',F16.6,'','',F16.6)') elemnt(labels(i)),
     .          idum,(f*coord(j,i),j=1,3)
                else
                write(ICH,'(i2,'','',a2,'','',i1,'','',
     .          F16.6,'','',F16.6,'','',F16.6)')iadex(i),
     .          elemnt(labels(i)),
     .          idum,(f*coord(j,i),j=1,3)
                endif
             ENDIF
             IF(MODE.EQ.3) THEN
                IONE=1
                WRITE(ICH,'(A2,5X,3(F14.6,I4))') 
     .          ELEMNT(LABELS(I)),
     .          COORD(1,I),IONE,
     .          COORD(2,I),IONE,
     .          COORD(3,I),IONE
             ENDIF
             IF(MODE.EQ.4) 
     .          WRITE(ICH,'(3F10.6,I5)') (F*COORD(J,I),J=1,3),
     .          LABELS(I)
         ENDIF
  240 CONTINUE

      if(ich.eq.1) then
      if(mode.le.1.or.mode.gt.4) write(1,'(''$end'')')
      if(mode.eq.2) write(1,'(''end'')')
      if(mode.eq.4) write(1,'('' $END'')')
      close (1)
      endif

      return
      end

      SUBROUTINE OUTI(ICH,MODE,COORD) 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'param'
      COMMON /GEOKST/ NATOMS,LABELS(NUMATM),
     .NA(NUMATM),NB(NUMATM),NC(NUMATM)
      COMMON /ELEMTS/ ELEMNT(107)
      CHARACTER*2 ELEMNT              
      DIMENSION COORD(3,NUMATM)
     
      F=1.0d0/0.529177260d0
      if(ich.eq.1) then
      open (unit=1,file='coord')
      rewind 1
      if(mode.le.1.or.mode.gt.4) write(1,'(''$coord'')')
      if(mode.eq.2) write(1,'(''geometry'')')
      endif
      
cmmg1 if a calculation of a numerical hessian is desired,
cmmg  the reference geometry will be saved on file refcoord

      if(ich.eq.11) then
      open (unit=11,file='refcoord')
      rewind 11
      if(mode.le.1.or.mode.gt.4) write(11,'(''$coord'')')
      if(mode.eq.2) write(11,'(''geometry'')')
      endif
cmmg2
cmmg
      
      DO I=1,NATOMS
         DO J=1,3
            C=COORD(J,I)
            IF(DABS(C).LT.1.0D-6) C=0.0D0
            COORD(J,I)=C
         ENDDO
      ENDDO
 
cmmg begin of additional zindo input
cmmg counting number of valence electrons (nvel),
cmmg basis sets (nbas?) etc. in order to complete
cmmg a zindo input

             IF(MODE.EQ.4) then
     		write(ich,'(a40)') '!number of atoms you fight with         '
		nbas=0
		nbass=0
		nbasp=0
		nbasd=0
		nbasf=0
		nbasg=0
		nbash=0
		nvel=0
		m=0
		do i=1,natoms
		   if((labels(i).le.2).and.(labels(i).gt.0)) then
		   nvel=nvel+labels(i)
		   m=m+1
		   nbass=nbass+1
		   else if((labels(i).ge.3).and.(labels(i).le.10)) then
      		   nvel=nvel+labels(i)-2
		   m=m+1
		   nbasp=nbasp+1
		   else if((labels(i).ge.11).and.(labels(i).le.18)) then
      		   nvel=nvel+labels(i)-10
		   m=m+1
		   nbasp=nbasp+1
		   else if((labels(i).ge.19).and.(labels(i).le.36)) then
      		   nvel=nvel+labels(i)-18
		   m=m+1
		   nbasd=nbasd+1
		   else if((labels(i).ge.37).and.(labels(i).le.54)) then
      		   nvel=nvel+labels(i)-36
		   m=m+1
		   nbasd=nbasd+1
		   else if((labels(i).ge.55).and.(labels(i).le.86)) then
      		   nvel=nvel+labels(i)-54
		   m=m+1
		   nbasf=nbasf+1
		   endif
                enddo
		write(ich,*) 'NAT',m
		write(ich,'(a40)') '!number of valence electrons             '
		write(ich,*) 'NEL',nvel
     		write(ich,'(a40)') '! basis set and ci info (latter ignored) '
		write(ich,'(a9,7i3)') 'DYNAL(1)=',nbaso,nbass,nbasp,nbasd,nbasf,
     .          nbasg,nbash
		write(ich,'(a40)') '!**** name of the outputfile ****        '
		write(ich,'(a40)') 'ONAME=zindo.out                          '
		write(ich,*) '  '
                write(ich,'('' $END'')')
		write(ich,*) '  '
                write(ich,'('' $DATAIN'')')
		write(ich,*) '  '
                endif 

cmmg end of additional zindo input

      IDUM=0
      DO 240 I=1,NATOMS
         IF(LABELS(I).LT.99) THEN
             IF(MODE.LE.1.or.mode.gt.4)
     .          WRITE(ICH,'(8X,3F14.5,5X,A2)') (F*COORD(J,I),J=1,3),
     .          ELEMNT(LABELS(I))
             IF(MODE.EQ.2)
     .          write(ICH,'(1x,a2,'','',i1,'','',
     .          F16.7,'','',F16.7,'','',F16.7)') elemnt(labels(i)),
     .          idum,(f*coord(j,i),j=1,3)
             IF(MODE.EQ.3)
     .          WRITE(ICH,'(8X,3F14.7,5X,A2)') (F*COORD(J,I),J=1,3),
     .          ELEMNT(LABELS(I))
             IF(MODE.EQ.4) 
     .          WRITE(ICH,'(3F10.6,I5)') (F*COORD(J,I),J=1,3),
     .          LABELS(I)
         ENDIF
  240 CONTINUE

      if(ich.eq.1) then
      if(mode.le.1.or.mode.gt.4) write(1,'(''$end'')')
      if(mode.eq.2) write(1,'(''end'')')
      if(mode.eq.4) write(1,'('' $END'')')
      close (1)
      endif

      return
      end
