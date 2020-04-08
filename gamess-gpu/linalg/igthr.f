      subroutine igthr ( NZ, Y, X, Indx )
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C $Id: iblas.m,v 1.13 2007-11-20 11:01:28 psh Exp $
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     igthr -- Gather from a vector using an index
C
C SYNOPSIS
      Integer NZ, INDX (NZ), Y(*), X(NZ)
C
C ARGUMENTS
C     NZ      Number of elements to be gathered (input)
C     Y       Source vector (input)
C     X       Destination vector (input)
C     Indx    Vector of indices of Y to put into X (input)
C
C DESCRIPTION
C     Performs the operation X(i) = Y( Indx(i) ) for i = 1, NZ.
C     Only elements of Y listed in Indx are accessed.
C
C NOTES
C     This is the integer analog to [DS]Gthr of the Dodson, Lewis, and
C     Grimes Sparse-BLAS1 proposal.
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C     LOCAL VARIABLES
C
      INTEGER             I
C
      IF  ( NZ .LE. 0 )  RETURN
C
      DO 10 I = 1, NZ
          X(I) = Y(INDX(I))
   10 CONTINUE
C
      RETURN
      END
