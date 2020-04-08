      subroutine isctr ( NZ, X, INDX, Y)
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C $Id: iblas.m,v 1.13 2007-11-20 11:01:28 psh Exp $
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     isctr -- Scatter elements of a vector into another
C
C SYNOPSIS
      Integer NZ, INDX (NZ), Y(*), X(NZ)
C
C ARGUMENTS
C     NZ      Number of elements to be gathered (input)
C     X       Source vector (input)
C     Indx    Vector of indices of Y to put elements of X into (input)
C     Y       Source vector (input)
C
C DESCRIPTION
C     Performs the operation Y( Indx(i) ) = X(i) for i = 1, NZ.
C     Only elements of Y listed in Indx are accessed.
C
C NOTES
C     This is the integer analog to [DS]Sctr of the Dodson, Lewis, and
C     Grimes Sparse-BLAS1 proposal.
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C     LOCAL VARIABLES
C
      INTEGER             I
C
      IF  ( NZ .LE. 0 )  RETURN
C
      DO 10 I = 1, NZ
         Y(INDX(I)) = X(I)
   10 CONTINUE
C
      RETURN
      END
