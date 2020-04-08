c 
c  $Author: jmht $
c  $Date: 2008-10-06 23:36:09 +0200 (Mon, 06 Oct 2008) $
c  $Locker:  $
c  $Revision: 5733 $
c  $Source: /c/qcg/cvs/psh/GAMESS-UK/m4/omitci.m,v $
c  $State: Exp $
c  
c ******************************************************
c ******************************************************
c             =   omitted routines (ci mode)   =
c ******************************************************
c ******************************************************
_IFN(convex)
      subroutine gamplt
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      call caserr('interactive graphics option unavailable')
      return
      end
      subroutine iniplt
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      return
      end
      subroutine stoplt
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      return
      end
      subroutine cntour
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      return
      end
      subroutine arch3d
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      return
      end
_ENDIF
_IFN(cray,convex,apollo,sgi,rs6000,sun,hp700,hpux11,dec,ksr)
      subroutine fulcin
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      call caserr('full-ci option unavailable')
      return
      end
      subroutine fci
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      call caserr('full-ci option unavailable')
      return
      end
_ENDIF
_IF(alliant,ipsc)
      subroutine tda
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      call caserr('tda option not available')
      return
      end
_ENDIF
_IF(ibm)
      subroutine ran
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      call caserr('random-number generator not available')
      return
      end
_ENDIF
_IF1(u)      subroutine tda
_IF1(u)      call caserr('tda option not available')
_IF1(u)      return
_IF1(u)      end
_IF1(u)c            =  omitted routines (mrdci)  =
_IF1(u)      subroutine amrdci
_IF1(u)      implicit REAL  (a-h,o-z)
_IF1(u)      entry tmrdci
_IF1(u)      entry tdata
_IF1(u)      entry smrdci
_IF1(u)      entry tabci
_IF1(u)      entry nmrdci
_IF1(u)      entry pmrdci
_IF1(u)      entry moment
_IF1(u)      entry dmrdci
_IF1(u)      entry closb3
_IF1(u)      call caserr('mrdci option not available')
_IF1(u)      return
      subroutine chfndr
      implicit REAL (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      entry raman
      entry mp2dd
      entry mp2eng
      entry hfpder
      entry anairr
      entry hfdder
      entry trndrv
      entry trnfkd
      entry hfhpol
      entry scfdd
      entry mpdder
      entry mptran
      entry chfidr
      entry indxsy
      entry emp23
      entry uhfmp2
      entry uhfmp3
      entry drcdrv
      entry mopin
_IFN(mopac)
      entry mopac
_ENDIF
      entry grmp23
      entry titandrv
      entry ccsdin
      entry rpdriv
      entry rpinit
      entry lrmain
      entry lrinit
       call caserr ('unavailable routine in ci-only code')
      return
      end
      entry uhfmp2
      entry uhfmp3
      entry drcdrv
      entry mopin
_IFN(mopac)
      entry mopac
_ENDIF
      entry grmp23
       call caserr (' unavailable routine in ci-only code')
      return
      end
