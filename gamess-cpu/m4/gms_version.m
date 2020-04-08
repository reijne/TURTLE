      subroutine gms_version(cdate,ctime,cname,cversion)
      character*10 cdate,cname,cversion
      character*5 ctime
      cdate="08-04-2020"
      ctime="10:48"
      cname="you"
      cversion=""
      return
      end
      subroutine getm4keys(m4keys)
      implicit none
      character*(*) m4keys
       m4keys=" opteron linux pclinux littleendian cio unix"
     &//" doublebackslash upck-equiv i8 serial ccpdft "
     &//" drf vdw sysmo dl-find nbo zora mrdci vb vdw "
     &//" sysmo full_build "
      return
      end
      subroutine wrtkeys(iwr)
      write(iwr,*)
     &"M4-Keys:  opteron linux pclinux littleendian cio unix "
      write(iwr,*)
     &"M4-Keys:  doublebackslash upck-equiv i8 serial ccpdft "
      write(iwr,*)
     &"M4-Keys:  drf vdw sysmo dl-find nbo zora mrdci vb vdw "
      write(iwr,*)
     &"M4-Keys:  sysmo full_build "
      return
      end
