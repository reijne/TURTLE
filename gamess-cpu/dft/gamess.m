c
c  implementation of DFT memory routines using GAMESS-UK
c  memory manager
c
      integer function push_memory_estimate()
      implicit none
      integer igmem_push_estimate
      push_memory_estimate = igmem_push_estimate()
      return
      end

      integer function pop_memory_estimate(icnt)
      implicit none
      integer icnt
      integer igmem_pop_estimate
      pop_memory_estimate = igmem_pop_estimate(icnt)
      return
      end

      integer function push_memory_count()
      implicit none
      integer igmem_push_usage
      push_memory_count = igmem_push_usage()
      return
      end

      integer function pop_memory_count(icnt)
      implicit none
      integer icnt
      integer igmem_pop_usage
      pop_memory_count = igmem_pop_usage(icnt)
      return
      end

c     subroutine init_memory_count
c     implicit none
c     call gmem_count_init
c     return
c     end

c     subroutine reset_memory_count
c     implicit none
c     call gmem_count_reset
c     return
c     end

      integer function max_memory_count()
      implicit none
      integer igmem_count_max
      max_memory_count = igmem_count_max()
      return
      end

      integer function current_memory_count()
      implicit none
      integer igmem_count_current
      current_memory_count = igmem_count_current()
      return
      end

      integer function incr_memory(length,type)
      implicit none
INCLUDE(common/dft_iofile)
      integer length
      integer lenwrd
      integer itemp
      integer ilen
      integer igmem_incr
      integer null_memory
      character*(*) type
      if(type .eq. 'd') then
c
         itemp = igmem_incr(length)
         incr_memory = itemp
         return 
c
      else if (type .eq. 'i') then
c 
         ilen = lenwrd()
c        itemp = igmem_incr((length-1)/ilen + 1)
c        the line above and line below are not the same in fortran77.
         incr_memory = igmem_incr((length+ilen-1)/ilen)
         return 
c
      else
         write(iwr,*)type
         incr_memory = -1
         call caserr('incr_memory: unknown data type')
         incr_memory = null_memory()
      endif
c
      end
      

      subroutine decr_memory(address,type)
INCLUDE(common/dft_iofile)
      integer ilen
      integer address
      character*(*) type
      integer lenwrd
      integer null_memory
c
      if(type .eq. 'd') then
c
         call gmem_decr(address)
c
      else if (type .eq. 'i') then
c
         ilen = lenwrd()
c        call gmem_decr( (address + (ilen-1)) /ilen)
         call gmem_decr(address)

      else
         write(6,*)type
         call caserr('decr_memory: unknown type')
      endif
c
      address = null_memory()
c
      end


      integer function incr_memory2(length,type,filename,subrname,
     &                              varid)
      implicit none
      integer length
      integer incr_memory
      character*(*) filename,subrname,varid
      character*(*) type
      incr_memory2=incr_memory(length,type)
      return
      end


      subroutine decr_memory2(address,type,filename,subrname,varid)
      implicit none
      integer address
      character*(*) type
      character*(*) filename,subrname,varid
      call decr_memory(address,type)
      end


      integer function memory_overhead()
      implicit none
      integer igmem_overhead
      memory_overhead = igmem_overhead()
      return
      end

      integer function null_memory()
      implicit none
      integer igmem_null
      null_memory = igmem_null()
      return
      end

      integer function allocate_memory(length,type)
      implicit none
INCLUDE(common/dft_iofile)
      integer length
      integer lenwrd
      integer itemp
      integer ilen
      integer igmem_alloc
      integer null_memory

      character*(*) type

      if(type .eq. 'd') then

         itemp = igmem_alloc(length)
         allocate_memory = itemp
         return 

      else if (type .eq. 'i') then

         ilen = lenwrd()
c        itemp = igmem_alloc((length-1)/ilen + 1)
c        the line above and line below are not the same in fortran77.
         itemp = igmem_alloc((length+ilen-1)/ilen)
         allocate_memory = itemp * ilen - (ilen -1)
         return 

      else
         write(6,*)type
         allocate_memory = -1
         call caserr('allocate_memory: unknown data type')
         allocate_memory = null_memory()
      endif

      end
      
      subroutine free_memory(address,type)
INCLUDE(common/dft_iofile)
      integer ilen
      integer address
      character*(*) type
      integer lenwrd
      integer null_memory

      if(type .eq. 'd') then

         call gmem_free(address)

      else if (type .eq. 'i') then

         ilen = lenwrd()
         call gmem_free( (address + (ilen-1)) /ilen)

      else
         write(6,*)type
         call caserr('free_memory: unknown type')
      endif

      address = null_memory()

      end

      integer function allocate_memory2(length,type,filename,subrname,
     &                                  varid)
      implicit none
INCLUDE(common/dft_iofile)
INCLUDE(../m4/common/gmempara)
      integer length
      integer lenwrd
      integer itemp
      integer ilen
      integer igmem_alloc_inf
      integer null_memory
      character*(*) filename,subrname,varid

      character*(*) type

      if(type .eq. 'd') then

         itemp = igmem_alloc_inf(length,filename,subrname,
     &                           varid,IGMEM_DEBUG)
         allocate_memory2 = itemp
         return 

      else if (type .eq. 'i') then

         ilen = lenwrd()
c        itemp = igmem_alloc_inf((length-1)/ilen + 1,filename,subrname,
c    &                           varid,IGMEM_DEBUG)
c        the line above and line below are not the same in fortran77.
         itemp = igmem_alloc_inf((length+ilen-1)/ilen,filename,subrname,
     &                           varid,IGMEM_DEBUG)
         allocate_memory2 = itemp * ilen - (ilen -1)
         return 

      else
         write(6,*)type
         allocate_memory2 = -1
         call caserr('allocate_memory: unknown data type')
         allocate_memory2 = null_memory()
      endif

      end
      
      subroutine free_memory2(address,type,filename,subrname,varid)
      implicit none
INCLUDE(common/dft_iofile)
      integer ilen
      integer address
      character*(*) type
      integer lenwrd
      integer null_memory
      character*(*) filename,subrname,varid

      if(type .eq. 'd') then

         call gmem_free_inf(address,filename,subrname,varid)

      else if (type .eq. 'i') then

         ilen = lenwrd()
         call gmem_free_inf( (address + (ilen-1)) /ilen,filename,
     &                      subrname,varid)

      else
         write(6,*)type
         call caserr('free_memory: unknown type')
      endif

      address = null_memory()

      end

      subroutine interface_gamess(rhf_sw,opshell_sw,ao_tag,iout)
C **********************************************************************
C *Description:							       *
C *Interface into gamess for ccp1 dft modules                          *
C **********************************************************************
      implicit none
C **********************************************************************
C *Declarations
C *
C *Parameters
INCLUDE(common/dft_parameters)
INCLUDE(../m4/common/sizes)
C *In variables
      logical rhf_sw
      integer ao_tag
      integer iout
INCLUDE(../m4/common/infoa)
INCLUDE(../m4/common/infob)
INCLUDE(../m4/common/nshel)
INCLUDE(../m4/common/runlab)
C *Out variables
INCLUDE(common/dft_basis)
INCLUDE(common/dft_basis_api)
INCLUDE(common/dft_mol_info)
INCLUDE(common/dft_module_comm)
      logical opshell_sw
C * Functions
      logical     opg_root
      integer     isubst
C *Local variables
      integer     necheck
      integer     atmnum(maxat)
      integer     latm,lgshl,lshl,kgshl,kshl
      logical     atomexist_sw, osame
      integer     atom_id(maxat),minshl,maxshl,ploc
c     character*8 z_at_id(max_atype)
c     integer     n_at_id,i_id, i, j
      integer     i_id, i, j
      integer     hybri
      integer     latmno, katmno
      integer     nlshl, nkshl, plloc, pkloc, nfunc
      integer     ierror
C *End declarations                                                    *
C **********************************************************************
c     if(opg_root())call ccpdft_banner(iout)
      necheck=mod(ne,2)
      if(necheck.ne.0.and.rhf_sw) opshell_sw = .true.

      ierror = BL_clear_basis_set(ao_tag)

      do i = 1,natoms
        atom_id(i) = 0
      enddo
C 
C Export AO basis set
c     n_at_id = 0
      do lgshl=1,nshell
        atomexist_sw=.false.
        latmno = isubst(zaname(katom(lgshl)))
        do kgshl=1,lgshl-1
           if(zaname(katom(kgshl)).eq.zaname(katom(lgshl))) then
              atomexist_sw = .true.
              atom_id(katom(lgshl)) = atom_id(katom(kgshl))
              goto 5
           endif
        enddo
 5      continue
c
        if(.not.atomexist_sw) then
c
c...       find all shells on centre l
c
           nlshl=lgshl
           do lshl=lgshl+1,nshell
             if(katom(lshl).eq.katom(lgshl)) nlshl=nlshl+1
           enddo
c
           if (lgshl.eq.1) then
              atomexist_sw = .false.
           else
              kgshl=1
c
c...          find all shells on centre k
c
 10           nkshl=kgshl
              do kshl=kgshl+1,nshell
                 if(katom(kshl).eq.katom(kgshl)) nkshl=nkshl+1
              enddo
c
c...          check if all exponents and coefficients are the same on 
c...          both centres (if the nuclear charges are different we
c...          will assume that the basis is different too).
c
              katmno = isubst(zaname(katom(kgshl)))
              if (nkshl-kgshl.eq.nlshl-lgshl.and.latmno.eq.katmno) then
                 osame = .true.
                 do i = 0,nkshl-kgshl
c
c...                check if the number of primitives is the same
c
                    osame = osame.and.(kng(lgshl+i).eq.kng(kgshl+i))
     +                           .and.(kmin(lgshl+i).eq.kmin(kgshl+i))
     +                           .and.(kmax(lgshl+i).eq.kmax(kgshl+i))
                    if (osame) then
                       nfunc = kng(lgshl+i)
                       plloc = kstart(lgshl+i)
                       pkloc = kstart(kgshl+i)
                       do j = 0, nfunc-1
                          osame = osame.and.(ex(plloc+j).eq.ex(pkloc+j))
                          osame = osame.and.(
     +                               cs(plloc+j).eq.cs(pkloc+j).or.
     +                               1.lt.kmin(lgshl+i) )
                          osame = osame.and.(
     +                               cp(plloc+j).eq.cp(pkloc+j).or.
     +                               2.gt.kmax(lgshl+i).or.
     +                               4.lt.kmin(lgshl+i) )
                          osame = osame.and.(
     +                               cd(plloc+j).eq.cd(pkloc+j).or.
     +                               5.gt.kmax(lgshl+i).or.
     +                               10.lt.kmin(lgshl+i) )
                          osame = osame.and.(
     +                               cf(plloc+j).eq.cf(pkloc+j).or.
     +                               11.gt.kmax(lgshl+i).or.
     +                               20.lt.kmin(lgshl+i) )
                          osame = osame.and.(
     +                               cg(plloc+j).eq.cg(pkloc+j).or.
     +                               21.gt.kmax(lgshl+i))
                       enddo
                    endif
                 enddo
              else
                 osame = .false.
              endif
              atomexist_sw = atomexist_sw.or.osame
              if (.not.atomexist_sw) then
                 kgshl=nkshl+1
                 if (kgshl.lt.lgshl) goto 10
              else
                 atom_id(katom(lgshl)) = atom_id(katom(kgshl))
              endif
           endif
c
c...       
           if(.not.atomexist_sw) then
c
c step 1 - create new atom tag
c            if (n_at_id.ge.max_atype) then
c               call caserr('Maximum number of atom types exceeded')
c            endif
c            n_at_id=n_at_id+1
             atom_id(katom(lgshl))=BL_create_atomtag(ao_tag,
     +                                               ian(katom(lgshl)))
c            z_at_id(n_at_id)=zaname(katom(kgshl))
             if (atom_id(katom(lgshl)).lt.0) then
                call caserr('Failed to create a new atom tag')
             endif
c
c step 2 - find out how many shells exist on this centre
c            nshl=lgshl
c            do lshl=lgshl+1,nshell
c               if(katom(lshl).eq.katom(lgshl)) nshl=nshl+1
c            enddo
c
c step3 - loop over shells on this centre and import the information
             minshl = lgshl
             maxshl = nlshl
             do lshl=minshl,maxshl
               ploc=kstart(lshl)
               hybri=ktype(lshl)
               if((ktype(lshl).eq.2).and.(kmax(lshl)-kmin(lshl).eq.3))
     &           hybri=1
               ierror=BL_import_shell(ao_tag,atom_id(katom(lshl)),
     &                              kng(lshl),ktype(lshl),hybri,
     &                              ex(ploc),cs(ploc),cp(ploc),cd(ploc),
     &                              cf(ploc),cg(ploc))
             enddo
           endif
        endif
      enddo
c
c     Assign basis functions to centres
c
c     do lgshl=1,nshell
c       do i_id=1,n_at_id
c          if (z_at_id(i_id).eq.zaname(katom(lgshl))) then
c             ierror = BL_assign_type(ao_tag,katom(lgshl),atom_id(i_id))
c             if(ierror .ne. 0) call caserr(
c    &             'basis set assignment failed for a.o. basis')
c          endif
c       enddo

      do i = 1, natoms
        ierror = BL_assign_type(ao_tag,i,atom_id(i))
        if(ierror .ne. 0) call caserr(
     &             'basis set assignment failed for a.o. basis')
      enddo
C 
C Assign basis functions to centres
c Currently only the atomic number is used
c
c     ierror=BL_assign_types_by_z(ao_tag)
c     if(ierror .ne. 0) call caserr(
c    &     'basis set assignment failed for a.o. basis')
C
C Check in basis
      call checkin_basis(ao_tag,.false.)
      bset_tags(ao_tag) = 1

      return
      end
c
c  $Author: hvd $ $Revision: 5774 $ $Date: 2008-12-05 00:26:07 +0100 (Fri, 05 Dec 2008) $
c  $Source: /c/qcg/cvs/psh/GAMESS-UK/dft/gamess.m,v $
c

      subroutine ver_dft_gamess(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/dft/gamess.m,v $
     +     "/
      data revision /
     +     "$Revision: 5774 $"
     +      /
      data date /
     +     "$Date: 2008-12-05 00:26:07 +0100 (Fri, 05 Dec 2008) $"
     +     /
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
