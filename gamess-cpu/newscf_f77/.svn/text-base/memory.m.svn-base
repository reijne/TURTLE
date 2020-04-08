c
c  Implementation of dynamic memory functionality
c  for newscf
c
c  In MA case, we can simply take all memory required
c  of the MA heap
c
c  In the absence of MA, we need a heap allocator
c  Alex Turner provided one
c
_IF(ma)

      subroutine memory_create()
c
c     Initialises heap
c
      implicit none
      integer i
      include 'matrix_ma.fh'
      logical opg_root
c
c     Set up handle array - all handles start of as -1
c
      do i=1,maxhandle
         hp(i)=-1
      enddo

      return
      end

      integer function stack_alloc(bytes)
      implicit none
INCLUDE(../m4/common/vcore)
      integer words, bytes
      integer igmem_alloc
      external igmem_alloc
      words = 1 + (bytes - 1) / 8
      stack_alloc = igmem_alloc(words)
c      write(6,*)'Address of stack alloc'
c      call chkadr(Q(stack_alloc))
      end

      integer function stack_free(addr)
      implicit none
      integer addr
      call gmem_free(addr)
      stack_free=0
      end


c
c handle_alloc allocates a heap memory segment
c and saves the relevant offset in the array hp
c

      integer function handle_alloc(bytes)

c
c     Allocates the memory and places
c     a pointer to it in hp, then
c     returns the index in hp of that pointer
c
      implicit none

INCLUDE(matrix_internal.fh)

      integer i,j,k,fmalloc,bytes, words
 
      logical ostat, opg_root
      character tag*30
c
INCLUDE(../m4/common/gmemdata)

_IF(64bitpointers)
      integer*8 iq64
      common/pointer64/iq64
_ENDIF

#include "mafdecls.fh"

      MA_INTEGER size_8, type_8, itag_8, iq_8

c
c     find the ammount of memory to allocate
c
      words = 1 + (bytes - 1) / 8

      do i=1,maxhandle
         if(hp(i).eq.-1) then

            write(tag,100)i
 100        format('gamess_',i3.3)

            type_8 = MT_DBL
            size_8 = words
c     MA tools wont allow zero allocation
            if(size_8.eq.0)size_8=1
            ostat = ma_alloc_get(
     &           type_8,
     &           size_8,
     &           tag,
     &           itag_8,
     &           iq_8)

            heap_tag(i) = itag_8



c      write(6,*)'Chek address after ma (memory.m)'
c      call chkadr(dbl_mb(iq_8))

c
c     store as the GAMESS-UK core
c	
_IF(64bitpointers)
c
c Adjust to get a small offset
c see computation of ivoff
c
            hp(i) = iq_8 - iq64 + 1
_ELSE
            hp(i) = iq_8 - iqoff
_ENDIF
c
            if (.not.ostat)then
               if(opg_root())then
                  call ma_summarize_allocated_blocks
                  write(6,*)'allocation failed:', words, 'words'
               endif
               call caserr('allocating heap ')
            endif
c
            if(opg_root() .and. ogmem_debug )
     &           write(6,*)'allocate ',tag(1:10),
     &           ' size=',words,' handle= ',
     &           itag_heap(numheap),'  gamess address=', 
     &           iq_heap(numheap)

      
            handle_alloc=i
            return
         endif
      enddo
c
c     Run out of handles
c
      call dlc_error(
     & 'Run out of handles for memory allocation','abort')
c for compiler
      handle_alloc=-1

      end


      integer function handle_free(handle)
c
c     Frees up memory and handle
c
      implicit none
INCLUDE(../m4/common/gmemdata)
      include 'matrix_ma.fh'
      integer handle

      MA_INTEGER itag_8
      logical ostat, opg_root
      external opg_root

#include "mafdecls.fh"


      if(handle.le.0.or.handle.gt.maxhandle)
     &call dlc_error('Non valid handle to handle_free','abort')


      itag_8 = heap_tag(handle)
      ostat = ma_free_heap(itag_8)

      if(opg_root() .and. ogmem_debug)
     &     write(6,*)'free heap memory handle= ',heap_tag(handle)
      if (.not.ostat)call caserr('problem freeing memory ')
      hp(handle)=-1
      handle_free = 0
      return
      end

_ELSE

c
c memory.f
c
c
c     Dymamic memory allocation method
c     Alexander J Turner
c     Junary 1997
c
c     Heap memory
c     ===========
c
c     Method works by having a linked list approach
c     Each memery allocation has a record of the
c     next allocation point, its own length 
c
c     [       record one        ]                 [      record two etc       ]
c      -----------------------------------------------------------------------
c     | next | length | data ...| < empty space > | next  | length | data ....|  
c      -----------------------------------------------------------------------
c
c     So a typical record will have the pointer to the start of the next 
c     record, the length of the data, and the data.  Then there may be some
c     free space before the next record.
c
c     As the list is tranversed from start to end any gaps are found by a
c     mismatch between the length and next point records. 
c
c     If the gap is long enough, a new record is inserted
c     its next record field is set to the value of the 
c     previous record's one and the previous records one 
c     is set to the start of the new record.
c
c     The last record of the list always has the MAXSIZE paramter (IE the size of
c     availible memory) as its next record.  If this is reached, and the gap between
c     the end of the last record and maxsize < length of the data to be inserted
c     no space was found on the heap and an exception is thrown.
c
c     Upon deletion of a record, the previous record must
c     have its next record field set to the removed records
c     next record field.  For simplicity The algorithm traverses
c     the linked list to find the previous record.
c
c     If the whole list is traversed before the pointer to the deleted
c     record is found, this means that the pointer is invalid and an
c     exception is thrown
c
c     The memory is in double precision.  Integer, logical and char allocation
c     is made into the double space
c
c     Stack memory
c     ============
c
c     This allocates from the top down.  The bottom of the stack is the top
c     of the heap, the keep track of one another via stackbottom and heaptop
c     variables.
c
c     Disk dumping
c     ============
c
c     Both the Heap and Stack and the structure of them both can be dumped and
c     retrived to/from disk.  The perpose of this is that the optimiser may
c     want two use a large ammount of memory, when it is not working, IE betweem
c     cycles when a function is being eveluated, all the optimiser memory can be
c     freed up and all data stored on disk, for retievel when it is neaded again.
c
c     Methods
c     =======
c      subroutine memory_create()
c         Initiate memory system
c
c      subroutine ffree(point)
c         Free heap memory with pointer 'point'
c
c      subroutine fcopy(in,out)
c         Copy heap memory from record in to record out 
c
c      subroutine finitialise(in,value)
c         Initialise a head memory record to the give value    
c
c      subroutine fcompress()
c         Defragment heap memory and reset handles
c
c      subroutine fmemory_error(size)
c         General error reprot routine
c
c      subroutine memory_dump(name)
c         Write head memory to disk
c
c      subroutine memory_undump(name)
c         Retrieve heap memory from disk
c
c      subroutine handle_free(handle)
c         Free heap memory pointed to by a given handle
c
c      subroutine handle_copy(in,out)
c         Copy a heap record from and area pointed to by one handle
c         to an area pointed to by another
c
c      subroutine handle_initialise(handle,value)
c         Initialise a section of heap memory pointed to by a handle to
c         the given value
c
c      subroutine handle_compress()
c         Calls fcompress()
c
c      subroutine handle_print()
c         Prints out memory based on handles
c
c      integer function fmalloc(size)
c         Allocate to heap, return a pointer
c
c      integer function handle_alloc(in)
c         Allocate to heap, return a handle
c
c      integer function stack_alloc(in)
c         Allocate to stack and return a pointer
c
c      integer function stack_free(point)
c         Free the stack pointed to by given pointer
c

      subroutine memory_create()
c
c     Initialises heap
c
      implicit none
      integer i
      include 'matrix.fh'
      logical opg_root

      if(opg_root())then
         write(6,'(/,a,/,/,5x,a,i14,/)') 
     &        'Setting up Heap','Size in bytes= ',SIZE_DOUBLE*MAXSIZE
      endif
c
c     Set up length and next record fields
c
      qm(1)=0
      qm(2)=2
c
c     Set up stack top and heap top
c
      stackbottom=maxsize
      heaptop =1
c
c     Set up handle array - all handles start of as -1
c
      do i=1,maxhandle
         hp(i)=-1
      enddo
c
c     Flag a clean heap
c
      heapclean=.true.
      return
      end
c
c-------------------------------------------------------------
c
      integer function fmalloc(size)
c
c     Creates a record or throws and error if
c     it cannot
c
c     It allocates size array elements
c
      implicit none
      include 'matrix.fh'
      integer size,i,j,k,next,length,pos
c      call timeSet()

      if(MEM_DEBUG) write(6,'(a,i8)') 
     &' Handle_alloc: ',size

c
c     First test the special case the first record in the
c     list is empty (next=0)
c
      if(qm(1).eq.0) then
         if(size+2.gt.stackbottom) call fmemory_error(size)
         qm(1)=stackbottom
         qm(2)=size+2
         fmalloc=3
         heaptop=1
         return
      endif
c
c     Traverse the linked list looking for a hole big enough
c
c     Or if heapclean, go straight to end of list
c

      if(heapclean) then
         pos=heaptop
      else
         pos=1
      endif
 100  continue
         length=qm(pos+1)
         next  =qm(pos) 
c
c        Is there big enough hole?
c
         if(next-pos-length.ge.size+2) then
            qm(pos)  =pos+length
            pos      =pos+length
            qm(pos)  =next
            qm(pos+1)=size+2 
            fmalloc=pos+2
c
c           Move the heaptop if allocation is at heaptop
            if(next.eq.stackbottom) heaptop=fmalloc-2
c            call timePrint('fmalloc')
            return
         endif
c
c        If there has been no bigenough hoels
c        and next=stackbottom, we have no room!
c
         if(next.eq.stackbottom) then
             call fmemory_error(size)
         endif
         if(next-pos.le.0.or.next.gt.stackbottom) then
            call handle_print()
            call dlc_error 
     &      ('Memory is damanged, check code','abort')
         endif
         pos=next
      goto 100
      end

c
c-------------------------------------------------------------
c
      subroutine ffree(point)
c
c     Frees up integer memory the pointer
c     to which is point
c
      implicit none
      include 'matrix.fh'
      integer pos,point,i,j,k,oldl
c      call timeSet()

c
c     Traverse the linked memory list to find the
c     record pointed to by point.
c
c     First reduce point by 2 so it points to the 
c     record next field
c
      point=point-2
      if(MEM_DEBUG) write(6,'(a,I8)') 
     &' Handle_free:    ',int(qm(point+1))-2
c
c     Special case, if the first record is to be removed
c     just set it to zero length (lenght=2 for next and lenght records!)
c
      if(point.eq.1) then
         qm(2)=2
         heaptop=1
         return
      endif
      pos=1
      oldl=qm(pos+1)
 100  continue
         if(qm(pos).eq.point) goto 200
         if(qm(pos).eq.stackbottom) 
     &      call dlc_error 
     &      ('Pointer for freeing memory invalid','abort')
            if(qm(pos)-pos.le.0.or.qm(pos).gt.stackbottom) 
     &      then
             call handle_print()
             call dlc_error 
     &       ('Memory is damanged, check code','abort')
            endif
         pos=qm(pos)
      goto 100
 200  continue
c
c     Free memory by moving the 'next' field of
c     the previous record to the next field of free record
c
      qm(pos)=qm(point)
c
c     Reset the top of the heap if necessary
c
      if(point.eq.heaptop) then
         heaptop=pos
      else
c
c        Flag a dirty heap
c
         heapclean=.false.
      endif

c      call timePrint('ffree')
      return
      end
c
c-------------------------------------------------------------
c
      subroutine fcopy(in,out)
c
c     Copies all of data in in to out
c
      implicit none
      include 'matrix.fh'
      integer pos,point,i,j,k,oldl,in,out

      if(MEM_DEBUG) write(6,'(a,I8)')
     &' Handle_copy in: ',int(qm(in-1))-2
      if(MEM_DEBUG) write(6,'(a,I8)')
     &' Handle_copy out:',int(qm(in-1))-2
c
c     error check lenghts
c
      if(qm(in-1).gt.qm(out-1))
     &   call dlc_error('Copy error in bigger than out','abort')
c
c     perform copy
c
      do i=0,qm(in-1)-3
        qm(out+i)=qm(in+i)
      enddo
      return
      end
c
c-------------------------------------------------------------
c
      subroutine finitialise(in,value)
c
c     Set memory  area pointed to by in to value
c
      implicit none
      include 'matrix.fh'
      integer pos,point,i,j,k,oldl,in,out
      double precision value

      if(MEM_DEBUG) write(6,'(a,I8,x,G13.5)')
     &' Handle_initialise : ',int(qm(in-1))-2,value
c
c     perform copy
c
      do i=0,qm(in-1)-3
        qm(in+i)=value
      enddo
      return
      end

c
c-------------------------------------------------------------
c
      subroutine fcompress()
c
c     Alexander J Turner - Jan 1998
c
c     Compresses heap to zero fragmentation
c
      implicit none
      include 'matrix.fh'
      integer i,j,k,pos,next,length,step,handle,gap

c      call timeSet()
      gap=0
 
      if(MEM_DEBUG) write(6,'(a)') ' Compressing heap'
      if(qm(1).eq.0) return
      if(qm(1).eq.stackbottom) return

      pos=1
100   continue
         next=qm(pos)
         length=qm(pos+1)
c
c        See if we have gone over the whole heap
c
         if(next.eq.stackbottom)  then
c
c           Yes - reset heaptop pointer,heapclean and return
c
            heapclean=.true.
            heaptop=pos
           if(MEM_DEBUG) 
     &        write(6,'(a,i12)')' Total saving ',step

c            call timePrint('fcompress')
            return
         endif
c
c        Test if there is empty space at the end of this record
c
         step=next-pos-length
         if(MEM_DEBUG) then
            if(step.ne.gap) then
               write(6,'(a,i12)')' Found gap of ',step-gap
               gap=step
            endif
         endif
         if(step.lt.0)then
            write(6,'(a,i12)')' Step =',step
            write(6,'(a,i12)')' Pos  =',pos 
            write(6,'(a,i12)')' Next =',next
            call dlc_error('Heap corrupted - check code','abort')
         endif
       
         if(step.gt.0) then
c
c           There is, so move the next record down to remove it
c           Use length record from next element in list to
c           give lentgh of loop required
c
c           Also scan handles for a match, and reset, error if
c           no match found
c
            j=pos+length
            k=next
            handle=k+2
            do i=0,qm(k+1)-1
              qm(j+i)=qm(k+i)
            enddo
            do i=1,maxhandle
               if(hp(i).eq.handle) then
                  hp(i)=j+2
                  goto 200
               endif
            enddo
            write(6,'(a,i9,/)')
     &      'The is a problem compressing memory at pointer ',handle
            call handle_print()
            call dlc_error
     &         ('Heap and Handles not aligned - check code','abort')

200         continue
c         
c           Reset next record of this list element
c
            qm(pos)=j
         endif
c
c           Iterate
c
         pos =qm(pos)
      goto 100
      end
            

c
c-------------------------------------------------------------
c

      subroutine fmemory_error(size)
c
c     Handle errors due to insufficiant space in heap
c
      implicit none
      include 'matrix.fh'
      integer size,next,length,pos,waist
      call dlc_error('Insufficiant space on heap','warn')
      write(6,'(a,i9,a)')
     &  'Failed asking for ',size*SIZE_DOUBLE,' bytes'
      call handle_print()
      call dlc_error('memory error','abort')
      end

c
c-------------------------------------------------------------
c

c
c     Routines to move the heap on an off disk
c
c     Also records internal data structures for matix object
c

      subroutine memory_dump(name)
c
c     Writes a dump of the heap into file 'name'
c
      implicit none
INCLUDE(matrix_internal.fh)
      character*(*) name
      integer i,j,k,pos

      OPEN(FILE=name,FORM='UNFORMATTED',STATUS='new',
     $    ACCESS='SEQUENTIAL',ERR=100,UNIT=9)
      write(6,'(a,/)')'Writing memory dumpfile'
      goto 200

 100  OPEN(FILE=name,FORM='UNFORMATTED',STATUS='old',
     $    ACCESS='SEQUENTIAL',ERR=2000,UNIT=9)
      write(6,'(a,/)')'Over-writing memory dumpfile'

 200  continue
      rewind(9)
c
c     Write out heap
c
      pos=1
      write(9) maxsize
      write(9) stackbottom
 300  continue
         write(9) qm(pos),qm(pos+1)
         write(9) (qm(pos+i),i=2,qm(pos+1)-1)
         if(qm(pos).eq.stackbottom) goto 400
         pos=qm(pos)
      goto 300
 400  continue
c
c     Write out stack
c
      write(9)(qm(i),i=stackbottom,MAXSIZE)
c
c     Write out handles and controls
c
      write(9)(hp(i),i=1,maxhandle)
      write(9)heapclean,heaptop
c
c     Write out matrix stuff
c
      write(9)(matrix_tag          (i),i=1,max_matrices)
      write(9)(matrix_type         (i),i=1,max_matrices)
      write(9)(matrix_handle       (i),i=1,max_matrices)
      write(9)(number_of_dimensions(i),i=1,max_matrices)
      do j = 1, max_dimensions
         write(9)(dimensions (j,i)     ,i=1,max_matrices)
      enddo
c
c      Close and return
c
      close(9)
      return
 2000 call dlc_error('Dump file error','abort')
      end
c
c-------------------------------------------------------------
c
      subroutine memory_undump(name)
c
c     Read a dump of the heap from file 'name'
c
      implicit none
INCLUDE(matrix_internal.fh)
      character*(*) name
      integer i,j,k,pos,fmax,stb

      write(6,'(a,/)')'Reading memory dumpfile'
      OPEN(FILE=name,FORM='UNFORMATTED',STATUS='old',
     $    ACCESS='SEQUENTIAL',ERR=100,UNIT=9)
      goto 200
 100  OPEN(FILE=name,FORM='UNFORMATTED',STATUS='new',
     $    ACCESS='SEQUENTIAL',ERR=2000,UNIT=9)

 200  continue
      rewind(9)
c
c     Read in heap
c
      pos=1
      read(9) fmax
      read(9) stackbottom
      if(fmax.gt.maxsize) 
     &  call dlc_error('Heap in file > present heap in memory','abort')
      if(fmax.lt.maxsize)
     &  call dlc_error('Heap in file mismatchs heap in memory','warn')
 300  continue
         read(9) qm(pos),qm(pos+1)
         read(9) (qm(pos+i),i=2,qm(pos+1)-1)
         if(qm(pos).eq.stackbottom) goto 400
         pos=qm(pos)
         if(pos+qm(pos+1).gt.stackbottom)
     &  call dlc_error
     &  ('Heap in file > present availible memory','abort')
      goto 300
 400  continue
c
c     Read in  stack
c
      read(9)(qm(i),i=stackbottom,MAXSIZE)
c
c     Read in handles and controls
c
      read(9)(hp(i),i=1,maxhandle)
      read(9)heapclean,heaptop
c
c     Read in matrix stuff
c
      read(9)(matrix_tag          (i),i=1,max_matrices)
      read(9)(matrix_type         (i),i=1,max_matrices)
      read(9)(matrix_handle       (i),i=1,max_matrices)
      read(9)(number_of_dimensions(i),i=1,max_matrices)
      do j = 1, max_dimensions
         read(9)(dimensions (j,i)     ,i=1,max_matrices)
      enddo
c
c      Close and return
c
      close(9)
      return
 2000 call dlc_error('Dump file error','abort')
      end

c-------------------------------------------------------------
c

c     HANDLES - Alexander J Turner - Jan 1998
c
c     Allocates handles in array hp to memory
c     allocated in array q
c
c     And frees this memory
c

c
c         ****** ALLOCATE *******
c
      integer function handle_alloc(in)
c
c     Allocates the memory and places
c     a pointer to it in hp, then
c     returns the index in hp of that pointer
c
      implicit none
      include 'matrix.fh'
      integer i,j,k,size,fmalloc,in
      external fmalloc

c
c     find the ammount of memory to allocate
c
      size=in
      size=int(real(size)/real(SIZE_DOUBLE))+1

      do i=1,maxhandle
         if(hp(i).eq.-1) then
c            if(.not.heapclean) call  fcompress()
            hp(i)=fmalloc(size)
            handle_alloc=i
            return
         endif
      enddo
c
c     Run out of handles
c
      call dlc_error('Run out of handles for memory allocation','abort')
c for compiler
      handle_alloc=-1
      end

c
c         ******   FREE   *******
c
      subroutine handle_free(handle)
c
c     Frees up memory and handle
c
      implicit none
      include 'matrix.fh'
      integer handle
      if(handle.le.0.or.handle.gt.maxhandle)
     &call dlc_error('Non valid handle to handle_free','abort')
      call ffree(hp(handle))
      hp(handle)=-1
      return
      end

c
c         ****** COPY *******
c
      subroutine handle_copy(in,out)
c
c     Copies one section of the handle memory 
c     to another
c
      implicit none
      include 'matrix.fh'
      integer i,j,k,size,fmalloc,in,out
      external fmalloc

      if(in.le.0.or.in.gt.maxhandle)
     &call dlc_error('Non valid handle to handle_free','abort')
      if(out.le.0.or.out.gt.maxhandle)
     &call dlc_error('Non valid handle to handle_free','abort')
c
c     Dereference and call fcopy to do the work
c
      call fcopy(hp(in),hp(out))
      return
      end
c
c         ****** INITIALISE *******
c
      subroutine handle_initialise(handle,value)
c
c     Set a chunk of handle memory to the supplied value
c
      implicit none
      include 'matrix.fh'
      integer handle
      double precision value
            if(handle.le.0.or.handle.gt.maxhandle)
     &call dlc_error('Non valid handle to handle_free','abort')
c
c     Dereference to finitialise to do work
c
      call finitialise(hp(handle),value)
      return
      end


c
c         ****** COMPRESS *******
c
      subroutine handle_compress()
      call fcompress()
      return
      end
c
c         ******  PRINT   *******
c
      subroutine handle_print()
c
c     Prints out the structure of the memory for debugging
c
      implicit none
      include 'matrix.fh'
      integer i,pos,length

      write(6,'(a,i9,/)')
     &'The stack bottom is at ',stackbottom
      write(6,'(a,i9,/)')
     &'The heap top is at     ',heaptop
      write(6,'(a,/)') 'Print out of handle memory'
      do i=1,maxhandle
         if(hp(i).ge.0) then
            write(6,'(5x,a,i10,a,i10)')
     &      'Handle:   ',i,' Pointer: ',hp(i)
_IFN(t3e,cray)
            call flushout()
_ENDIF
            write(6,'(5x,a,i10,a,i10)')'->Length: ',
     &      int(qm(hp(i)-1)),' Next:    ',int(qm(hp(i)-2))
_IFN(t3e,cray)
            call flushout()
_ENDIF
         endif
      enddo
      write(6,'(/,/)')
      return
      end

c
c-------------------------------------------------------------
c

c     STACK - Alexander J Turner - Jan 1998
c

c
c         ****** ALLOCATE *******
c

      integer function stack_alloc(in)
      implicit none
      include 'matrix.fh'
      integer i,j,k,size,room,in
c
c
c     Set size to the correct double pricion records
c     add one for round off error and one for the length
c     record
c
      size=in
      size=int(real(size)/real(SIZE_DOUBLE))+2
c
c     Compute the space between the top of the heap and
c     the bottom of the stack
c
      room=stackbottom-(heaptop+qm(heaptop+1))-1

      if(MEM_DEBUG)write(6,'(a,i12,a,i12,a,i12)')
     &'Stack_alloc: STB ',
     &stackbottom,' Room  ',room,' Length ',size
c
c    If not enough room - error
c
      if(size.gt.room) 
     &     call dlc_error('Insufficiant room on stack','abort')
c
c     return present bottom of stack as start of new record
c
      stack_alloc=stackbottom-size+2
c
c     Move bottom of stack down
c
      stackbottom=stackbottom-size
c
c     Record the length of this section
c
      qm(stackbottom+1)=size
      qm(heaptop)=stackbottom
      return
      end
      
c
c         ******   FREE   *******
c

      integer function stack_free(point)
      implicit none
      include 'matrix.fh'
      integer i,j,k,point,length    

c
c     Get length of current record
c
      length=int(qm(stackbottom+1))

      if(MEM_DEBUG)write(6,'(a,i12,a,i12,a,i12)')
     &' Stack_free: STB ',
     &stackbottom,' Point ',point,' Length ',length
c
c     Check that pointer points to this record' top
c
      if(point.ne.stackbottom+2) 
     &   call dlc_error('Non stack memory call to stack_free','abort')
c
c     Move up stack bottom
c
      stackbottom=stackbottom+length 
      qm(heaptop)=stackbottom

      stack_free = 0

      return
      end
_ENDIF
