      subroutine init_nxtval_mpi2

      use mpi_global_counter

      implicit none

      call set_tasks

      end

      subroutine finalize_nxtval_mpi2

      use mpi_global_counter

      implicit none

      call end_tasks

      end

      integer function nxtval_mpi2( flag )

      use mpi_global_counter

      implicit none

      integer flag

      integer value

      if( flag < 0 ) then
         call reset_tasks
         value = 0
      else
         call get_task( value )
      end if

      nxtval_mpi2 = value

      end

      
