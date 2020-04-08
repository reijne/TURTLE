! Wrapper routines for MPI
Module mpi_p_to_p_data

  Use newscf_numbers

  Implicit None

  Private

  Save

  Real( wp ), Dimension( : ), Pointer, Public :: send_ptr
  Real( wp ), Dimension( : ), Pointer, Public :: recv_ptr

  Integer, Public :: recv_request = 1
  Integer, Public :: outstanding_recv_request

End Module mpi_p_to_p_data

Module mpi_comm_data

  Implicit None

  Include 'mpif.h'

  Integer :: comm_base = mpi_comm_world

End Module mpi_comm_data

Subroutine mpi_init( error )

  Use newscf_numbers

  Implicit None

  Integer, Intent( Out ) :: error

  error = 0

End Subroutine mpi_init

Subroutine mpi_finalize( error )

  Use newscf_numbers

  Implicit None

  Integer, Intent( Out ) :: error

  error = 0

End Subroutine mpi_finalize

Subroutine mpi_abort( comm, code, error )


  Use newscf_numbers

  Implicit None

  Integer, Intent( In    ) :: comm
  Integer, Intent( In    ) :: code
  Integer, Intent(   Out ) :: error

  Stop

  error = 0

End Subroutine mpi_abort

Subroutine mpi_isend

  Write( *, * ) 'isend!!!!!!!!!!!!!!!!!!!!'

End Subroutine mpi_isend

Subroutine mpi_recv

  Write( *, * ) 'recv!!!!!!!!!!!!!!!!!!!!'

End Subroutine mpi_recv

! Communicator routines

Subroutine mpi_comm_size( comm, size, error )

  Implicit None

  Integer, Intent( In    ) :: comm
  Integer, Intent(   Out ) :: size
  Integer, Intent(   Out ) :: error

  size  = 1
  error = 0

End Subroutine mpi_comm_size

Subroutine mpi_comm_rank( comm, rank, error )

  Implicit None

  Integer, Intent( In    ) :: comm
  Integer, Intent(   Out ) :: rank
  Integer, Intent(   Out ) :: error

  rank  = 0
  error = 0

End Subroutine mpi_comm_rank

Subroutine mpi_comm_dup( comm1, comm2, error )

  Implicit None

  Integer, Intent( In    ) :: comm1
  Integer, Intent(   Out ) :: comm2
  Integer, Intent(   Out ) :: error

  comm2 = comm1
  error = 0

End Subroutine mpi_comm_dup

Subroutine mpi_comm_split( comm1, colour, key, comm2, error )

  Use mpi_comm_data

  Implicit None

  Integer, Intent( In    ) :: comm1
  Integer, Intent( In    ) :: colour
  Integer, Intent( In    ) :: key
  Integer, Intent(   Out ) :: comm2
  Integer, Intent(   Out ) :: error

  comm_base = comm_base + 1
  comm2 = comm_base

  error = 0

End Subroutine mpi_comm_split

Subroutine mpi_comm_free( comm, error )

  Integer, Intent( In    ) :: comm
  Integer, Intent(   Out ) :: error

  error = 0

End Subroutine mpi_comm_free

! Collectives

Subroutine mpi_barrier( comm, error )

  Integer, Intent( In    ) :: comm
  Integer, Intent(   Out ) :: error

  error = 0

End Subroutine mpi_barrier

Subroutine mpi_bcast( buffer, n, datatype, root, comm, error )

  Integer,                   Intent( In    ) :: n
  Integer, Dimension( 1:n ), Intent( InOut ) :: buffer
  Integer,                   Intent( In    ) :: datatype
  Integer,                   Intent( In    ) :: root
  Integer,                   Intent( In    ) :: comm
  Integer,                   Intent(   Out ) :: error

  error = 0

End Subroutine mpi_bcast

Subroutine mpi_allreduce( sendbuf, recvbuf, n, datatype, op, comm, error )

  Implicit None

  Include 'mpif.h'

  Integer,                   Intent( In    ) :: n
  Integer, Dimension( 1:n ), Intent( In    ) :: sendbuf
  Integer, Dimension( 1:n ), Intent(   Out ) :: recvbuf
  Integer,                   Intent( In    ) :: datatype
  Integer,                   Intent( In    ) :: op
  Integer,                   Intent( In    ) :: comm
  Integer,                   Intent(   Out ) :: error

  Select Case( datatype )

  Case( MPI_INTEGER )
     Call mpi_allred_int( n, sendbuf, recvbuf )
     
  Case( MPI_DOUBLE_PRECISION )
     Call mpi_allred_dp ( n, sendbuf, recvbuf )

  Case Default
     Call dlc_error( 'Unknown data type in mpi_allreduce wrapper', 'abort' )

  End Select

  error = 0

End Subroutine mpi_allreduce

Subroutine mpi_allred_int( n, a, b )

  Implicit None

  Integer,                   Intent( In    ) :: n
  Integer, Dimension( 1:n ), Intent( In    ) :: a
  Integer, Dimension( 1:n ), Intent(   Out ) :: b

  b = a

End Subroutine mpi_allred_int

Subroutine mpi_allred_dp( n, a, b )

  Use newscf_numbers
  
  Implicit None

  Integer   ,                   Intent( In    ) :: n
  Real( wp ), Dimension( 1:n ), Intent( In    ) :: a
  Real( wp ), Dimension( 1:n ), Intent(   Out ) :: b

  b = a

End Subroutine mpi_allred_dp

! MPI point to point

Subroutine mpi_send( a, n, datatype, dest, tag, comm, error )

  Use newscf_numbers
  Use mpi_p_to_p_data

  Implicit None

  Include 'mpif.h'

  Integer                             , Intent( In    ) :: n
  Real( wp ), Dimension( 1:n ), Target, Intent( In    ) :: a
  Integer                             , Intent( In    ) :: datatype
  Integer                             , Intent( In    ) :: dest
  Integer                             , Intent( In    ) :: tag
  Integer                             , Intent( In    ) :: comm
  Integer                             , Intent(   Out ) :: error

  Select Case( datatype )

  Case( MPI_DOUBLE_PRECISION )
     send_ptr => a

  Case Default
     Call dlc_error( 'Unknown data type in mpi_send wrapper', 'abort' )

  End Select

  error = 0

End Subroutine mpi_send

Subroutine mpi_irecv( a, n, datatype, source, tag, comm, request, error )

  Use newscf_numbers
  Use mpi_p_to_p_data

  Implicit None

  Include 'mpif.h'

  Integer                             , Intent( In    ) :: n
  Real( wp ), Dimension( 1:n ), Target, Intent( In    ) :: a
  Integer                             , Intent( In    ) :: datatype
  Integer                             , Intent( In    ) :: source
  Integer                             , Intent( In    ) :: tag
  Integer                             , Intent( In    ) :: comm
  Integer                             , Intent(   Out ) :: request
  Integer                             , Intent(   Out ) :: error

  Select Case( datatype )

  Case( MPI_DOUBLE_PRECISION )
     recv_ptr => a

  Case Default
     Call dlc_error( 'Unknown data type in mpi_irecv wrapper', 'abort' )

  End Select

  request = recv_request
  recv_request = recv_request + 1

  outstanding_recv_request = request

  error = 0

End Subroutine mpi_irecv

Subroutine mpi_wait( request, status, error )

  Use mpi_p_to_p_data

  Implicit None

  Include 'mpif.h'

  Integer                                , Intent( In    ) :: request
  Integer, Dimension( 1:MPI_STATUS_SIZE ), Intent(   Out ) :: status
  Integer                                , Intent(   Out ) :: error

  If( request /= outstanding_recv_request ) Then
     Call dlc_error( 'Internal error in point to point comms wrapper', 'abort' )
  End If

  If( Size( send_ptr ) /= Size( recv_ptr ) ) Then
     Call dlc_error( 'Internal error in point to point comms wrapper', 'abort' )
  End If

  recv_ptr = send_ptr

  status = 0

  error = 0

End Subroutine mpi_wait
