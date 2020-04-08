/*

   Plug and play type wrapper for the PeIGS library routine pspev 
  so that it can be used directly in GAMESS-UK

   The major difficulty is that GAMESS gives the Hamiltonian matrix
  in a packed format which is incompatible with the packed format
  wanted by the PeIGS routine PSPEV. The format conversion is done
  mostly in a_to_a_alias and is described in some detail in that routine.
  Most of the rest is to do with setting up the dynamic memory required
  by PeIGS and subsequent re-replication of the evecs after the diag.

*/
 
#include <stdlib.h>
#include <string.h>
#ifdef PROFILE_DIAG
#include <time.h>
#endif

#ifdef CRAYXX
#define gms_pdiag_ GMS_PDIAG
#define ipg_nnodes_ IPG_NNODES
#define ipg_nodeid_ IPG_NODEID
#define reconstitute_evecs_ RECONSTITUTE_EVECS
#endif



/* HACK - These should come from the Makefile eventually */
/*      - Another ugly hack. Taking the explicit declaration of EXT_INT
          out would break every i8 build. So instead I am going to 
          suppress EXT_INT if the makefile provides a STD_INT. This
          way all old builds keeping working. The only proper way out
          of this mess is to provide either EXT_INT or STD_INT from
          the makefile in all cases and all builds. 
          Huub van Dam, July 2006 */
/*      - The moment of truth has arrived. It turns out that one cannot 
          rely on STD_INT being specified when required. So this leaves only
          the nuclear option: you have to specify what you want STD_INT or
          EXT_INT otherwise the code will not compile!
          Huub van Dam, June 2009 */

#if !defined(STD_INT)
#if !defined(EXT_INT)
 You have to specify which integer type you want, either:
 - STD_INT - normally for int
 - EXT_INT - normally for long
 Make sure the specification is consistent with the PeIGS library.
 Otherwise there will be no executable... 
#endif
#endif
#define STD_DBL

#include "macdecls.h"

/* A PeIGS header. This typedefs appropriate datatypes to ease
   the passing of integers and reals between Fortran and C 
   PS Hack - MA_DEFINES_TYPES suppresses redefinitions from
   macdecls.h
*/ 

#define MA_DEFINES_TYPES

#include "globalp.c.h"

/* Do some things need initializing */
static short int initialized = 0;

/* Data about who I am and how many mates I have */
static Integer me, nodes;

/* The prototypes */

/* The main wrapper function */
void    gms_pdiag_( Integer *, DoublePrecision *, DoublePrecision *, 
               DoublePrecision * , Integer *);

/* This generates the arrays that tell the diag how the matrix and
   its evecs are to be distributed across the processors */
static Integer generate_map( Integer, Integer * );

/* One of the major problems in all this is that GAMESS provides us
   with a packed UPPER triangle, while PeIGS wants a
   packed LOWER triangle. This routine solves this in
   a slightly sneaky way, and also avoids the need for any scratch
   space */
static void    a_to_a_alias( Integer, Integer *, 
			     DoublePrecision *, DoublePrecision ** );

/* Get the evec array in a form suitable for PeIGS */
static void    evec_to_evec_alias( Integer, Integer *, 
			           DoublePrecision *, 
				   DoublePrecision ** );

/* PeIGS generates the eval/evec pairs in a slightly clumsy
   ordering. Reorganize things so that GAMESS doesn't know anything
   has changed */
static void    reorganize_evals( Integer, Integer *, Integer *,
				 DoublePrecision *, Integer *ierr );

/* GAMESS wants its evecs replicated, regenerate them. This is actually
   a Fortran subroutine. Why ? To call the internal global sumer
   in GAMMES would require that nightmare scenario, passing characters
   between C and Fortran. Hence a Fortran wrapper routine */
extern void    reconstitute_evecs_( Integer *, DoublePrecision *);

/* At certain points need to turn vectors upside down */

static void    upside_downify_vector( Integer, DoublePrecision* );

void gms_pdiag_( Integer *order, DoublePrecision *a, 
            DoublePrecision *evecs, DoublePrecision *evals,
	    Integer *ierr)
{

  DoublePrecision **a_alias, **evecs_alias;
  DoublePrecision *real_work1, **real_work2;
  Integer         *map_a, *map_evecs;
  Integer         n;
  Integer         i_hold;
  Integer         *int_work;
  Integer         size_int_work, size_real_work1, size_real_work2;
  Integer         method = 1;
  Integer         info;
  Integer         i;
#ifdef PROFILE_DIAG
  float           freq = 1.0 / CLOCKS_PER_SEC;
  Integer         t_start, t_end;
  float           t_init, t_diag, t_reorg, t_permute, t_replicate;
#endif

  Integer  i_int_work, h_int_work;
  Integer  i_real_work1, h_real_work1;
  Integer  i_real_work2, h_real_work2;
  Integer  i_map_a, h_map_a;
  Integer  i_map_evecs, h_map_evecs;
  Integer  i_a_alias, h_a_alias;
  Integer  i_evecs_alias, h_evecs_alias;
  Integer  ptr_to_long;
  
  extern Integer         ipg_nodeid_( void ), ipg_nnodes_( void );

#ifdef PROFILE_DIAG
  t_start = clock();
#endif

  *ierr = 0;

/* Get the matrix order dereferenced to make life a bit easier */
  n = *order;

  /* printf("In Peigs interface int len = %d Order = %d\n",
	 sizeof(Integer), n); */

/* If necessary initalize */
  if( !initialized ){
    initialized = 1;
    me    = ipg_nodeid_();
    nodes = ipg_nnodes_();
  }

  /* Set up the mapping arrays */

  if(!MA_alloc_get(MT_F_INT, n, "map_a", &h_map_a, &i_map_a)){*ierr=-100; return; }
  if(!MA_get_pointer(h_map_a,&map_a)){*ierr=-2; return;}

  if(!MA_alloc_get(MT_F_INT, n, "map_evecs", &h_map_evecs, &i_map_evecs)){*ierr=-101; return; }
  if(!MA_get_pointer(h_map_evecs,&map_evecs)){*ierr=-2; return;}

  i_hold = generate_map( n, map_a );

/* The evecs will have the same mapping as the original matrix */

  for( i = 0; i < n; i++ ){

    map_evecs[ i ] = map_a[ i ];
  }


  ptr_to_long = sizeof( DoublePrecision * ) / sizeof ( long );


/* Set up the arrays of pointers to the appropriate bits of
   the matrix and evec arrays, whilst also changing the
   data packing format of a */

  if(!MA_alloc_get(MT_LONGINT, ptr_to_long*n, "a_alias", &h_a_alias, &i_a_alias)){*ierr=-102; return;}
  if(!MA_get_pointer(h_a_alias,&a_alias)){*ierr=-2; return;}

  a_to_a_alias( n, map_a, a, a_alias );

  /* There seems to be no pointer allocation in the MA tools, so we need
     to know how many characters to allow for a pointer type */

  if(!MA_alloc_get(MT_LONGINT, ptr_to_long*n, "evecs_alias", &h_evecs_alias, &i_evecs_alias)){*ierr=-103; return;}

  if(!MA_get_pointer(h_evecs_alias,&evecs_alias)){*ierr=-2; return;}

  evec_to_evec_alias( n, map_evecs, evecs, evecs_alias );
  
/* Use the PeIGS utility to find out how much scratch space
   it will need */

  if(!MA_alloc_get(MT_F_INT, 3*n, "int_work", &h_int_work, &i_int_work)){*ierr=-104; return; }
  if(!MA_get_pointer(h_int_work,&int_work)){*ierr=-2; return;}

  memreq_( &method,  &n, map_a, map_evecs, map_a,
	   &size_int_work, &size_real_work1, &size_real_work2,
	   int_work );



  if(!MA_free_heap(h_int_work)){*ierr=-3; return;}

/* And set up that scratch space */

  if(!MA_alloc_get(MT_F_INT, size_int_work, "int_work", &h_int_work, &i_int_work)){*ierr=-105; return; }
  if(!MA_get_pointer(h_int_work,&int_work)){*ierr=-2; return;}

  if(!MA_alloc_get(MT_DBL, size_real_work1, "real_work1", &h_real_work1, &i_real_work1)){*ierr=-106; return; }
  if(!MA_get_pointer(h_real_work1,&real_work1)){*ierr=-2; return;}

  if(!MA_alloc_get(MT_LONGINT, ptr_to_long*size_real_work2, "real_work2", &h_real_work2, &i_real_work2)){*ierr=-107; return; }
  if(!MA_get_pointer(h_real_work2,&real_work2)){*ierr=-2; return;}

/* Need to zero the evecs array for the reconstitution stage */

  for( i = 0; i < n * n; i++ ){
    evecs[ i ] = ( DoublePrecision ) 0.0;
  }

/* At last, the diag */

#ifdef PROFILE_DIAG
  t_end   = clock();
  t_init  = freq * ( t_end - t_start );
  t_start = clock();
#endif

  pdspev( &n, a_alias, map_a, evecs_alias, map_evecs, evals,
	  int_work  , &size_int_work,
	  real_work2, &size_real_work2,
	  real_work1, &size_real_work1, &info );

#ifdef PROFILE_DIAG
  t_end   = clock();
  t_diag  = freq * ( t_end - t_start );
  t_start = clock();
#endif

/* Clear up the workspace arrays */

  if(!MA_free_heap(h_int_work)){*ierr=-3; return;}
  if(!MA_free_heap(h_real_work1)){*ierr=-3; return;}
  if(!MA_free_heap(h_real_work2)){*ierr=-3; return;}

/* Tidy up the evals  - This may in fact not be necessary, I'll have
   to think it through a bit more. However it is very quick so no major
   damage is done */

  reorganize_evals( n, map_a, map_evecs, evals, ierr );  

  if(*ierr != 0)return;

/* Don't need the maps anymore */

  if(!MA_free_heap(h_map_a)){*ierr=-3; return;}
  if(!MA_free_heap(h_map_evecs)){*ierr=-3; return;}

#ifdef PROFILE_DIAG
  t_end   = clock();
  t_reorg = freq * ( t_end - t_start );
  t_start = clock();
#endif

/* Permuting the evecs back to their proper form is so easy
   can't be bothered to write a seperate function */

  for( i = 0; i < i_hold; i++ ){
    upside_downify_vector( n, evecs_alias[ i ] );
  }

#ifdef PROFILE_DIAG
  t_end     = clock();
  t_permute = freq * ( t_end - t_start );
  t_start   = clock();
#endif

/* And make the evecs replicated. Remember this is a Fortran routine
   so need to pass pointers */

  reconstitute_evecs_( &n, evecs ); 

  if(!MA_free_heap(h_a_alias)){*ierr=-3; return;}
  if(!MA_free_heap(h_evecs_alias)){*ierr=-3; return;}

#ifdef PROFILE_DIAG
  t_end       = clock();
  t_replicate = freq * ( t_end - t_start );

  if( me == 0 ){
    printf( "Node 0 internal profiling for the parallel diag: \n" );
    printf( "Processor times are in seconds.\n" );
    printf( " Initialization     : %8.4f\n", t_init  );
    printf( " Diagonalization    : %8.4f\n", t_diag  );
    printf( " Eval reorganization: %8.4f\n", t_reorg );
    printf( " Evec permutation   : %8.4f\n", t_permute );
    printf( " Evec replication   : %8.4f\n", t_replicate );
  }
#endif

}

/* Map generation function */

static Integer generate_map( Integer n, Integer *map )
{
  Integer location;
  Integer i_hold;
  Integer i;

  i_hold = 0;

/* Simple cyclic distribution for present, might be worth investigating
   if other forms are better though */

  for( i = 0; i < n; i++ ){
    location = i % nodes;
    map[ i ] = location;
    if( location == me ){
      i_hold++;
    }
  }

  return i_hold;

}

/* Set up the arrays of pointers */

static void evec_to_evec_alias( Integer n, Integer *map, 
			        DoublePrecision *evecs,
			        DoublePrecision **evecs_alias )
{

  Integer my_next;
  Integer i;

  my_next = 0;

/* Whack through the memory allocated for the Fortran array. Look
   at each column, and if I am meant to own it add it to the array
   of pointers. */

  for( i = 0; i < n; i++ ){
    if( map[ i ] == me ){
      evecs_alias[ my_next ] = evecs;
      my_next++;
    }
    evecs += n;
  }

}


/* eval order rationalization */

static void reorganize_evals( Integer n, Integer *map_a, 
			      Integer *map_evecs, 
			      DoublePrecision *evals, Integer *ierr )
{

  DoublePrecision *sorted_evals;
  Integer         *search_starts;
  Integer          target;
  Integer          i, j;

  Integer          h_search_starts, i_search_starts;
  Integer          h_sorted_evals, i_sorted_evals;

/* This will hold the re-ordered evals */

  if(!MA_alloc_get(MT_DBL, n, "sorted_evals", &h_sorted_evals, &i_sorted_evals)){*ierr=-108; return; }
  if(!MA_get_pointer(h_sorted_evals,&sorted_evals)){*ierr=-2; return;}

/* search_starts stops us searching though the wole array of evals.
   If the last eval this node had an eigenvector for was number m,
   it need only search the array from location m + 1 onward */
  
  if(!MA_alloc_get(MT_F_INT, nodes, "search_starts", &h_search_starts, &i_search_starts)){*ierr=-109; return; }
  if(!MA_get_pointer(h_search_starts,&search_starts)){*ierr=-2; return;}

  for( i = 0; i < nodes; i++ ){
    search_starts[ i ] = 0;
  }

  for( i = 0; i < n; i++ ){

/* For the i'th eval/evec pair find out on which node the corresponding
   column of the original matrix was upon. This is the column of the
   replicated evec matrix where the i'th evec will end up */
    target = map_a[ i ];

/* Now search for the eval which corresponds to this evec */
    for( j = search_starts[ target ]; j < n; j++ ){
      if( map_evecs[ j ] == target ){
	break;
      }
    }

/* Update the search_starts and store the eval in its proper place */
    search_starts[ target ] = j + 1;
    sorted_evals[ i ] = evals[ j ];
  }

/* Copy the sorted evals back */

  evals = memcpy( evals, sorted_evals, n * sizeof( DoublePrecision ) );

/* Tidy up */

  if(!MA_free_heap(h_sorted_evals)){*ierr=-3; return;}
  if(!MA_free_heap(h_search_starts)){*ierr=-3; return;}

}

/* Convert the replicated, packed upper triangle usable by Fortran
   to a distributed, packed lower triangel better suited to C */ 

static void a_to_a_alias( Integer n, Integer *map_a, 
			  DoublePrecision *a, DoublePrecision **a_alias )
{
  Integer length, my_next;
  Integer i;

/* Now this is where things get sneaky ! The data conversions occurs
   conceptually in two stages:

   1) Transfer the relevant columns of the fortran matrix to the C form 
   IN REVERSE ORDER. ( When I say transfer I really  mean do a whole load 
   of pointer assignments, but coming from a fortran background I find 
   thinking in terms of real data transfers easier. ) 

   2) Turn the column upside down 

   Now you scream Ian, you thick divvy, we havent got the same matrix
   as before so all the evals and evecs will be wrong so it won't work
   and it's all your fault.' True it ain't the same matrix, BUT
   1) It's only a similarity transform away so the evals will be correct
   2) I know what the similarity transform is so at the end I can use it
      to generate the correct evecs from those found by the diag
   3) The matrix that does the similarity transform is extremely simple
      and applying it is quick, you simply turn the matrix upside down
      again 

   Also note the whole data conversion is ( almost ) totally parallel */

/*  So first move to the start of the last column of the Fortran matrix */

  a += n * ( n - 1 ) / 2;

/* Length hold the length of the current column, my_next keeps tabs on
   where the next bit of data I will look after will be stored ( see note
   on pointers above ) */

  length  = n;
  my_next = 0;

/* Look at each column in turn, storing data about it if I need it */

  for( i = 0; i < n; i++ ){
    if( map_a[ i ] == me ){
      a_alias[ my_next ] = a;
/* Turn the column upside down, it might be worth inlining this routine */
      upside_downify_vector( length, a_alias[ my_next ] );
      my_next++;
    }
    length--;
    a -= length;
  }

}

/* Function to turn a vector upside down. This probably could be improved 
   markedly by some loop unrolling, the data locality is abysmal. Could
   also consider calls to DSWAP with negative strides */

static void upside_downify_vector( Integer n, DoublePrecision *vector )
{
  DoublePrecision swap_temp;
  Integer         i;

  for( i = 0; i < n / 2; i++ ){
    swap_temp           = vector[ i ];
    vector[ i         ] = vector[ n - i - 1 ];
    vector[ n - i - 1 ] = swap_temp;
  }
  
}

char source[]="$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/peigs_interface.m,v $";
char revision[]="$Revision: 6131 $";
char date[]="$Date: 2010-05-17 13:17:36 +0200 (Mon, 17 May 2010) $";

void ver_peigs_interface_( char *string1, int *length1, char *string2, int *length2,char
*string3, int *length3)
{
  int l, start;
  l = strlen(source);
  start = 8;
  l-=start;
  if(l > *length1)l = *length1;
  strncpy(string1,source+start,l);

  l = strlen(revision);
  start = 10;
  l-=start;
  if(l > *length2)l = *length2;
  strncpy(string2,revision+start,l);

  l = strlen(date);
  start = 6;
  l-=start;
  if(l > *length3)l = *length3;
  strncpy(string3,date+start,l);
}
