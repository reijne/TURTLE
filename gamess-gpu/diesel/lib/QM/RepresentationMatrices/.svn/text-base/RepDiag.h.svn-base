//***********************************************************************
//
//	Name:			RepDiag.h
//
//	Description:	stores information for an efficient calculation
//					of P=5 diagonal representation matrices
//
//					
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			08.04.1997
//
//
//
//
//
//***********************************************************************

#ifndef __RepDiag_h
#define __RepDiag_h

#include "../../../config.h"


#include "RepresentationMatrix.h"
#include "../MO/Iterators/MOListIterator.h"
#include "../Configuration/Configuration.h"
#include "HMatElements.h"


template <class RepMatType, class VectorType>
class RepDiag {
public:
	RepDiag(
		const HMatElements<RepMatType, VectorType>	* const *repMatsP5,
		Configuration<MOType> conf,
		INT minOpenShells,
		INT maxOpenShells);
	~RepDiag();

	RepMatType	getDiag(INT) const;

	//------------------------------------------------------------------------

	void	addExt(
		Configuration<MOType> conf,
		INT totalOpenShells,
		INT	moExtPos1,
		INT moExtPos2);

	void	addExtFirst(MOListIterator moiter);

	void	addExtLast(MOListIterator moiter);

	//------------------------------------------------------------------------


private:
RepMatType	*diag;								// complete diagonal

RepMatType	*DiagIntPart;						// diagonal values (to be
												// multiplied with unity matrix)
												// including internal MOs only

RepMatType	*DiagExtPart;						// pointer to array of
												// diagonal (including P=5
												// representation matrices)

RepMatType	*DiagExtPartLast;					// (last step in rs-tupel)

Configuration<MOType> intConf;					// internal part of configurations


Configuration<MOType> restConf;					// configurations without
												// running MOs

Configuration<MOType> restConfLast;				// (last step in rs-tupel)

const HMatElements<RepMatType, VectorType>	* const *repMatsP5;		// pointer to array of P=5
												// representation matrices for
												// a different number of open shells
												// (used in calculation of 
												// diagonal elements)

INT	dim;										// dimension of currently
												// contained diagonal
												
INT	nOpen;										// open shells to be added
												// by addExt(MOListIterator)

const HMatElements<RepMatType, VectorType>	*repMatP5;		// actually used P=5 RepMat-object

												
RepMatType	*pP5;								// pointer to actually used
												// P=5 matrix elements

RepMatType	*pP5Last;							// (last step in rs-tupel)
};


template <class RepMatType, class VectorType>
inline
RepMatType	RepDiag<RepMatType, VectorType>::getDiag(INT i) const
{	return	diag[i];	}




template <class RepMatType, class VectorType>
inline
void	RepDiag<RepMatType, VectorType>::addExtFirst(MOListIterator moiter)
{
RepMatType	diagI = 0;

	if ( moiter.j>0 && moiter.mode!=MOListIterator::noIndex ) 
		diagI = repMatP5->getCase14DiagAppendOpen(
			restConf,
			moiter.j
			);

	for ( INT i=0 ; i<dim ; i++ )
		DiagExtPartLast[i] = DiagExtPart[i] + diagI;

	if ( moiter.j>0 && moiter.mode!=MOListIterator::noIndex ) 
		repMatP5->getCase14UDiagRest(
			DiagExtPartLast,
			pP5,
			restConf,
			moiter.j
			);
	restConfLast = restConf;
//	cout << "restConf=" << restConf;
	if ( moiter.j>0 && moiter.mode!=MOListIterator::noIndex ) 
		restConfLast.insertOpenMO(moiter.j);
//	cout << "  restConfLast=" << restConfLast;
}

template <class RepMatType, class VectorType>
inline
void	RepDiag<RepMatType, VectorType>::addExtLast(MOListIterator moiter)
{
RepMatType	diagI = 0;

	switch ( nOpen ) {
	case 0:
		if ( moiter.i>0 && moiter.mode!=MOListIterator::noIndex ) 
		{
			diagI = repMatP5->getCase14DiagAppendClosed(
				restConf,
				moiter.i
				);
		}
		for ( INT i=0 ; i<dim ; i++ )
			diag[i] = DiagExtPart[i] + diagI;
		break;
		
	case 1:	
		if ( moiter.i>0 && moiter.mode!=MOListIterator::noIndex ) 
			diagI = repMatP5->getCase14DiagAppendOpen(
				restConf,
				moiter.i
				);

		for ( INT i=0 ; i<dim ; i++ )
			diag[i] = DiagExtPart[i] + diagI;

		if (  moiter.mode!=MOListIterator::noIndex ) 
			repMatP5->getCase14UDiagRest(
				diag,
				pP5,
				restConf,
				moiter.i
				);
		break;
		
	case 2:	
		if ( moiter.i>0 && moiter.mode!=MOListIterator::noIndex ) 
			diagI = repMatP5->getCase14DiagAppendOpen(
				restConfLast,
				moiter.i
				);

		for ( INT i=0 ; i<dim ; i++ )
			diag[i] = DiagExtPartLast[i] + diagI;

		if (  moiter.mode!=MOListIterator::noIndex ) 
			repMatP5->getCase14UDiagRest(
				diag,
				pP5Last,
				restConfLast,
				moiter.i
				);
		break;
	}
}


#endif

