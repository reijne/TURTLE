//***********************************************************************
//
//	Name:			MatchingMOListIterator.h
//
//	Description:	iterators for MO lists
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			26.03.1997
//
//
//
//
//
//***********************************************************************


#ifndef __MatchingMOListIterator_h
#define __MatchingMOListIterator_h

#include "../../../../config.h"


#include "MOList.h"
#include "../../Configuration/ConfigurationSAFNr.h"
#include "../../Configuration/DiffConf.h"
#include "MOListIterator.h"

#include "../../../Container/SLList.h"


struct	MatchingMOList {
MOListIterator *moiter;					// pointer to MO iterators
DiffConf<GeneralizedMO>	*dcGen;			// pointer to gen. DiffConfs
Configuration<MOType>	*extNotRunning;	// pointer to non running external MOs 
};




class MatchingMOListIterator {
public:
	MatchingMOListIterator(
		ConfigurationSAFNr<MOType> base,
		Configuration<MOType> toMatch,
		const MOList *molist,
		INT minOpenShells);

	~MatchingMOListIterator();

	void	next();
	INT	isEnd() const;
	
	const DiffConf<GeneralizedMO>	&	getGenDiffConf() const;
	const MOListIterator &	getMOListIterator() const;
	const Configuration<MOType> &	getExtNotRunning() const;

private:
DiffConf<GeneralizedMO>		dcGenBase;
const MOList *molist;
INT	excite;							// order of excitation between internal parts
SLList<MatchingMOList>	*matchList;	// list of to MO iterators
Pix	SLIter;							// list iterator
};


//	problem with GNU C/C++ Compiler v 2.7.2:
//	private variable "SLIter" is not used correctly
/*
inline
const DiffConf<GeneralizedMO> &	MatchingMOListIterator::getGenDiffConf() const
{	printf("iter = %lx\n", SLIter);

	printf("%lx\n", matchList(SLIter).dcGen);
	return *matchList(SLIter).dcGen;	}


inline
const MOListIterator &	MatchingMOListIterator::getMOListIterator() const
{	return *matchList(SLIter).moiter;	}

inline
const Configuration<MOType> &	MatchingMOListIterator::getExtNotRunning() const
{	return *matchList(SLIter).extNotRunning;	}

inline
void	MatchingMOListIterator::next()
{	matchList.next(SLIter);	}


inline
INT	MatchingMOListIterator::isEnd() const
{	printf("iter = %lx\n", SLIter);

	printf("%lx\n", matchList(SLIter).dcGen);
	return SLIter==0;	}
*/

#endif
