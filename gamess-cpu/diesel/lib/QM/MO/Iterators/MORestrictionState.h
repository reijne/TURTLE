//***********************************************************************
//
//	Name:			MORestrictionState.h
//
//	Description:	restrictions for MO occupation patterns
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			24.03.1998
//
//
//
//
//
//***********************************************************************


#ifndef __MORestrictionState_h
#define __MORestrictionState_h

#include "../../../../config.h"

#include "../../../Container/BitSet.h"

#include "../MOType.h"
#include "MORestriction.h"

template <class TMOType> class	Configuration;


class MORestrictionState {
public:
	MORestrictionState(const MORestriction &restriction,
		const Configuration<MOType> &conf);


	~MORestrictionState();


	//	attention: the corresponding configurations must be disjunct
	MORestrictionState &	operator |= (const MORestrictionState &);



	INT	check() const;

	INT isImpossible() const;

private:
	void Conf2BitSet(
		const Configuration<MOType> &conf,
		BitSet & singleMOs,
		BitSet & doubleMOs) const;


const MORestriction	*moRestriction;		//	pointer to applied MO restriction
INT					*occupation;		//	actual occupation pattern
INT					impossible;			//	flag if restriction can not be fullfilled
};


inline
INT	MORestrictionState::check() const
{
	if ( impossible )
		return 0;
	if ( !moRestriction )
		return 1;
	for ( INT i=0 ; i<moRestriction->nRestrictions ; i++ )
	{
//		cout << moRestriction->restrictions[i].op << "===" << occupation[i] << " " << moRestriction->restrictions[i].occupation << endl;
		switch ( moRestriction->restrictions[i].op ) {
		case 0:
			if ( !(occupation[i]<moRestriction->restrictions[i].occupation) )
				return 0;
			break;
		
		case 1:
			if ( !(occupation[i]<=moRestriction->restrictions[i].occupation) )
				return 0;
			break;
		
		case 2:
			if ( !(occupation[i]==moRestriction->restrictions[i].occupation) )
				return 0;
			break;
		
		case 3:
			if ( !(occupation[i]>moRestriction->restrictions[i].occupation) )
				return 0;
			break;
		
		case 4:
			if ( !(occupation[i]>=moRestriction->restrictions[i].occupation) )
				return 0;
			break;
		
		case 5:
			if ( !(occupation[i]!=moRestriction->restrictions[i].occupation) )
				return 0;
			break;
		}
	}
	return 1;
}


inline
INT MORestrictionState::isImpossible() const
{	return	impossible;	}


#endif
