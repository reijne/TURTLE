//***********************************************************************
//
//	Name:			MORestrictionState.cc
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




#include "MORestrictionState.h"


#include "../../Configuration/Configuration.h"



MORestrictionState::MORestrictionState(const MORestriction &_moRestriction,
		const Configuration<MOType> &conf)
{
	moRestriction = &_moRestriction;
	if ( moRestriction )
	{
	BitSet	doubleMOs, singleMOs;

		occupation = new INT[moRestriction->nRestrictions];

		Conf2BitSet(conf, singleMOs, doubleMOs);

		impossible = 0;
		for ( INT i=0 ; i<moRestriction->nRestrictions ; i++ )
		{
			occupation[i] = 
				2*(doubleMOs & moRestriction->restrictions[i].mask).count() +
				(singleMOs & moRestriction->restrictions[i].mask).count();

	//		cout << i << " " << occupation[i] << endl;

			if ( moRestriction->restrictions[i].op==0 && 
				occupation[i]>=moRestriction->restrictions[i].occupation )
				impossible = 1;
			else
			if ( moRestriction->restrictions[i].op<=2 && 
				occupation[i]>moRestriction->restrictions[i].occupation )
				impossible = 1;
		}
	}
	else
	{
		impossible = 0;
		occupation = NULL;
	}
}




MORestrictionState::~MORestrictionState()
{
	if ( occupation )
		delete occupation;
}




void	MORestrictionState::Conf2BitSet(
	const Configuration<MOType> &conf,
	BitSet & singleMOs,
	BitSet & doubleMOs) const
{
	for ( INT i=0 ; i<conf.getNumberOfOpenShells() ; i++ )
		singleMOs.set(conf.getOpenShell(i)-1);

	for ( INT i=0 ; i<conf.getNumberOfClosedShells() ; i++ )
		doubleMOs.set(conf.getClosedShell(i)-1);
}




MORestrictionState &	MORestrictionState::operator |= (
	const MORestrictionState &m)
{
	impossible |= m.impossible;
	if ( moRestriction!=m.moRestriction || moRestriction==NULL || impossible )
		return *this;
	
	for ( INT i=0 ; i<moRestriction->nRestrictions ; i++ )
	{
		occupation[i] += m.occupation[i];
//		cout << occupation[i] << endl;
		
		if ( moRestriction->restrictions[i].op==0 && 
			occupation[i]>=moRestriction->restrictions[i].occupation )
			impossible = 1;
		else
		if ( moRestriction->restrictions[i].op<=2 && 
			occupation[i]>moRestriction->restrictions[i].occupation )
			impossible = 1;
	}

	return *this;
}




