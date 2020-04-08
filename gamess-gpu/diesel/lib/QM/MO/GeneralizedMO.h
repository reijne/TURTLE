//***********************************************************************
//
//	Name:			GeneralizedMO.h
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			14.08.1996
//
//
//
//
//
//***********************************************************************


#ifndef __GENERALIZEDMO_H
#define __GENERALIZEDMO_H

#include "../../../config.h"


#include "MOType.h"
#include "../Symmetry/IrRep.h"
#include "MRMOs.h"




//***********************************************************************
// The following is not implemented by the use of inheritance
// and virtual functions because a very special behaviour of
// the comparison operators is needed.



class GeneralizedMO {
public:
enum Type { Number, IrRep_Space };
	
//------------------------------------------------------------------

	GeneralizedMO();
	GeneralizedMO(MOType);
	GeneralizedMO(MOType, Type, INT sortFlag = 0);
	GeneralizedMO(MOType, INT, INT sortFlag = 0);
	GeneralizedMO(
		IrRep _irrep,
		INT _internal,
		MOType mo,
		INT signature,
		INT sortFlag);
//	virtual ~GeneralizedMO() {}


	
	// Copy-Konstruktor
//	Configuration(const Configuration<TMOType> &);
	
	// Zuweisungs-Operator
//	Configuration & operator = (const Configuration<TMOType> &);

	// type conversion
//	GeneralizedMO(const MOType &);

//------------------------------------------------------------------

	void	setMRMOs(MRMOs *);
	MRMOs *	getMRMOs() const;

//------------------------------------------------------------------
	
	void	set(MOType mo, Type type);
	void	set(MOType mo);
	
	void	setSignature(MOType signature);

//------------------------------------------------------------------

	void	setIrRep_Space();

	Type	getType() const;
	void	setType(Type);

	INT	getSortFlag() const;
	
//------------------------------------------------------------------

	MOType	getMONumber() const;
	INT	getSignature() const;

	IrRep	getIrRep() const;

	INT	isInternal() const;
	INT	isExternal() const;

//------------------------------------------------------------------

	INT	operator == (const GeneralizedMO &) const;
	INT	operator != (const GeneralizedMO &) const;
	INT	operator <(const GeneralizedMO &) const;
	INT	operator <= (const GeneralizedMO &) const;
	INT	operator > (const GeneralizedMO &) const;
	INT	operator >= (const GeneralizedMO &) const;

//------------------------------------------------------------------

	friend ostream& operator<<(ostream & s, const GeneralizedMO &);

//------------------------------------------------------------------

static MRMOs	*mrmos;		// information about symmetry and INT/ext status

protected:
Type	type;				// type flag: class contains MO either by
							// explicit number or by irrep and INT/ext status

MOType	mo;					// MO number
INT	signature;				// signature used to address integral indices
IrRep	irrep;				// irreducibel representation
INT	internal;				// flag if internal MO
INT	sortFlag;				// ordering in case of comparison Number<-->Irrep
							// 0: ==> external>internal
							// 1: ==> external<internal
							// 2: ==> signature contains MO number, 
							//        comparison based upon MO number
};





inline
GeneralizedMO::GeneralizedMO()
{	type = Number;	}


inline
void	GeneralizedMO::setMRMOs(MRMOs * _mrmos)
{	mrmos = _mrmos;	}


inline
MRMOs *	GeneralizedMO::getMRMOs() const
{	return mrmos;	}


/*inline
GeneralizedMO::GeneralizedMO(const MOType & _mo)
{
	type = Number;
	mo = _mo;
}
*/
//------------------------------------------------------------------

inline
GeneralizedMO::Type	GeneralizedMO::getType() const
{	return type;	}

inline
void	GeneralizedMO::setType(GeneralizedMO::Type _type)
{	type = _type;	}

inline
INT	GeneralizedMO::getSortFlag() const
{	return	sortFlag;	}

//------------------------------------------------------------------

inline
void	GeneralizedMO::set(MOType _mo)
{
	if ( type==Number )
		mo = _mo;
	else
		signature = _mo;
		
	irrep = mrmos->getIrRep(_mo);
	internal = mrmos->isInternal(_mo);
}

inline
void	GeneralizedMO::set(MOType mo, Type _type)
{
	type = _type;
	set(mo);
}

inline
void	GeneralizedMO::setIrRep_Space()
{
	if ( type==Number )
	{	irrep = mrmos->getIrRep(mo);
		internal = mrmos->isInternal(mo);
		type = IrRep_Space;
	}
}

inline
GeneralizedMO::GeneralizedMO(MOType i)
{	set(i, Number);	}

inline
GeneralizedMO::GeneralizedMO(MOType i, Type type, INT _sortFlag)
{
	set(i, type);
	sortFlag = _sortFlag;
}

inline
GeneralizedMO::GeneralizedMO(MOType i, INT j, INT _sortFlag)
{	
	if ( j<0 )
	{
		set(i, Number);
		signature = 0;
	}
	else
	{
		type = IrRep_Space;
		set(i);
		signature = j;
//	for back compatibility:
		mo = j;
	}
	sortFlag = _sortFlag;
}

inline
GeneralizedMO::GeneralizedMO(
	IrRep _irrep,
	INT _internal,
	MOType _mo,
	INT _signature,
	INT _sortFlag)
{
	type = IrRep_Space;
	mo = _mo;
	signature = _signature;
	sortFlag = _sortFlag;
	internal = _internal;
	irrep = _irrep;
}

inline
void	GeneralizedMO::setSignature(MOType _signature)
{	signature = _signature;	}

//------------------------------------------------------------------

inline
MOType	GeneralizedMO::getMONumber() const
{	return mo;	}

inline
INT	GeneralizedMO::getSignature() const
{	return signature;	}

inline
IrRep	GeneralizedMO::getIrRep() const
{	return irrep;	}

inline
INT	GeneralizedMO::isInternal() const
{	return internal;	}

inline
INT	GeneralizedMO::isExternal() const
{	return !internal;	}

//------------------------------------------------------------------

inline
INT	GeneralizedMO::operator == (const GeneralizedMO & mo2) const
{
	if ( (type!=IrRep_Space || sortFlag==2) &&
		 (mo2.getType()!=IrRep_Space || mo2.getSortFlag()==2) )
//	if ( type==mo2.getType() && type==Number)
		return mo==mo2.getMONumber();
	return 0;
}

inline
INT	GeneralizedMO::operator != (const GeneralizedMO & mo2) const
{
	if ( (type!=IrRep_Space || sortFlag==2) &&
		 (mo2.getType()!=IrRep_Space || mo2.getSortFlag()==2) )
//	if ( type==mo2.getType() && type==Number)
		return mo!=mo2.getMONumber();
	return 1;
}

inline
INT	GeneralizedMO::operator < (const GeneralizedMO & mo2) const
{
	if ( (type!=IrRep_Space || sortFlag==2) &&
		 (mo2.getType()!=IrRep_Space || mo2.getSortFlag()==2) )
//	if ( type==mo2.getType() && type==Number)
		return mo<mo2.getMONumber();

	if ( irrep<mo2.getIrRep() )
		return 1;
	if ( irrep>mo2.getIrRep() )
		return 0;

	if ( type==mo2.getType() )
		return mo<mo2.getMONumber();

	if ( internal>mo2.isInternal() )
		return 1;

	if ( internal<mo2.isInternal() )
		return 0;
		
/*	if ( sortFlag )
		return type==IrRep_Space;
	else
		return type!=IrRep_Space;
*/
//??????????????????????????????????
	return  ( type==IrRep_Space && sortFlag );
}

inline
INT	GeneralizedMO::operator > (const GeneralizedMO & mo2) const
{
	if ( (type!=IrRep_Space || sortFlag==2) &&
		 (mo2.getType()!=IrRep_Space || mo2.getSortFlag()==2) )
//	if ( type==mo2.getType() && type==Number)
		return mo>mo2.getMONumber();

	if ( irrep>mo2.getIrRep() )
		return 1;
	if ( irrep<mo2.getIrRep() )
		return 0;

	if ( type==mo2.getType() )
		return mo>mo2.getMONumber();


	if ( internal<mo2.isInternal() )
		return 1;

	if ( internal>mo2.isInternal() )
		return 0;
		
/*	if ( sortFlag )
		return type!=IrRep_Space;
	else
		return type==IrRep_Space;
*/
//??????????????????????????????????
	return  !( type==IrRep_Space && sortFlag );
}

inline
INT	GeneralizedMO::operator <= (const GeneralizedMO & mo2) const
{	return !((*this)>mo2);	}

inline
INT	GeneralizedMO::operator >= (const GeneralizedMO & mo2) const
{	return !((*this)<mo2);	}




#endif
