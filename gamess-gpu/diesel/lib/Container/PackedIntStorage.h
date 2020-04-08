//***********************************************************************
//
//	Name:			PackedIntStorage.h
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			19.10.1998
//
//
//
//
//
//***********************************************************************


#ifndef __PackedIntStorage_h
#define __PackedIntStorage_h

#include "../../config.h"



class PackedIntStorage {
public:
	PackedIntStorage(INT n, INT ldBits);
	~PackedIntStorage();
	
	void	set(INT i, unsigned INT v);
	unsigned INT	get(INT i) const;
	
	
	
private:
unsigned INT	*p;
unsigned INT	mask;
INT	binDiv;
INT	binMod;
INT	ldBits;
};


#include <iostream>

inline
void	PackedIntStorage::set(INT i, unsigned INT v)
{
	v &= mask;

INT	j = (i & binMod) << ldBits;
INT	k = i >> binDiv;
	p[k] &= ~(mask << j);
	p[k] |= (v << j);
}


inline
unsigned INT	PackedIntStorage::get(INT i) const
{
	return (p[i >> binDiv] >> ((i & binMod) << ldBits)) & mask;
}
 




#endif
