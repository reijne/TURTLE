//***********************************************************************
//
//	Name:			Permutable.h
//
//	Description:	abstract bases class for permutable objects
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			05.07.1996
//
//
//
//
//
//***********************************************************************


#ifndef __PERMUTABLE_H
#define __PERMUTABLE_H

#include "../../../../config.h"

class Permutator;

class Permutable {
public:
	virtual void permute(Permutator &) = 0;

	void	antisymmetrize();
	void	symmetrize();
	
private:
	virtual void DoSymmetrize(INT Anti) = 0;	
		
};


inline
void	Permutable::antisymmetrize()
{	DoSymmetrize(1);	}

inline
void	Permutable::symmetrize()
{	DoSymmetrize(0);	}



#endif
