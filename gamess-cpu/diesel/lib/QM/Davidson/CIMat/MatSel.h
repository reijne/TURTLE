//***********************************************************************
//
//	Name:			MatSel.h
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			01. Aug 1998
//
//***********************************************************************

#ifndef __MatSel_h
#define __MatSel_h

#include <iostream>

#include "../../../../config.h"

/* FD class ostream; */

class MatSel {
public:
	MatSel(INT dim);
	~MatSel();
	
	INT	getDim() const;

	INT & operator [] (INT);
	INT  operator [] (INT) const;

	INT		del(INT where);
	void	ins(INT where, INT what);
	void	neutral();
	
	
        friend std::ostream & operator << (std::ostream &, const MatSel &);
	
	

protected:
INT	dim;
INT	*p;
};


inline
INT	MatSel::getDim() const
{	return	dim;	}


inline
INT &	MatSel::operator [] (INT i)
{	return p[i];	}

inline
INT	MatSel::operator [] (INT i) const
{	return p[i];	}

#endif
