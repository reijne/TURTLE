//***********************************************************************
//
//	Name:			MatSelSort.h
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




#ifndef __MatSelSort_h
#define __MatSelSort_h

#include "../../../../config.h" 

#include "MatSel.h"

//FD class ostream; 
using std::ostream;

class MatSelSort : public MatSel {
public:
	MatSelSort(INT dim);
	~MatSelSort();



	double	getValue(INT i) const;
	void	setValue(INT i, double v);
	void	sort();
	


	friend ostream & operator << (ostream &, const MatSelSort &);

struct	TRecord {
	INT	n;
	double	v;
	};

private:


double	*e;
};


inline
double	MatSelSort::getValue(INT i) const
{	return e[i];	};

inline
void	MatSelSort::setValue(INT i, double v)
{	e[i] = v;	};




#endif

