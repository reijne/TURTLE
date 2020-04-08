//***********************************************************************
//
//	Name:			Permutator.h
//
//	Description:	implements permutator objects
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

#ifndef __PERMUTATOR_H
#define __PERMUTATOR_H

#include "../../../../config.h"


#include <iostream>


class Permutator  {
public:
	Permutator(INT Ordnung = 0);
	~Permutator();
	
	
	Permutator & operator = (const Permutator & s);
	Permutator(const Permutator &);
	

	INT	getOrdnung() const;
	void	setOrdnung(INT n) ;
	void	setC(INT _c);
	
	INT	getC() const;
	INT	*getP() const;
	INT	getStart() const;
	void	setStart(INT n);


	Permutator &	shorten();
	Permutator &	expand(); 
	INT	getParity();


	friend INT operator == (const Permutator & a, const Permutator & b);
	friend INT operator != (const Permutator & a, const Permutator & b);
	
	Permutator & operator ++ (INT);
	Permutator & operator -- (INT);

	void	setIndex(INT n);
	INT	getIndex();

	LONG_INT	resetIndex();
	INT	nextIndex();	
	
	void	setTransposition(INT, INT);
	
	INT & operator[] (INT i) const;

	INT	check();

	Permutator Inv() const;

	friend Permutator operator * (const Permutator &, const Permutator &);
	friend Permutator operator / (const Permutator & p1, 
		const Permutator & p2);
	

	friend Permutator operator * (const INT &, const Permutator &);
	friend Permutator operator * (const Permutator & s, const INT & c);

	friend ostream& operator << (ostream& s, const Permutator & y);
	friend istream& operator >> (istream& s, Permutator & y);


	void	permute(void * pv, INT size);

private:
	void	Long2Short(INT *pi);
	void	Short2Long(INT *pi);
	void	inc(INT *pi);
	void	dec(INT *pi);

INT	Ordnung;	// Ordnung
INT	c;			// Vorfaktor
INT	start;		// Abbildung von (start  ...  start+Ordnung)
INT	*p;			//           auf (p[0]   ...  p[Ordnung-1])
INT	*sp;		// komprimierte Darstellung fuer Iterator "nextIndex()"
LONG_INT	nfak;	// restliche Iterationen
};

inline
INT	Permutator::getOrdnung() const
{	return Ordnung;	}	

inline
void	Permutator::setOrdnung(INT n) 
{	if ( Ordnung )
		delete p;
	Ordnung = n;
	if ( Ordnung )
		p = new INT[Ordnung];
}	

inline
void	Permutator::setC(INT _c)
{	c = _c;	}
	
inline
INT	Permutator::getC() const
{	return c;	}

inline
INT	*Permutator::getP() const
{	return p;	}


inline
INT	Permutator::getStart() const
{	return start;	}

inline
void	Permutator::setStart(INT n)
{	start = n;	}


inline
INT operator != (const Permutator & a, const Permutator & b)
{	return !(a==b);	}
	
	
	
inline
INT & Permutator::operator[] (INT i) const
{	return p[i];	}


inline
Permutator operator / (const Permutator & p1, const Permutator & p2)
{	return p1*p2.Inv();	}
	

inline
Permutator operator * (const Permutator & s, const INT & c)
{	return c*s;	}



#endif
