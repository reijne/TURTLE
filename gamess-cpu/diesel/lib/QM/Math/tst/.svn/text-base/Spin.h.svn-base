//***********************************************************************
//
//	Name:			Spin.h
//
//	Description:	implements spin algebra
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

#ifndef __SPIN_H
#define __SPIN_H

#include "../../../../config.h"

#include "Permutable.h"
#include "../../../Container/Complex.h"
#include <iostream>
#include <string>


typedef char SpinBaseType;


class Permutator;


class Spin : public Permutable {
public:
	Spin(INT nsum = 0);	
	virtual ~Spin();	
	
	Spin & operator = (const Spin & s);
	Spin(const Spin &);
	

	INT	getNumberOfTerms() const;
	INT	getNumberOfElectrons() const;
	INT	getNumberOfShells(INT Summand = 0) const;

	void	setNumberOfShells(INT shells) const;
	void	setNumberOfShells(INT Summand, INT shells) const;
	
	SpinBaseType	*getP(INT Summand = 0) const;
	

	void	setSpin(LONG_INT l);
	void	setSpin(INT Summand, LONG_INT l);
	void	setSpinBinary(INT n, SpinBaseType c);


	LONG_INT	getSpin();
	
	LONG_INT	getSpin(INT Summand);

	LONG_INT	getSpinCompact();	
	LONG_INT	getSpinCompact(INT Summand);

	void	setSpinShell(INT Summand, INT shell, SpinBaseType c) const;	
	void	setSpinShell(INT shell, SpinBaseType c) const;	

	SpinBaseType	getSpinShell(INT shell) const;	
	SpinBaseType	getSpinShell(INT Summand, INT shell) const;	
	
	void	setC(Complex _c);
	void	setC(INT Summand, Complex _c);	
	Complex	getC(INT i) const;	


	friend INT operator == (Spin, Spin);
	friend INT operator != (Spin & a, Spin & b);
	
	INT	isMultipleOf(Spin a, Complex & e);
	INT	isMultipleOf(Spin & a);	
	
	Spin & operator += (const Spin &);
	
	friend Spin operator + (const Spin &, const Spin &);
	friend Spin operator - (const Spin &, const Spin &);
	friend Spin operator * (const Spin &, const Spin &);

	friend Spin operator * (const Complex &, const Spin &);
	friend Spin operator * (const Spin & s, const Complex & c);
	friend Spin operator * (const double &, const Spin &);
	friend Spin operator * (const Spin & s, const double & c);



//	Spin-Operatoren
	Spin 	Sx();
	Spin 	Sy();
	Spin 	Sz();
	Spin 	SS();		// S^2
	Spin 	Sp();		// S+
	Spin 	Sm();		// S-

	Complex	Integrate();

	virtual void permute(Permutator &);


	Spin 	simplify();	// Ziehe Terme mit gleichem Spin zusammen
						// und ellimiere Null-Terme


	friend ostream& operator << (ostream& s, const Spin & y);
	friend istream& operator >> (istream& s, Spin & y);
                

private:
	void	init(INT);
	void	kill();

	Spin 	Sxy(INT);
	Spin 	Spm(INT);
	INT	compareSpin(INT n1, SpinBaseType *p1, INT n2, SpinBaseType *p2);

	virtual void DoSymmetrize(INT Anti);
	
INT	nsum;				// Anzahl Summanden
INT	electrons;			// Anzahl Teilchen
INT	*n;					// Anzahl Schalen
Complex	*c;				// Vorfaktoren
SpinBaseType	**p;	// Zeiger auf Spin-Array
};

Spin	operator - (Spin & s);
Spin alpha(INT n);
Spin beta(INT n);


inline
Spin::Spin(INT nsum = 0)
{	init(nsum);	}
	
inline
Spin::~Spin()
{	kill();	}
	
	
	
inline
INT	Spin::getNumberOfTerms() const
{	return nsum;	}	
	
		
inline
INT	Spin::getNumberOfShells(INT Summand = 0) const
{	return	n[Summand];	}
	
inline
void	Spin::setNumberOfShells(INT Summand, INT shells) const
{	if ( n[Summand] )
		delete p[Summand];
	p[Summand] = new SpinBaseType[shells];
	memset(p[Summand], 0, shells*sizeof(SpinBaseType));
	n[Summand] = shells;
}

inline
void	Spin::setNumberOfShells(INT shells) const
{	setNumberOfShells(0, shells);	}
	

inline
SpinBaseType	*Spin::getP(INT Summand = 0) const
{	return p[Summand];	}
	

inline
void	Spin::setSpin(LONG_INT l)
{	setSpin(0, l);	}
	

inline
LONG_INT	Spin::getSpin()
{	return getSpin(0);	}
	

inline
LONG_INT	Spin::getSpinCompact()
{	return getSpinCompact(0);	}
	
inline
void	Spin::setSpinShell(INT Summand, INT shell, SpinBaseType c) const
{	p[Summand][shell] = c;	}
	
inline
void	Spin::setSpinShell(INT shell, SpinBaseType c) const
{	p[0][shell] = c;	}
	

inline
SpinBaseType	Spin::getSpinShell(INT shell) const
{	return p[0][shell];	}
	
inline
SpinBaseType	Spin::getSpinShell(INT Summand, INT shell) const
{	return p[Summand][shell];	}
	
	
inline
void	Spin::setC(Complex _c)
{	c[0] = _c;	}

inline
void	Spin::setC(INT Summand, Complex _c)
{	c[Summand] = _c;	}
	
inline
Complex	Spin::getC(INT i) const
{	return c[i];	}
	


inline
INT operator != (Spin & a, Spin & b)
{	return !(a==b);	}

inline
INT	Spin::isMultipleOf(Spin & a)
{	Complex c;
	return isMultipleOf(a, c);
}
	
	
	
inline
Spin	operator * (const Spin & s, const Complex & c)
{	return c*s;	}

inline
Spin	operator * (const Spin & s, const double & c)
{	return c*s;	}




inline
Spin 	Spin::Sx()
{	return Sxy(0);	}

inline
Spin 	Spin::Sy()
{	return Sxy(1);	}
	
inline
Spin 	Spin::Sp()		// S+
{	return Spm(0);	}

inline
Spin 	Spin::Sm()		// S-
{	return Spm(1);	}




#endif
