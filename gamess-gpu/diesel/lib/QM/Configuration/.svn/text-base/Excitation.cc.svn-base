//***********************************************************************
//
//	Name:	Excitation.cc
//
//	Description:	general representation of excitations
//
//	Author:	Michael Hanrath
//
//	Version:	0.0
//
//	Date:	11.06.1996
//
//
//
//
//
//***********************************************************************

#include "Excitation.h"
#include "Configuration.h"

#include <stdlib.h>

Excitation::Excitation(INT Order, MOType mo ...)
{
	order = Order;

va_list	ap;
	
	pFrom[0] = mo;
	va_start(ap, mo);
	for ( INT i=1 ; i<order ; i++ )
// `MOType' is promoted to `INT' when passed through '...'
// (so you should pass `INT' not `MOType' to `va_arg')
		pFrom[i] = va_arg(ap, INT);
	for ( INT i=0 ; i<order ; i++ )
// `MOType' is promoted to `INT' when passed through '...'
// (so you should pass `INT' not `MOType' to `va_arg')
		pTo[i] = va_arg(ap, INT);
	va_end(ap);


}


Excitation::Excitation(const Excitation & exc)
{
	order = exc.getOrder();
	memcpy(pFrom, exc.getFromP(), order*sizeof(MOType));
	memcpy(pTo, exc.getToP(), order*sizeof(MOType));
}
	
Excitation & Excitation::operator = (const Excitation & exc)
{
	order = exc.getOrder();
	memcpy(pFrom, exc.getFromP(), order*sizeof(MOType));
	memcpy(pTo, exc.getToP(), order*sizeof(MOType));
	return *this;
}

//--------------------------------------------------------------------------


/*
static int cmp(const void *p1, const void *p2)
{
	if ( *((MOType *) p1) > *((MOType *) p2) )
		return 1;
	if ( *((MOType *) p1) < *((MOType *) p2) )
		return -1;
	return 0;
}
*/

void	Excitation::calcExcitation(const Configuration<MOType> & a,
	const Configuration<MOType> & b)
{
/*
struct	TSortArray {
	MOType	mo;
	INT			m;
	};
	
INT	ng = 	a.getNumberOfOpenShells() + 
			a.getNumberOfClosedShells() + 
			b.getNumberOfOpenShells() +
			b.getNumberOfClosedShells();
			
TSortArray	*pl, *p = new TSortArray [ng];


	pl = p;
	for ( INT i=0 ; i<a.getNumberOfOpenShells() ; i++ )
	{	pl->mo = a.getOpenShell(i);
		pl++->m = 1;
	}

	for ( INT i=0 ; i<a.getNumberOfClosedShells() ; i++ )
	{	pl->mo = a.getClosedShell(i);
		pl++->m = 2;
	}

	for ( INT i=0 ; i<b.getNumberOfOpenShells() ; i++ )
	{	pl->mo = b.getOpenShell(i);
		pl++->m = -1;
	}

	for ( INT i=0 ; i<b.getNumberOfClosedShells() ; i++ )
	{	pl->mo = b.getClosedShell(i);
		pl++->m = -2;
	}
	
	qsort(p, ng, sizeof(TSortArray), cmp);
	
		
INT	mo, n;
INT	nf = 0;
INT	nt = 0;

	for ( INT i=0 ; i<ng ;  )
	{	mo = p[i].mo;
		n = 0;
		while ( i<ng && p[i].mo==mo )
			n += p[i++].m;
		
		switch ( n ) {
		case -2:
			pTo[nt++] = mo;
			pTo[nt++] = mo;
			break;

		case -1:
			pTo[nt++] = mo;
			break;

		case 1:
			pFrom[nf++] = mo;
			break;

		case 2:
			pFrom[nf++] = mo;
			pFrom[nf++] = mo;
			break;
		}
	}
	order = nf;
	
	delete p;

*/


INT	iao, iac, ibo, ibc;
INT	fao, fac, fbo, fbc;

//	cout << "a, b=" << a << ", " << b << endl;




	iao = iac = ibo = ibc = 0;
	fao = a.getNumberOfOpenShells()>0;
	fac = a.getNumberOfClosedShells()>0;
	fbo = b.getNumberOfOpenShells()>0;
	fbc = b.getNumberOfClosedShells()>0;

INT	n, mo;	
INT	nf = 0;
INT	nt = 0;

	while ( fao || fac || fbo || fbc ) 
	{	mo = 1000000;
		if ( fao && a.getOpenShell(iao)<mo )
			mo = a.getOpenShell(iao);

		if ( fac && a.getClosedShell(iac)<mo )
			mo = a.getClosedShell(iac);

		if ( fbo && b.getOpenShell(ibo)<mo )
			mo = b.getOpenShell(ibo);

		if ( fbc && b.getClosedShell(ibc)<mo )
			mo = b.getClosedShell(ibc);
		
		n = 0;

		if ( fao && mo==a.getOpenShell(iao) )
		{	n++;
			if ( ++iao==a.getNumberOfOpenShells() )
				fao = 0;
		}
	
		if ( fac && mo==a.getClosedShell(iac) )
		{	n+=2;
			if ( ++iac==a.getNumberOfClosedShells() )
				fac = 0;
		}
	
		if ( fbo && mo==b.getOpenShell(ibo) )
		{	n--;
			if ( ++ibo==b.getNumberOfOpenShells() )
				fbo = 0;
		}
	
		if ( fbc && mo==b.getClosedShell(ibc) )
		{	n-=2;
			if ( ++ibc==b.getNumberOfClosedShells() )
				fbc = 0;
		}
	
		switch ( n ) {
		case -2:
			pTo[nt++] = mo;
			pTo[nt++] = mo;
			break;

		case -1:
			pTo[nt++] = mo;
			break;

		case 1:
			pFrom[nf++] = mo;
			break;

		case 2:
			pFrom[nf++] = mo;
			pFrom[nf++] = mo;
			break;
		}
	}
	order = nf;		

	

/*
Configuration	da, db;
INT	i, j;
	i = j = 0;
	while ( 1 ) 
	{	if ( i==a.getNumberOfClosedShells() )
		{	for ( ; j<b.getNumberOfClosedShells() ; j++ )
				da.appendClosedShell(b.getClosedShell(j));
			break;
		}
		if ( j==b.getNumberOfClosedShells() )
		{	for ( ; i<a.getNumberOfClosedShells() ; i++ )
				db.appendClosedShell(a.getClosedShell(i));
			break;
		}
		if ( a.getClosedShell(i)==b.getClosedShell(j) )
		{	i++;
			j++;
		}
		else
		if ( a.getClosedShell(i)<b.getClosedShell(j) )
		{	db.appendClosedShell(a.getClosedShell(i));
			i++;
		}
		else
		{	da.appendClosedShell(b.getClosedShell(j));
			j++;
		}
	}

	i = j = 0;
	while ( 1 ) 
	{	if ( i==a.getNumberOfOpenShells() )
		{	for ( ; j<b.getNumberOfOpenShells() ; j++ )
				da.appendOpenShell(b.getOpenShell(j));
			break;
		}
		if ( j==b.getNumberOfOpenShells() )
		{	for ( ; i<a.getNumberOfOpenShells() ; i++ )
				db.appendOpenShell(a.getOpenShell(i));
			break;
		}
		if ( a.getOpenShell(i)==b.getOpenShell(j) )
		{	i++;
			j++;
		}
		else
		if ( a.getOpenShell(i)<b.getOpenShell(j) )
		{	db.appendOpenShell(a.getOpenShell(i));
			i++;
		}
		else
		{	da.appendOpenShell(b.getOpenShell(j));
			j++;
		}
	}
	
	.
	.
	.
*/

}

//--------------------------------------------------------------------------


INT	operator == (const Excitation & exc1, const Excitation & exc2)
{
	if ( exc1.getOrder() != exc2.getOrder() )
		return 0;
	for ( INT i=0 ; i<exc1.getOrder() ; i++ )
		if ( exc1.getFrom(i) != exc2.getFrom(i) )
			return 0;
	for ( INT i=0 ; i<exc1.getOrder() ; i++ )
		if ( exc1.getTo(i) != exc2.getTo(i) )
			return 0;
	return 1;

}

INT	operator != (Excitation & exc1, Excitation & exc2)
{	return !(exc1==exc2);	}


//--------------------------------------------------------------------------


ostream& operator<<(ostream& s, const Excitation & exc)
{
	switch ( exc.getOutputMode() ) {
	case MathObject::TerminalASCII:
		for ( INT i=0 ; i<exc.getOrder() ; i++ )
			s << " " << exc.getFrom(i);
		s << " --> ";
		for ( INT i=0 ; i<exc.getOrder() ; i++ )
			s << " " << exc.getTo(i);
		break;
		
	case MathObject::TeX:
		break;
	}
	return s;
}



