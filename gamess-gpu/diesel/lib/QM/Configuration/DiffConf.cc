//***********************************************************************
//
//	Name:	DiffConf.cc
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

#include "DiffConf.h"
#include "Configuration.h"

#include <stdlib.h>




/*
template <class TMOType>
DiffConf<TMOType>::DiffConf(Configuration<TMOType> & _from, Configuration<TMOType> & _to)
{
	from = _from;
	to = _to;
}
*/


/*
DiffConf<GeneralizedMO>::DiffConf(const DiffConf<MOType> & exc)
{
	same = exc.getSame();

	from = exc.getFrom();
	memcpy(posFrom, exc.getPosFromP(), from.getNumberOfOpenShells()*sizeof(INT));
	openShellsFrom = exc.getOpenShellsFrom();
		
	to = exc.getTo();
	memcpy(posTo, exc.getPosToP(), to.getNumberOfOpenShells()*sizeof(INT));
	openShellsTo = exc.getOpenShellsTo();
}
*/
	
/*
DiffConf<MOType>::operator DiffConf<GeneralizedMO> ()
{
DiffConf<GeneralizedMO>	diff;
	diff.same = same;

	diff.from = from;
	memcpy(diff.posFrom, getPosFromP(), to.getNumberOfOpenShells()*sizeof(INT));
	diff.openShellsFrom = getOpenShellsFrom();
		
	diff.to = to;
	memcpy(diff.posTo, getPosToP(), to.getNumberOfOpenShells()*sizeof(INT));
	diff.openShellsTo = getOpenShellsTo();

	return diff;
}
*/

template <>
DiffConf<MOType>::DiffConf(
	const DiffConf<GeneralizedMO> & exc,
	const MOListIterator & moiter)
{
	same = Configuration<MOType>(exc.getSame(), moiter);

	from = Configuration<MOType>(exc.getFrom(), moiter);
	memcpy(posFrom, exc.getPosFromP(), from.getNumberOfOpenShells()*sizeof(INT));
	openShellsFrom = exc.getOpenShellsFrom();
		
	to = Configuration<MOType>(exc.getTo(), moiter);
	memcpy(posTo, exc.getPosToP(), to.getNumberOfOpenShells()*sizeof(INT));
	openShellsTo = exc.getOpenShellsTo();
}
	
//--------------------------------------------------------------------------



template <class TMOType>
void	DiffConf<TMOType>::calcDiffConf(const Configuration<TMOType> & a, const Configuration<TMOType> & b,
	INT append)
{
INT	i, j;

//	cout << "a, b=" << a << ", " << b << endl;
//	cout << "DiffConf this: " << *this << endl;

	if ( !append )
	{	same.clear();
		from.clear();
		to.clear();
		openShellsFrom = a.getNumberOfOpenShells();
		openShellsTo = b.getNumberOfOpenShells();

		i = j = 0;
		while ( 1 ) 
		{	if ( i==a.getNumberOfClosedShells() )
			{	for ( ; j<b.getNumberOfClosedShells() ; j++ )
					to.appendClosedShell(b.getClosedShell(j));
				break;
			}
			if ( j==b.getNumberOfClosedShells() )
			{	for ( ; i<a.getNumberOfClosedShells() ; i++ )
					from.appendClosedShell(a.getClosedShell(i));
				break;
			}
			if ( a.getClosedShell(i)==b.getClosedShell(j) )
			{	same.appendClosedShell(a.getClosedShell(i));
				i++;
				j++;
			}
			else
			if ( a.getClosedShell(i)<b.getClosedShell(j) )
			{	from.appendClosedShell(a.getClosedShell(i));
				i++;
			}
			else
			{	to.appendClosedShell(b.getClosedShell(j));
				j++;
			}
		}

		i = j = 0;
		while ( 1 ) 
		{	if ( i==a.getNumberOfOpenShells() )
			{	for ( ; j<b.getNumberOfOpenShells() ; j++ )
				{	posTo[to.getNumberOfOpenShells()] = j;
					to.appendOpenShell(b.getOpenShell(j));
				}
				break;
			}
			if ( j==b.getNumberOfOpenShells() )
			{	for ( ; i<a.getNumberOfOpenShells() ; i++ )
				{	posFrom[from.getNumberOfOpenShells()] = i;
					from.appendOpenShell(a.getOpenShell(i));
				}
				break;
			}
			if ( a.getOpenShell(i)==b.getOpenShell(j) )
			{	same.appendOpenShell(a.getOpenShell(i));
				i++;
				j++;
			}
			else
			if ( a.getOpenShell(i)<b.getOpenShell(j) )
			{	posFrom[from.getNumberOfOpenShells()] = i;
				from.appendOpenShell(a.getOpenShell(i));
				i++;
			}
			else
			{	posTo[to.getNumberOfOpenShells()] = j;
				to.appendOpenShell(b.getOpenShell(j));
				j++;
			}
		}

	}
	else
	{	openShellsFrom += a.getNumberOfOpenShells();
		openShellsTo += b.getNumberOfOpenShells();

		i = j = 0;
		while ( 1 ) 
		{	if ( i==a.getNumberOfClosedShells() )
			{	for ( ; j<b.getNumberOfClosedShells() ; j++ )
					to.insertClosedMO(b.getClosedShell(j));
				break;
			}
			if ( j==b.getNumberOfClosedShells() )
			{	for ( ; i<a.getNumberOfClosedShells() ; i++ )
					from.insertClosedMO(a.getClosedShell(i));
				break;
			}
			if ( a.getClosedShell(i)==b.getClosedShell(j) )
			{	same.insertClosedMO(a.getClosedShell(i));
				i++;
				j++;
			}
			else
			if ( a.getClosedShell(i)<b.getClosedShell(j) )
			{	from.insertClosedMO(a.getClosedShell(i));
				i++;
			}
			else
			{	to.insertClosedMO(a.getClosedShell(i));
				j++;
			}
		}

		i = j = 0;
		while ( 1 ) 
		{	if ( i==a.getNumberOfOpenShells() )
			{	for ( ; j<b.getNumberOfOpenShells() ; j++ )
				{	posTo[to.getNumberOfOpenShells()] = j;
					to.insertOpenMO(b.getOpenShell(j));
				}
				break;
			}
			if ( j==b.getNumberOfOpenShells() )
			{	for ( ; i<a.getNumberOfOpenShells() ; i++ )
				{	posFrom[from.getNumberOfOpenShells()] = i;
					from.insertOpenMO(a.getOpenShell(i));
				}
				break;
			}
			if ( a.getOpenShell(i)==b.getOpenShell(j) )
			{	same.insertOpenMO(a.getOpenShell(i));
				i++;
				j++;
			}
			else
			if ( a.getOpenShell(i)<b.getOpenShell(j) )
			{	posFrom[from.getNumberOfOpenShells()] = i;
				from.insertOpenMO(a.getOpenShell(i));
				i++;
			}
			else
			{	posTo[to.getNumberOfOpenShells()] = j;
				to.insertOpenMO(b.getOpenShell(j));
				j++;
			}
		}
	}

}

template <class TMOType>
void	DiffConf<TMOType>::addExternal(const Configuration<TMOType> & a,
	const Configuration<TMOType> & b)
{
//	cout << a << " | " << b << endl;

//	cout << "From1: " << openShellsFrom << ", To: " << openShellsTo << endl;
	
//	cout << "DiffConf0 this: " << *this << "          " << endl;
	calcDiffConf(a, b, 1);
//	cout << "From2: " << openShellsFrom << ", To: " << openShellsTo << endl;
//	cout << "DiffConf1 this: " << *this << "          " << endl;
	
	updatePos(from, posFrom);
	updatePos(to, posTo);

//	cout << "DiffConf2 this: " << *this << "          " << endl;
//	cout << "ok." << endl;
}



template <class TMOType>
INT	DiffConf<TMOType>::addExternalTo(INT nOpen, INT nClosed, const TMOType *p)
{
INT	order = 0;
	openShellsTo += nOpen;
	for ( INT i=0 ; i<nOpen ; i++ , p++ )
	{
		for ( INT j=0 ; j<from.getNumberOfClosedShells() ; j++ )
			if ( *p==from.getClosedShell(j) )
				order--;

		for ( INT j=0 ; j<from.getNumberOfOpenShells() ; j++ )
			if ( *p==from.getOpenShell(j) )
			{
				from.deleteSingleMO(*p);
				same.insertOpenMO(*p);
				goto breakOpen;
			}
		to.insertOpenMO(*p);
		order++;
			
breakOpen:
		;
	}
	
	for ( INT i=0 ; i<nClosed ; i++ , p++ )
	{
		for ( INT j=0 ; j<from.getNumberOfClosedShells() ; j++ )
			if ( *p==from.getClosedShell(j) )
			{
				from.deleteDoubleMO(*p);
				same.insertClosedMO(*p);
				goto breakClosed;
			}
		for ( INT j=0 ; j<from.getNumberOfOpenShells() ; j++ )
			if ( *p==from.getOpenShell(j) )
			{
				to.insertClosedMO(*p);
				order++;
				goto breakClosed;
			}
		to.insertClosedMO(*p);
		order += 2;
			
breakClosed:
		;
	}
	
	updatePos(from, posFrom);
	updatePos(to, posTo);
	
	return order;
}



template <class TMOType>
void	DiffConf<TMOType>::updatePos(Configuration<TMOType> & conf, INT *ipos)
{
INT	add = 0;
INT	j = 0;
//	cout << "same: " << same << endl;
//	cout << "conf: " << conf << endl;
//	cout << conf.getNumberOfOpenShells() << " " << same.getNumberOfOpenShells() << endl;
	for ( INT i=0 ; i<conf.getNumberOfOpenShells() ; i++  )
	{	//cout << "A: i=" << i << ", j=" << j << ", " <<
		//	conf.getOpenShell(i) << ", " << same.getOpenShell(j) << endl;
		while ( j<same.getNumberOfOpenShells() &&
			conf.getOpenShell(i)>same.getOpenShell(j) )
		{//	cout << "i=" << i << ", j=" << j << ", " <<
			//	conf.getOpenShell(i) << ", " << same.getOpenShell(j) << endl;
			j++;
		}
		ipos[i] = j + add++;
	}
}


/*
void	DiffConf<GeneralizedMO>::addExternal(const Configuration<MOType> & a,
	const Configuration<MOType> & b)
{
INT	fromOpen = from.getNumberOfOpenShells();
INT	toOpen = to.getNumberOfOpenShells();
INT	fromClosed = from.getNumberOfClosedShells();
INT	toClosed = to.getNumberOfClosedShells();

//	cout << a << " | " << b << endl;

//	cout << "From1: " << openShellsFrom << ", To: " << openShellsTo << endl;
	
//	Umwandlung !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Configuration<GeneralizedMO>	agen = a, bgen = b;
	agen.setIrRep_Space();
	bgen.setIrRep_Space();
	calcDiffConf(agen, bgen, 1);
	
//	cout << "From2: " << openShellsFrom << ", To: " << openShellsTo << endl;
//	cout << "DiffConf1 this: " << *this << "          " << endl;
	
	from.sortMerge(fromOpen, fromClosed, posFrom);
	to.sortMerge(toOpen, toClosed, posTo);

//	cout << "DiffConf2 this: " << *this << "          " << endl;
//	cout << "ok." << endl;
}
*/


//--------------------------------------------------------------------------

//*********************************************
// extend to posFrom and posTo
template <class TMOType>
INT	operator == (const DiffConf<TMOType> & exc1, const DiffConf<TMOType> & exc2)
{
	return	(exc1.getSame()==exc2.getSame()) &&
			(exc1.getFrom()==exc2.getFrom()) &&
			(exc1.getTo()==exc2.getTo());
}

template <class TMOType>
INT	operator != (DiffConf<TMOType> & exc1, DiffConf<TMOType> & exc2)
{
	return	(exc1.getSame()!=exc2.getSame()) ||
			(exc1.getFrom()!=exc2.getFrom()) ||
			(exc1.getTo()!=exc2.getTo());
}


//--------------------------------------------------------------------------

template <>
Configuration<MOType>	DiffConf<GeneralizedMO>::simplify()
{
Configuration<GeneralizedMO>	v(to & from);
	to -= v;
	from -= v;

//	cout << "vA=" << v << endl;
	for ( INT i=0 ; i<v.getNumberOfOpenShells() ; i++ )
	{
	GeneralizedMO mo = v.getOpenShell(i);
		if ( mo.getType()==GeneralizedMO::IrRep_Space &&	
			 mo.getSortFlag()==2 )
		{
			mo.setType(GeneralizedMO::Number);
			v.setOpenShell(i, mo);
		}
	}

	for ( INT i=0 ; i<v.getNumberOfClosedShells() ; i++ )
	{
	GeneralizedMO mo = v.getClosedShell(i);
		if ( mo.getType()==GeneralizedMO::IrRep_Space &&	
			 mo.getSortFlag()==2 )
		{
			mo.setType(GeneralizedMO::Number);
			v.setClosedShell(i, mo);
		}
	}
//	cout << "vB=" << v << endl;
	same += v;
Configuration<MOType>	vv(v.getNumberOfOpenShells(), v.getNumberOfClosedShells());
	for ( INT i=0 ; i<v.getNumberOfOpenShells() ; i++ )
		vv.setOpenShell(i, v.getOpenShell(i).getMONumber());
	for ( INT i=0 ; i<v.getNumberOfClosedShells() ; i++ )
		vv.setClosedShell(i, v.getClosedShell(i).getMONumber());
	return vv;
}


template <>
Configuration<MOType>	DiffConf<MOType>::simplify()
{
Configuration<MOType>	v(to & from);
	to -= v;
	from -= v;
	same += v;

	updatePos(from, posFrom);
	updatePos(to, posTo);
	return v;
}

//--------------------------------------------------------------------------

template <class TMOType>
ostream& operator<<(ostream& s, const DiffConf<TMOType> & exc)
{
	switch ( exc.getOutputMode() ) {
	case MathObject::TerminalASCII:
		s << "[" << exc.getSame() << endl;
		s << "(" << exc.getFrom() << ") --> " << "(" << exc.getTo() << ")";
		s << endl;
		s << " ( ";
		for ( INT i=0 ; i<exc.getFrom().getNumberOfOpenShells() ; i++ )
			s << " " << exc.getPosFrom(i);
		for ( INT i=0 ; i<exc.getFrom().getNumberOfClosedShells() ; i++ )
			s << "  ";
		s << ")     ( ";
		for ( INT i=0 ; i<exc.getTo().getNumberOfOpenShells() ; i++ )
			s << " " << exc.getPosTo(i);
		for ( INT i=0 ; i<exc.getTo().getNumberOfClosedShells() ; i++ )
			s << "  ";
		s << ")]";
		s << endl;
		break;
		
	case MathObject::TeX:
		s << "[" << exc.getSame() << endl;
		s << "(" << exc.getFrom() << ",";
		for ( INT i=0 ; i<exc.getFrom().getNumberOfOpenShells() ; i++ )
			s << " " << exc.getPosFrom(i);
		s << ") --> ";

		s << "(" << exc.getTo() << ",";
		for ( INT i=0 ; i<exc.getTo().getNumberOfOpenShells() ; i++ )
			s << " " << exc.getPosTo(i);
		s << ")]" ;
		break;
	}
	return s;
}



template class DiffConf<MOType>;
template class DiffConf<GeneralizedMO>;



//**********************************************************************
//	explicit instantiation for g++ 2.7.2 
//	("template <class TMOType> void calcInteraction(...)" does not work)
void	DiffConfInstanciateTemplates()
{
#define TEMPLATE_INSTANCIATE \
	{\
	DiffConf<T> t;\
		cout << t;\
	}

#define T MOType
TEMPLATE_INSTANCIATE
#undef T	
	
#define T GeneralizedMO
TEMPLATE_INSTANCIATE
#undef T	
	
#undef TEMPLATE_INSTANCIATE
}


