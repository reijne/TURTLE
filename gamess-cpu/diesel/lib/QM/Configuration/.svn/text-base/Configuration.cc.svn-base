//***********************************************************************
//
//	Name:	Configuration.cc
//
//	Description:	implements configuration handling
//					based on second quantization
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			21.06.1996
//
//
//
//
//
//***********************************************************************




#include "Configuration.h"
#include "Excitation.h"

#include "../IO/Array.h"

#include "../MO/MOMapping.h"

#include "../IO/IntegerStream.h"

#include "../MO/Iterators/MOListIterator.h"

using namespace std;

static int	cmp(const void *p1, const void *p2)
{
	if ( *((const MOType *) p1)<*((const MOType *) p2) )
		return -1;
	if ( *((const MOType *) p1)>*((const MOType *) p2) )
		return 1;
	return 0;
}



template <>
void	Configuration<MOType>::sort()
{
	qsort(pOpen, openShells, sizeof(MOType), cmp);
	qsort(pClosed, closedShells, sizeof(MOType), cmp);
}


template <>
void	Configuration<GeneralizedMO>::sort()
{
}


template <class TMOType>
Configuration<TMOType>::Configuration(istream &s)
{
/*
INT	pintArray[MAXELECTRONS];
Array<INT>	intArray(pintArray, MAXELECTRONS);
//	s >> intArray;	// skip last <lf> (!)
	s >> intArray;
//	cout << "intArray:" << endl;
//	cout << intArray << endl;
	if ( intArray.getNumber() )
	{
		openShells = intArray[0];
		closedShells = intArray.getNumber()-(openShells+1);
	}
	else
		openShells = closedShells = 0;
		
	for ( INT i=0 ; i<openShells ; i++ )
		pOpen[i] = moMapping.getContinuous(intArray[i+1]);
	for ( INT i=0 ; i<closedShells ; i++ )
		pClosed[i]= moMapping.getContinuous(intArray[i+1+openShells]);
*/

	s >> openShells;
IntegerStream	is(' ');
	s >> is;
Pix	pix = is.first();
	closedShells = is.length()-getNumberOfOpenShells();
	for ( INT i=0 ; i<getNumberOfOpenShells() ; i++ )
	{
		pOpen[i] = moMapping.getContinuous(is(pix));
		is.next(pix);
	}
	for ( INT i=0 ; i<closedShells ; i++ )
	{
		pClosed[i] = moMapping.getContinuous(is(pix));
		is.next(pix);
	}


	sort();
}
 

/*
template <class TMOType>
Configuration<TMOType>::Configuration(INT NumberOfOpenShells, INT NumberOfClosedShells, 
	TMOType mo ...)
{
	openShells = NumberOfOpenShells;
	closedShells = NumberOfClosedShells;	

va_list	ap;


	if ( openShells )
		pOpen[0] = mo;
	else
		pClosed[0] = mo;

	va_start(ap, mo);
	for ( INT i=1 ; i<openShells ; i++ )
		pOpen[i] = va_arg(ap, TMOType);
	for ( INT i=(openShells ? 0 : 1 ) ; i<closedShells ; i++ )
		pClosed[i] = va_arg(ap, TMOType);
	va_end(ap);


}
*/


template <>
Configuration<MOType>::Configuration(INT NumberOfOpenShells, INT NumberOfClosedShells,
	INT mo1)
{
	openShells = NumberOfOpenShells;
	closedShells = NumberOfClosedShells;

	if ( openShells )
		pOpen[0] = mo1;
	else
		pClosed[0] = mo1;
}


template <>
Configuration<GeneralizedMO>::Configuration(INT NumberOfOpenShells, INT NumberOfClosedShells,
        GeneralizedMO mo1)
{
	openShells = NumberOfOpenShells;
	closedShells = NumberOfClosedShells;

	if ( openShells )
		pOpen[0] = mo1;
	else
		pClosed[0] = mo1;
}


template <>
Configuration<MOType>::Configuration(INT NumberOfOpenShells, INT NumberOfClosedShells,
	INT mo1, INT mo2)
{
	openShells = NumberOfOpenShells;
	closedShells = NumberOfClosedShells;


	if ( openShells )
		pOpen[0] = mo1;
	else
		pClosed[0] = mo1;

	for ( INT i=1 ; i<openShells ; i++ )
		pOpen[i] = mo2;
	for ( INT i=(openShells ? 0 : 1 ) ; i<closedShells ; i++ )
		pClosed[i] = mo2;
}


template <>
Configuration<GeneralizedMO>::Configuration(INT NumberOfOpenShells, INT NumberOfClosedShells,
        GeneralizedMO mo1, GeneralizedMO mo2)
{
	openShells = NumberOfOpenShells;
	closedShells = NumberOfClosedShells;


	if ( openShells )
		pOpen[0] = mo1;
	else
		pClosed[0] = mo1;

	for ( INT i=1 ; i<openShells ; i++ )
		pOpen[i] = mo2;
	for ( INT i=(openShells ? 0 : 1 ) ; i<closedShells ; i++ )
		pClosed[i] = mo2;
}

/*
Configuration<MOType>::Configuration(INT NumberOfOpenShells, INT NumberOfClosedShells,
	INT mo, ...)
{
	openShells = NumberOfOpenShells;
	closedShells = NumberOfClosedShells;

va_list ap;

	if ( openShells )
		pOpen[0] = mo;
	else
		pClosed[0] = mo;

	va_start(ap, mo);
	for ( INT i=1 ; i<openShells ; i++ )
		pOpen[i] = va_arg(ap, INT);
	for ( INT i=(openShells ? 0 : 1 ) ; i<closedShells ; i++ )
		pClosed[i] = va_arg(ap, INT);
	va_end(ap);
}


Configuration<GeneralizedMO>::Configuration(INT NumberOfOpenShells, INT NumberOfClosedShells,
        GeneralizedMO mo, ...)
{
	openShells = NumberOfOpenShells;
	closedShells = NumberOfClosedShells;

va_list	ap;


	if ( openShells )
		pOpen[0] = mo;
	else
		pClosed[0] = mo;

	va_start(ap, mo);
	for ( INT i=1 ; i<openShells ; i++ )
		pOpen[i] = va_arg(ap, GeneralizedMO);
	for ( INT i=(openShells ? 0 : 1 ) ; i<closedShells ; i++ )
		pClosed[i] = va_arg(ap, GeneralizedMO);
	va_end(ap);


}

*/



/*
Configuration<GeneralizedMO>::Configuration(const Configuration<MOType> & conf)
{
	openShells = conf.getNumberOfOpenShells();
	closedShells = conf.getNumberOfClosedShells();
	for ( INT i=0 ; i<openShells ; i++ )
		pOpen[i] =  conf.getOpenShell(i);
	for ( INT i=0 ; i<closedShells ; i++ )
		pClosed[i] =  conf.getClosedShell(i);
}
*/


template <>
Configuration<MOType>::operator Configuration<GeneralizedMO> ()
{
Configuration<GeneralizedMO>	conf(openShells, closedShells);
	for ( INT i=0 ; i<openShells ; i++ )
		conf.setOpenShell(i, pOpen[i]);
	for ( INT i=0 ; i<closedShells ; i++ )
		conf.setClosedShell(i, pClosed[i]);
	return conf;
}


template <>
Configuration<MOType>::Configuration(
	const Configuration<GeneralizedMO> & conf,
		const MOListIterator &moiter)
{
	openShells = conf.getNumberOfOpenShells();
	closedShells = conf.getNumberOfClosedShells();

	for ( INT i=0 ; i<openShells ; i++ )
		if ( conf.getOpenShell(i).getType() == GeneralizedMO::Number ||
			moiter.mode == MOListIterator::noIndex )
			pOpen[i] =  conf.getOpenShell(i).getMONumber();
		else
			if ( conf.getOpenShell(i).getSignature()==0 )
				pOpen[i] = moiter.i;
			else
				pOpen[i] = moiter.j;

	for ( INT i=0 ; i<closedShells ; i++ )
		if ( conf.getClosedShell(i).getType() == GeneralizedMO::Number ||
			moiter.mode == MOListIterator::noIndex )
			pClosed[i] =  conf.getClosedShell(i).getMONumber();
		else
			if ( conf.getClosedShell(i).getSignature()==0 )
				pClosed[i] = moiter.i;
			else
				pClosed[i] = moiter.j;
}


template <>
IrRep	Configuration<MOType>::calcIrRep(const MOIrReps & moirreps) const
{
IrRep	prod = 0;
	for ( INT i=0 ; i<openShells ; i++ )
		prod = moirreps.getProd(prod, moirreps.getIrRep(pOpen[i]));
	return prod;
}


template <>
void Configuration<GeneralizedMO>::setIrRep_Space()
{
	for ( INT i=0 ; i<openShells ; i++ )
		pOpen[i].setIrRep_Space();
	for ( INT i=0 ; i<closedShells ; i++ )
		pClosed[i].setIrRep_Space();
}

//--------------------------------------------------------------------------

template <class TMOType>
void	Configuration<TMOType>::calcInteraction(
	const Configuration<TMOType> & a,
	const Configuration<TMOType> & b, 
	Configuration<TMOType> & da,
	Configuration<TMOType> & db,
	INT *iposa, INT *iposb)
{
INT	i, j;

//	cout << "a, b=" << a << ", " << b << endl;

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
			{	iposa[da.getNumberOfOpenShells()] = j;
				da.appendOpenShell(b.getOpenShell(j));
			}
			break;
		}
		if ( j==b.getNumberOfOpenShells() )
		{	for ( ; i<a.getNumberOfOpenShells() ; i++ )
			{	iposb[db.getNumberOfOpenShells()] = i;
				db.appendOpenShell(a.getOpenShell(i));
			}
			break;
		}
		if ( a.getOpenShell(i)==b.getOpenShell(j) )
		{	i++;
			j++;
		}
		else
		if ( a.getOpenShell(i)<b.getOpenShell(j) )
		{	iposb[db.getNumberOfOpenShells()] = i;
			db.appendOpenShell(a.getOpenShell(i));
			i++;
		}
		else
		{	iposa[da.getNumberOfOpenShells()] = j;
			da.appendOpenShell(b.getOpenShell(j));
			j++;
		}
	}



//	cout << da << " ! " << db << endl;
}


template <class TMOType>
INT	Configuration<TMOType>::calcExcitationOrder(
	const Configuration<TMOType> & a,
	const Configuration<TMOType> & b)
{
INT	iao, iac, ibo, ibc;
INT	fao, fac, fbo, fbc;

//	cout << "a, b=" << a << ", " << b << endl;


//
//	differences in numbers of electrons in configurations are
//	treated as an excitation!
//



//	to symmetrize behavior with respect to configuration with
//	different numbers of electrons:
INT	sign = a.getNumberOfElectrons()>=b.getNumberOfElectrons() ? 1 : -1;


	iao = iac = ibo = ibc = 0;
	fao = a.getNumberOfOpenShells()>0;
	fac = a.getNumberOfClosedShells()>0;
	fbo = b.getNumberOfOpenShells()>0;
	fbc = b.getNumberOfClosedShells()>0;

INT	n;
TMOType	mo;
INT	nn = 0;

	while ( fao || fac || fbo || fbc ) 
	{	mo = 30000;
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
		{	n += sign;
			if ( ++iao==a.getNumberOfOpenShells() )
				fao = 0;
		}
	
		if ( fac && mo==a.getClosedShell(iac) )
		{	n += 2*sign;
			if ( ++iac==a.getNumberOfClosedShells() )
				fac = 0;
		}
	
		if ( fbo && mo==b.getOpenShell(ibo) )
		{	n -= sign;
			if ( ++ibo==b.getNumberOfOpenShells() )
				fbo = 0;
		}
	
		if ( fbc && mo==b.getClosedShell(ibc) )
		{	n -= 2*sign;
			if ( ++ibc==b.getNumberOfClosedShells() )
				fbc = 0;
		}

		if ( n>0 )
			nn += n;

//	short cut disabled 
//	(needed by MatchingMOListIterator)
//		if ( nn>2 )
//			return 3;
	}
//	cout << "nn=" << nn << endl;
	return nn;
}


/*
#ifndef __LINUX
INT	calcExcitationOrder(
	const Configuration<MOType> & a,
	const Configuration<MOType> & b)
{
INT	iao, iac, ibo, ibc;
INT	fao, fac, fbo, fbc;

//	cout << "a, b=" << a << ", " << b << endl;


//
//	differences in numbers of electrons in configurations are
//	treated as an excitation!
//



//	to symmetrize behavior with respect to configuration with
//	different numbers of electrons:
INT	sign = a.getNumberOfElectrons()>=b.getNumberOfElectrons() ? 1 : -1;



	iao = iac = ibo = ibc = 0;
	fao = a.getNumberOfOpenShells()>0;
	fac = a.getNumberOfClosedShells()>0;
	fbo = b.getNumberOfOpenShells()>0;
	fbc = b.getNumberOfClosedShells()>0;

INT	n;
INT	mo;
INT	nn = 0;

	while ( fao || fac || fbo || fbc ) 
	{	mo = 100000;
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
		{	n += sign;
			if ( ++iao==a.getNumberOfOpenShells() )
				fao = 0;
		}
	
		if ( fac && mo==a.getClosedShell(iac) )
		{	n += 2*sign;
			if ( ++iac==a.getNumberOfClosedShells() )
				fac = 0;
		}
	
		if ( fbo && mo==b.getOpenShell(ibo) )
		{	n -= sign;
			if ( ++ibo==b.getNumberOfOpenShells() )
				fbo = 0;
		}
	
		if ( fbc && mo==b.getClosedShell(ibc) )
		{	n -= 2*sign;
			if ( ++ibc==b.getNumberOfClosedShells() )
				fbc = 0;
		}


		if ( n>0 )
			nn += n;

//	short cut disabled 
//	(needed by MatchingMOListIterator)
//		if ( nn>2 )
//			return 3;
	}
	return nn;
}

#endif
*/


template <class TMOType>
INT	Configuration<TMOType>::calcExcitationOrderFast(
	const Configuration<TMOType> & a,
	const Configuration<TMOType> & b,
	INT maxOrder)
{
INT	order = 0;

	{
	INT	ko = 0;
	INT	kc = 0;
		for ( INT i=0 ; i<a.getNumberOfClosedShells() ; i++ )
		{
//			cout << i << " c " << a.getClosedShell(i) << ", order=" << order << endl;
			if ( order>maxOrder )
				return maxOrder+1;
			while ( ko<b.getNumberOfOpenShells() && 
				b.getOpenShell(ko)<a.getClosedShell(i) )
				ko++;
			if ( ko<b.getNumberOfOpenShells() && b.getOpenShell(ko)==a.getClosedShell(i) )
			{
				order++;
				continue;
			}

			while ( kc<b.getNumberOfClosedShells() && 
				b.getClosedShell(kc)<a.getClosedShell(i) )
				kc++;
			if ( kc<b.getNumberOfClosedShells() && b.getClosedShell(kc)==a.getClosedShell(i) )
				continue;

			order += 2;
		}
	}
//	cout << "::::  " << order << endl;
	{
	INT	ko = 0;
	INT	kc = 0;
		for ( INT i=0 ; i<a.getNumberOfOpenShells() ; i++ )
		{
//			cout << i << " o " << a.getOpenShell(i) << ", order=" << order << endl;
			if ( order>maxOrder )
				return maxOrder+1;
			while ( ko<b.getNumberOfOpenShells() && 
				b.getOpenShell(ko)<a.getOpenShell(i) )
				ko++;
			if ( ko<b.getNumberOfOpenShells() && b.getOpenShell(ko)==a.getOpenShell(i) )
				continue;

			while ( kc<b.getNumberOfClosedShells() && 
				b.getClosedShell(kc)<a.getOpenShell(i) )
				kc++;
			if ( kc<b.getNumberOfClosedShells() && b.getClosedShell(kc)==a.getOpenShell(i) )
				continue;
			
			order++;
		}
	}
	return order;
}


//--------------------------------------------------------------------------

template <class TMOType>
void	Configuration<TMOType>::create(const MOListIterator & molist)
{
	if ( molist.mode==MOListIterator::noIndex )
		return;
		
	if ( molist.mode==MOListIterator::oneIndex )
	{
		if ( molist.i==0 )
			return;
		if ( molist.isSquare() )
			insertMO(closedShells, pClosed, molist.i);
		else
			insertMO(openShells, pOpen, molist.i);
	}
	else
	{
		if ( molist.i!=0 )
			insertMO(openShells, pOpen, molist.i);
		if ( molist.j!=0 )
			insertMO(openShells, pOpen, molist.j);
	}
}


//--------------------------------------------------------------------------

template <>
void	Configuration<MOType>::writeToStream(ostream & s) const
{
	((Configuration<MOType> *) this)->sort();
	switch ( getOutputMode() ) {
	case MathObject::TerminalASCII:
	{
		s.setf(ios::dec);
		{
		IntegerStream	is(' ');
			s << getNumberOfOpenShells();
			for ( INT i=0 ; i<getNumberOfOpenShells() ; i++ )
			{
			INT h = moMapping.getReal(getOpenShell(i));
				is.append(h);
			}
			s << is;
		}
		{
		IntegerStream	is(' ');
			for ( INT i=0 ; i<getNumberOfClosedShells() ; i++ )
			{
			INT	h = moMapping.getReal(getClosedShell(i));
				is.append(h);
			}
			s << is;
		}
	}
		break;
		
	case MathObject::TeX:
		break;
	}
}

template <>
void	Configuration<GeneralizedMO>::writeToStream(ostream & s) const
{
	((Configuration<GeneralizedMO> *) this)->sort();
	switch ( getOutputMode() ) {
	case MathObject::TerminalASCII:
		s.setf(ios::dec);
		s << getNumberOfOpenShells();
		for ( INT i=0 ; i<getNumberOfOpenShells() ; i++ )
			s << " " << moMapping.getReal(getOpenShell(i));
		for ( INT i=0 ; i<getNumberOfClosedShells() ; i++ )
			s << " " << moMapping.getReal(getClosedShell(i));
		break;
		
	case MathObject::TeX:
		break;
	}
}


template <class TMOType>
ostream& operator<<(ostream& s, const Configuration<TMOType> & conf)
{
	conf.writeToStream(s);
	return s;
}


template <>
istream& operator >> (istream & s, Configuration<MOType> &conf)
{
/*
INT	pintArray[MAXELECTRONS];
Array<INT>	intArray(pintArray, MAXELECTRONS);
	s >> intArray;
	conf.openShells = intArray[0];
	conf.closedShells = intArray.getNumber()-(conf.getNumberOfOpenShells()+1);
	for ( INT i=0 ; i<conf.getNumberOfOpenShells() ; i++ )
		conf.pOpen[i] = moMapping.getContinuous(intArray[i+1]);
	for ( INT i=0 ; i<conf.closedShells ; i++ )
		conf.pClosed[i]= moMapping.getContinuous(intArray[i+1+conf.getNumberOfOpenShells()]);
	conf.sort();
	return s;
*/

	s >> conf.openShells;
IntegerStream	is(' ');
	s >> is;
Pix	pix = is.first();
	conf.closedShells = is.length()-conf.getNumberOfOpenShells();
	for ( INT i=0 ; i<conf.getNumberOfOpenShells() ; i++ )
	{
		conf.pOpen[i] = moMapping.getContinuous(is(pix));
		is.next(pix);
	}
	for ( INT i=0 ; i<conf.closedShells ; i++ )
	{
		conf.pClosed[i] = moMapping.getContinuous(is(pix));
		is.next(pix);
	}
	conf.sort();
	return s;
}


//--------------------------------------------------------------------------


template <class TMOType>
void	Configuration<TMOType>::sortMerge(INT openExtStart, INT closedExtStart)
{
TMOType	buf[MAXELECTRONS];
INT	ext;

//	cout << *this << endl;

	ext = openShells-openExtStart;
//	cout << openShells << " " << openExtStart << endl;
	memcpy(buf, pOpen+openExtStart, ext*sizeof(TMOType));
	for ( INT i=openExtStart-1 ; ext ; )
	{	if ( pOpen[i]>buf[ext-1] && i>=0 )
			pOpen[i+ext] = pOpen[i--];
		else
			pOpen[i+ext] = buf[ext-- - 1];
	}
	
	ext = closedShells-closedExtStart;
//	cout << closedShells << " " << closedExtStart << endl;
	memcpy(buf, pClosed+closedExtStart, ext*sizeof(TMOType));
	for ( INT i=closedExtStart-1 ; ext ; )
	{	if ( pClosed[i]>buf[ext-1] && i>=0 )
			pClosed[i+ext] = pClosed[i--];
		else
			pClosed[i+ext] = buf[ext-- - 1];
	}
}




template <class TMOType>
void	Configuration<TMOType>::sortMerge(INT openExtStart, INT closedExtStart,
	INT *ipos)
{
TMOType	buf[MAXELECTRONS];
INT	ext;

	cout << openExtStart << endl;
	if ( openExtStart>0 )
	{	ext = openShells-openExtStart;
		cout << "ext= " << ext << endl;
		memcpy(buf, pOpen+openExtStart, ext*sizeof(TMOType));
		for ( INT i=openExtStart-1 ; ext ; )
		{	if ( pOpen[i]>buf[ext-1] && i>=0 )
			{	pOpen[i+ext] = pOpen[i];
				ipos[i+ext] = ipos[i] + ext;
				cout << "ipos[i]" << ipos[i] << 
					", ipos[i+ext]"  << ipos[i+ext] << endl;
				i--;
			}
			else
			{	pOpen[i+ext] = buf[ext - 1];
				ipos[i+ext] = ((i>=0) ? ipos[i] : -1) + ext;
				ext--;
			}
		}
	}
	
	if ( closedExtStart>0 )
	{	ext = closedShells-closedExtStart;
		memcpy(buf, pClosed+closedExtStart, ext*sizeof(TMOType));
		for ( INT i=closedExtStart-1 ; ext ; )
		{	if ( pClosed[i]>buf[ext-1] && i>=0 )
				pClosed[i+ext] = pClosed[i--];
			else
				pClosed[i+ext] = buf[ext-- - 1];
		}
	}
}


//--------------------------------------------------------------------------

template <class TMOType>
INT	Configuration<TMOType>::insertMO(INT & shells, TMOType *p, TMOType mo)
{
INT	i, j;
	for ( i=0 ; i<shells ; i++ )
		if ( mo<=p[i] )
			break;

	if ( i<shells && mo==p[i] )
		return 1;

	j = i;
	for ( i=shells ; i>j ; i-- )
		p[i] = p[i-1];
	p[j] = mo;
	shells++;
	return 0;
}

template <class TMOType>
void	Configuration<TMOType>::insertOpenMO(TMOType mo)
{
INT	i, j;
	for ( i=0 ; i<openShells ; i++ )
		if ( mo<=pOpen[i] )
			break;

	if ( i<openShells && mo==pOpen[i] )
		return;

	j = i;
	for ( i=openShells ; i>j ; i-- )
		pOpen[i] = pOpen[i-1];
	pOpen[j] = mo;
	openShells++;
}

template <class TMOType>
void	Configuration<TMOType>::insertClosedMO(TMOType mo)
{
INT	i, j;
	for ( i=0 ; i<closedShells ; i++ )
		if ( mo<=pClosed[i] )
			break;

	if ( i<closedShells && mo==pClosed[i] )
		return;

	j = i;
	for ( i=closedShells ; i>j ; i-- )
		pClosed[i] = pClosed[i-1];
	pClosed[j] = mo;
	closedShells++;
}


template <class TMOType>
INT	Configuration<TMOType>::deleteMO(INT & shells, TMOType *p, TMOType mo)
{
INT	i;
	for ( i=0 ; i<shells ; i++ )
		if ( p[i]==mo )
			break;
	if ( i==shells )
		return 1;
	shells--;

	for ( ; i<shells ; i++ )
		p[i] = p[i+1];
	return 0;

}


template <class TMOType>
INT	Configuration<TMOType>::deleteSingleMO(TMOType mo)
{	return deleteMO(openShells, pOpen, mo);	}

template <class TMOType>
INT	Configuration<TMOType>::deleteDoubleMO(TMOType mo)
{	return deleteMO(closedShells, pClosed, mo);	}


template <class TMOType>
INT	Configuration<TMOType>::annihilate(TMOType mo)
{
	if ( deleteMO(openShells, pOpen, mo) )
	{	if ( deleteMO(closedShells, pClosed, mo) )
			return 1;
		return insertMO(openShells, pOpen, mo);
	}
	else
		return 0;
}



template <class TMOType>
INT	Configuration<TMOType>::create(TMOType mo)
{
	for ( INT i=0 ; i<closedShells ; i++ )
		if ( pClosed[i]==mo )
			return 1;

	if ( insertMO(openShells, pOpen, mo) )
	{	if ( insertMO(closedShells, pClosed, mo) )
			return 1;
		return deleteMO(openShells, pOpen, mo);
	}
	else
		return 0;
}


template <class TMOType>
INT	Configuration<TMOType>::excite(TMOType from, TMOType to)
{
	if ( from==to )
		return 0;

	if ( annihilate(from) )
		return 1;
	return create(to);
}




template <class TMOType>
INT	Configuration<TMOType>::excite(const Excitation & excitation)
{
	for ( INT i=0 ; i<excitation.getOrder() ; ++i )
		if ( excite(excitation.getFrom(i), excitation.getTo(i)) )
			return 1;
	return 0;
}


template <class TMOType>
INT	Configuration<TMOType>::check() const
{
	for ( INT i=1 ; i<openShells ; i++ )
		if ( pOpen[i-1]>=pOpen[i] )
			return 1;
	for ( INT i=1 ; i<closedShells ; i++ )
		if ( pClosed[i-1]>=pClosed[i] )
			return 1;
			
	for ( INT i=0 ; i<openShells ; i++ )
		for ( INT j=0 ; j<closedShells ; j++ )
			if ( pOpen[i] == pClosed[j] )
				return 1;
	return 0;
}


template <>
INT	Configuration<GeneralizedMO>::getNthExtPosOpen(INT n) const
{
	for ( INT i=0 ; i<openShells ; i++ )
		if ( pOpen[i].getType()==GeneralizedMO::IrRep_Space &&
                    pOpen[i].getSignature()==n )
			return i;
	return -1;
}


//--------------------------------------------------------------------------

template <class TMOType>
Configuration<TMOType>	Configuration<TMOType>::operator - 
	(const Configuration<TMOType> & conf2)
{
Configuration<TMOType>	conf;


	for ( INT i=0, j=0 ; i<getNumberOfOpenShells() ; i++ )
	{	while ( j<conf2.getNumberOfOpenShells() && 
			conf2.getOpenShell(j)<getOpenShell(i) )
				j++;
		if ( conf2.getOpenShell(j)!=getOpenShell(i) )
			conf.appendOpenShell(getOpenShell(i));
	}

	for ( INT i=0, j=0 ; i<getNumberOfClosedShells() ; i++ )
	{	while ( j<conf2.getNumberOfClosedShells() && 
			conf2.getClosedShell(j)<getClosedShell(i) )
				j++;
		if ( conf2.getClosedShell(j)!=getClosedShell(i) )
			conf.appendClosedShell(getClosedShell(i));
	}
	return conf;
}





//--------------------------------------------------------------------------


template <class TMOType>
INT	operator == (const Configuration<TMOType> & conf1, const Configuration<TMOType> & conf2)
{
	if ( conf1.getNumberOfOpenShells() != conf2.getNumberOfOpenShells() )
		return 0;
	if ( conf1.getNumberOfClosedShells() != conf2.getNumberOfClosedShells() )
		return 0;
	for ( INT i=0 ; i<conf1.getNumberOfOpenShells() ; i++ )
		if ( conf1.getOpenShell(i) != conf2.getOpenShell(i) )
			return 0;
	for ( INT i=0 ; i<conf1.getNumberOfClosedShells() ; i++ )
		if ( conf1.getClosedShell(i) != conf2.getClosedShell(i) )
			return 0;
	return 1;
}


INT operator == (const Configuration<MOType> & conf1, const Configuration<MOType> & conf2)
{
	if ( conf1.getNumberOfOpenShells() != conf2.getNumberOfOpenShells() )
		return 0;
	if ( conf1.getNumberOfClosedShells() != conf2.getNumberOfClosedShells() )
		return 0;
	for ( INT i=0 ; i<conf1.getNumberOfOpenShells() ; i++ )
		if ( conf1.getOpenShell(i) != conf2.getOpenShell(i) )
			return 0;
	for ( INT i=0 ; i<conf1.getNumberOfClosedShells() ; i++ )
		if ( conf1.getClosedShell(i) != conf2.getClosedShell(i) )
			return 0;
	return 1;
}




INT	operator <= (const Configuration<MOType> & conf1, const Configuration<MOType> & conf2)
{
	if ( conf1.getNumberOfOpenShells() > conf2.getNumberOfOpenShells() )
		return 0;
	if ( conf1.getNumberOfOpenShells() < conf2.getNumberOfOpenShells() )
		return 1;
	if ( conf1.getNumberOfClosedShells() > conf2.getNumberOfClosedShells() )
		return 0;
	if ( conf1.getNumberOfClosedShells() < conf2.getNumberOfClosedShells() )
		return 1;
	for ( INT i=0 ; i<conf1.getNumberOfOpenShells() ; i++ )
	{
		if ( conf1.getOpenShell(i) > conf2.getOpenShell(i) )
			return 0;
		if ( conf1.getOpenShell(i) < conf2.getOpenShell(i) )
			return 1;
	}
	for ( INT i=0 ; i<conf1.getNumberOfClosedShells() ; i++ )
	{
		if ( conf1.getClosedShell(i) > conf2.getClosedShell(i) )
			return 0;
		if ( conf1.getClosedShell(i) < conf2.getClosedShell(i) )
			return 1;
	}
	return 1;
}

INT	operator < (const Configuration<MOType> & conf1, const Configuration<MOType> & conf2)
{
	if ( conf1.getNumberOfOpenShells() > conf2.getNumberOfOpenShells() )
		return 0;
	if ( conf1.getNumberOfOpenShells() < conf2.getNumberOfOpenShells() )
		return 1;
	if ( conf1.getNumberOfClosedShells() > conf2.getNumberOfClosedShells() )
		return 0;
	if ( conf1.getNumberOfClosedShells() < conf2.getNumberOfClosedShells() )
		return 1;
	for ( INT i=0 ; i<conf1.getNumberOfOpenShells() ; i++ )
	{
		if ( conf1.getOpenShell(i) > conf2.getOpenShell(i) )
			return 0;
		if ( conf1.getOpenShell(i) < conf2.getOpenShell(i) )
			return 1;
	}
	for ( INT i=0 ; i<conf1.getNumberOfClosedShells() ; i++ )
	{
		if ( conf1.getClosedShell(i) > conf2.getClosedShell(i) )
			return 0;
		if ( conf1.getClosedShell(i) < conf2.getClosedShell(i) )
			return 1;
	}
	return 0;
}


//--------------------------------------------------------------------------

template <class TMOType>
Configuration<TMOType>	&Configuration<TMOType>::operator &=
	(const Configuration<TMOType> & conf2)
{
	*this = *this & conf2;
	return *this;
}

template <class TMOType>
Configuration<TMOType> Configuration<TMOType>::
	operator & (const Configuration & b)
{
INT	i, j;
Configuration<TMOType>	conf;
//	cout << "a, b=" << a << ", " << b << endl;

	i = j = 0;
	while ( 1 ) 
	{	if ( i==getNumberOfClosedShells() )
			break;

		if ( j==b.getNumberOfClosedShells() )
			break;

		if ( getClosedShell(i)==b.getClosedShell(j) )
		{
			conf.appendClosedShell(getClosedShell(i));
			i++;
			j++;
		}
		else
		if ( getClosedShell(i)<b.getClosedShell(j) )
			i++;
		else
			j++;
	}

	i = j = 0;
	while ( 1 ) 
	{
		if ( i==getNumberOfOpenShells() )
			break;

		if ( j==b.getNumberOfOpenShells() )
			break;

		if ( getOpenShell(i)==b.getOpenShell(j) )
		{
			conf.appendOpenShell(getOpenShell(i));
			i++;
			j++;
		}
		else
		if ( getOpenShell(i)<b.getOpenShell(j) )
			i++;
		else
			j++;
	}
	return conf;
}

template <class TMOType>
Configuration<TMOType> & Configuration<TMOType>::
	operator += (const Configuration &conf)
{
	for ( INT i=0 ; i<conf.getNumberOfOpenShells() ; i++ )
//		insertOpenMO(conf.getOpenShell(i));
		if ( create(conf.getOpenShell(i)) )
		{
			clear();
			return *this;
		}

	for ( INT i=0 ; i<conf.getNumberOfClosedShells() ; i++ )
//		insertClosedMO(conf.getClosedShell(i));
		if ( create(conf.getClosedShell(i)) )
		{
			clear();
			return *this;
		}
		else
		if ( create(conf.getClosedShell(i)) )
		{
			clear();
			return *this;
		}

	return *this;
}

template <class TMOType>
Configuration<TMOType> & Configuration<TMOType>::
	operator -= (const Configuration &conf)
{
//	cout << "-=" << conf << ",   " << *this;
	for ( INT i=0 ; i<conf.getNumberOfOpenShells() ; i++ )
//		deleteSingleMO(conf.getOpenShell(i));
		if ( annihilate(conf.getOpenShell(i)) )
		{
			clear();
			return *this;
		}

	for ( INT i=0 ; i<conf.getNumberOfClosedShells() ; i++ )
//		deleteDoubleMO(conf.getClosedShell(i));
		if ( annihilate(conf.getClosedShell(i)) )
		{
			clear();
			return *this;
		}
		else
		if ( annihilate(conf.getClosedShell(i)) )
		{
			clear();
			return *this;
		}
		
//	cout << ",    " << *this << endl;
	return *this;
}


template <>
void	Configuration<MOType>::split(const MRMOs *mrmos,
		Configuration<MOType> &internal, Configuration<MOType> &external) const
{
	internal = Configuration<MOType>();
	external = Configuration<MOType>();

	for ( INT i=0 ; i<openShells ; i++ )
		if ( mrmos->isInternal(pOpen[i]) )
			internal.appendOpenShell(pOpen[i]);
		else
			external.appendOpenShell(pOpen[i]);

	for ( INT i=0 ; i<closedShells ; i++ )
		if ( mrmos->isInternal(pClosed[i]) )
			internal.appendClosedShell(pClosed[i]);
		else
			external.appendClosedShell(pClosed[i]);
}

		
template <>
INT	Configuration<MOType>::getNumberOfExternals(const MRMOs *mrmos) const
{
INT	n = 0;
	for ( INT i=0 ; i<openShells ; i++ )
		if ( !mrmos->isInternal(pOpen[i]) )
			n++;
	for ( INT i=0 ; i<closedShells ; i++ )
		if ( !mrmos->isInternal(pClosed[i]) )
			n+=2;
	return n;
}


template class Configuration<MOType>;
template class Configuration<GeneralizedMO>;


template ostream& operator << (ostream& s, const
    Configuration<GeneralizedMO> & v);
template ostream& operator << (ostream& s, const
    Configuration<MOType> & v);

INT operator==(Configuration<MOType> const &,
	Configuration<MOType> const &);


/*
template <>
 INT	calcExcitationOrder(const Configuration<MOType> &,
	const Configuration<MOType> &);
*/

/*template <> INT
	calcExcitationOrder(const Configuration<GeneralizedMO> &,
	const Configuration<GeneralizedMO> &);
*/



//**********************************************************************
//	explicit instantiation for g++ 2.7.2 
//	("template <class TMOType> void calcInteraction(...)" does not work)
void	ConfigurationInstantiateTemplates()
{
	{
	Configuration<MOType> a, b, c, d;
	INT	ia[10], ib[10];

		a == b;
		Configuration<MOType>::calcInteraction(a, b, c, d, ia, ib);
		Configuration<MOType>::calcExcitationOrder(a, b);
		a = b - c;
		cout << a;
		cin >> a;
	}

	{
	Configuration<GeneralizedMO> a, b, c, d;
	INT	ia[10], ib[10];

		Configuration<GeneralizedMO>::calcInteraction(a, b, c, d, ia, ib);
                a = b - c;
		cout << a;
	}
}

#include "../../Container/AVLSet.cc"
template class AVLSet<Configuration<MOType> >;

#include "../../Container/Set.cc"
template class Set<Configuration<MOType> >;

