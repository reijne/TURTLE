//***********************************************************************
//
//	Name:			Permutator.cc
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

#include "Permutator.h"
#include <stdlib.h>
#include <string>




Permutator::Permutator(INT _Ordnung)
{
	Ordnung = _Ordnung;
	if ( Ordnung )
	{	start = 0;
		p = new INT[Ordnung];
	}
	c = 1;
	sp = NULL;
}

Permutator::~Permutator()
{
	if ( Ordnung )
		delete p;
	if ( sp )
		delete sp;	
	Ordnung = 0;
}



Permutator	&Permutator::operator = (const Permutator & s)
{
	if ( Ordnung != s.getOrdnung() )
	{	if ( Ordnung ) 
			delete p;
		Ordnung = s.getOrdnung();
		if ( Ordnung )
			p = new INT[Ordnung];
	}

	start = s.getStart();
	memcpy(p, s.getP(), Ordnung * sizeof(INT));
	c = s.getC();
	sp = NULL;
	return *this;
}


Permutator::Permutator(const Permutator & s)
{
	Ordnung = s.getOrdnung();
	c = s.getC();
	start = s.getStart();
	sp = NULL;
	p = new INT[Ordnung];
	memcpy(p, s.getP(), Ordnung * sizeof(INT));
}


INT	operator == (const Permutator & a, const Permutator & b)
{
	if ( a.getOrdnung()!=b.getOrdnung() )
		return 0;
	if ( a.getC()!=b.getC() )
		return 0;
	return 1;	
}


Permutator operator * (const Permutator &b, const Permutator &a)
{
INT	uu, oo;
	if ( a.getStart()<b.getStart() )
		uu = a.getStart();
	else
		uu = b.getStart();

	if ( a.getStart()+a.getOrdnung()>b.getStart()+b.getOrdnung() )
		oo = a.getStart()+a.getOrdnung();
	else
		oo = b.getStart()+b.getOrdnung();

Permutator	p(oo-uu);
	if ( a.getStart()<b.getStart() )
		p.setStart(a.getStart());
	else
		p.setStart(b.getStart());

INT	pi[oo-uu];

	for ( INT i=uu ; i<a.getStart() ; i++ )
		p[i-uu] = i;
	memcpy(p.getP()+a.getStart()-uu, 
		a.getP(), a.getOrdnung()*sizeof(INT));
	for ( INT i=a.getStart()+a.getOrdnung() ; i<oo ; i++ )
		p[i-uu] = i;
	memcpy(pi, p.getP(), (oo-uu)*sizeof(INT));

	for ( INT i=0 ; i<b.getOrdnung() ; i++ )
		p[b[i]-uu] = pi[i+b.getStart()-uu]; 
	
	
	p.setC(a.getC()*b.getC());
	return p;
}


Permutator operator * (const INT &a, const Permutator &b)
{
Permutator	p(b.getOrdnung());
	p = b;
	p.setC(a*b.getC());
	return p;
}



void	Permutator::setTransposition(INT a1, INT a2)
{
	for ( INT i=0 ; i<Ordnung ; i++ )
		p[i] = i;
	p[a1] = a2;
	p[a2] = a1;
}


INT	Permutator::getParity()
{
INT	pi[Ordnung];
	for ( INT i=0 ; i<Ordnung ; i++ )
		pi[i] = p[i]-start;
INT	par = 1;
	for ( INT i=0 ; i<Ordnung-1 ; i++ )
	{	if ( pi[i]==i )
			continue;
		for ( INT j=i+1 ; j<Ordnung ; j++ )
			if ( pi[j]==i )
			{	for ( INT k=j ; k>i ; k-- )
					pi[k] = pi[k-1];			
				if ( (j-i) & 1 )
					par = -par;
				break;
			}	
	}
	return par;
}


Permutator Permutator::Inv() const
{
Permutator	s(Ordnung);

	for ( INT i=0 ; i<Ordnung ; i++ )
		s[p[i]-start] = i+start; 
	s.setC(c);
	s.setStart(start);
	return s;
}



Permutator & Permutator::shorten()
{
INT	u, o;
	for ( u=0 ; u<Ordnung ; u++ )
		if ( p[u]!=u+start )
			break;

	for ( o=Ordnung-1 ; o>u ; o-- )
		if ( p[o]!=o+start )
			break;
			
INT	*pi = NULL;
	if ( o-u+1>0 )
	{	pi = new INT[o-u+1];
		memcpy(pi, p+u, (o-u+1)*sizeof(INT));
	}	
	delete p;
	p = pi;
	Ordnung = o-u+1;
	start = u;	 		
	return *this;
}


Permutator & Permutator::expand()
{
	if ( !start )
		return *this;
		
INT	*pi = new INT[start+Ordnung];

	for ( INT i=0 ; i<start ; i++ )
		pi[i] = i;
	memcpy(pi+start, p, Ordnung*sizeof(INT));
	delete p;
	p = pi;
	Ordnung += start;

	return *this;
}

static int Permutator_compare_ints(const void *_a, const void *_b)
{
const INT *a = (const INT *) _a;
const INT *b = (const INT *) _b;

	return (*a > *b);
}


INT	Permutator::check()
{
INT	*pi = new INT[Ordnung];

	memcpy(pi, p, Ordnung*sizeof(INT));
	qsort(pi, Ordnung, sizeof(INT), Permutator_compare_ints);
INT	flag = 1;
	for ( INT i=0 ; i<Ordnung ; i++ )
		if ( pi[i]!=i+start )
		{	flag = 0;
			break;
		}
	delete pi;
	return flag;
}


void	Permutator::permute(void * pv, INT size)
{
char	*pc = new char[Ordnung*size];
	memcpy(pc, pv, Ordnung*size);
	for ( INT i=0 ; i<Ordnung ; i++ )
		memcpy(((char *) pv)+p[i]*size, pc+(i+start)*size, size);
	delete pc;
}	



void	Permutator::Long2Short(INT *pi)
{
INT	pos[Ordnung];
	for ( INT i=0 ; i<Ordnung ; i++ )
		pos[i] = i;

INT	k;
	for	( INT i=0 ; i<Ordnung-1 ; i++ )
	{	k=pi[i];
		if ( i )
			pi[i] = pos[k];
		for ( INT j=k+1 ; j<Ordnung ; j++ )
			pos[j]--;
	}
}

void	Permutator::Short2Long(INT *pi)
{
INT	pos[Ordnung];
	for ( INT i=0 ; i<Ordnung ; i++ )
		pos[i] = i;

INT	k;
	for	( INT i=0 ; i<Ordnung-1 ; i++ )
	{	k=pi[i];
		if ( i )
			pi[i] = pos[k];
		memmove(pos+k, pos+k+1, (Ordnung-k-i-1)*sizeof(INT));	
	}
	pi[Ordnung-1] = pos[0];
}


void	Permutator::inc(INT *pi)
{
	for ( INT i=Ordnung-1-1 , j=2 ; i>=0 ; i-- , j++ )
	{	pi[i]++;
		if ( pi[i]>=j )
			pi[i] = 0;
		else
			break;		 
	}
}


void	Permutator::dec(INT *pi)
{
	for ( INT i=Ordnung-1-1 , j=1 ; i>=0 ; i-- , j++ )
	{	pi[i]--;
		if ( pi[i]<0 )
			pi[i] = j;
		else
			break;		 
	}
}


Permutator & Permutator::operator ++ (INT)
{
	expand();
	Long2Short(p);
	inc(p);
	Short2Long(p);
	return *this;
}


Permutator & Permutator::operator -- (INT)
{
	expand();
	Long2Short(p);
	dec(p);	
	Short2Long(p);
	return *this;
}


void	Permutator::setIndex(INT n)
{
LONG_INT	fak = 1;

	for ( INT i=Ordnung-2 , j=2 ; i>=0 ; i-- , j++ )
	{	p[i] = n % j;
		fak *= j;
		n /= fak;
	}
	Short2Long(p);
}


INT	Permutator::getIndex()
{
INT	ind = 0;
LONG_INT	fak = 1;

	expand();
INT	pi[Ordnung];
	memcpy(pi, p, Ordnung*sizeof(INT));
	Long2Short(pi);
	for ( INT i=Ordnung-2 ; i>=0 ; i-- )
	{	ind += fak * pi[i];
		fak *= (Ordnung-i);
	}
	return ind;
}
	

LONG_INT	Permutator::resetIndex()
{
	if ( sp )
		delete sp;
	sp = new INT[Ordnung-1];
	memset(sp, 0, (Ordnung-1)*sizeof(INT));
	nfak = 1;
	for ( INT i=1 ; i<=Ordnung ; i++ )
		nfak *= i;
	setIndex(0);
	return nfak;
}


INT	Permutator::nextIndex()
{
 	inc(sp);
	nfak--;
	memcpy(p, sp, (Ordnung-1)*sizeof(INT));
	Short2Long(p);
	if ( nfak==0 )
	{	delete sp;
		sp = NULL;
		return 0;
	}
	return 1;
}


ostream&	operator << (ostream & s, const Permutator & y)
{
	s << "   (";
	for ( INT i=0 ; i<y.getOrdnung() ; i++ )
	{	if ( i )
			s << " ";
		s << i+y.getStart();
	}
	s << ")" << endl;

	if ( y.getC()>=0 )
		s << " ";
	s << y.getC();
	s << "*(";
	for ( INT i=0 ; i<y.getOrdnung() ; i++ )
	{	if ( i )
			s << " ";
		s << y[i];
	}
	s << ")";

	return s;
}


/*
istream&	operator >> (istream& s, Permutator & y)
{
	return s;
}
	
*/
