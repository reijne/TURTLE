//***********************************************************************
//
//	Name:			Spin.cc
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

using namespace std;

#include "Spin.h"
#include "Permutator.h"

#include <iostream>
#include <sstream>


#define MAX(a,b)	((a)>(b) ? (a) : (b))


void	Spin::init(INT _nsum)
{
	nsum = _nsum;
	n = new INT[nsum];
	c = new Complex[nsum];
	p = new SpinBaseType * [nsum];
	for ( INT i=0 ; i<nsum ; i++ )
	{	c[i] = 1;
		n[i] = 0;
	}
}



void	Spin::kill()
{
	if ( nsum )
	{	for ( INT i=0 ; i<nsum ; i++ )
			if ( n[i] )
				delete p[i];
		delete p;
		delete n;
		delete[] c;
	}
	nsum = 0;
}




Spin	&Spin::operator = (const Spin & s)
{
	if ( nsum != s.getNumberOfTerms() )
	{	kill();
		init(s.getNumberOfTerms());
	}

	for ( INT i=0 ; i<nsum ; i++ )
	{	if ( n[i] != s.getNumberOfShells(i) )
		{	if ( n[i] )
				delete p[i];
			n[i] = s.getNumberOfShells(i);
			p[i] = new SpinBaseType[n[i]];
		}
		c[i] = s.getC(i);
		memcpy(p[i], s.getP(i), n[i] * sizeof(SpinBaseType));
	}
	return *this;
}


Spin::Spin(const Spin & s)
{
	init(s.getNumberOfTerms());
	for ( INT i=0 ; i<nsum ; i++ )
	{	n[i] = s.getNumberOfShells(i);
		p[i] = new SpinBaseType[n[i]];
		c[i] = s.getC(i);
		memcpy(p[i], s.getP(i), n[i] * sizeof(SpinBaseType));
	}
}

INT	Spin::getNumberOfElectrons() const
{
INT	ne = 0;
	for ( INT i=0 ; i<n[0] ; i++ )
	{	ne += p[0][i] & 1;
		ne += (p[0][i] & 2) >> 1;
	}
	return ne;	
}


Spin	alpha(INT n)
{
Spin	s(1);
	s.setSpinBinary(n, 1);
	return s;
}


Spin	beta(INT n)
{
Spin	s(1);
	s.setSpinBinary(n, 2);
	return s;
}


INT operator == (Spin a, Spin b)
{
	a.simplify();
	b.simplify();
	if ( a.getNumberOfTerms() != b.getNumberOfTerms() )
		return 0;

INT found;
	for ( INT i=0 ; i<a.getNumberOfTerms() ; i++ )
	{	found = 0;
		for ( INT j=0 ; j<b.getNumberOfTerms() ; j++ )
		{	if ( a.compareSpin(	a.getNumberOfShells(i), a.getP(i), 
								b.getNumberOfShells(j), b.getP(j)) )
			{	if ( a.getC(i) != b.getC(j) )
					return 0;
				found = 1;
				continue;
			}
		}
		if ( !found )
			return 0;
	}
	return 1;
}


INT	Spin::isMultipleOf(Spin a, Complex & e)
{
	a.simplify();
	simplify();

	if ( !a.getNumberOfTerms() )
	{	e = 0;
		return 1;
	}

	if ( a.getNumberOfTerms() != nsum )
		return 0;

INT found;
Complex	ee;
const double eps = 1e-6;
	for ( INT i=0 ; i<a.getNumberOfTerms() ; i++ )
	{	found = 0;
		for ( INT j=0 ; j<nsum ; j++ )
		{	if ( compareSpin(	a.getNumberOfShells(i), a.getP(i), 
								n[j], p[j]) )
			{	ee = a.getC(i)/c[j];
				if ( i )
				{	if ( abs(e-ee)>eps )
						return 0;
				}
				else
					e = ee;
				found = 1;
				continue;
			}
		}
		if ( !found )
			return 0;
	}
	return 1;
}


Spin	& Spin::operator += (const Spin &a)
{
INT	j;

	for ( j=nsum-1 ; j>=0 ; j-- )
		if ( n[j] )
			break;
	j++;

	if ( nsum-j<a.getNumberOfTerms() )
	{	j = nsum;
		nsum += a.getNumberOfTerms();
	
INT	*nn = new INT[nsum];
		memcpy(nn, n, j*sizeof(INT));
		memset(nn+j, 0, a.getNumberOfTerms()*sizeof(INT));
		delete n;
		n = nn;

SpinBaseType	**pn = new SpinBaseType * [nsum];
		memcpy(pn, p, j*sizeof(SpinBaseType * ));
		delete p;
		p = pn;

Complex	*nc = new Complex[nsum];
		memcpy(nc, c, j*sizeof(Complex));
		delete[] c;
		c = nc;
	}

	for ( INT i=0 ; i<a.getNumberOfTerms() ; i++ , j++ )
	{	setNumberOfShells(j, a.getNumberOfShells(i));
		c[j]= a.getC(i);
		memcpy(p[j], a.getP(i), n[j] * sizeof(SpinBaseType));
	}
	return *this;
}


Spin	operator + (const Spin &a, const Spin &b)
{
Spin	s(a.getNumberOfTerms()+b.getNumberOfTerms());
INT	j = 0;
	for ( INT i=0 ; i<a.getNumberOfTerms() ; i++ , j++ )
	{	s.setNumberOfShells(j, a.getNumberOfShells(i));
		s.c[j]= a.getC(i);
		memcpy(s.p[j], a.getP(i), s.n[j] * sizeof(SpinBaseType));
	}
	for ( INT i=0 ; i<b.getNumberOfTerms() ; i++ , j++ )
	{	s.setNumberOfShells(j, b.getNumberOfShells(i));
		s.c[j]= b.getC(i);
		memcpy(s.p[j], b.getP(i), s.n[j] * sizeof(SpinBaseType));
	}
	return s;
}


Spin	operator - (const Spin &a, const Spin &b)
{
Spin	s(a.getNumberOfTerms()+b.getNumberOfTerms());
INT	j = 0;
	for ( INT i=0 ; i<a.getNumberOfTerms() ; i++ , j++ )
	{	s.setNumberOfShells(j, a.getNumberOfShells(i));
		s.c[j]= a.getC(i);
		memcpy(s.p[j], a.getP(i), s.n[j] * sizeof(SpinBaseType));
	}
	for ( INT i=0 ; i<b.getNumberOfTerms() ; i++ , j++ )
	{	s.setNumberOfShells(j, b.getNumberOfShells(i));
		s.c[j]= -b.getC(i);
		memcpy(s.p[j], b.getP(i), s.n[j] * sizeof(SpinBaseType));
	}
	return s;
}


Spin	operator * (const Spin &a, const Spin &b)
{
Spin	s(a.getNumberOfTerms()*b.getNumberOfTerms());
INT	k = 0;
	for ( INT i=0 ; i<a.getNumberOfTerms() ; i++ )
		for ( INT j=0 ; j<b.getNumberOfTerms() ; j++ , k++ )
		{	s.setNumberOfShells(k, 
				MAX(a.getNumberOfShells(i), b.getNumberOfShells(j)));
			s.c[k] = a.getC(i) * b.getC(j);
			if ( a.getNumberOfShells(i) > b.getNumberOfShells(j) )
			{	memcpy(s.p[k], a.getP(i), s.n[k] * sizeof(SpinBaseType));
				for ( INT l=0 ; l<b.getNumberOfShells(j) ; l++ )
					s.p[k][l] |= b.getP(j)[l];
			}
			else
			{	memcpy(s.p[k], b.getP(j), s.n[k] * sizeof(SpinBaseType));
				for ( INT l=0 ; l<a.getNumberOfShells(i) ; l++ )
					s.p[k][l] |= a.getP(i)[l];
			}
		}
	return s;
}


Spin	operator * (const Complex & c, const Spin & a)
{
Spin	s(a.getNumberOfTerms());

	for ( INT i=0 ; i<a.getNumberOfTerms() ; i++ )
	{	s.setNumberOfShells(i, a.getNumberOfShells(i));
		s.c[i]= c*a.getC(i);
		memcpy(s.p[i], a.getP(i), s.n[i] * sizeof(SpinBaseType));
	}
	return s;
}


Spin	operator * (const double & c, const Spin & a)
{
Spin	s(a.getNumberOfTerms());

	for ( INT i=0 ; i<a.getNumberOfTerms() ; i++ )
	{	s.setNumberOfShells(i, a.getNumberOfShells(i));
		s.c[i]= c*a.getC(i);
		memcpy(s.p[i], a.getP(i), s.n[i] * sizeof(SpinBaseType));
	}
	return s;
}


Spin	operator - (Spin & s)
{
Spin	r(s.getNumberOfTerms());
	
	for ( INT i=0 ; i<s.getNumberOfTerms() ; i++ )
	{	r.setNumberOfShells(i, s.getNumberOfShells(i));
		r.setC(i, -s.getC(i));
		memcpy(r.getP(i), s.getP(i), 
			s.getNumberOfShells(i) * sizeof(SpinBaseType));
	}
	return r;
}


void	Spin::setSpinBinary(INT _n, SpinBaseType _c)
{
	if ( n[0]!=_n )
	{	if ( n[0] )
			delete p[0];
		n[0] = _n;
		p[0] = new SpinBaseType[n[0]];
		c[0] = 1;
	}
	for ( INT i=0 ; i<_n-1 ; i++ )
		p[0][i] = 0;
	p[0][_n-1] = _c;
}


void	Spin::setSpin (INT Summand, LONG_INT l)
{
INT	i;

	if ( Summand >= nsum )
		return;
	
	for ( i=15 ; i >= 0 ; i-- )
	{	if ( l & (3 << (i*2)) )
			break; 
	}

	if ( n[Summand] != i+1 )
	{	if ( n[Summand] )
			delete p[Summand];
		n[Summand] = i + 1;
		p[Summand] = new SpinBaseType[i+1];
	}
	c[Summand] = 1;
	for ( ; i >= 0 ; i-- )
		p[Summand][i] = (l >> (i*2)) & 3;
}

LONG_INT	Spin::getSpin (INT Summand)
{
	if ( Summand >= nsum )
		return 0;

INT	l = 0;	
	for ( INT i = 0 ; i<n[Summand] ; i++ )
		l |= ((p[Summand][i]==2) ? 2 : 1) << (2*i);
	return l;
}


LONG_INT	Spin::getSpinCompact (INT Summand)
{
	if ( Summand >= nsum )
		return 0;

INT	l = 0;	
	for ( INT i = 0 ; i<n[Summand] ; i++ )
		l |= (p[Summand][i]==2) << i;
	return l;
}



Spin	Spin::Sz()
{
Complex	a;
Spin	s(nsum);

	for ( INT j=0 ; j<nsum ; j++ )
	{	a = 0;
		s.setNumberOfShells(j, n[j]);
		for ( INT i=0 ; i<n[j] ; i++ )
		{	a += 0.5 * ((p[j][i] & 1)!=0);
			a -= 0.5 * ((p[j][i] & 2)!=0);
			s.setSpinShell(j, i, getSpinShell(j, i));
		}
		s.setC(j, a*c[j]);
	}
	return s;
}


Spin	Spin::Sxy(INT isY)
{
Complex	a;
INT	nn = 0;
static INT swap[] = { 0, 2, 1, 3 };

	for ( INT j=0 ; j<nsum ; j++ )
		nn += n[j];

Spin	s(nn);


	nn = 0;
	for ( INT j=0 ; j<nsum ; j++ )
		for ( INT k=0 ; k<n[j] ; k++ )
		{	s.setNumberOfShells(nn, n[j]);
			for ( INT i=0 ; i<n[j] ; i++ )
			{	if ( i==k )
				{	s.setSpinShell(nn, i, swap[getSpinShell(j, i)]);
					a = 0.5 * ((p[j][i] & 1)!=0);
					a += (isY ? -0.5 : 0.5) * ((p[j][i] & 2)!=0);
					if ( isY )
						a *= Complex(0, 1);
				}
				else
					s.setSpinShell(nn, i, getSpinShell(j, i));
			}
			s.setC(nn, a*c[j]);
			nn++;
		}
	return s;
}



Spin	Spin::Spm(INT isM)
{
Complex	a;
INT	nn = 0;
static INT table[2][4] = 
	{	{ 0, 0, 1, 3 },
		{ 0, 2, 0, 3 }
	};

	for ( INT j=0 ; j<nsum ; j++ )
		nn += n[j];

Spin	s(nn);
INT	h;

	nn = 0;
	for ( INT j=0 ; j<nsum ; j++ )
		for ( INT k=0 ; k<n[j] ; k++ )
		{	s.setNumberOfShells(nn, n[j]);
			a = 1;
			for ( INT i=0 ; i<n[j] ; i++ )
			{	if ( i==k )
				{	s.setSpinShell(nn, i, (h=table[isM][getSpinShell(j, i)]));
					a = (h>0);
				}
				else
					s.setSpinShell(nn, i, getSpinShell(j, i));
			}
			s.setC(nn, a*c[j]);
			nn++;
		}

	return s;
}


Spin	Spin::SS()
{
Spin	s = Sm().Sp() - Sz() + Sz().Sz();

	return s;
}



Complex	Spin::Integrate()
{
Complex	_c = 0;
INT j;
	for ( INT i=0 ; i<nsum ; i++ )
	{	for ( j=0 ; j<n[i] ; j++ )
			if ( p[i][j]==3 )
				break;
		if ( p[i][j]!=3 )
			_c += c[i];
	}
	return _c;
}




void	Spin::permute(Permutator &P)
{
Spin	s(nsum);
INT	maxn = 0;

	for ( INT i=0 ; i<nsum ; i++ )
	{	maxn = MAX(maxn, n[i]);
		if ( n[i] )
			s.setNumberOfShells(i, n[i]);
	}

INT	lorb[maxn*2];
char	lspin[maxn*2];	
INT	e;
SpinBaseType	*sp;

	for ( INT i=0 ; i<nsum ; i++ )
	{	e = 0;
		sp = s.getP(i);
		for ( INT j=0 ; j<n[i] ; j++ )
		{	if ( p[i][j] & 1 )
			{	lorb[e] = j;
				lspin[e] = 1;
				e++;
			}
			if ( p[i][j] & 2 )
			{	lorb[e] = j;
				lspin[e] = 2;
				e++;
			}
		}
		P.permute(lspin, sizeof(char));
		for ( INT j=0 ; j<e ; j++ )
			sp[lorb[j]] |= lspin[j];
		s.setC(i, this->getC(i));
	}
	*this = s;
}



void	Spin::DoSymmetrize(INT Anti)
{
Permutator	P(getNumberOfElectrons());

LONG_INT	N = P.resetIndex();
Spin	*s = new Spin(N*nsum);
Spin	h;
	do
	{	h = *this;
		h.permute(P);
		if ( Anti )
			*s += P.getParity()*h;
		else
			*s += h;
	}		
	while ( P.nextIndex() );
	
	*this = (1/sqrt(N)) * *s;
	delete s;
}


INT	Spin::compareSpin(INT n1, SpinBaseType *p1, INT n2, SpinBaseType *p2)
{
INT	i;
INT	nmin = (n1>n2) ? n2 : n1;

	for ( i=0 ; i<nmin ; i++ )
		if ( p1[i]!=p2[i] )
			return 0;

	if ( n1>n2 )
	{	for ( ; i<n1 ; i++ )
			if ( p1[i] )
				return 0;
	}
	else
	{	for ( ; i<n2 ; i++ )
			if ( p2[i] )
				return 0;
	}

	return 1;
}



Spin Spin::simplify()
{
	for ( INT i=0 ; i<nsum ; i++ )
	{	if ( abs(c[i])<1e-10 )
			continue;
		for ( INT j=i+1 ; j<nsum ; j++ )
		{	if ( compareSpin(n[i], p[i], n[j], p[j]) )
			{	c[i] += c[j];
				c[j] = 0;
			}
		}
	}

INT	j = 0;
INT	nullen = 0;

	for ( INT i=0 ; i<nsum ; i++ )
	{	if ( abs(c[i])>1e-10 )
		{	if ( i>j )
			{	c[j] = c[i];
				c[i] = 0;
				n[j] = n[i];
				p[j] = p[i];
			}
			j++;
		}
		else
		{	nullen++;
			delete p[i];
		}
	}
	nsum -= nullen;
	return *this;
}






ostream&	operator << (ostream & s, const Spin & y)
{
INT	i, j;
ostringstream	*dummy;


/*	for ( j=0 ; j<y.getNumberOfTerms() ; j++ )
	{	s << "\n[";
		dummy = new ostringstream;
		*dummy << y.getC(j) << " ";
		for ( i=0 ; i<dummy->tellp() ; i++ )
			s << " ";
		delete dummy;
		for ( i=0 ; i<y.getNumberOfShells(j) ; i++ )
		{	if ( y.getSpinShell(j, i) & 1 )
				s << " ";
			if ( y.getSpinShell(j, i) & 2 )
				s << "_";
		}
		s << "]\n[" << y.getC(j) << " ";
		for ( i=0 ; i<y.getNumberOfShells(j) ; i++ )
		{	if ( y.getSpinShell(j, i) & 1 )
				s << ((char) ('a'+i));
			if ( y.getSpinShell(j, i) & 2 )
				s << ((char) ('a'+i));
		}
		s << "]";
		if ( j == y.getNumberOfTerms()-1 )
			s << "\n";
		else
			s << " + ";
	}
*/

	s << "\n";
	if ( !y.getNumberOfTerms() )
	{	s << "0\n";
		return s;
	}

	for ( j=0 ; j<y.getNumberOfTerms() ; j++ )
	{	dummy = new ostringstream;
		*dummy << y.getC(j) << "*";
		for ( i=0 ; i<dummy->tellp() ; i++ )
			s << " ";
		delete dummy;
		for ( i=0 ; i<y.getNumberOfShells(j) ; i++ )
		{	if ( y.getSpinShell(j, i) & 1 )
				s << " ";
			if ( y.getSpinShell(j, i) & 2 )
				s << "_";
		}
		s << "   ";
	}
	s << "\n";
	for ( j=0 ; j<y.getNumberOfTerms() ; j++ )
	{	s << y.getC(j) << "*";
		for ( i=0 ; i<y.getNumberOfShells(j) ; i++ )
		{	if ( y.getSpinShell(j, i) & 1 )
				s << ((char) ('a'+i));
			if ( y.getSpinShell(j, i) & 2 )
				s << ((char) ('a'+i));
		}
		if ( j == y.getNumberOfTerms()-1 )
			s << "\n";
		else
			s << " + ";
	}
	return s;
}

/*
istream&	operator >> (istream& s, Spin & y)
{
	return s;
}
*/
