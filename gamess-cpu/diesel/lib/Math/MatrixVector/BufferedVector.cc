//***********************************************************************
//
//	Name:			BufferedVector.cc
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			09. Nov 1998
//
//***********************************************************************

#include "BufferedVector.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#include <string>	
#include <stdlib.h>
#include <math.h>

#include "../../Container/DiskBuffer.h"
#include "../../Container/TempDir.h"

using namespace std;

template <class T>
INT	BufferedVector<T>::id = 0;



template <class T>
BufferedVector<T>::BufferedVector(INT n, INT _bufSize) :
	Vector<T>(n, (T*) 0)
{
INT	hi = n>_bufSize ? _bufSize : n;
	shift = (INT) ceil(log((double)hi)/log((double)2));
	bufSize = 1 << shift;
	mask = bufSize - 1;
	dirty = 0;
//	cout << shift << " " << bufSize
//		<< " " << mask << endl;
//	cout << bufSize << ", N=" << N << endl;

	createFile();
}

template <class T>
BufferedVector<T>::BufferedVector(const DiskBuffer &buf, INT i, INT _bufSize) :
	Vector<T>(buf.getObjectSize()/sizeof(T), (T*) 0)
{
INT	hi = buf.getObjectSize()/(INT)sizeof(T)>_bufSize ? _bufSize : buf.getObjectSize()/sizeof(T);

	shift = (INT) ceil(log((double)hi)/log((double)2));
	bufSize = 1 << shift;
	mask = bufSize - 1;
	dirty = 0;
	createFile();


T	*h = new T[this->N];
	buf.get(i, h);
	put(h);
	delete h;
}

template <class T>
BufferedVector<T>::BufferedVector(const BufferedVector<T> &b)
{
	bufSize = b.bufSize;
	mask = b.mask;
	shift = b.shift;
	this->N = b.N;
	blockInBuffer = -1;
	dirty = 0;
	createFile();

	((BufferedVector<T> &) b).flush();
	lseek(b.fd, 0, SEEK_SET);
	lseek(fd, 0, SEEK_SET);
	for ( INT i=0 ; i<this->N/bufSize ; i++ )
	{
		read(b.fd, this->p, bufSize*sizeof(T));
		write(fd, this->p, bufSize*sizeof(T));
	}
	read(b.fd, this->p, (this->N%bufSize)*sizeof(T));
	write(fd, this->p, (this->N%bufSize)*sizeof(T));
}


template <class T>
BufferedVector<T> & BufferedVector<T>::operator = (const BufferedVector<T> &b)
{
	if ( this==&b )
		return *this;
	bufSize = b.bufSize;
	mask = b.mask;
	shift = b.shift;
	this->N = b.N;
	blockInBuffer = -1;
	dirty = 0;
	
	if ( this->p )
		delete this->p;
	this->p = new T[bufSize];
	
	((BufferedVector<T> &) b).flush();
	lseek(b.fd, 0, SEEK_SET);
	lseek(fd, 0, SEEK_SET);
	for ( INT i=0 ; i<this->N/bufSize ; i++ )
	{
		read(b.fd, this->p, bufSize*sizeof(T));
		write(fd, this->p, bufSize*sizeof(T));
	}
	read(b.fd, this->p, (this->N%bufSize)*sizeof(T));
	write(fd, this->p, (this->N%bufSize)*sizeof(T));
	return *this;
}


template <class T>
BufferedVector<T>::~BufferedVector()
{
	close(fd);
	unlink(FileName);
}



template <class T>
void	BufferedVector<T>::createFile()
{
	sprintf(FileName, "%s/Vector.%d.%d", TempDir.chars(), getpid(), id);
	fd = open(FileName, O_RDWR | O_CREAT | O_TRUNC,
		S_IRUSR | S_IRGRP | S_IROTH | S_IWUSR);
	
//	cout << bufSize << ", N=" << N << endl;
	this->p = new T[bufSize];
	clear();
	
	id++;
}


template <class T>
void	BufferedVector<T>::put(const T *pp)
{
	blockInBuffer = -1;
	dirty = 0;
	lseek(fd, 0, SEEK_SET);
	for ( INT i=0 ; i<this->N/bufSize ; i++ , pp += bufSize )
		write(fd, pp, bufSize*sizeof(T));
	write(fd, pp, (this->N%bufSize)*sizeof(T));
}


template <class T>
void	BufferedVector<T>::get(T *pp) const
{
	((BufferedVector<T> *) this)->flush();
	lseek(fd, 0, SEEK_SET);
	for ( INT i=0 ; i<this->N/bufSize ; i++ , pp += bufSize )
		read(fd, pp, bufSize*sizeof(T));
	read(fd, pp, (this->N%bufSize)*sizeof(T));
}


//---------------------------------------------------------------

template <class T>
void	BufferedVector<T>::fillBuf(INT i)
{
	if ( dirty )
		flush();

	blockInBuffer = i >> shift;
//	cout << this << " R: " << blockInBuffer << endl;
	lseek(fd, blockInBuffer*bufSize*sizeof(T), SEEK_SET);
	read(fd, this->p, bufSize*sizeof(T));
}

template <class T>
void	BufferedVector<T>::flush()
{
	if ( dirty )
	{
//		cout << this << " W: " << blockInBuffer << endl;
		lseek(fd, blockInBuffer*bufSize*sizeof(T), SEEK_SET);
		write(fd, this->p, bufSize*sizeof(T));
		dirty = 0;
	}
}

//---------------------------------------------------------------

BufferedVector<Complex> operator+(const BufferedVector<Complex> &v1, const BufferedVector<double> &v2)
{
INT	i;
	if ( v1.getDim() != v2.getDim() )
	{	v1.Fehler(1);
BufferedVector<Complex> 	v;
		return v;
	}
BufferedVector<Complex>	v(v1.getDim());
	for ( i=0 ; i<v1.getDim() ; i++ )
		v[i] = v1[i] + v2[i];
	return v;
}

BufferedVector<Complex> operator-(const BufferedVector<Complex> &v1, const BufferedVector<double> &v2)
{
INT	i;
	if ( v1.getDim() != v2.getDim() )
	{	v1.Fehler(1);
BufferedVector<Complex> 	v;
		return v;
	}
BufferedVector<Complex>	v(v1.getDim());
	for ( i=0 ; i<v1.getDim() ; i++ )
		v[i] = v1[i] - v2[i];
	return v;
}

Complex operator*(const BufferedVector<Complex> &v1, const BufferedVector<double> &v2)
{
INT	i;
Complex	sum = 0;
	if ( v1.getDim() != v2.getDim() )
	{	v1.Fehler(1);
		return sum;
	}
	for ( i=0 ; i<v1.getDim() ; i++ )
		sum += v1[i] * v2[i];
	return sum;
}

Complex operator^(const BufferedVector<Complex> &v1, const BufferedVector<double> &v2)
{
INT	i, n;
Complex	sum = 0;
	if ( v1.getDim() != v2.getDim() )
	{	v1.Fehler(1);
		return sum;
	}
	n = v1.getDim();
	for ( i=0 ; i<n ; i++ )
		sum += v1[n-i-1] * v2[i];
	return sum;
}




BufferedVector<Complex> operator+(const BufferedVector<double> &v1, const BufferedVector<Complex> &v2)
{	return	v2+v1;	};

BufferedVector<Complex> operator-(const BufferedVector<double> &v1, const BufferedVector<Complex> &v2)
{	return	-(v2-v1);	};

Complex operator*(const BufferedVector<double> &v1, const BufferedVector<Complex> &v2)
{	return	v2*v1;	};

Complex operator^(const BufferedVector<double> &v1, const BufferedVector<Complex> &v2)			// diskrete Faltung
{	return	v2^v1;	};




//---------------------------------------------------------------






//*********************************************************************
//*********************************************************************
//*********************************************************************
//*********************************************************************
//*********************************************************************
//*********************************************************************





template <class T> BufferedVector<T> operator*(T a, const BufferedVector<T> &v1)
{
INT	i;
BufferedVector<T>	v(v1.getDim());
	for ( i=0 ; i<v1.getDim() ; i++ )
		v[i] = a*v1[i];
	return v;
}


template <class T> 
BufferedVector<T> BufferedVector<T>::operator*(T a) const
{
INT	i;
BufferedVector<T>	v(getDim());
	for ( i=0 ; i<this->N ; i++ )
		v[i] = a*v[i];
	return v;
}


template <class T> BufferedVector<T> operator+(const BufferedVector<T> &v1, const BufferedVector<T> &v2)
{
INT	i;
	if ( v1.getDim() != v2.getDim() )
	{	v1.Fehler(1);
BufferedVector<T> 	v;
		return v;
	}
BufferedVector<T>	v(v1.getDim());
	for ( i=0 ; i<v1.getDim() ; i++ )
		v[i] = v1[i] + v2[i];
	return v;
}

template <class T> BufferedVector<T> operator-(const BufferedVector<T> &v1, const BufferedVector<T> &v2)
{
INT	i;
	if ( v1.getDim() != v2.getDim() )
	{	v1.Fehler(1);
BufferedVector<T> 	v;
		return v;
	}
BufferedVector<T>	v(v1.getDim());
	for ( i=0 ; i<v1.getDim() ; i++ )
		v[i] = v1[i] - v2[i];
	return v;
}

template <class T> T operator*(const BufferedVector<T> &v1, const BufferedVector<T> &v2)
{
INT     i;
T       sum = 0;
        if ( v1.getDim() != v2.getDim() )
        {       v1.Fehler(1);
                return sum;
        }
        for ( i=0 ; i<v1.getDim() ; i++ )
                sum += v1[i] * v2[i];
        return sum;
}

double operator*(const BufferedVector<double> &v1, const BufferedVector<double> &v2)
{
INT	i;
double	sum = 0;
	if ( v1.getDim() != v2.getDim() )
	{	v1.Fehler(1);
		return sum;
	}
	for ( i=0 ; i<v1.getDim() ; i++ )
		sum += v1[i] * v2[i];
	return sum;
}

float operator*(const BufferedVector<float> &v1, const BufferedVector<float> &v2)
{
INT	i;
double	sum = 0;
	if ( v1.getDim() != v2.getDim() )
	{	v1.Fehler(1);
		return sum;
	}
	for ( i=0 ; i<v1.getDim() ; i++ )
		sum += v1[i] * v2[i];
	return sum;
}

template <class T> T operator^(const BufferedVector<T> &v1, const BufferedVector<T> &v2)
{
INT	i, n;
T 	sum = 0;
	if ( v1.getDim() != v2.getDim() )
	{	v1.Fehler(1);
		return sum;
	}
	n = v1.getDim();
	for ( i=0 ; i<n ; i++ )
		sum += v1[n-i-1] * v2[i];
	return sum;
}




template <class T>	T BufferedVector<T>::Sum()
{
INT	i;
T	sum = 0;
	for ( i=0 ; i<this->N ; i++ )
		sum += (*this)[i];
	return sum;
}



template <class T> 	void BufferedVector<T>::Fehler(INT FehlerNr) const
{
	puts("Fehler: ");
	switch ( FehlerNr )
	{	case 1	:	puts("Falsche Dimension");
					break;
		case 2	:	puts("Falsche Dimension (!=3) bei BufferedVectorprodukt");
					break;
	}
}


template <class T> 
INT BufferedVector<T>::operator==(const BufferedVector<T> & v) const
{
INT	i;
	if ( getDim() != v.getDim() )
	{	Fehler(1);
		return 0;
	}
	for ( i=0 ; i<getDim() ; i++ )
		if ( (*this)[i] != v[i] )
			return 0;
	return 1;
}


template <class T> 
INT BufferedVector<T>::operator!=(const BufferedVector<T> & v) const
{
INT	i;
	if ( getDim() != v.getDim() )
	{	Fehler(1);
		return 0;
	}
	for ( i=0 ; i<getDim() ; i++ )
		if ( (*this)[i] != v[i] )
			return 1;
	return 0;
}



template <class T> BufferedVector<T> & BufferedVector<T>::operator+=(const BufferedVector<T> &v1)
{
INT	i;
	if ( this->getDim() != v1.getDim() )
	{	v1.Fehler(1);
		return *this;
	}
	for ( i=0 ; i<v1.getDim() ; i++ )
		(*this)[i] += v1[i];
	return *this;
}

template <class T> BufferedVector<T> & BufferedVector<T>::operator-=(const BufferedVector<T> &v1)
{
INT	i;
	if ( this->getDim() != v1.getDim() )
	{	v1.Fehler(1);
		return *this;
	}
	for ( i=0 ; i<v1.getDim() ; i++ )
		(*this)[i] -= v1[i];
	return *this;
}




template <class T> BufferedVector<T> BufferedVector<T>::operator+()
{
	return *this;
}

template <class T> BufferedVector<T> BufferedVector<T>::operator-()
{
INT	i;
BufferedVector<T>	v(this->getDim());
	for ( i=0 ; i< this->getDim() ; i++ )
		v[i] = - (*this)[i];

	return v;
}



template <class T> BufferedVector<T> & BufferedVector<T>::operator*=(T a)
{
	for ( INT i=0 ; i<getDim() ; i++ )
		(*this)[i] *= a;
	return *this;
}

template <class T> BufferedVector<T> & BufferedVector<T>::operator/=(T a)
{
	for ( INT i=0 ; i<getDim() ; i++ )
		(*this)[i] /= a;
	return *this;
}


template <class T> void BufferedVector<T>::normalize(double Order)
{
	*this /= (T)getNorm(Order);
}

template <class T> void BufferedVector<T>::normalize1()
{
	*this /= (T)getNorm1();
}

template <class T> void BufferedVector<T>::normalize2()
{
	*this /= (T)getNorm2();
}

template <class T> void BufferedVector<T>::normalizeInf()
{
	*this /= (T)getInfNorm();
}



template <class T> double BufferedVector<T>::getNorm(double Order) const
{
double	Norm = 0;

	for ( INT i=0 ; i<this->N ; i++ )
		Norm += pow((double) abs((*this)[i]), Order);
			
	return pow(Norm, 1/Order);
}


	


template <class T> double BufferedVector<T>::getNorm1() const
{
double	Norm = 0;

	for ( INT i=0 ; i<this->N ; i++ )
		Norm += abs((*this)[i]);
			
	return Norm;
}


template <class T> double BufferedVector<T>::getNorm2() const
{
double	Norm = 0;

	for ( INT i=0 ; i<this->N ; i++ )
	{
	T	h = (*this)[i];
		Norm += abs(h*h);
	}
			
	return sqrt(Norm);
}


template <class T> double BufferedVector<T>::getNorm2Sqr() const
{
double	Norm = 0;

	for ( INT i=0 ; i<this->N ; i++ )
	{
	T	h = (*this)[i];
		Norm += abs(h*h);
	}
			
	return Norm;
}


template <class T> double BufferedVector<T>::getInfNorm() const
{
double	Norm = 0;
double	n;
	for ( INT i=0 ; i<this->N ; i++ )
		if ( (n=abs((*this)[i]))>Norm )
			Norm = n;
			
	return Norm;
}


template <class T>
void	BufferedVector<T>::clear()
{
//	cout << bufSize << ", N=" << N << endl;
	memset(this->p, 0, bufSize*sizeof(T));
	
	lseek(fd, 0, SEEK_SET);
	for ( INT i=0 ; i<this->N/bufSize ; i++ )
		write(fd, this->p, bufSize*sizeof(T));
	write(fd, this->p, (this->N%bufSize)*sizeof(T));
}


template <class T> ostream& operator<<(ostream& s, const BufferedVector<T> & v)
{
INT	i;
	s << "[ ";
	for ( i=0 ; i<v.getDim() ; i++ )
	{	s << v[i] << ' ';
	}
	s << ']';
	return s;
}


//*************************************************************************
//*************************************************************************
//*************************************************************************
//*************************************************************************


template class BufferedVector<INT>;
template class BufferedVector<float>;
template class BufferedVector<double>;
template class BufferedVector<Complex>;


#define INSTANTIATE_TEMPLATE_FUNCTIONS \
template class BufferedVector<T> operator*(T, const BufferedVector<T> &);\
template class BufferedVector<T> operator+(const BufferedVector<T> &, const BufferedVector<T> &);\
template class BufferedVector<T> operator-(const BufferedVector<T> &, const BufferedVector<T> &);\
template ostream& operator << (ostream& s, const BufferedVector<T> & v); \
T operator*(const BufferedVector<T> &, const BufferedVector<T> &);\
T operator^(const BufferedVector<T> &, const BufferedVector<T> &);



// Macke im GNU C++ 2.7.0:
// Instantiierung der Nicht-Memberfunctions, die integralen Typ
// zurückgeben, funktioniert nicht! Daher folgende
// "Zwangs-Instantiierung":

static void	dummy()
{
Complex	h;
BufferedVector<INT>	vi;
	vi!=vi;
	vi==vi;
	h = vi*vi;
	h = vi^vi;
	
BufferedVector<float>	vf;
	vf!=vf;
	vf==vf;
	h = vf*vf;
	h = vf^vf;
	
BufferedVector<double>	vd;
	vd!=vd;
	vd==vd;
	h = vd*vd;
	h = vd^vd;
	
BufferedVector<Complex>	vc;
	vc!=vc;
	vc==vc;
	h = vc*vc;
	h = vc^vc;
	
}



//template class INT operator!=(const BufferedVector<INT> &, const BufferedVector<INT> &);

#define T INT
INSTANTIATE_TEMPLATE_FUNCTIONS
#undef T

#define T float
INSTANTIATE_TEMPLATE_FUNCTIONS
#undef T

#define T double
INSTANTIATE_TEMPLATE_FUNCTIONS
#undef T

#define T Complex
INSTANTIATE_TEMPLATE_FUNCTIONS
#undef T

