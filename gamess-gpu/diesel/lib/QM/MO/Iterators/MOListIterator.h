//***********************************************************************
//
//	Name:			MOListIterator.h
//
//	Description:	iterators for MO lists
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			26.03.1997
//
//
//
//
//
//***********************************************************************


#ifndef __MOListIterator_h
#define __MOListIterator_h

#include "../../../../config.h"


#include "MOList.h"

using std::cout;
using std::endl;

class MOListIterator {
public:
enum Mode { noIndex, oneIndex, twoIndex, lowerTriangular};

	MOListIterator();
	MOListIterator(MOType start, MOType end, INT square, Mode mode = oneIndex);
	MOListIterator(MOType istart, MOType iend, MOType jstart, MOType jend);
	
	~MOListIterator() {}
	
		
	INT	next();	// returns true if not end
	
	INT	getN() const;	

	//----------------------------------------------------------------	

	INT isNewBlock() const;
	INT isSquare() const;

	//----------------------------------------------------------------	

	MOListIterator & operator += (const MOListIterator &);


	//----------------------------------------------------------------	

	friend ostream& operator<<(ostream & s, const MOListIterator &);

	//----------------------------------------------------------------

Mode	mode;	//	mode of loop
MOType	i;			//	loop variable i
MOType	j;			//	loop variable j
MOType	si;			//	start value i
MOType	ei;			//	end value i
MOType	sj;			//	start value j
MOType	ej;			//	end value j

private:
INT	newBlock;	//	flag if first loop variable has changed
				//	(only for "twoIndex" and "lowerTriangular" modes)

INT	square;		//	flag if running index is a doubly occupied MO
};



inline
MOListIterator::MOListIterator()
{
	mode = noIndex;
	i = j = -1;
	newBlock = 0;
}


inline
MOListIterator::MOListIterator(MOType start, MOType end, INT _square, Mode _mode)
{
	mode = _mode;
	i = si = start;
	square = _square;
	if ( mode==lowerTriangular )
	{
		j = i+1;
		ej = end;
		ei = j-1;
		newBlock = 1;
	}
	else
	{
		ei = end;
		newBlock = 0;
	}	
}

inline
MOListIterator::MOListIterator(MOType istart, MOType iend, MOType jstart, MOType jend)
{
	mode = twoIndex;
	i = si = istart;
	j = sj = jstart;
	ei = iend;
	ej = jend;
	newBlock = 1;
	square = 0;
}



inline
INT	MOListIterator::next()
{
	if ( mode == noIndex )
		return 0;
		
	if ( ++i <= ei )
	{
		newBlock = 0;
		return 1;
	}

	if ( mode == oneIndex )
		return 0;

	newBlock = 1;
	i = si;
	if ( mode==lowerTriangular )
		ei = j;

	if ( ++j <= ej )
		return 1;

	return 0;
}

inline
INT	MOListIterator::getN() const
{
	switch ( mode ) {
	case noIndex:
		return 1;

	case oneIndex:
		return ei-si+1;

	case twoIndex:
		return (ei-si+1)*(ej-sj+1);

	case lowerTriangular:
		return (ej-si)*(ej-si+1)/2;
	}
	return 0;
}


inline
MOListIterator & MOListIterator::operator += (const MOListIterator & iter)
{
	switch ( mode ) {
	case noIndex:
		return *this = iter;
		
	case oneIndex:
		switch ( iter.mode ) {
		case noIndex:
			return *this;

		case oneIndex:
			mode = twoIndex;
			j = sj = iter.si;
			ej = iter.ei;
			newBlock = 1;
			square = 0;
			return *this;

		default:
			cout << "Error: in MOListIterator::operator +=" << endl;
			exit(1);
		}
	
	default:
		cout << "Error: in MOListIterator::operator +=" << endl;
		exit(1);
	}
	return *this;
}

inline
INT	MOListIterator::isNewBlock() const
{	return newBlock;	}

inline
INT	MOListIterator::isSquare() const
{	return square;	}

#endif

