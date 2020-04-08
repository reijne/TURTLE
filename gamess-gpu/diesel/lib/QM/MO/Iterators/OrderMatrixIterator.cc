//***********************************************************************
//
//	Name:			OrderMatrixIterator.cc
//
//	Description:	iterator for order relation matrices
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			22.07.1997
//
//
//
//
//
//***********************************************************************

#include "OrderMatrixIterator.h"

#include <string>
#include <iomanip>


OrderMatrixIterator::OrderMatrixIterator(INT _rows, INT _columns)
{
	rows = _rows;
	cols = _columns;
	dim = rows*cols;
	if ( dim )
	{
		p = new INT[rows*cols];
		row0 = new INT[rows];
		memset(row0, 0, rows*sizeof(INT));
		col0 = new INT[cols];
		memset(col0, 0, cols*sizeof(INT));
	}
}






OrderMatrixIterator::~OrderMatrixIterator()
{
	if ( dim )
	{
		delete p;
		delete row0;
		delete col0;
	}
}



void	OrderMatrixIterator::first()
{
	if ( dim )
	{
		i = 0;
		p[i] = -2;
		end = 0;
		row0[i/cols] += p[i]==0;
		col0[i%cols] += p[i]==0;
		next();
	}
	else
		end = 1;
}


void	OrderMatrixIterator::next()
{
INT	flag = 0;

	if ( end )
	{
		end = 2;
		return;
	}

	// that's really magic...
	for ( ; ; )
	{
		row0[i/cols] -= p[i]==0;
		col0[i%cols] -= p[i]==0;
		p[i]++;
		row0[i/cols] += p[i]==0;
		col0[i%cols] += p[i]==0;
		if ( row0[i/cols]<=1 && col0[i%cols]<=1 )
		{
			if ( i<dim-1 )
			{
				i++;
				// succeding row must not be smaller
				p[i] = (i>=cols) ? p[i-cols]-1 : -2;
				row0[i/cols] += p[i]==0;
				col0[i%cols] += p[i]==0;
				continue;
			}

			flag = (i==dim-1);

			while ( i>=0 && p[i]==1 )
				i--;

			if ( i<0 )
			{
				end = 1;
				break;
			}
		}

		// succeding column must not be greater
		if ( i % cols && p[i]>=p[i-1] )
		{
			row0[i/cols] -= p[i]==0;
			col0[i%cols] -= p[i]==0;
			i--;
		}

		if ( i<0 )
		{
			end = 1;
			break;
		}


		if ( flag )
			break;
	}

}


INT	OrderMatrixIterator::calcMOListIterator(
	MOType start, MOType end,
	INT tupelClosed,
	const MOType *baseMOsInIrRep,
	MOListIterator & molist,
	Configuration<MOType> &extNotRunning,
	DiffConf<GeneralizedMO> &dcGen,
	INT &sig,
	IrRep irrep) const
{
	// !!!!
	// restriction: max. two externals on right side only!
MOType* s = new MOType[cols];
MOType* e = new MOType[cols];
INT	icols = 0;
	for ( INT c=0 ; c<cols ; c++ )
	{
		s[icols] = rows ? baseMOsInIrRep[rows-1] + 1 : start;
		e[icols] = end;

	INT	cmp = 1;
	INT	iopen = (c>=tupelClosed);
	INT iclosed = !iopen;
		for ( INT r=0 ; r<rows ; r++ )
		{
			cmp = operator()(r, c);
			if ( !cmp )
			{
				dcGen.addExternal(
					Configuration<GeneralizedMO>(),
					Configuration<GeneralizedMO>(iopen, iclosed,
						GeneralizedMO(baseMOsInIrRep[r])));
				extNotRunning += Configuration<MOType>(iopen, iclosed,
					baseMOsInIrRep[r]);
				break;
			}
			if ( cmp==1 )
			{
				s[icols] = r ? baseMOsInIrRep[r-1] + 1 : start;
				e[icols] = baseMOsInIrRep[r] - 1;
				break;
			}
		}
		if ( cmp!=0 )
		{
		INT	add = ( icols>0 && s[icols-1]==s[icols] );
			dcGen.addExternal(
				Configuration<GeneralizedMO>(),
				Configuration<GeneralizedMO>(iopen, iclosed,
					GeneralizedMO(irrep, 0, 
						s[icols] + add, sig++, 2)));

			icols++;
		}
	}

//	cout << "cols = " << cols << endl;
//	cout << "icols = " << icols << endl;

	for ( INT i=0 ; i<icols ; i++ )
		if ( s[i]>e[i] )
                {
                        delete[] e;
                        delete[] s;
			return 1;
                 }

	switch ( icols ) {
	case 0:
                delete[] e;
                delete[] s;
		return 0;
		
	case 1:
		molist += MOListIterator(s[0], e[0], (tupelClosed>0));
                delete[] e;
                delete[] s;
		return 0;
	
	case 2:
		if ( s[0]==s[1] && e[0]==e[1] )
		{
			if ( s[0]<e[0] )
				molist += MOListIterator(s[0], e[0], 0, MOListIterator::lowerTriangular);
			else 
                        {
                                delete[] e;
                                delete[] s;
				return 1;
                        }
		}
		else
			molist += MOListIterator(s[0], e[0], s[1], e[1]);
                delete[] e;
                delete[] s;
        	return 0;
		
	default:
		cout << "Error: OrderMatrixIterator::calcMOListIterator" << endl;
		exit(1);
	}
	delete[] e;
	delete[] s;
	return -1;
}


		

ostream& operator<<(ostream & s, const OrderMatrixIterator & order)
{
INT	*pp = order.p;
	if ( order.dim )
	{
		for ( INT i=0 ; i<order.rows ; i++ )
		{
			s << "(";
			for ( INT j=0 ; j<order.cols ; j++ )
				s << std::setw(3) << *pp++;
			s << ")" << endl;
		}
			s << "{";
		for ( INT j=0 ; j<order.cols ; j++ )
			s << std::setw(3) << order.col0[j];
			s << "}" << endl;
	}
	else
		s << "( free )" << endl;
	s << endl;
	return s;
}


