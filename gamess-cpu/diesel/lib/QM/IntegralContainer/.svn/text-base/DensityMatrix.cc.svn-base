//***********************************************************************
//
//	Name:			DensityMatrix.cc
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			23.05.1998
//
//
//
//
//
//***********************************************************************

#include "DensityMatrix.h"

#include "../../Math/FortranNumeric/FortranEigenProblems.h"

#include "MOTrafo.h"

#include <iostream>
#include <string>
#include <iomanip>

using namespace std;

DensityMatrix::DensityMatrix() :
	MOIrReps()
{
	diagOcc = 0;
	nFrozen = NULL;
	nDeleted = NULL;
	inIrRepStored = NULL;
	p = NULL;
}



DensityMatrix::DensityMatrix(istream &is) :
	MOIrReps(is)
{
	is >> diagOcc;

	nFrozen = new INT [IrReps];
	nDeleted = new INT [IrReps];

	for ( INT i=0 ; i<IrReps ; i++ )
		is >> nFrozen[i];
	for ( INT i=0 ; i<IrReps ; i++ )
		is >> nDeleted[i];

	inIrRepStored = new INT [IrReps];
	for ( INT i=0 ; i<IrReps ; i++ )
		inIrRepStored[i] = inIrRep[i] - nFrozen[i] - nDeleted[i];
	
INT* flag = new INT[IrReps*IrReps];

	for ( INT i=0 ; i<IrReps*IrReps ; i++ )
		is >> flag[i];

	p = new Type *[IrReps*IrReps];
	for ( INT i=0 ; i<IrReps*IrReps ; i++ )
		if ( flag[i] )
		{
//		INT	iSym = i/IrReps;
//		INT jSym = i%IrReps;
//		cout << iSym << " " << jSym << endl;
		INT n = inIrRepStored[i/IrReps]*inIrRepStored[i%IrReps];
//			cout << "n=" << n << endl;
			p[i] = new Type[n];
			for ( INT j=0 ; j<n ; j++ )
			{
				is >> p[i][j];
//				cout << i << " " << j << " " << p[i][j] << endl;
			}
		}
		else
			p[i] = NULL;

//	for ( INT i=0 ; i<IrReps ; i++ )
//		inIrRep[i] += nFrozen[i] + nDeleted[i];

	maxMO = 0;
	init();	

	delete flag;
}



DensityMatrix::DensityMatrix(
	const MOIrReps &moirreps, const INT *_nFrozen, const INT *_nDeleted,
	const Type *_p, INT transition) : MOIrReps(moirreps)
{
INT* flag = new INT[IrReps*IrReps];
	memset(flag, 0, IrReps*IrReps*sizeof(INT));


	for ( INT i=0 ; i<maxMO*maxMO ; i++ )
		flag[getIrRep(i/maxMO+1)*IrReps + getIrRep(i%maxMO+1)] |= (_p[i] != 0);


	p = new Type *[IrReps*IrReps];
	for ( INT i=0 ; i<IrReps*IrReps ; i++ )
		if ( flag[i] )
		{
		INT n = inIrRep[i/IrReps]*inIrRep[i%IrReps];
			p[i] = new Type[n];
			for ( INT j=0 ; j<n ; j++ )
			{
			INT	iMO = j / inIrRep[i%IrReps] + getStartMO(i/IrReps) - 1;
			INT	jMO = j % inIrRep[i%IrReps] + getStartMO(i%IrReps) - 1;

				p[i][j] = _p[maxMO*iMO + jMO];
			}
		}
		else
			p[i] = NULL;

		
	nFrozen = new INT[IrReps];
	memcpy(nFrozen, _nFrozen, IrReps*sizeof(INT));
	nDeleted = new INT[IrReps];
	memcpy(nDeleted, _nDeleted, IrReps*sizeof(INT));
	inIrRepStored = new INT [IrReps];
	memcpy(inIrRepStored, inIrRep, IrReps*sizeof(INT));

	for ( INT i=0 ; i<IrReps ; i++ )
		inIrRep[i] += nFrozen[i] + nDeleted[i];

	delete MOSymmetry;
	delete MOinIrRep;
	delete IrRepStart;

	diagOcc = transition ? 0 : 2;

	maxMO = 0;
	init();	

	delete flag;
}


DensityMatrix::~DensityMatrix()
{
	if ( p )
	{
		for ( INT i=0 ; i<IrReps*IrReps ; i++ )
			if ( p[i] )
				delete p[i];
		delete p;
	}
	if ( nFrozen )
		delete nFrozen;
	if ( nDeleted )
		delete nDeleted;
	if ( inIrRepStored )
		delete inIrRepStored;
		
}


DensityMatrix::DensityMatrix(const DensityMatrix &d) :
	MOIrReps(d)
{
	p = new Type *[IrReps*IrReps];
	for ( INT i=0 ; i<IrReps*IrReps ; i++ )
		if ( d.p[i] )
		{
		INT n = inIrRep[i/IrReps]*inIrRep[i%IrReps];
			p[i] = new Type[n];
			for ( INT j=0 ; j<n ; j++ )
				p[i][j] = d.p[i][j];
		}
		else
			p[i] = NULL;

		
	nFrozen = new INT[IrReps];
	memcpy(nFrozen, d.nFrozen, IrReps*sizeof(INT));
	nDeleted = new INT[IrReps];
	memcpy(nDeleted, d.nDeleted, IrReps*sizeof(INT));
	inIrRepStored = new INT [IrReps];
	memcpy(inIrRepStored, d.inIrRepStored, IrReps*sizeof(INT));


	diagOcc = d.diagOcc;

	maxMO = 0;
	init();	
}


DensityMatrix & DensityMatrix::operator = (const DensityMatrix &d)
{
	((MOIrReps &) *this) = d;
	if ( p )
	{
		for ( INT i=0 ; i<IrReps*IrReps ; i++ )
			if ( p[i] )
				delete p[i];
		delete p;
	}

	if ( nFrozen )
		delete nFrozen;
	if ( nDeleted )
		delete nDeleted;
	if ( inIrRepStored )
		delete inIrRepStored;

	
	p = new Type *[IrReps*IrReps];
	for ( INT i=0 ; i<IrReps*IrReps ; i++ )
		if ( d.p[i] )
		{
		INT n = inIrRep[i/IrReps]*inIrRep[i%IrReps];
			p[i] = new Type[n];
			for ( INT j=0 ; j<n ; j++ )
				p[i][j] = d.p[i][j];
		}
		else
			p[i] = NULL;

		
	nFrozen = new INT[IrReps];
	memcpy(nFrozen, d.nFrozen, IrReps*sizeof(INT));
	nDeleted = new INT[IrReps];
	memcpy(nDeleted, d.nDeleted, IrReps*sizeof(INT));
	inIrRepStored = new INT [IrReps];
	memcpy(inIrRepStored, d.inIrRepStored, IrReps*sizeof(INT));


	diagOcc = d.diagOcc;

	maxMO = 0;
	init();	

	return *this;
}

void	DensityMatrix::calcNaturalOrbitals(
	Matrix<DensityMatrix::Type> & trafo,
	Vector<DensityMatrix::Type> & occupationNumbers) const
{
	trafo = Matrix<Type>(maxMO, maxMO);
	occupationNumbers = Vector<Type>(maxMO);

	for ( IrRep irrep=0 ; irrep<getNumberOfIrReps() ; irrep++ )
	{
	INT	n = getInIrRep(irrep);
		if ( !n )
			continue;
	
Type	*buf = new Type[(n+1)*n/2];
	INT	k = 0;
		for ( MOType i=0 ; i<n ; i++ )
			for ( MOType j=0 ; j<=i ; j++ )
				buf[k++] = (*this)(irrep, irrep, i+1, j+1);

	Type	*lambda = new Type[n];
	Type	*ev = new Type[n*n];

		hqrii1(&n, &n, buf, lambda, ev, &n);
	MOType	offset = getStartMO(irrep);
	
		
		k = 0;
		for ( MOType i=0 ; i<n ; i++ )
		{
			for ( MOType j=0 ; j<n ; j++ )
				trafo(offset+(n-i-1)-1, offset+j-1) = ev[k++];
			occupationNumbers[offset+(n-i-1)-1] = lambda[i];
		}

		delete ev;
		delete lambda;
		delete buf;
	}
}

Matrix<DensityMatrix::Type>	DensityMatrix::calcAtomic(MOTrafo & MOtrafo) const
{
Matrix<Type>	trafo(MOtrafo);
	trafo.transpose();

//	Inv(trafo);
//	*this = trafo * *this * Transpose(trafo);

Matrix<Type>	D(maxMO, maxMO);
	for  ( INT i=0 ; i<maxMO ; ++i )
		for ( INT j=0 ; j<maxMO ; ++j )
			D(i, j) = (*this)(i+1,j+1);
	
//	cout << "D=" << endl;
//	cout << D << endl;
	trafo.inverse();
//	cout << "Inverse:" << endl;
//	cout << trafo << endl;
	D = trafo*D;
	trafo.transpose();
	D = D*trafo;

	return D;
}



DensityMatrix & DensityMatrix::operator += (const DensityMatrix &d)
{
	if ( !p )
		*this = d;
	else
	{
		for ( INT i=0 ; i<IrReps*IrReps ; i++ )
			if ( d.p[i] )
			{
			INT n = inIrRepStored[i/IrReps]*inIrRepStored[i%IrReps];
				for ( INT j=0 ; j<n ; j++ )
				{
					p[i][j] += d.p[i][j];
				}
			}
		diagOcc += d.diagOcc;
	}

	return *this;
}


DensityMatrix operator * (DensityMatrix::Type a, const DensityMatrix &d)
{
DensityMatrix	mat(d);

	for ( INT i=0 ; i<d.IrReps*d.IrReps ; i++ )
		if ( d.p[i] )
		{
		INT n = d.inIrRepStored[i/d.IrReps]*d.inIrRepStored[i%d.IrReps];
			for ( INT j=0 ; j<n ; j++ )
				mat.p[i][j] = a*d.p[i][j];
		}

	mat.diagOcc = a*d.diagOcc;
	return mat;
}

 


ostream & DensityMatrix::writeToStream(ostream &os)
{
	MOIrReps::writeToStream(os);

	os << diagOcc << endl;
	for ( INT i=0 ; i<IrReps ; i++ )
		os << nFrozen[i] << " ";
	os << endl;

	for ( INT i=0 ; i<IrReps ; i++ )
		os << nDeleted[i] << " ";
	os << endl;

	for ( INT i=0 ; i<IrReps*IrReps ; i++ )
		os << (p[i] ? 1 : 0) << " ";
	os << endl;

	for ( INT i=0 ; i<IrReps*IrReps ; i++ )
		if ( p[i] )
		{
			for ( INT j=0 ; j<inIrRepStored[i/IrReps]*inIrRepStored[i%IrReps] ; j++ )
				os << setprecision(14) << p[i][j] << " ";
			os << endl;
		}

	return os;
}


ostream & operator << (ostream &o, const DensityMatrix &d)
{
	for ( INT i=0 ; i<d.maxMO ; i++ )
	{
		for ( INT j=0 ; j<d.maxMO ; j++ )
			if ( fabs(d(i+1, j+1))>0 )
			cout << setw(4) << i+1 << setw(4) << j+1 
				<< setprecision(8) << setw(16) << d(i+1,j+1) << endl;
		cout << endl;
	}
	cout << "===================" << endl;
	cout << "===================" << endl;
	cout << "===================" << endl;
	for ( INT i=0 ; i<d.IrReps*d.IrReps ; i++ )
		if ( d.p[i] )
		{
		INT r = i/d.IrReps;
		INT c = i%d.IrReps;
		
			for ( INT j=0 ; j<d.inIrRepStored[r] ; j++ )
				for ( INT k=0 ; k<d.inIrRepStored[c] ; k++ )
					cout << setw(4) << r << setw(4) << c
						<< setw(4) << j+1 << setw(4) << k+1
						<< setprecision(8) << setw(16) << d(r, c, j+1, k+1) << endl;
		}
/*	for ( INT i=0 ; i<d.maxMO ; i++ )
	{
		for ( INT j=0 ; j<d.maxMO ; j++ )
		{
			cout << setprecision(8) << setw(16) << d(i+1,j+1);
		}
		cout << endl;
	}
*/	return o;
}
