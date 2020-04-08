//***********************************************************************
//
//	Name:			RepresentationMatrices.cc
//
//	Description:	stores and handles the 
//					representation matrices of the 
//					symmetric group approach (SGA)
//
//					
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			23.08.1996
//
//
//
//
//
//***********************************************************************

#include "RepresentationMatrices.h"

#include "../Configuration/TableCase.h"

#include "../Math/SpinEigenFunctionDegeneration.h"
#include "RepDiag.h"

#include "../../Math/etc/Histogram.h"


#include <stdlib.h>


template <class MatrixType>
RepresentationMatrixFortranInterface<MatrixType> *
	RepresentationMatrices<MatrixType>::repMatFInt = NULL;



template <class MatrixType>
RepresentationMatrices<MatrixType>::RepresentationMatrices()
{
	CbMat = NULL;
	ExP3Mat = NULL;
}

template <class MatrixType>
RepresentationMatrices<MatrixType>::RepresentationMatrices(const RepresentationMatrices &r)
{
//	printf("+++A %x\n", this);
//	fflush(stdout);
	occupiedMemory = r.getOccupiedMemory();
	matDim.rows = r.getNumberOfRows();
	matDim.cols = r.getNumberOfColumns();
	P = r.getP();
	P3Mats = r.getNumberOfP3Mats();
	CbMat = new RepresentationMatrix<MatrixType>(*r.getCbMat());
	ExP3Mat = NULL;
	if ( P==3 || P==1 )
		ExP3Mat = new RepresentationMatrix<MatrixType>(*r.getExP3Mat());
}




template <class MatrixType>
RepresentationMatrices<MatrixType>::RepresentationMatrices(TableKey key)
{
//	printf("+++A %x\n", this);
//	fflush(stdout);
//	printf("Im Constructor!!!\n");
//	fflush(stdout);

TableCase<MOType>	tc = key.getTableCase();

	CbMat = NULL;
	ExP3Mat = NULL;
	occupiedMemory = 0;


//	cout << key << endl;	
//	cout << key.getKey() << endl;

//	cout << "table case: " << tc << endl;
//	printf("%x\n", this);
//	fflush(stdout);

	P = tc.getP();
//	cout << "P=" << P << endl;
//	cout << "table case: " << tc << endl;
//		return;

	if ( tc.getdK()>=0 )
	{//	cout << "dK>=0: " << tc.getNumberOfMoreOpenShells() - 2*tc.getdK() << " " <<
	//		tc.getNumberOfMoreOpenShells() << endl;
	
		matDim.rows = repMatFInt->getSpinEigs(tc.getNumberOfMoreOpenShells() - 2*tc.getdK());
		matDim.cols = repMatFInt->getSpinEigs(tc.getNumberOfMoreOpenShells());
	}
	else
	{//	cout << "dK<0: " << tc.getNumberOfOpenShells() << " " <<
	//		tc.getNumberOfMoreOpenShells() + 2*tc.getdK() << endl;

		matDim.rows = repMatFInt->getSpinEigs(tc.getNumberOfMoreOpenShells());
		matDim.cols = repMatFInt->getSpinEigs(tc.getNumberOfMoreOpenShells() + 2*tc.getdK());
	}
//	cout << "cols: " << matDim.cols << ", rows: " << matDim.rows << endl;
//	fflush(stdout);
//	return;
	
INT	NumberOfElements = 0;

	switch ( P ) {
	case 1:		// double excitation: coulomb and exchange RepMat
		NumberOfElements = 2 * matDim.cols * matDim.rows;
		break;

	case 2:		// double excitation: only coulomb RepMat
	case 4:
		NumberOfElements = matDim.cols * matDim.rows;
		break;
		
	case 3:
		NumberOfElements = matDim.cols * matDim.rows;
		
		//	number of P=3 matrices depends on dK:
		//	abs(dK)=0: i!=a		==> more open shells - 1  +  1 UMat
		//	abs(dK)=1: i!=a, c	==> more open shells - 2  +  1 UMat
		P3Mats = tc.getNumberOfMoreOpenShells()-abs(tc.getdK()) - 1;
		NumberOfElements *= 1 + P3Mats;
		break;
		
	case 5:
		P3Mats = tc.getNumberOfMoreOpenShells()*(tc.getNumberOfMoreOpenShells()-1)/2;
		NumberOfElements = P3Mats * matDim.cols * (matDim.cols+1)/2;
		break;
	}
	

//	if ( tc.getP()!=3 )
	{
MatrixType	*buf = new MatrixType[NumberOfElements];
	repMatFInt->calc(
			buf, buf + matDim.rows*matDim.cols, 
			tc.getdK(), tc.getP(), tc.getR(), tc.getqR(), tc.getqL(),
			matDim.rows, matDim.cols,
			tc.getNumberOfMoreOpenShells());

		calcDense(buf);
		
		delete[] buf;
//		cout << "Hallo 3" << endl;
		
			
//		if ( tc.getP()!=5 )
//			cout << *this << endl;
	}

	if ( CbMat )
		occupiedMemory += CbMat->getOccupiedMemory();
	if ( ExP3Mat )
		occupiedMemory += ExP3Mat->getOccupiedMemory();

	occupiedMemory += sizeof(RepresentationMatrices);

//	printf("Ende Constructor!!!\n");
//	fflush(stdout);
}


template <class MatrixType>
RepresentationMatrices<MatrixType>::~RepresentationMatrices()
{
//	printf("+++D %x\n", this);
//	fflush(stdout);
	if ( CbMat )
		delete CbMat;
	if ( ExP3Mat )
		delete ExP3Mat;
//	printf("~~~\n", this);
//	fflush(stdout);
}



template <class MatrixType>
void	RepresentationMatrices<MatrixType>::calcDense(MatrixType *buf)
{
//	cout << "Hallo A" << endl;
//	cout << "RepresentationMatrices this = " << this << endl;
//	cout << "cols: " << matDim.cols << ", rows: " << matDim.rows << endl;
	switch ( P ) {
	case 1:
		CbMat = new RepresentationMatrix<MatrixType>(buf, matDim.rows, matDim.cols);
		ExP3Mat = new RepresentationMatrix<MatrixType>(
			buf + matDim.rows*matDim.cols, matDim.rows, matDim.cols);
			
			
		if ( CbMat->getDense()==RepresentationMatrix<MatrixType>::sparse &&
			ExP3Mat->getDense()==RepresentationMatrix<MatrixType>::sparse &&
			CbMat->getNSparse()==ExP3Mat->getNSparse() )
		{
		INT	same = 1;
			for ( INT i=0 ; i<CbMat->getNSparse() ; i++ )
				if ( CbMat->getMPosP()[i].r != ExP3Mat->getMPosP()[i].r ||
					CbMat->getMPosP()[i].c != ExP3Mat->getMPosP()[i].c )
				{	same = 0;
					break;
				}
			if ( same )
			{	CbMat->setDense(RepresentationMatrix<MatrixType>::sameSparse);
				ExP3Mat->setDense(RepresentationMatrix<MatrixType>::sameSparse);
			}
		}
//		CbDense = checkType(buf, Elements); 
		
		break;
		
	case 2:
	case 4:
		CbMat = new RepresentationMatrix<MatrixType>(buf, matDim.rows, matDim.cols);
		break;
		
	case 3:
//		cout << "P3Mats = " << P3Mats << endl;
		CbMat = new RepresentationMatrix<MatrixType>(buf, matDim.rows, matDim.cols, 0);
//		cout << "P=3 CbMat=" << CbMat << endl;
		ExP3Mat = new RepresentationMatrix<MatrixType>(P3Mats, 
			buf + matDim.rows*matDim.cols, matDim.rows, matDim.cols,
			RepresentationMatrix<MatrixType>::full);
		break;

	case 5:
		CbMat = new RepresentationMatrix<MatrixType>(P3Mats, 
			buf, matDim.rows, matDim.cols, RepresentationMatrix<MatrixType>::lowerTriangle);
		break;
	}
//	cout << "Hallo B" << endl;
}


template <class MatrixType>
void	RepresentationMatrices<MatrixType>::printMats(ostream & s) const
{

	s << *CbMat << endl;

	switch ( P ) {
	case 1:
		s << *ExP3Mat << endl;
		break;
	
	case 3:
		s << endl;
		const MatrixType	*p = ExP3Mat->getP();
		for ( INT i=0 ; i<matDim.rows ; i++ )
		{	for ( INT k=0 ; k<P3Mats ; k++ )
			{	s << "(";
				for ( INT j=0 ; j<matDim.cols ; j++ )
				{	s << p[i*matDim.cols*P3Mats + j*P3Mats + k];
					if ( j<matDim.cols-1 )
						s << " ";
				}
				s << ")\t";
			}
			s << endl;
		}
//						*p++ = pEMat[k*kml*kmj + i*kmj + j];
		break;
	
	}
}



//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------


template class RepresentationMatrices<float>;
template class RepresentationMatrices<double>;
