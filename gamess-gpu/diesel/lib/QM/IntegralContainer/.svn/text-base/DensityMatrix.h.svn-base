#include "../../../config.h"
//***********************************************************************
//
//	Name:			DensityMatrix.h
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


#ifndef __DensityMatrix_h
#define __DensityMatrix_h

#include "../MO/MOIrReps.h"


#include "../../Math/MatrixVector/Matrix.h"


/* FD class istream; */
/* FD class ostream; */
class  MOTrafo;

class DensityMatrix : public MOIrReps {
public:
typedef double Type;

	DensityMatrix();
	DensityMatrix(istream &);
	DensityMatrix(const MOIrReps &, const INT *nFrozen, const INT *nDeleted,
		const Type *, INT transition);
	~DensityMatrix();
	
	DensityMatrix(const DensityMatrix &);
	DensityMatrix & operator = (const DensityMatrix &);


	Type operator () (INT IrRepi, INT IrRepj, INT i, INT j) const;
	Type operator () (MOType i, MOType j) const;

	
	void	calcNaturalOrbitals(
		Matrix<DensityMatrix::Type> & trafo,
		Vector<DensityMatrix::Type> & OccupationNumbers) const;
	
	Matrix<DensityMatrix::Type>	calcAtomic(MOTrafo & trafo) const;


	DensityMatrix & operator += (const DensityMatrix &);
	
	friend	DensityMatrix operator * (DensityMatrix::Type, const DensityMatrix &);


	
	ostream & writeToStream(ostream &);
	friend ostream & operator << (ostream &, const DensityMatrix &);

private:


INT	*nFrozen;				// number of frozen orbitals in irrep
INT	*nDeleted;				// number of deleted orbitals in irrep
INT	*inIrRepStored;			// number of orbitals without frozen and deleted ones

Type	**p;				// pointer to density matrix blocks
Type	diagOcc;			// if not transition matrix =2, =0 otherwise
};



inline
DensityMatrix::Type DensityMatrix::operator () (
	INT IrRepi, INT IrRepj, INT i, INT j) const
{
//	cout << "----  " << IrRepi << " " << i << " " << IrRepj << " " << j << endl;
	if ( (i<=nFrozen[IrRepi] || j<=nFrozen[IrRepj]) )
	{
		if ( IrRepi==IrRepj && i==j )
			return diagOcc;
		else
			return 0;
	}

	if ( (i>inIrRepStored[IrRepi]+nFrozen[IrRepi] || 
		j>inIrRepStored[IrRepj] + nFrozen[IrRepj]) )
		return 0;

//	cout << "::::::::: " << IrRepi << " " << IrReps << " " << IrRepj << endl;
//	cout << "::::::::: " << IrRepi*IrReps + IrRepj << endl;
//	cout << "::::::::: " << p[IrRepi*IrReps + IrRepj] << endl;
//	cout << ":::::::::A " << inIrRepStored[IrRepi] << " " <<
//		(i-nFrozen[IrRepi]-1) << " " << (j-nFrozen[IrRepj]-1) << endl;
	if ( p[IrRepi*IrReps + IrRepj] )
		return p[IrRepi*IrReps + IrRepj][inIrRepStored[IrRepj]*(i-nFrozen[IrRepi]-1) +
			(j-nFrozen[IrRepj]-1)];
	else
		return 0;
}


inline
DensityMatrix::Type DensityMatrix::operator () (
	MOType i, MOType j) const
{
//	cout << "+++  " << i << " " << j << endl;
	return operator () (
		getIrRep(i), getIrRep(j), 
		getNumberInIrRep(i)+1, getNumberInIrRep(j)+1);
	
}



#endif
