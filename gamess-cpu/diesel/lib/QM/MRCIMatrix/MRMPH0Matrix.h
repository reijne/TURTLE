//***********************************************************************
//
//	Name:			MRMPH0Matrix.h
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			16.09.1998
//
//
//
//
//
//***********************************************************************


#ifndef __MRMPH0Matrix_h
#define __MRMPH0Matrix_h

#include "../../../config.h"


#include "../Davidson/SLEMatrix.h"
#include "CICalculation.h"
#include "../MRTree/Diag/NExternalsDiag.h"
#include "../RepresentationMatrices/MRMPH0MatElements.h"

class NExternalsDiag;

class MOAccess;
template <class VectorType> class BufferedVector;


#include "../../Container/AVLSet.h"
#include "../Configuration/ConfigurationStart.h"


template <class MatrixType, class VectorType>
class MRMPH0Matrix :
	public CICalculation<MatrixType, VectorType>,
	public SLEMatrix<MatrixType, VectorType> {
public:

enum ProjectionMode { P_no0, P_0Complement };

	MRMPH0Matrix(
		const CICalculation<MatrixType, VectorType> &ciCalc,
		const NExternalsDiag & Psi0,
		const NExternalsDiag & imhomogenity,
		const BufferedVector<VectorType> &Psi0EV,
		const BufferedVector<VectorType> &inhomogenityEV,
		const MRFockMatrix<VectorType> *mrFockMatrix,
		ProjectionMode projectionMode = P_no0
		);
	~MRMPH0Matrix();


	const AVLSet<ConfigurationStart<MOType> > & getInternal(INT) const;


	void	printMat() const;

	//----------------------------------------------------------------	
	//	virtual functions from SLEMatrix

	// calculate resiual vector
	//          x  = x - b
	void	calcResidual(BufferedVector<VectorType> &x) const;


	// perform multiplication with invers of diagonal
	//                  1
	//          x  = -------  x 
	//           i      A      i
	//                   ii
	void	multInvDiag(BufferedVector<VectorType> &x) const;


	// perform multiplication on n vectors
	// y = A*x
	void	mult(const BufferedVector<VectorType> &x, BufferedVector<VectorType> &y) const;


	MatrixType	calcMP3(const BufferedVector<VectorType> &x) const;

	//----------------------------------------------------------------	
	//	virtual functions from IterationMatrix

	// get diagonal of matrix
	void	getDiagonal(MatrixType *) const {}


	// get total matrix (attention: probably VERY large)
	void	getTotalMatrix(MatrixType *) const {}
	void	getTotalMatrix(
		MatrixStorage<MatrixType, VectorType> &) const{}
        


private:
	INT	getConfNr(const Configuration<MOType> &) const;
	
	void	init(const MRFockMatrix<VectorType> *);
	void	buildInternalExcitations();
	void	createInteractionLists();
	void	calcE0();
	void	calcE1();
	void	calcDim();
	void	calcAlpha();
	void	clearPsi0(BufferedVector<VectorType> &) const;
	
	void	multFock(
		const BufferedVector<VectorType> &x, 
		BufferedVector<VectorType> &y) const;
	void	multP_no0(
		const BufferedVector<VectorType> &x, 
		BufferedVector<VectorType> &y) const;
	void	multP_0Complement(
		const BufferedVector<VectorType> &x, 
		BufferedVector<VectorType> &y) const;





const NExternalsDiag & Psi0;				// configuration tree to be excluded
											// with reference configurations

const NExternalsDiag & inhomogenity;		// configuration tree to be used as
											// inhomogenity


const BufferedVector<VectorType>	&Psi0EV;// Psi0 Eigenvectors

const BufferedVector<VectorType>	&inhomogenityEV;// inhomgenity Eigenvectors



VarSizeReadOnlyCache<TableKey, HMatElements<MatrixType, VectorType> >	*Hcache;
VarSizeReadOnlyCache<TableKey, MRMPH0MatElements<MatrixType, VectorType> >	*H0cache;

INT	totalConfs;								// number of total MR-CI configurations
IrRep	irrep;								// irrep of total wave function

static const INT	maxExc = 2;						// max excitation level
AVLSet<ConfigurationStart<MOType> > 		// internal configuration rests
	internal[maxExc+1];						// with start addresses

Pix	****interactionList;					// list of pixes, indices:
											//   1. number of holes (left side)
											//   2. index of configuration (left side)
											//   3. number of holes (right side)
											//   4. index of configuration (right side)
											// the lists are organized as an upper
											// triangular, with diagonal

IrRep	****interactionSymList;				// list of IrReps, structure as above

INT	***nInteractions;						// list of number of entries in
											// 4th dimension of interactionList

MOAccess	*moAccess[maxExc+1];			// index list from external tupel
											// to address with 1 and 2 externals


MatrixType	E0;								// E0 = < Psi0 | H0 | Psi0 >
MatrixType	E1;								// E1 = < Psi0 | H1 | Psi0 >


											//            n
											//             sel
											//           -----
BufferedVector<VectorType>	*alpha;			//           \      |i>  
											//	alpha  =  >    a    <Phi | F | Phi >
											//       i   /      j       j         i
											//           -----
											//            j=1

//BufferedVector<VectorType>	&xSpread;		// selected wave function coefficients
											// expanded to total MRCI space


ProjectionMode	projectionMode;
};




template <class MatrixType, class VectorType>
inline
const AVLSet<ConfigurationStart<MOType> > & MRMPH0Matrix<MatrixType, VectorType>::
	getInternal(INT i) const
{	return internal[i];	}



#endif
