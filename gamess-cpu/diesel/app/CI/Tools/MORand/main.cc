#include "../../../../lib/QM/IntegralContainer/DensityMatrix.h"

#include "../../../../lib/QM/IntegralContainer/MOTrafo.h"

#include "../../../../lib/Math/FortranNumeric/FortranEigenProblems.h"

#include <stdlib.h>
#include "../../../../lib/Container/String.h"
#include <fstream>
#include <iomanip>
#include <vector>

#include "../../../../config.h"
#include "../../../../VersionDate.h"

using namespace std;

int	main(INT argc, char **argv)
{
	cerr << "morand (Part of DIESEL-MR-CI), " << VERSION << ", " << DATE << endl << endl;
	
	if ( argc!=2  )
	{
		cerr << endl;
		cerr << "purpose:" << endl;
		cerr << "MO random transformation" << endl; 
		cerr << endl;
		cerr << "usage: natorb Occupationthreshold <InputOrbitals >OutputOrbitals" << endl;
		cerr << endl;
		cerr << endl;
		exit(1);
	}

typedef double Type;

Type	T = atof(argv[1]);

	
//		cout << T << endl;

MOTrafo	moTrafo(cin);
const INT n = moTrafo.getMaxMO();

Matrix<double>	mat(n, n);
	mat.setOne();
Vector<double>	v = moTrafo.getOccupationNumbers();

time_t	_time;
	srand(time(&_time));

	for ( INT irrep=0 ; irrep<moTrafo.getNumberOfIrReps() ; ++irrep )
	{
	vector<INT>	a;

		for ( INT i=moTrafo.getStartMO(irrep) ; i<=moTrafo.getEndMO(irrep) ; i++ )
			if ( v[i-1] < T )
				a.push_back(i-1);
			else
				cerr << i << endl;

	INT m = a.size();
	Matrix<double>	anti(m, m);
		anti.setFill(0);
	

		for ( INT i=0 ; i<m ; ++i )
			for ( INT j=0 ; j<i ; ++j )
			{
			double a = 2*M_PI*static_cast<double>(rand())/RAND_MAX;
				anti(i, j) = a;
				anti(j, i) = -a;
			}

		{
		INT ideg = 6;
		double t = 1;
		INT lwsp = 4*m*m+ideg+1;
		double	*wsp = new double[lwsp];
		
		INT iexph;
		INT	ns;
		INT iflag;
		INT ipiv[m];
		INT ldh = m;
		
//			cout << anti << endl;
		
			FORTRAN_LINKAGE(dgpadm)(
				&ideg, &m, &t, anti.getP(), &ldh, wsp, &lwsp, ipiv, &iexph, &ns, &iflag);

			
		Matrix<double>	T(m, m);
			for ( INT i=0 ; i<m ; ++i )
				for ( INT j=0 ; j<m ; ++j )
					T(i, j) = wsp[iexph+i*m+j-1];
			delete wsp;
					
		Matrix<double>	T1(T);
			T1.transpose();
			T1 *= T;
//			cout << T << endl;
//			cout << T1 << endl;

			for ( INT i=0 ; i<m ; ++i )
				for ( INT j=0 ; j<m ; ++j )
					mat(a[i], a[j]) = T(i, j);
				

			
		}
	}

/*	cout << densityMatrix << endl;
	cout << endl;
	cout << endl;
	cout << endl;
*/
//cout << "::::::::::" << endl;
//	cout << mat << endl;
//	cout << v << endl;

	Matrix<double>	m(mat);
	m.transpose();
	m*= mat;
//	cout << "E=" << endl;
//	cout << m << endl;

/*	cout << moTrafo << endl;
	cout << endl;
	cout << endl;
	cout << endl;
*/
	moTrafo.transform(mat);
//	cout << moTrafo << endl;

	moTrafo.writeToStream(cout);
}
