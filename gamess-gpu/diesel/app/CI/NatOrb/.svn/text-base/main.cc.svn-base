#include "../../../lib/QM/IntegralContainer/DensityMatrix.h"

#include "../../../lib/QM/IntegralContainer/MOTrafo.h"

#include <stdlib.h>
#include "../../../lib/Container/String.h"
#include <fstream>

#include <iomanip>

#include "../../../config.h"
#include "../../../VersionDate.h"

using namespace std;

int	main(INT argc, char **argv)
{
	cerr << "natorb (Part of DIESEL-MR-CI), " << VERSION << ", " << DATE << endl << endl;
	
	if ( argc<2  )
	{
		cerr << endl;
		cerr << "purpose:" << endl;
		cerr << "the program calculates natural orbitals based on density matrices" << endl; 
		cerr << endl;
		cerr << "usage: natorb DensityMatrix ... [-weight w1 ... ] <InputOrbitals >OutputOrbitals" << endl;
		cerr << endl;
		cerr << endl;
		exit(1);
	}

typedef double Type;

Type	*weight;
INT	N = 0;
	{
	INT i;
		for ( i=1 ; i<argc ; i++ )
		{
			if ( !strcmp(argv[i], "-weight") )
				break;
		}

		if ( i<argc )
		{
			N = i-1;
			if ( argc-i-1 != N )
			{
				cerr << "number of weights wrong." << endl;
				exit(1);
			}
			weight = new Type[N];
			for ( i++ ; i<argc ; i++ )
				weight[i-N-2] = atof(argv[i]);
		}
		else
		{
			N = argc-1;
			weight = new Type[N];
			for ( i=0 ; i<N ; i++ )
				weight[i] = 1;
		}
	}

	{
	Type	sum = 0;
		for ( INT i=0 ; i<N ; i++ )
			sum += weight[i];
		for ( INT i=0 ; i<N ; i++ )
			weight[i] /= sum;
	}


DensityMatrix   densityMatrix;

	for ( INT i=0 ; i<N ; i++ )
	{
	ifstream	fDens(argv[i+1]);
		if ( !fDens )
		{
			cerr << "no file \"" << argv[i+1] << "\"" << endl;
			exit(1);
		}
		DensityMatrix	dens(fDens);
//		cout << dens;
		densityMatrix += weight[i]*dens;
//		cout << densityMatrix;
	}

	delete weight;
	
Matrix<double>	mat;
Vector<double>	v;

	densityMatrix.calcNaturalOrbitals(mat, v);

/*	cout << densityMatrix << endl;
	cout << endl;
	cout << endl;
	cout << endl;
*/
//	cout << mat << endl;
//	cout << v << endl;


	MOTrafo	moTrafo(densityMatrix, cin);

/*	cout << moTrafo << endl;
	cout << endl;
	cout << endl;
	cout << endl;
*/
	moTrafo.transform(mat);
	moTrafo.setOccupationNumbers(v);
//	cout << moTrafo << endl;

	moTrafo.writeToStream(cout);
}
