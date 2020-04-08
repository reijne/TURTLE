using namespace std;

#include "../../../../lib/QM/IntegralContainer/DensityMatrix.h"

#include <unistd.h>
#include <fstream>
#include "../../../../lib/Container/String.h"
#include <iomanip>

#include "../../../../config.h"


int	main(INT argc, char **argv)
{
ifstream	fDens(argv[1]);

	if ( !fDens )
	{
		cerr << "no file \"" << argv[1] << "\"" << endl;
		exit(1);
	}
DensityMatrix	densityMatrix(fDens);


	cout << "densityMatrix" << endl;
	cout << densityMatrix << endl;



	return 0;
}
