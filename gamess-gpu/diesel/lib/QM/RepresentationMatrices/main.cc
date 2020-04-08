#include "../../../config.h"

using namespace std;

#include <iostream>

#include "RepresentationMatrices.h"
#include "RepresentationMatrix.h"
#include "../Configuration/Configuration.h"
#include "../Configuration/TableCase.h"
#include "../Configuration/DiffConf.h"


int	main()
{
RepresentationMatrix<double>	r();



BinomialCoefficient	binom(MAXOPENSHELLS);
RepresentationMatrixFortranInterface<double>	repMatFInt(2);

TableCase<MOType>	tc(&binom);

/*	{
		tc.setNr(5, 25, 2, 2);

	TableKey	tk(tc);
	RepresentationMatrices	rmc(tk);


		cout << tc << endl << rmc << endl;
	}
*/
	{
	Configuration<INT>	conf1;
		conf1.create(1);
//		conf1.create(2);
//		conf1.create(2);
		conf1.create(3);
		conf1.create(5);
	Configuration<INT>	conf2;
		conf2.create(1);
		conf2.create(1);
//		conf2.create(2);
//		conf2.create(2);
		conf2.create(4);

		cout << conf1 << endl;
		cout << conf2 << endl;
		tc.calc(conf1, conf2);
		cout << tc << endl;

	TableKey	tk(tc);
	RepresentationMatrices<double>	rmc(tk);


		cout << rmc << endl;
		cout << "--------------------------------" << endl << endl;
	}
	
/*	for ( INT i=1 ; i<7 ; ++i )
	for ( INT j=1 ; j<7 ; ++j )
	{
	Configuration<INT>	conf1;
		conf1.create(1);
		conf1.create(2);
		conf1.create(2);
		conf1.create(5);
		if ( conf1.create(i) )
			continue;
	Configuration<INT>	conf2;
		conf2.create(1);
		conf2.create(1);
		conf2.create(2);
		conf2.create(4);
		if ( conf2.create(j) )
			continue;
		cout << conf1 << endl;
		cout << conf2 << endl;
		tc.calc(conf1, conf2);
		cout << tc << endl;

	TableKey	tk(tc);
	RepresentationMatrices<double>	rmc(tk);


		cout << rmc << endl;
		cout << "--------------------------------" << endl << endl;
	}
*/	

}
