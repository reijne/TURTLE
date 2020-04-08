#include "../../../config.h"

using namespace std;

#include "QCProgOutputFormat.h"

#include <fstream>
#include <iostream>

#include "InternalCoord.h"

#include "Molecule.h"


int	main()
{

ifstream	f("/tmp1/hanrath/Inputs/intcoord");
QCProgOutputFormat	qcFormat(f);

	cout << qcFormat.getIdentifier() << endl;

	f.seekg(0);
	f.clear();
Molecule	molecule(f, qcFormat.getPackage());
	cout << molecule << endl;



}
