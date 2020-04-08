#include "../../../../config.h"

using namespace std;

#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>

#include "../../../../lib/QM/Configuration/Configuration.h"

#include "../../../../lib/Math/MatrixVector/Vector.h"
#include "../../../../lib/Container/DiskBuffer.h"

#include "../../../../lib/Container/GrepAwk.h"

#include <string>

#include <vector>

#include "../../../../VersionDate.h"

typedef vector<Configuration<MOType> > Confs;

Confs read(const char *fname)
{
ifstream	f(fname);
GrepAwk	ga(f);

Confs	A;

	do {
	const INT max = 10000;
	char	buf[max];
	stringstream	sbuf(buf, stringstream::in | stringstream::out);
	string	s = ga.getLine().chars();
		for ( INT i=1 ; i<=ga.getNumberOfWords() ; i++ )
		{
			if ( ga.getWord(i)=="#" )
				break;
			sbuf << ga.getWord(i) << " ";
		}
		sbuf << endl;
	Configuration<MOType>	conf(sbuf);

		A.push_back(conf);			
		ga++;
	} while ( !ga.illegal() && ga.getNumberOfWords()>2 );
	
	return A;
}







int	main(INT argc, char **argv)
{
	cout << "diffConf (Part of DIESEL-MR-CI), " << VERSION << ", " << DATE 
		<< endl << endl;

	if ( argc!=3 )
	{
		cerr << endl;
		cerr << "purpose:" << endl;
		cerr << "the program calculates excitations of to configuration sets." << endl; 
		cerr << endl;
		cerr << "usage: confStat confsA confsB" << endl;
		cerr << endl;
		exit(1);
	}




vector<Configuration<MOType> >	A = read(argv[1]);
vector<Configuration<MOType> >	B = read(argv[2]);


	for ( unsigned INT i=0 ; i<A.size() ; i++ )
		for ( unsigned INT j=0 ; j<B.size() ; j++ )
		{
			cout << i << " " << j << " " << 
				Configuration<MOType>::calcExcitationOrder(A[i], B[j]) << endl;
		}
	return 0;
}
