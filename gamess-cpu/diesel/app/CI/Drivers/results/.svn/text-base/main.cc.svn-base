using namespace std;

#include "DieselResults.h"

#include "../../../../config.h"

#include "../../../../lib/Container/GrepAwk.h"

#include <iomanip>
#include <iostream>
#include <sstream>


int main(int argc, char **argv)
{
	argc--;

String* arg = new String[argc];
Number<double>	Thresh;
INT	Root = 0;
INT	Cr = -1;
INT	Cc = -1;
INT	noHeaders = 0;
INT	noWave = 0;
INT 	m3=3;
INT 	m0=0;
	for ( INT i=0 ; i<argc ; i++ )
		arg[i] = argv[i+1];
		
	for ( INT i=0 ; i<argc ; i++ )
	{
		if ( arg[i]=="-h" )
			noHeaders = 1;

		if ( arg[i]=="-w" )
			noWave = 1;

		if ( arg[i].at(m0, m3)=="-T=" )
		{
		double	h;
			sscanf(arg[i].chars(), "-T=%lf", &h);
			Thresh = h;
		}

		if ( arg[i].at(m0, m3)=="-R=" )
			sscanf(arg[i].chars(), "-R=%d", &Root);

		if ( arg[i].at(m0, m3)=="-C=" )
			sscanf(arg[i].chars(), "-C=%d,%d", &Cr, &Cc);

	}



String	MRPTThresh;

	if ( Cc==-1 )
	{
	DieselResults	dr(cout, Thresh*1e3, Root-1, noWave, noHeaders, MRPTThresh);
		cout << dr;
	}
	else
	{
	char	buf[1000000];
	stringstream	s(buf,stringstream::in | stringstream::out);

	DieselResults	dr(s, Thresh*1e3, Root-1, noWave, noHeaders, MRPTThresh);
	
		s << dr;
	GrepAwk	ga(s);
		ga.tail();
		ga -= 2;
		ga += Cr-1;
		
		cout << ga.getWord(Cc) << endl;
	
	
	}

delete[] arg;

}


