#include "../../../../VersionDate.h"
#include "../../../../lib/QM/IntegralContainer/PropInts.h"
#include "../../../../lib/QM/IntegralContainer/MOTrafo.h"
#include "../../../../lib/QM/IntegralContainer/DensityMatrix.h"

#include <unistd.h>

#include <fstream>

#include <string>

#include "../../../../lib/Container/String.h"
#include <iomanip>

#include "../../../../config.h"

using namespace std;

int	main(int argc, char **argv)
{
	cerr << "prop (Part of DIESEL-MR-CI), " << VERSION << ", " << DATE << endl << endl;
	if ( argc<5 )
	{
		cerr << "purpose:" << endl;
		cerr << "calculate properties from density matrices." << endl << endl;
		cerr << "usage:" << endl;
		cerr << "prop Operator OrbitalPath IntegralPath DensityPath ..." << endl;
		cerr << "       Operator     = { MLTPL1, KINETIC, ONEHAM, ANGMOM }" << endl;
		cerr << "       OrbitalPath  = path to MOLCAS orbitals file (e.g. \"SCFORB\", \"RASORB\")" << endl;
		cerr << "       IntegralPath = path to MOLCAS ONEINT File" << endl;
		cerr << "       DensityPath  = path to CI density matrix(ces) (generated with \"dens\")" << endl;
		cerr << endl;
		return 1;
	}

INT	nOps = 4;
char	*OpNames[] = { "MLTPL1", "KINETIC", "ONEHAM", "ANGMOM" };
bool	anti[] = { false, false, false, true };

PropInts::Operator op = (PropInts::Operator) -1;
	{
	
	String	s(argv[1]);
//	String	s("ANGMOM");
//	String	s("MLTPL1");
		s.upcase();
		for ( INT i=0 ; i<nOps ; i++ )
		{
		String	s1(OpNames[i]);

			if ( s==s1 )
			{
				op = (PropInts::Operator) i;
				break;
			}

		}
		if ( op==-1 )
		{
			cout << "unknown operator " << argv[1] << endl;
			return 1;
		}
	}
	
	if ( strcmp(argv[3], "ONEINT") )
		if ( symlink(argv[3], "ONEINT") )
		{
//			cerr << "error linking \"" << argv[3] << "\" to \"ONEINT\"." << endl;
//			exit(1);
		}



PropInts	*propInts[3];

	for ( INT axes=0 ; axes<3 ; axes++ )
	{
	ifstream	fOrbs(argv[2]);
		if ( !fOrbs )
		{
			cerr << "no file \"" << argv[2] << "\"" << endl;
			exit(1);
		}
	
//	PropInts	propInts(op, atoi(argv[3]), "ONEINT");
		propInts[axes] = new PropInts(op, axes+1, "ONEINT");
//		MOTrafo	moTrafo(*propInts[axes], fOrbs);

//		cout << "PropInts(untransformed)" << endl;
//		cout << *propInts[axes] << endl;


//		cout << "moTrafo" << endl;
//		cout << moTrafo << endl;

//		propInts[axes]->transform(moTrafo);

//		cout << "PropInts(transformed)" << endl;
//		cout << *propInts[axes] << endl;
	}
	

	for ( INT i=0 ; i<argc-4 ; i++ )
	{
	ifstream	fDens(argv[i+4]);

		if ( !fDens )
		{
			cerr << "no file \"" << argv[i+4] << "\"" << endl;
			exit(1);
		}
//	cout << "   OK" << endl;
	DensityMatrix	densityMatrix(fDens);


//		cout << "densityMatrix" << endl;
//		cout << densityMatrix << endl;
	INT	i1, i2, r1, r2;
	double	dummy;
        double oscstr[3];
		sscanf(strstr(argv[i+4], "Density.dat.")+strlen("Density.dat."),
			 "I%dR%d_I%dR%d.%lf\n", 
			&i1, &r1, &i2, &r2, &dummy);
//		cout << i1 << " " << r1 << " " << i2 << " " << r2 << "  ";

		for ( INT axes=0 ; axes<3 ; axes++ )
		{
		double	v = propInts[axes]->multDensity(densityMatrix, anti[op]);

			cout.precision(12);
			cout << setw(20) << v;
                        oscstr[axes]=v;
		}
		cout << endl;
		cout << "To get osc. strengths, multiply with DeltaE!!" <<endl;
		for ( INT axes=0 ; axes<3 ; axes++ )
		{
			oscstr[axes]=2.0/3.0*oscstr[axes]*oscstr[axes];
			cout.precision(12);
			cout << setw(20) << oscstr[axes];
		}
		cout << endl;
	}


	if ( strcmp(argv[3], "ONEINT") )
		unlink("ONEINT");

	return 0;
}
