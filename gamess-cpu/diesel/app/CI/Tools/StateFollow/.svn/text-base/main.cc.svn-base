using namespace std;

#include "../../../../lib/QM/IntegralContainer/PropInts.h"
#include "../../../../lib/QM/IntegralContainer/MOTrafo.h"
#include "../../../../lib/QM/IntegralContainer/DensityMatrix.h"

#include <unistd.h>

#include <fstream>

#include <string>
#include <dirent.h>

#include <iomanip>

#include "../../../../config.h"
#include "../../../../VersionDate.h"

#include <map>

struct IR {
	IR(INT irrep, INT root) :
		irrep(irrep), root(root) {}

	INT irrep;
	INT root;

	bool operator < (IR const & a) const
	{
		if ( irrep<a.irrep )
			return true;
		if ( irrep>a.irrep )
			return false;
		if ( root<a.root )
			return true;
		return false;
	}
};




int	main(INT argc, char **argv)
{
	cerr << "stateFollow (Part of DIESEL-MR-CI), " << VERSION << ", " << DATE << endl << endl;
	if ( argc<4 )
	{
		cerr << "purpose:" << endl;
		cerr << "follow states by AO transformed density matrices." << endl << endl;
		cerr << "usage:" << endl;
		cerr << "stateFollow Threshold CIPath Dir1 ..." << endl;
		cerr << "       OrbitalPath  = path to MOLCAS orbitals file (e.g. \"SCFORB\", \"RASORB\")" << endl;
		cerr << "       IntegralPath = path to MOLCAS ONEINT File" << endl;
		cerr << "       DensityPath  = path to CI density matrix(ces) (generated with \"dens\")" << endl;
		cerr << endl;
		return 1;
	}


typedef  map<IR, Matrix<double> > IMap;
map<string, IMap >	data;


double	Thresh = atof(argv[1]);
	
	for ( INT i=3 ; i<argc ; i++ )
	{
	string	path = argv[i];

		cout << "scanning directory \"" << path << "\"" << endl;
	
	string orbs = path+"/INPORB";
	ifstream	fOrbs(orbs.c_str());
		if ( !fOrbs )
		{
			cerr << "no file \"" << orbs << "\"" << endl;
			exit(1);
		}

		unlink("ONEINT");
	string oneInt = path+"/ONEINT";
		symlink(oneInt.c_str(), "ONEINT");
	PropInts	propInts((PropInts::Operator) 0, 1, "ONEINT");
		unlink("ONEINT");

	MOTrafo	moTrafo(propInts, fOrbs);


	//		cout << "moTrafo" << endl;
	//		cout << moTrafo << endl;


	string ciPath = path+"/"+argv[2]+"/Densities";
	IMap	_data;

	DIR	*dir = opendir(ciPath.c_str());
		if ( !dir )
		{
			cout << "no such directory \"" << ciPath << "\"" << endl;
			exit(1);
		}
	struct dirent	*entry;

		while ( (entry=readdir(dir)) )
		{
			if ( !strncmp("Density.dat.", entry->d_name, 12) )
			{
			INT	i1, i2, r1, r2;
			double	T;
				sscanf(strstr(entry->d_name, "Density.dat.")+strlen("Density.dat."),
					 "I%dR%d_I%dR%d.%lf\n", 
					&i1, &r1, &i2, &r2, &T);
				if ( T==Thresh && i1==i2 && r1==r2 )
				{
//					cout << entry->d_name << endl;
	
					string	f = ciPath + "/" + entry->d_name;
					ifstream	fDens(f.c_str());
	
					if ( !fDens )
					{
						cerr << "no file \"" << f << "\"" << endl;
						exit(1);
					}
//					cout << "   OK" << endl;
					DensityMatrix	densityMatrix(fDens);


//					cout << "densityMatrix" << endl;
//					cout << densityMatrix << endl;


					_data[IR(i1, r1)] = Matrix<double>(densityMatrix.calcAtomic(moTrafo));
//					cout << D << endl;
				}

			}
		}
		closedir(dir);
	
		data[path] = _data;



	}

	for ( INT i=4 ; i<argc ; i++ )
	{
	string	from = argv[i-1];
	string	to = argv[i];
	
		cout << endl << endl << endl << from << " --> " << to << endl << endl;
	map<IR, map<IR, double> >	sim;
	INT	oirrep = 0;
		for ( IMap::const_iterator j=data[from].begin() ; j!=data[from].end() ; ++j )
		{
			if ( j->first.irrep!=oirrep )
				cout << endl;
			oirrep = j->first.irrep;
			for ( IMap::const_iterator k=data[to].begin() ; k!=data[to].end() ; ++k )
			{
				if ( j->first.irrep!=k->first.irrep )
					continue;
			Matrix<double>	M = k->second;
				M -= j->second;
				for ( INT _i=0 ; _i<M.getRows() ; ++_i )
					for ( INT _j=0 ; _j<M.getCols() ; ++_j )
						if ( _i!=_j )
							M(_i ,_j) = 0;
//				cout << j->second << endl << endl;
//				cout << k->second << endl << endl;
//				cout << M << endl << endl;
				cout << setw(10) << M.getNorm2() << "\t";
//				exit(0);
				sim[j->first][k->first] = M.getNorm2();
				
			}
			cout << endl;
		}
			
	
	}


	return 0;
}
