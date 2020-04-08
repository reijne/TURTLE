#include <iostream>
#include <sstream>
#include <fstream>

using namespace std;

#include "../../../../lib/QM/MRTree/Set/NExternalsSet.h"
#include "../../../../lib/QM/MO/MOMapping.h"
#include "../../../../lib/QM/Configuration/ConfigurationSet.h"
#include "../../Selector/MRConfInput.h"

#include "../../../../lib/Container/GrepAwk.h"

#include <ctype.h>

#include "../../../../config.h"
#include "../../../../VersionDate.h"

int	main(INT argc, char **argv)
{
enum	Mode { ModeUndefined, merge, intersect, difference } mode = ModeUndefined;
enum	Format { FormatUndefined, plain, tree, selector, selectorPlain } format = FormatUndefined;

	cerr << "setops (Part of DIESEL-MR-CI), " << VERSION << ", " << DATE 
		<< endl << endl;;
	
	if ( argc<3 )
	{
		cerr << endl;
		cerr << "purpose:" << endl;
		cerr << "the program calculates (m)erge, (i)ntersection or (d)ifference sets" << endl; 
		cerr << "of the given configurations and writes them to the standard output." << endl;
		cerr << endl;
		cerr << "usage: setops -{c|s|o|t} -{m|i|d} file1 file2 ..." << endl;
		cerr << endl;
		cerr << "options:" << endl;
		cerr << endl;
		cerr << "  -c: configurations in plain format" << endl;
		cerr << "  -s: configurations to be read from selector input (with interpretation)" << endl;
		cerr << "  -o: configurations to be read from selector input (without interpretation)" << endl;
		cerr << "  -t: configurations in tree format" << endl;
		cerr << endl;
		cerr << "  -m: perform merge" << endl;
		cerr << "  -i: perform intersection" << endl;
		cerr << "  -d: perform difference" << endl;
		cerr << endl;
		exit(1);
	}
	
	for ( INT i=1 ; i<argc ; i++ )
	{
		if ( argv[i][0]=='-' )
		{
			argv[i][0] = 0;
			switch ( argv[i][1] ) {
			case 'm':
				mode = merge;
				break;
				
			case 'i':
				mode = intersect;
				break;
			
			case 'd':
				mode = difference;
				break;
				
			case 'c':
				format = plain;
				break;
				
			case 't':
				format = tree;
				break;

			case 's':
				format = selector;
				break;
				
			case 'o':
				format = selectorPlain;
				break;
				
			default:
				cerr << "unknown option " << argv[i] << endl;
				exit(1);
			}
		}
	}
	
	if ( mode == difference || mode == intersect )
	{
		cerr << "sorry, not yet implemented." << endl;
		exit(1);
	}
	
	if ( mode == ModeUndefined )
	{
		cerr << "operation mode missing." << endl;
		exit(1);
	}
	if ( format == FormatUndefined )
	{
		cerr << "format missing." << endl;
		exit(1);
	}
	
	switch ( format ) {
	case selector:
	{	
	ConfigurationSet	merged;
	INT	first = 1;
	MRConfInput	*mrinp = NULL;
		
	stringstream	o;
//FD	ostream &	coutOrg = cout;
//FD		cout = o;

		for ( INT i=1 ; i<argc ; i++ )
		{
			if ( argv[i][0]==0 )
				continue;


			cerr << "reading " << argv[i] << "..." << flush;

	ifstream	in(argv[i]);
	ConfigurationSet	*toMerge;

		MRConfInput	t;
		
			in >> t;
			if ( first )
				mrinp = new MRConfInput(t);
				
			toMerge = new ConfigurationSet;
			Pix j = t.getRefConfSet().first();
			while ( j )
			{
				toMerge->add(t.getRefConfSet()(j));
				t.getRefConfSet().next(j);
			}

			cerr << " o.k." << endl;

			switch ( mode ) {
			case merge:
				merged |= *toMerge;
				break;

			case intersect:
				merged &= *toMerge;
				break;

			case difference:
				merged -= *toMerge;
				break;

			default:
				break;
			}
			delete toMerge;
			first = 0;
		}
		cerr << "writing selector input file with merged reference configurations" << endl;
		cerr << "to standard output..." << flush;
//FD		cout = coutOrg;
		if ( mrinp )
		{
			mrinp->getRefConfSet() = merged;
			mrinp->getPTRefConfSet() = merged;
			cout << *mrinp << endl;
		}
		cerr << " o.k." << endl;
		break;
	}	

	case plain:
	{	
	ConfigurationSet	merged;
		
		for ( INT i=1 ; i<argc ; i++ )
		{
			if ( argv[i][0]==0 )
				continue;


			cerr << "reading " << argv[i] << "..." << flush;

	ifstream	in(argv[i]);
	ConfigurationSet	*toMerge;

			toMerge = new ConfigurationSet(in);

			cerr << " o.k." << endl;

			switch ( mode ) {
			case merge:
				merged |= *toMerge;
				break;

			case intersect:
				merged &= *toMerge;
				break;

			case difference:
				merged -= *toMerge;
				break;

			default:
				break;
			}
			delete toMerge;
		}
		cerr << "writing merged configurations to standard output..." << flush;
		merged.writeToStream(cout);
		cerr << " o.k." << endl;
		break;
	}	

	case selectorPlain:
	{	
	ConfigurationSet	merged;
	INT	first = 1;		
		for ( INT i=1 ; i<argc ; i++ )
		{
			if ( argv[i][0]==0 )
				continue;


			cerr << "reading " << argv[i] << "..." << flush;

	ifstream	in(argv[i]);
	GrepAwk	ga(in);
	
		
	
	
	ConfigurationSet	*toMerge;
		toMerge = new ConfigurationSet;

		while ( !ga.illegal() )
		{
			if ( ga.getWord(1)!=String("RefConfs") )
			{
				if ( first )
					cout << ga.getLine() << endl;
			}
			else
			{
				while ( !ga.illegal() && ga.getWord(1)!="}" )
				{
				Configuration<MOType>	conf;
					if ( isdigit(ga.getWord(1)[0]) )
					{
					char	buf[10000];
					stringstream	s(buf, stringstream::in | stringstream::out);
						s << ga.getLine() << endl;
					char	*p = strchr(buf, '#');
						if ( p )
							*p = 0;
						s >> conf;
						toMerge->add(conf);
					}
					ga++;
				}
				break;
			}
			ga++;
		}


			cerr << " o.k." << endl;

			switch ( mode ) {
			case merge:
				merged |= *toMerge;
				break;

			case intersect:
				merged &= *toMerge;
				break;

			case difference:
				merged -= *toMerge;
				break;

			default:
				break;
			}
			delete toMerge;
			first = 0;
		}
		cerr << "writing merged configurations to standard output..." << flush;
		cout << "RefConfs = {" << endl;
		merged.writeToStream(cout);
		cout << "}" << endl;
		cerr << " o.k." << endl;
		break;
	}	

	case tree:
	{
	NExternalsSet	*merged = NULL;
	INT	first = 1;
	MOMapping	*moMap = NULL;

		for ( INT i=1 ; i<argc ; i++ )
		{
			if ( argv[i][0]==0 )
				continue;

			cerr << "reading " << argv[i] << "..." << flush;

	ifstream	in(argv[i]);
			if ( first )
			{
				moMap = new MOMapping(in);
				merged = new NExternalsSet(in);
			}
			else
			{
			MOMapping	moMapping(in);
			NExternalsSet	toMerge(in);
				cerr << " o.k." << endl;

				if ( !merged->checkSameRefs(toMerge) )
				{
					cerr << "error:" << endl;
					cerr << "configuration tree in file " << argv[i] << endl;
					cerr << "is based on different references." << endl;
					exit(1);
				}
				cerr << "working..." << flush;

				switch ( mode ) {
				case merge:
					*merged |= toMerge;
					break;

				case intersect:
					*merged &= toMerge;
					break;

				case difference:
					*merged -= toMerge;
					break;

				default:
					break;
				}
			}
			cerr << " o.k." << endl;
			first = 0;
		}
		cerr << "writing merged configurations to standard output..." << flush;
		moMap->writeToStream(cout);
		delete moMap;
		merged->writeToStream(cout);
		delete merged;
		cerr << " o.k." << endl;
		break;
	}	
	}
	return 0;
}
