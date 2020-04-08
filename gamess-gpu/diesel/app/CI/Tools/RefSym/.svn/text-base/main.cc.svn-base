#include "../../../../config.h"

using namespace std;

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

#include "../../../../lib/Container/GrepAwk.h"

#include <vector>
#include <string>

#include "../../../../lib/QM/Configuration/Configuration.h"

#include "../../../../lib/QM/Symmetry/DegenPointgroups/Occupation.h"
#include "../../../../lib/QM/Symmetry/DegenPointgroups/ConfSet.h"
#include "../../../../lib/QM/Symmetry/DegenPointgroups/FullConf.h"
#include "../../../../lib/QM/Symmetry/DegenPointgroups/OrbitalEquivalence.h"


int	main(INT argc, char **argv)
{

OrbitalEquivalence	equiv;

	{
	GrepAwk	ga(cin);
		cout << "MO equivalences:" << endl;
		
		while ( !ga.illegal() )
		{
			cout << ga.getLine() << endl;
		vector<INT>	v;
			for ( INT i=0 ; i<ga.getNumberOfWords() ; i++ )
				v.push_back(atoi(ga.getWord(i+1)));
			equiv.push_back(v);
			ga++;
		}
	}






Fort31File	f31("fort.31");
MOIrReps	moirreps(f31);

IrRep	irrep = 0;

	irrep = atoi(argv[1]);

set<FullConf>	fcs;

ConfSet	orgRef;
	orgRef.clear();

	for ( INT ii=1 ; ii<argc ; ii++ )
	{
	ifstream	selInAll(argv[ii]);
		if ( !selInAll )
		{
			cerr << "no such file \"" << argv[ii] << "\"" << endl;
			exit(1);
		}

	GrepAwk	ga(selInAll);

		ga.grep("RefConfs");
		ga++;
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
//			irrep = conf.calcIrRep(moirreps);
		
			orgRef.insert(conf);			
		FullConf	fc(equiv(conf));
			fcs.insert(fc);
			ga++;
		} while ( ga.getNumberOfWords()>2 );
	}



	cout << endl << endl;
	cout << endl << endl;
	cout << "original reference configuration set (all irreps merged):" << endl;
	cout << orgRef << endl;
	cout << endl << endl;
	cout << endl << endl;
	cout << endl << endl;


	cout << "====================================" << endl;


ConfSet	newRef;
	newRef.clear();
	
	for ( set<FullConf>::const_iterator i = fcs.begin() ; i!=fcs.end() ; ++i )
	{
		cout << *i << ":" << endl; 
	ConfSet	confSet = equiv(*i);
	INT	inSym[8];
	INT count = 0;
		memset(inSym, 0, 8*sizeof(INT));
		for ( ConfSet::const_iterator j = confSet.begin() ; j!=confSet.end() ; ++j )
		{
			cout << (*j).calcIrRep(moirreps);
			if ( (*j).calcIrRep(moirreps)==irrep )
				cout << " +++ ";
			else
				cout << "     ";

		INT ii = orgRef.getInd(*j);
			if ( ii>=0 )
			{
				if ( (*j).calcIrRep(moirreps) == irrep )
					count++;
				cout << " (" << setw(3) << ii+1 << ") ";
			}
			else
				cout << "  " << "   " << "  ";
			

				
			cout << "\t";
			inSym[(*j).calcIrRep(moirreps)]++;
			cout << *j << endl;
		}
		cout << count << " / " << inSym[irrep] << endl;
		if ( 1.0*count/inSym[irrep]>=0*0.25 )
		{
			cout << "using " << inSym[irrep] << " configurations" << endl;
			for ( ConfSet::const_iterator j = confSet.begin() ; j!=confSet.end() ; ++j )
				if ( (*j).calcIrRep(moirreps)==irrep )
					newRef.insert(*j);
		}
		else
			cout << "discarding " << inSym[irrep] << " configurations" << endl;
	
		cout << endl;
		cout << endl;
		cout << endl;
	}
	
	cout << endl << endl;
	cout << endl << endl;
	cout << "====================================" << endl;
	cout << endl << endl;
	cout << "new reference configuration set:" << endl;
	cout << newRef << endl;
	cout << endl << endl;
	cout << endl << endl;
	cout << endl << endl;
	cout << endl;
	cout << "configurations added:" << endl;
	cout << newRef-orgRef << endl;

	cout << endl << endl;
ConfSet	discRef(orgRef-newRef);
	cout << "configurations discarded:" << endl;
	for ( ConfSet::const_iterator j = discRef.begin() ; j!=discRef.end() ; ++j )
	{
		if ( (*j).calcIrRep(moirreps) == irrep )
			cout << *j << endl;
	}
	
}








