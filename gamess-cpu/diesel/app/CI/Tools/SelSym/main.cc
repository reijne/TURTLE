using namespace std;

#include <iostream>
#include <sstream>
#include <fstream>

#include "../../../../lib/QM/MRTree/Diag/NExternalsDiag.h"
#include "../../../../lib/QM/MRTree/Set/NExternalsSet.h"
#include "../../../../lib/QM/MO/MOMapping.h"
#include "../../../../lib/QM/MO/MOIrReps.h"
#include "../../../../lib/QM/MRTree/MRTreeIterator.h"
#include "../../../../lib/QM/IO/Fortran/Fort31File.h"

#include "../../../../lib/QM/Configuration/Excitation.h"

#include "../../../../lib/QM/MO/Iterators/MOEquivalence.h"
#include "../../../../lib/QM/Symmetry/DegenPointgroups/Occupation.h"
#include "../../../../lib/QM/Symmetry/DegenPointgroups/ConfSet.h"
#include "../../../../lib/QM/Symmetry/DegenPointgroups/FullConf.h"
#include "../../../../lib/QM/Symmetry/DegenPointgroups/OrbitalEquivalence.h"

#include "../../../../config.h"
#include "../../../../VersionDate.h"

int	main(INT argc, char **argv)
{
	cout << "symTree (Part of DIESEL-MR-CI), " << VERSION << ", " << DATE 
		<< endl << endl;

	if ( argc<3 || argc>5 ) 
	{
		cerr << endl;
		cerr << "purpose:" << endl;
		cerr << "the program symmetrizes the configurations in a tree (\"ConfTree.dat\"-format)." << endl; 
		cerr << endl;
		cerr << "usage: selsym \"ConfTree.dat\"-file [-h] [-e  \"equivalence\"] [\"fort.31\"-file]" << endl;
		cerr << endl;
		exit(1);
	}

ifstream	ConfIn(argv[1]);
MOMapping	dummyMoMapping(ConfIn);


bool	homo = false;
MOEquivalence	*equiv = 0;

	for ( INT i=2 ; i<argc ; i++ )
	{
		if ( !strcmp(argv[i], "-h") )
		{
			homo = true;
			continue;
		}
		if ( !strcmp(argv[i], "-e") )
		{
			equiv = new MOEquivalence(argv[i+1]);
			cout << "MO equivalence:" << *equiv << endl << endl;
			continue;
		}
	}
	


Fort31File	f31("fort.31");
//Fort31File	f31(argc == 3 ? argv[2] : "fort.31");
MOIrReps	moirreps(f31);


ConfigurationSet	sel;
ConfigurationSet	refs;

INT NumberOfElectrons = 0;
INT Multiplicity = 0;
INT NumberOfRoots = 0;

MRMOs	mrmos;
	{
		cout << "reading tree..." << flush;
	NExternalsDiag	tree(ConfIn, &moirreps);
//		cout << mrmos << endl;

		mrmos = MRMOs(*tree.getMRMOs());

		NumberOfElectrons = tree.getNumberOfElectrons();
		Multiplicity =  tree.getMultiplicity();
		NumberOfRoots = tree.getNumberOfRoots();
		cout << "OK" << endl;
		cout << "converting Tree..." << flush;
		for ( MRTreeIterator i = tree.firstInTree() ;
				!tree.isLastInTree(i) ; tree.nextInTree(i))
		{
		Configuration<MOType> conf = tree.getConfigurationSAFNr(i);
			sel.add(conf);
			if ( tree.getConfigurationSAFNr(i).isReference() )
				refs.add(conf);
		}
		cout << "OK" << endl;
	}
	
INT	irrep = sel(sel.first()).calcIrRep(moirreps);
	

INT	n = sel.length();
	cout << n << " configurations." << endl;
	
	if ( homo )
	{
		cout << "analyzing..." << endl;


	set<Excitation>	excs;

		{
		Pix i = refs.first();
			while ( i )
			{
			Configuration<MOType>	ref(refs(i));

				cout << "ref " << ref << "...";
			Pix j = sel.first();
				while ( j )
				{
					if ( Configuration<MOType>::calcExcitationOrderFast(
							ref, sel(j))<=2 )
					{		
					Excitation	exc;
						exc.calcExcitation(ref, sel(j));
						excs.insert(exc);

	//					cout << ref << "\t" << sel(j) << "\t" << exc << endl;
					}

					sel.next(j);
				}
	//			cout << endl << "---------------" << endl << endl;
				refs.next(i);
				cout << endl;
			}
		}


		cout << endl << endl << endl;
		cout << excs.size() << " excitations." << endl;
	//	for ( set<Excitation>::iterator i=excs.begin() ; i!=excs.end() ; ++i )
	//		cout << *i << endl;


		sel.clear();


		cout << "applying..." << endl;

		{
		Pix i = refs.first();
			while ( i )
			{
				cout << "ref " << refs(i) << "...";
				for ( set<Excitation>::iterator j=excs.begin() ; j!=excs.end() ; ++j )
				{
				Configuration<MOType>	conf(refs(i));
	//				cout << conf << "\t" << *j << "\t";
					if ( !conf.excite(*j) )
					{
						sel.add(conf);
	//					cout << conf;
					}
	//				cout << endl;
				}
				cout << endl;
				refs.next(i);
			}
		}
		cout << endl;
		cout << "homogeneous excitations:" << endl;
		cout << sel.length()-n << " configurations added." << endl;
		cout << sel.length() << " configurations total." << endl;
		cout << endl << endl;
	}
	

	if ( equiv )
	{
	INT nn = sel.length();


	OrbitalEquivalence	orbEquiv(*equiv);
	set<FullConf>	sfc;
	
		{
		Pix	i = sel.first();
			while ( i )
			{
				sfc.insert(orbEquiv(sel(i)));
				sel.next(i);
			}
		}
		
		cout << "full sym confs:" << sfc.size() << endl;
		sel.clear();
		
	INT ii = 0 ;
		for ( set<FullConf>::const_iterator i=sfc.begin() ; i!=sfc.end() ; ++i )
		{
		ConfSet	cs(orbEquiv(*i));
			cout << ++ii << endl << *i << endl << cs << endl << endl;
			cs.projectOnIrRep(moirreps, irrep);
			for ( ConfSet::const_iterator j=cs.begin() ; j!=cs.end() ; ++j )
			{
			Configuration<MOType>	conf(*j);
				sel.add(conf);
			}
		}
		




		cout << endl;
		cout << "symmetrized selection:" << endl;
		cout << sel.length()-nn << " configurations added." << endl;
		cout << sel.length() << " configurations total." << endl;
		cout << endl << endl;
	}



	cout << endl;
	cout << "total:" << endl;
	cout << sel.length()-n << " configurations added." << endl;
	cout << endl << endl;

	



	
//	cout << refs << endl;

String	of = String(argv[1])+String(".selsym");
ofstream	ConfOut(of);
	{
		cout << "converting Tree..." << flush;
	NExternalsSet	tree
		(&mrmos,  NumberOfElectrons, Multiplicity, NumberOfRoots, refs);
//		tree.add(internal);
	Pix	i = sel.first();
		while ( i )
		{
			tree.add(sel(i));
			sel.next(i);
		}
		cout << "OK" << endl;
		cout << "writing symmetrized tree to " << of << "..." << flush;
		dummyMoMapping.writeToStream(ConfOut);
		tree.cutTreeLevel2();
		tree.writeToStream(ConfOut);
		cout << "OK" << endl;
	}

	if ( equiv )
		delete equiv;		

	return 0;
}
