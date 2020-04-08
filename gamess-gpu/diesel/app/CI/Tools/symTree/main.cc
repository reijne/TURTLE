#include <iostream>
#include <sstream>
#include <fstream>

using namespace std;

#include "../../../../lib/QM/MRTree/Diag/NExternalsDiag.h"
#include "../../../../lib/QM/MRTree/Set/NExternalsSet.h"
#include "../../../../lib/QM/MO/MOMapping.h"
#include "../../../../lib/QM/MO/MOIrReps.h"
#include "../../../../lib/QM/MRTree/MRTreeIterator.h"
#include "../../../../lib/QM/IO/Fortran/Fort31File.h"

#include "../../../../lib/QM/MO/Iterators/MOEquivalence.h"

#include "../../../../config.h"
#include "../../../../VersionDate.h"

int	main(INT argc, char **argv)
{
	cout << "symTree (Part of DIESEL-MR-CI), " << VERSION << ", " << DATE 
		<< endl << endl;

	if ( argc<2 || argc>4 )
	{
		cerr << endl;
		cerr << "purpose:" << endl;
		cerr << "the program symmetrizes the configurations in a tree (\"ConfTree.dat\"-format)." << endl; 
		cerr << endl;
		cerr << "usage: symTree \"ConfTree.dat\"-file \"equivalence\" [\"fort.31\"-file]" << endl;
		cerr << endl;
		exit(1);
	}
	
ifstream	ConfIn(argv[1]);
MOMapping	dummyMoMapping(ConfIn);

MOEquivalence	equiv(argv[2]);
	cout << "MO equivalence:" << equiv << endl << endl;

Fort31File	f31(argc == 4 ? argv[3] : "fort.31");
MOIrReps	moirreps(f31);

ConfigurationSet	set;
ConfigurationSet	ref;

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
			set.add(conf);
			if ( tree.getConfigurationSAFNr(i).isReference() )
				ref.add(conf);
		}
		cout << "OK" << endl;
	}
	

INT	n = set.length();
	cout << n << " configurations." << endl;
	cout << "symmetrizing..." << flush;
	equiv.symmetrize(set, &mrmos);
	cout << "OK" << endl;
	cout << set.length()-n << " configurations added." << endl;
	cout << set.length() << " configurations total." << endl;
	
//	cout << ref << endl;

String	of = String(argv[1])+String(".sym");
ofstream	ConfOut(of);
	{
		cout << "converting Tree..." << flush;
	NExternalsSet	tree
		(&mrmos,  NumberOfElectrons, Multiplicity, NumberOfRoots, ref);
//		tree.add(internal);
	Pix	i = set.first();
		while ( i )
		{
			tree.add(set(i));
			set.next(i);
		}
		cout << "OK" << endl;
		cout << "writing symmetrized tree to " << of << "..." << flush;
		dummyMoMapping.writeToStream(ConfOut);
		tree.cutTreeLevel2();
		tree.writeToStream(ConfOut);
		cout << "OK" << endl;
	}
		
	return 0;
}
