#include "MRConfInput.h"
#include "Selector.h"

#include "../../../lib/QM/MRTree/Sel/NExternalsSel.h"
#include "../../../lib/QM/MRTree/Diag/NExternalsDiag.h"
#include "../../../lib/QM/MRTree/Sel/InternalConfsSel.h"

#include "../../../lib/QM/MO/Iterators/MOList.h"

#include "../../../lib/QM/RepresentationMatrices/RepresentationMatrices.h"

#include "../../../lib/QM/MO/MOMapping.h"
#include "../../../lib/QM/MO/Iterators/MOEquivalence.h"
#include "../../../lib/QM/IO/TimeTicks.h"

#include "../../common/StartUp.h"

#include "../../../lib/Container/String.h"
#include <sstream>
#include <fstream>
#include <stdlib.h>

#include <iomanip>

#include "EnlargeReferenceSpace.h"

#include "../../../lib/QM/Cache/VarSizeReadOnlyCache.h"

#include "../../../lib/QM/MRTree/EnergyMap/PTSum.h"

#include "../../../lib/QM/MRTree/Set/InternalConfsSet.h"
#include "../../../lib/QM/MRTree/Set/TupelStructureSet.h"
#include "../../../lib/QM/MRTree/Set/extMOsSet.h"
#include "../../../lib/QM/MRTree/Base/extEntry.h"

#include "../../../lib/QM/IO/Verbosity.h"

#include "../../../Configured.h"
#include "Compiled.h"
#include "../../common/Banner.h"

#include "../../../config.h"
#include "../../../VersionDate.h"

using namespace std;

void	MakeBanner()
{
const INT w = 80;

	MakeBannerTop(w, "S e l e c t o r");
	center(VERSION, w);
	center(DATE, w);
	MakeBannerBottom(w);
}






void	excitationStatistics(ConfigurationSet refs, const NExternalsSel *tree)
{
const INT maxRef = refs.length();
const INT maxExc = 100;
INT	exc[maxRef][maxExc];

	memset(exc, 0, maxRef*maxExc*sizeof(INT));

	for ( MRTreeIterator i = tree->firstInTree() ;
			!tree->isLastInTree(i) ; tree->nextInTree(i))
	{
	Pix	iRef = refs.first();
	INT	j = 0;
		while ( iRef )
		{
			cout << "A" << endl;
		Configuration<MOType> conf = tree->getConfigurationSAFNr(i);
			cout << "A" << endl;
			cout << conf << ":::" << refs(iRef) << endl;
			cout << "j=" << j << ", order =" << Configuration<MOType>::calcExcitationOrder(conf, refs(iRef)) << endl;
			exc[j][Configuration<MOType>::calcExcitationOrder(conf, refs(iRef))]++;
			cout << "A" << endl;
			refs.next(iRef);
			j++;
		}
		cout << "OK mriter " << i.n << endl;
	}
		cout << "-----OK" << endl;

INT i = 0;
	for ( i=maxExc-1 ; i>=0 ; i-- )
	{
	INT	breaked = 0;
		for ( INT j=0 ; j<maxRef ; j++ )
			if ( (breaked=exc[j][i]) )
				break;
		if ( breaked )
			break;
	}

Pix	iRef = refs.first();
INT	k = 0;
	for ( INT j=0 ; j<=i ; j++ )
		cout << setw(10) << j;
	cout << "    configuration" << endl;
	for ( INT j=0 ; j<=i ; j++ )
		cout << "----------";
	cout << "----------------------------------------------------------------" << endl;
	while ( iRef )
	{
		for ( INT j=0 ; j<=i ; j++ )
			cout << setw(10) << exc[k][j];
		cout << "    " << refs(iRef) << endl;
		refs.next(iRef);
		k++;
	}
}





TimeTicks	globalTime;

int main(int argc, char **argv)
{
	globalTime.start();

	MakeBanner();


	StartUp();

time_t	T;
	srand(time(&T)*getpid());
			
MRConfInput	mrinp;

	cout.precision(6);
	
	cin >> mrinp;
	cout << mrinp << endl;


RootEnergies	re;
	re.setStorePTEnergy(mrinp.getStorePTEnergy());
	re.setStorePTCoef(mrinp.getStorePTCoef());
	


INT	recalc = 0;

	if ( argc>1 )
	{
		if ( !strcmp(argv[1], "-r") )
		{
			recalc = 1;
			cout << "recalculation mode" << endl;
		}
		else
		{
			cerr << "unkown option " << argv[1] << endl;
			exit(1);
		}
	}

	if ( mrinp.getActiveReferenceThreshold()!=0 )
	{
	Selector<double, double>	selector(mrinp, NULL);

	ConfigurationSet	confs;
	PTReferenceSpace<double, double>	PTReference(mrinp, selector);
		confs = PTReference.getRefSet();
		
		if ( mrinp.getMOEquivalence() )
		{
			cout << "symmetrizing reference space according to MO equivalences..." << endl;
		INT	n = confs.length();
		ConfigurationSet	confsOrg(confs);
			((MOEquivalence *) mrinp.getMOEquivalence())->symmetrize(confs, mrinp.getMRMOs());
			cout << confs.length()-n << " configurations added:" << endl;
		ConfigurationSet	diff(confs);
			diff -= confsOrg;
			cout << diff << endl;
			cout << endl;
		}


		cout << endl;
		cout << endl;
		cout << "selected reference configuration (active space references above "
			<< mrinp.getActiveReferenceThreshold() << ")" << endl;
		cout << confs << endl;
		cout << endl;
		cout << endl;
		
		exit(0);
		mrinp.getRefConfSet() = mrinp.getPTRefConfSet() = confs;
		mrinp.initInternalExternal();
	}

PTSum	ptSum(
		mrinp.getNumberOfSelectionThresholds(), mrinp.getSelectionThresholdP(),
		mrinp.getEstimationMode()==EnergyMap::EpsteinNesbet,
		mrinp.getNumberOfRoots());


	for ( INT t=0 ; t<mrinp.getNumberOfSelectionThresholds() ; t++ )
	{
	NExternalsSet	*preSelConfs = NULL;
	NExternalsDiag	*diag = NULL;
		if ( recalc )
		{
		ifstream	in("ConfTree.dat."+String(mrinp.getSelectionThresholdString(t)));
			if ( !in )
			{
				cerr << "no such file ConfTree.dat." << mrinp.getSelectionThresholdString(t) << endl;
				exit(1);
			}
		MOMapping	moMapping(in);
			cout << "reading configuration tree " << 
				mrinp.getSelectionThresholdString(t) << "..." << flush;
			preSelConfs = new NExternalsSet(in);
			{
			ifstream	in("ConfTree.dat."+String(mrinp.getSelectionThresholdString(t)));
			MOMapping	moMapping(in);
				diag = new NExternalsDiag(in, mrinp.getMRMOs());
			}
			cout << "o.k." << endl;
			preSelConfs->init(0, mrinp.getMRMOs()->getMOSymmetryP());
			
		}
			
	Selector<double, double>	selector(mrinp, preSelConfs);

	PTReferenceSpace<double, double>	*PTReference = 
		new PTReferenceSpace<double, double>(mrinp, selector);

	NExternalsSel	*genConfs = selector.getGenConfs();
		genConfs->setNumberOfRoots(mrinp.getNumberOfRoots());
		for ( INT i=0 ; i<mrinp.getNumberOfRoots() ; i++ )
			genConfs->setRootNumber(i, mrinp.getRootNumber(i));



		selector.useAnnihilators();

		for ( INT i=0 ; i<=mrinp.getExcitationLevel() ; i++ )
			cout << "intern-" << i << ": " << (*genConfs)[i]->getNumberOfElements() << endl;


		cout << "**********************************************" << endl;

		cout << setiosflags(ios::fixed);


		selector.useCreators(PTReference,
		mrinp.getSelectionThreshold(mrinp.getNumberOfSelectionThresholds()-1));



		if ( verbosity.isActive(Verbosity::CacheStatistics) )
			cout << endl << endl << *selector.getCache() << endl << endl;

		cout << flush;

		if ( preSelConfs )
		{
			ptSum.setConfSAF(t, mrinp.getSelectionThreshold(t),
				diag->getNumberOfLeaves()-diag->getNumberOfReferences(),
				diag->getNumberOfTotalSpinAdaptedFunctions()-
				diag->getNumberOfRefConfSpinAdaptedFunctions());
			for ( INT r=0 ; r<mrinp.getNumberOfRoots() ; r++ )
				ptSum.setE(t, r, selector.getPreESum(r));
		}
		else
		{
			selector.printResults();
			cout << endl;
			cout << endl;


		/*	cout << "selected configurations due to nearly degenerated energy denominator:" << endl;
			cout << endl;
		INT	ndegen = 0;
			for ( MRTreeIterator j = genConfs->firstInTree() ;
					!genConfs->isLastInTree(j) ; genConfs->nextInTree(j) )
			{
				ConfigurationSAFNr<MOType> main = genConfs->getConfigurationSAFNr(j);
				if ( main.getSelectionFlags() & 1)
				{
					cout << main << endl;
					ndegen++;
				}
			}
			cout << endl;
			cout << ndegen << " configurations selected due to near degeneration of denominator. " << endl;
			cout << endl;
			cout << endl;
		*/

			for ( INT i=mrinp.getNumberOfSelectionThresholds()-1 ; i>=0  ; i-- )
			{
				cout << "cutting tree... " << flush;
			TimeTicks	ticks;
				ticks.start();
				selector.getGenConfs()->cutTreeLevel2();
				ticks.stop();
				cout << ticks << endl;	

				cout << selector.getGenConfs()->getNumberOfLeaves()
					<< " configurations in tree." << endl;


		/*		cout << endl << endl;
				cout << "-----------------------------------------------------------------------------" << endl;
				cout << "excitation statistics: " << endl;

				excitationStatistics(mrinp.getRefConfSet(), selector.getGenConfs());


				cout << "-----------------------------------------------------------------------------" << endl;
				cout << endl << endl;
		*/

				ticks.start();

			char	filename[1000];
				sprintf(filename, "ConfTree.dat.%s", mrinp.getSelectionThresholdString(i));

				cout << "writing to " << filename << "..." << flush;

			ofstream	ConfOut(filename);
				mrinp.getMOMappingNoPTRef().writeToStream(ConfOut);

				selector.getGenConfs()->writeToStream(ConfOut);
				ticks.stop();
				cout << ticks << "  ";	
				cout << "OK." << endl << endl;

				if ( i )
				{
					cout << endl << endl << endl;
					cout << "performing selection on threshold " 
						<< mrinp.getSelectionThreshold(i-1) << "..." << flush;
					ticks.start();
					selector.getGenConfs()->threshCut(mrinp.getSelectionThreshold(i-1), 
						mrinp.getEstimationMode()==EnergyMap::EpsteinNesbet);
					ticks.stop();
					cout << ticks << endl;	
				}
			}
			exit(0);
		}
		cout << flush;
//		delete preSelConfs;
//		delete diag;
	}	
	cout << "=================================================================" << endl;
	cout << endl;
	cout << "                      R e s u l t s" << endl;	
	cout << endl;
	cout << "=================================================================" << endl;
	cout << endl;
	cout << endl;
	cout << ptSum << endl << endl;
}


