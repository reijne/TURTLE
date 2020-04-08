#include "../../../../lib/QM/MRTree/Diag/NExternalsDiag.h"

#include "../../../../lib/QM/Davidson/Davidson.h"

#include "../../../../lib/QM/MRTree/MRTreeIterator.h"
#include "../../../../lib/QM/MO/MOMapping.h"

#include <stdlib.h>
#include "../../../../lib/Container/String.h"
#include <fstream>

#include <iomanip>

#include "../../../../lib/Math/etc/Histogram.h"

#include "../../../../lib/QM/Configuration/ConfigurationSet.h"

#include "../../Selector/MRConfInput.h"

#include "../../../../config.h"
#include "../../../../VersionDate.h"

using namespace std;

void	printDimHist(Histogram<double> ** hist, INT n)
{
	cout << "<\t";
	for ( INT j=0 ; j<n ; j++ )
		cout << hist[j]->getNumberOfHitsBelow() << "\t";
	cout << endl;
	for ( INT i=0 ; i<hist[0]->getSteps() ; i++ )
	{
		cout << (hist[0]->getIntervalMin(i)+hist[0]->getIntervalMax(i))/2 << "\t";
		for ( INT j=0 ; j<n ; j++ )
			cout << hist[j]->getFrequencyByIntervalNo(i) << "\t";
		cout << endl;
	}
	cout << "<\t";
	for ( INT j=0 ; j<n ; j++ )
		cout << hist[j]->getNumberOfHitsAbove() << "\t";
	cout << endl;
}





int	main()
{
	cout << "cmpPt (Part of DIESEL-MR-CI), " << VERSION << ", " << DATE
		<< endl << endl;
			
MRConfInput	mrinp;


	cin >> mrinp;
	cout << mrinp << endl;

Fort31File	f31("fort.31", Fort31RecordFormatNew);
MOIrReps	moirreps(f31);

	cout << "reading configurations..." << endl;
	
ifstream	ConfIn("ConfTree.dat");

	moMapping = MOMapping(ConfIn);
NExternalsDiag	MRDCIConfs(ConfIn, &moirreps);

INT	SAFs = MRDCIConfs.getNumberOfTotalSpinAdaptedFunctions();
//INT	RefSAFs = MRDCIConfs.getNumberOfRefConfSpinAdaptedFunctions();

	cout.precision(10);


typedef double MatrixType;
typedef double VectorType;

Vector<VectorType>	*ev = new Vector<VectorType>(SAFs);
DiskBuffer	*evBuf = new DiskBuffer(SAFs*sizeof(VectorType),
	"Eigenvectors.dat", 0);



ConfigurationSet	RefSet = mrinp.getRefConfSet();
ConfigurationSet	PTRefSet = mrinp.getPTRefConfSet();





const INT maxOpen = MAXSGATABLEOPENSHELLS;
Histogram<double>	*histOpenShells[maxOpen];

	for ( INT i=0 ; i<=maxOpen ; i++ )
		histOpenShells[i] = new Histogram<double>(-1, +5, Histogram<double>::Linear, 60);


const INT specs = 20;
Histogram<double>	histRel(-1, +5, Histogram<double>::Linear, 60);
Histogram<double>	*histRelAmount[specs];
Histogram<double>** histRelRefRel   = new Histogram<double>*[RefSet.length()+1];
Histogram<double>** histRelPTRefRel = new Histogram<double>*[PTRefSet.length()+1];




	for ( INT i=0 ; i<specs ; i++ )
		histRelAmount[i] = new Histogram<double>(-1, +5, Histogram<double>::Linear, 60);

	for ( INT i=0 ; i<=RefSet.length() ; i++ )
		histRelRefRel[i] = new Histogram<double>(-1, +5, Histogram<double>::Linear, 60);

	for ( INT i=0 ; i<=PTRefSet.length() ; i++ )
		histRelPTRefRel[i] = new Histogram<double>(-1, +5, Histogram<double>::Linear, 60);


Histogram<double>	histAbs(-1e-5, 1e-5, Histogram<double>::Linear, 60);
Histogram<double>	*histAbsAmount[specs];
Histogram<double>** histAbsRefRel   = new Histogram<double>*[RefSet.length()+1];
Histogram<double>** histAbsPTRefRel = new Histogram<double>*[PTRefSet.length()+1];

	for ( INT i=0 ; i<specs ; i++ )
		histAbsAmount[i] = new Histogram<double>(-1e-5, 1e-5, Histogram<double>::Linear, 60);

	for ( INT i=0 ; i<=RefSet.length() ; i++ )
		histAbsRefRel[i] = new Histogram<double>(-1e-5, 1e-5, Histogram<double>::Linear, 60);

	for ( INT i=0 ; i<=PTRefSet.length() ; i++ )
		histAbsPTRefRel[i] = new Histogram<double>(-1e-5, 1e-5, Histogram<double>::Linear, 60);


INT	roots = 1;
	for ( INT i=0 ; i<roots ; i++ )
	{
		evBuf->get(i, ev->getP());
INT	safNr = 0;
		for ( MRTreeIterator j = MRDCIConfs.firstInTree() ;
				!MRDCIConfs.isLastInTree(j) ; MRDCIConfs.nextInTree(j) )
		{
			ConfigurationSAFNr<MOType> main = MRDCIConfs.getConfigurationSAFNr(j);
			if ( !main.isReference() )
			{
				RootEnergies	re = main.getEnergy();	
			double	c2_CI = 0;
			double	c2_PT = 0;
			double	h, habs, hrel;
				for ( INT k=0 ; k<main.getSAFInc() ; k++ )
				{
					h = (*ev)[safNr+k];
					c2_CI += h*h;

					h = re.getCICoef(i, k);
					c2_PT += h*h;
				}

			INT	ind = (INT) floor(log(c2_CI)/log(4.0)+23);
				if ( ind<0 )
					ind = 0;
				if ( ind>=specs )
					ind = specs-1;
					
				habs = c2_PT-c2_CI;
				histAbs.count(habs);
				histAbsAmount[ind]->count(habs);
				hrel = habs / c2_CI;
				histRel.count(hrel);
				histRelAmount[ind]->count(hrel);
				
			Pix	pix = RefSet.first();
			INT	interact = 0;
				while ( pix )
				{
					interact += (Configuration<MOType>::calcExcitationOrder(
						main, RefSet(pix))<=2);
					RefSet.next(pix);
				}
				histAbsRefRel[interact]->count(habs);
				histRelRefRel[interact]->count(hrel);


				pix = PTRefSet.first();
				interact = 0;
				while ( pix )
				{
					interact += (Configuration<MOType>::calcExcitationOrder(
						main, PTRefSet(pix))<=2);
					PTRefSet.next(pix);
				}
				histAbsPTRefRel[interact]->count(habs);
				histRelPTRefRel[interact]->count(hrel);
				
				histOpenShells[main.getNumberOfOpenShells()]->count(hrel);


				


/*				cout << setw(20) << c2_CI <<
					setw(20) << (*ev)[safNr+0] << ": " 
						<< main << endl;
				for ( INT k=1 ; k<main.getSAFInc() ; k++ )
				cout << setw(20) << " " <<
					setw(20) << (*ev)[safNr+k] << endl;
*/			}

			
			safNr += main.getSAFInc();
		}
	}	
	delete ev;
	delete evBuf;
	
//	histRel.printPrettyASCII(80);


	cout << endl << endl << endl << endl;
	cout << "rel. Err.\t total\t" << endl;
	for ( INT i=0 ; i<histRel.getSteps() ; i++ )
		cout << ((INT) floor((histRel.getIntervalMin(i)+histRel.getIntervalMax(i))/2*100+0.5)) << "%\t"
			<< histRel.getFrequencyByIntervalNo(i) << endl;
			
	cout << endl << endl << endl << endl;
	cout << "rel. Err.\t";
	for ( INT j=0 ; j<specs ; j++ )
		cout << pow(4.0, j-23.0) << "\t";
	cout << endl;
	printDimHist(histRelAmount, specs);

	cout << endl << endl << endl << endl;
	cout << "rel. Err.\t";
	for ( INT j=0 ; j<=RefSet.length() ; j++ )
		cout << j << "\t";
	cout << endl;
	printDimHist(histRelRefRel, RefSet.length());
			
	cout << endl << endl << endl << endl;
	cout << "rel. Err.\t";
	for ( INT j=0 ; j<=PTRefSet.length() ; j++ )
		cout << j << "\t";
	cout << endl;
	printDimHist(histRelPTRefRel, PTRefSet.length());
			
			

	cout << endl << endl << endl << endl;
	cout << "abs. Err.\t total\t" << endl;
	for ( INT i=0 ; i<histAbs.getSteps() ; i++ )
		cout << (histAbs.getIntervalMin(i)+histAbs.getIntervalMax(i)) << "%\t"
			<< histAbs.getFrequencyByIntervalNo(i) << endl;
			
	cout << endl << endl << endl << endl;
	cout << "abs. Err.\t";
	for ( INT j=0 ; j<specs ; j++ )
		cout << pow(4.0, j-23.0) << "\t";
	cout << endl;
	printDimHist(histAbsAmount, specs);



	cout << endl << endl << endl << endl;
	cout << "abs. Err.\t";
	for ( INT j=0 ; j<=RefSet.length() ; j++ )
		cout << j << "\t";
	cout << endl;
	printDimHist(histAbsRefRel, RefSet.length());


			
	cout << endl << endl << endl << endl;
	cout << "abs. Err.\t";
	for ( INT j=0 ; j<=PTRefSet.length() ; j++ )
		cout << j << "\t";
	cout << endl;
	printDimHist(histAbsPTRefRel, PTRefSet.length());
			
			
	cout << endl << endl << endl << endl;
	cout << "rel. Err.\t";
	for ( INT j=0 ; j<maxOpen ; j++ )
		cout << j << "\t";
	cout << endl;
	printDimHist(histOpenShells, maxOpen);



	for ( INT i=0 ; i<specs ; i++ )
	{
		delete histAbsAmount[i];
		delete histRelAmount[i];
	}

	for ( INT i=0 ; i<=RefSet.length() ; i++ )
	{
		delete histAbsRefRel[i];
		delete histRelRefRel[i];
	}

	for ( INT i=0 ; i<=PTRefSet.length() ; i++ )
	{
		delete histAbsPTRefRel[i];
		delete histRelPTRefRel[i];
	}

	for ( INT i=0 ; i<=maxOpen ; i++ )
		delete histOpenShells[i];

        delete[] histAbsPTRefRel;
        delete[] histAbsRefRel;
        delete[] histRelPTRefRel;
        delete[] histRelRefRel;
	
}
