#include "RefSelInput.h"

#include "../../../lib/QM/MRTree/Diag/NExternalsDiag.h"
#include "../../../lib/QM/MRTree/Diag/InternalConfsDiag.h"

#include "../../../lib/QM/MRCIMatrix/MRCIMatrix.h"

#include "../../../lib/QM/Davidson/DavidsonCI.h"

#include "../../../lib/QM/MRTree/MRTreeIterator.h"
#include "../../../lib/QM/MO/MOMapping.h"

#include "../../../lib/QM/Configuration/ConfigurationSet.h"
#include "../../../lib/Math/etc/Histogram.h"

#include "../../../lib/QM/Symmetry/DegenPointgroups/RefSelMode.h"

#include "../../common/StartUp.h"

#include <string>
#include <sstream>
#include <fstream>
#include <stdlib.h>
#include <iomanip>
#include <vector>
#include <algorithm>

#include "../../../lib/QM/IO/Verbosity.h"

#include "../../../lib/QM/IO/TimeTicks.h"


#include "../../../Configured.h"
#include "Compiled.h"
#include "../../common/Banner.h"

#include "../../../config.h"
#include "../../../VersionDate.h"

using namespace std;

typedef double MatrixType;
typedef double VectorType;

void	MakeBanner()
{
const INT w = 80;

	MakeBannerTop(w, "R e f e r e n c e   S e l e c t i o n");
	center(VERSION, w);
	center(DATE, w);
	MakeBannerBottom(w);
}




struct ConfCISqr {
	ConfCISqr(Configuration<MOType> conf, double ci2) :
		conf(conf), ci2(ci2) {}

	bool operator < (const ConfCISqr c) const
	{	return ci2>c.ci2;	}

	Configuration<MOType>	conf;
	double	ci2;
};







#include "../../../lib/QM/Symmetry/DegenPointgroups/Occupation.h"
#include "../../../lib/QM/Symmetry/DegenPointgroups/ConfSet.h"
#include "../../../lib/QM/Symmetry/DegenPointgroups/FullConf.h"
#include "../../../lib/QM/Symmetry/DegenPointgroups/FullConfC2.h"
#include "../../../lib/QM/Symmetry/DegenPointgroups/OrbitalEquivalence.h"


struct FullConfCISqr {
	FullConfCISqr(FullConf fc, double ci2) :
		fc(fc), ci2(ci2) {}

	bool operator < (const FullConfCISqr c) const
	{	return ci2>c.ci2;	}

FullConf	fc;
double	ci2;
};





class FullConfSet : 
	public set<FullConfC2> {
public:

	FullConfSet(const MOEquivalence &moequiv, INT irreps, INT roots) :
		orbEquiv(moequiv), irreps(irreps), roots(roots) {}
	
	void	add(Configuration<MOType> conf, INT irrep, INT root, double _c2)
	{
	FullConf	fc(orbEquiv(conf));
	set<FullConfC2>::iterator i = find(fc);
		if ( i==end() )
			i = insert(FullConfC2(fc, irreps, roots)).first;
		
		const_cast<FullConfC2 &>(*i).c2[irrep][root] += _c2;
	}
	
	friend ostream & operator << (ostream &s, const FullConfSet &f);


private:
OrbitalEquivalence	orbEquiv;
INT	irreps;
INT	roots;
};



ostream & operator << (ostream &s, const FullConfSet &f)
{
	for ( FullConfSet::const_iterator i = f.begin() ; i!=f.end() ; ++i )
		s << *i << endl;

return s;
}


static INT ConfEcitationSortMode;

class ConfExcitation :
	public Configuration<MOType>,
	public vector<INT> {
public:
	ConfExcitation(Configuration<MOType> const & conf) :
		Configuration<MOType>(conf), vector<INT>(conf.getNumberOfElectrons())
	{}

	bool operator < (ConfExcitation const & a) const
	{
		if ( getNumberOfElectrons()<a.getNumberOfElectrons() )
			return true;
		if ( getNumberOfElectrons()>a.getNumberOfElectrons() )
			return false;
		switch ( ConfEcitationSortMode ) {
		case 0:
			for ( vector<INT>::const_iterator i=begin() ; i!=end() ; ++i )
			{
				if ( *i<a[i-begin()] )
					return true;
				if ( *i>a[i-begin()] )
					return false;
			}
			break;
		
		case 1:
			for ( vector<INT>::const_reverse_iterator i=rbegin() ; i!=rend() ; ++i )
			{
				if ( *i<a[rend()-i-1] )
					return true;
				if ( *i>a[rend()-i-1] )
					return false;
			}
			break;

		case 2:
			{
			INT sa = 0;
			INT sb = 0;
				for ( vector<INT>::const_iterator i=begin() ; i!=end() ; ++i )
				{
					sa += *i*(i-begin());
					sb += a[i-begin()]*(i-begin());
				}
				if ( sa<sb )
					return true;
				if ( sb<sa )
					return false;
			}
			break;
		}
		return 	static_cast<Configuration<MOType> const &>(*this) <
			static_cast<Configuration<MOType> const &>(a);
	}


};

ostream & operator << (ostream & s, ConfExcitation const & ce)
{
	s << static_cast<Configuration<MOType> const &>(ce) << ": ";
	for ( vector<INT>::const_iterator i=ce.begin() ; i!=ce.end() ; ++i )
		s << setw(4) << *i;
	return s;
}






TimeTicks	globalTime;


int main(INT argc, char **argv)
{
	globalTime.start();
	
	MakeBanner();


	StartUp();

bool	noDir = false;
INT	maxRoot = 100;
	if ( argc==2 )
	{
		noDir = true;
		maxRoot = atoi(argv[1]);
	}

			
RefSelInput	mrconf;

	cout.precision(6);
	
	cin >> mrconf;
	cout << mrconf << endl;
	
	
Fort31File	f31(mrconf.getMOIntegralFile());
MOIrReps	moirreps(f31);



SLList<INT>	irreps = mrconf.getIrReps();

INT	maxRoots = 0;
	{
	Pix	iIrrep = irreps.first();
		while ( iIrrep )
		{
		char	buf[100];
			if ( noDir )
				sprintf(buf, "Eigenvectors.dat.ref");
			else
				sprintf(buf, "%d/Eigenvectors.dat.ref", irreps(iIrrep));
				
		DiskBuffer	evBuf(buf, DiskBuffer::noTempDir);

			if ( evBuf.getNumberOfObjects()>maxRoots )
				maxRoots = evBuf.getNumberOfObjects();
			irreps.next(iIrrep);
		}
	}
	if ( maxRoots>maxRoot )
		maxRoots = maxRoot;
	

FullConfSet	fcs(*mrconf.getMOEquivalence(), moirreps.getNumberOfIrReps(), maxRoots);


vector<INT>	nRoots(8);



Pix	iIrrep = irreps.first();
	while ( iIrrep )
	{
	char	buf[10];
	INT irrep = irreps(iIrrep);
		sprintf(buf, "%d", irrep);
		if ( !noDir )
			chdir(buf);


		cout << "reading configurations..." << endl;

	ifstream	ConfIn("ConfTree.dat.ref");
		{
			MOMapping	dummy(ConfIn);
		}


	NExternalsDiag	MRDCIConfs(ConfIn, &moirreps);


	INT	SAFs = MRDCIConfs.getNumberOfTotalSpinAdaptedFunctions();




		cout.precision(10);



	Vector<VectorType>	*ev = new Vector<VectorType>(SAFs);
	DiskBuffer	*evBuf = new DiskBuffer(SAFs*sizeof(VectorType),
		"Eigenvectors.dat.ref", DiskBuffer::noTempDir);
		
		nRoots[irrep] = min(maxRoots, evBuf->getNumberOfObjects());

	ConfigurationSet	newRefs;
	ConfigurationSet	cs;





	//vector<Histogram<double> >	hists;

		for ( INT i=0 ; i<nRoots[irrep] ; i++ )
		{
		EnergyType	RefSum = 0;
		Histogram<double>	hist(1e-4, 0.5, 1);


			if ( verbosity.isActive(Verbosity::WaveFunction) )
			{
				cout << "root # " << i+1 << endl;
				cout << "                ci^2           CSF coef.   status   ";
				for ( INT j=0 ; j<nRoots[irrep] ; j++ )
					cout << "  PT dE root #" << j+1 << "  ";
				cout << "configuration" << endl;
				cout << "----------------------------------------------------";
				for ( INT j=0 ; j<nRoots[irrep] ; j++ )
					cout << "------------------";
				cout << "-----------------------------------------" << endl;
			}


			cout.setf(ios::fixed);
			cout.precision(10);

			evBuf->get(i, ev->getP());
	INT	safNr = 0;
	INT	EigenvectorSAF = 0;
			for ( MRTreeIterator j = MRDCIConfs.firstInTree() ;
					!MRDCIConfs.isLastInTree(j) ; MRDCIConfs.nextInTree(j) )
			{
				ConfigurationSAFNr<MOType> main = MRDCIConfs.getConfigurationSAFNr(j);
			double	c2 = 0;
			double	h;
				for ( INT k=0 ; k<main.getSAFInc() ; k++ )
				{
					h = (*ev)[safNr+k];
					c2 += h*h;
				}
				hist.count(c2);
				fcs.add(main, irreps(iIrrep), i, c2);

				if ( main.isReference() )
				{
					RefSum += c2;
					EigenvectorSAF += main.getSAFInc();
				}

				if ( verbosity.isActive(Verbosity::WaveFunction) )
				{
	//				if ( c2>mrconf.getRefThreshold() || main.isReference() )
	//				changed: 16.01.1998
					if ( c2/main.getSAFInc()>0.01 || main.isReference() )
					{
						cout << setw(20) << c2 <<
							setw(20) << (*ev)[safNr+0] << ": " 
								<< main;
						if ( main.isReference() )
						{
						InternalConfsDiag *internal = MRDCIConfs.operator [] (0);
						ContainerIterator iteri = internal->first();
						INT	n = 1;
							while ( (!(*internal)[iteri]->isReference()) || !(main == *(*internal)[iteri]) )
							{
								if ( (*internal)[iteri]->isReference() )
									n++;
								internal->next(iteri);
							}
							cout << "\t\t# " << n;
						}
						cout << endl;
						for ( INT k=1 ; k<main.getSAFInc() ; k++ )
						cout << setw(20) << " " <<
							setw(20) << (*ev)[safNr+k] << endl;
					}
				}

	//			if ( c2>mrconf.getRefThreshold() )
	//			changed: 16.01.1998
				if ( c2/main.getSAFInc()>mrconf.getReferenceThreshold() )
						newRefs.add(main);


				safNr += main.getSAFInc();
			}
			if ( verbosity.isActive(Verbosity::WaveFunction) )
			{
				cout << endl;
				cout << "ci^2 of references: " << RefSum << endl;
				cout << endl;
				cout << "===========================================" << endl;
				cout << endl;
			}

			hist.printPrettyASCII();


		}
		chdir("..");
		irreps.next(iIrrep);
		delete ev;
		delete evBuf;
	}

	
//	cout << fcs << endl;
	

set<FullConf>	refs;
	iIrrep = irreps.first();
	while ( iIrrep )
	{
		cout << "======================================================" << endl;
		cout << endl;
		cout << "irrep " << irreps(iIrrep) << endl;
		for ( INT i=0 ; i<nRoots[irreps(iIrrep)] ; i++ )
		{
			cout << "root #" << i+1 << ":" << endl;
		vector<FullConfCISqr>	v;


			for ( FullConfSet::const_iterator j=fcs.begin() ; j!=fcs.end() ; ++j )
				v.push_back(FullConfCISqr(*j, j->c2[irreps(iIrrep)][i]));

			sort(v.begin(), v.end());

		double	sum = 0;

			switch ( mrconf.getRefSelMode() ) {
			case RefSelMode::SumThresh:
				for ( vector<FullConfCISqr>::const_iterator j=v.begin() ; j!=v.end() ; ++j )
				{
					sum += j->ci2;
				FullConf	fc(j->fc);
					refs.insert(fc);
					cout << fc << "\t" << j->ci2 << "\t" << sum << endl;


					if ( sum>mrconf.getReferenceThreshold() )
						break;
				}
				break;

			case RefSelMode::ConfThresh:
				for ( vector<FullConfCISqr>::const_iterator j=v.begin() ; j!=v.end() ; ++j )
					if ( j->ci2>=mrconf.getReferenceThreshold() )
					{
						sum += j->ci2;
					FullConf	fc(j->fc);
						refs.insert(fc);
						cout << fc << "\t" << j->ci2 << "\t" << sum << endl;
					}
				break;
			}
			cout << endl << endl;
		}
		irreps.next(iIrrep);
	}



OrbitalEquivalence	orbEquiv(*mrconf.getMOEquivalence());

vector<INT>	activeOcc(orbEquiv.size(), 0);


ConfSet	total;
	total.clear();
	iIrrep = irreps.first();
	while ( iIrrep )
	{
		cout << "irrep " << irreps(iIrrep) << endl;
		cout << "configurations summing up to " << mrconf.getReferenceThreshold() << ":" << endl;
		cout << "{" << endl;		
	char	buf[10];
		sprintf(buf, "%d", irreps(iIrrep));
		chdir(buf);
	ofstream	of("refs.out");
		for ( set<FullConf>::const_iterator i=refs.begin() ; i!=refs.end() ; ++i )
		{
			for ( FullConf::const_iterator j=i->begin() ; j!=i->end() ; ++j )
			{
			INT	_j = j-i->begin();
				if ( activeOcc[_j]!=-1 )
				{
					if ( activeOcc[_j]==0 )
						activeOcc[_j] = *j;
					else
					if ( activeOcc[_j]!=*j )
						activeOcc[_j] = -1;
				}
			}
//			cout << *i << endl;
//			cout << orbEquiv(*i) << endl;
		ConfSet	cs(orbEquiv(*i));
			cs.projectOnIrRep(*mrconf.getMRMOs(), irreps(iIrrep));
			total += cs;
			cout << cs;
			cout << "#" <<endl;
			of << cs;
		}
		cout << "}" << endl;
		cout << endl;


		chdir("..");
		irreps.next(iIrrep);
	}

	{
	ofstream	of("refs.out");
		of << total;
	}

	{
	ofstream	of("active.out");
		for ( OrbitalEquivalence::const_iterator i=orbEquiv.begin() ; i!=orbEquiv.end() ; ++i )
		{
			if ( activeOcc[i-orbEquiv.begin()]==-1 )
				for ( vector<INT>::const_iterator j=i->begin() ; j!=i->end() ; ++j )
					of << *j << " ";
		}
		of << endl;
	}



	{
	ofstream	of("refExc.out");


		for ( ConfEcitationSortMode=0 ; ConfEcitationSortMode<=2 ; ++ConfEcitationSortMode )
		{
		vector<ConfExcitation>	confExc;
			for ( set<Configuration<MOType> >::const_iterator i=total.begin() ; i!=total.end() ; ++i )
			{
			ConfExcitation	ce(*i);
				for ( set<Configuration<MOType> >::const_iterator j=total.begin() ; j!=total.end() ; ++j )
				++ce[Configuration<MOType>::calcExcitationOrder(ce, *j)];
				confExc.push_back(ce);
			}
			sort(confExc.begin(), confExc.end());
			for ( vector<ConfExcitation>::const_iterator i=confExc.begin() ; i!=confExc.end() ; ++i )
				of << *i << endl;
			of << endl << char(12);
		}
	}

}


