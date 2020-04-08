//***********************************************************************
//
//	Name:			TwoIndexIntegralContainer.cc
//
//	Description:	root of tree containing four index integrals
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			03.09.1996
//
//
//
//
//
//***********************************************************************

#include "TwoIndexIntegralContainer.h"

#include "../IO/Fortran/FortranFileIO.h"
#include "../IO/Fortran/Fort31FirstRecord.h"
#include "../IO/Fortran/Fort31SecondRecord.h"

#include "../IO/TimeTicks.h"

#include "../MO/MOMapping.h"
#include "../IO/Verbosity.h"

using namespace std;

TwoIndexIntegralContainer::TwoIndexIntegralContainer(const MRMOs & _mrmos) 
{
	mrmos = &_mrmos;
	n = mrmos->getMaxMO()*(mrmos->getMaxMO()+1)/2;
	p = new OneElectronIntegralType[n];
	tab = new INT[mrmos->getMaxMO()];
	for ( INT i=0 ; i<mrmos->getMaxMO() ; i++ )
		tab[i] = i*(i+1)/2;
}


TwoIndexIntegralContainer::~TwoIndexIntegralContainer()
{
	delete[] p;
	delete[] tab;
}


double	TwoIndexIntegralContainer::loadIntegrals(Fort31File _f31)
{
FortranFileIO	fort31(_f31.name);
OneElectronIntegralIndex	ind;
INT	integrals = 0;

	switch ( _f31.format ) {
	case GAMESSC1:
		{
			fort31.nextRecord();
		double	core = 0;
			fort31.read(&core);
			
		INT	recSize = mrmos->getMaxMO()*(mrmos->getMaxMO()+1)/2;
		OneElectronIntegralType* buf = new OneElectronIntegralType[recSize];

		OneElectronIntegralType	*p = buf;

			fort31.read(buf, recSize*sizeof(OneElectronIntegralType));


			for ( INT i=0 ; i<mrmos->getMaxMO() ; i++ )
			{	
				ind.setI(moMapping.getContinuous(i+1));
				for ( INT j=0 ; j<=i ; j++ )
				{
					if ( mrmos->getProd(mrmos->getIrRep(i+1), mrmos->getIrRep(j+1)) )
						continue;
					ind.setJ(moMapping.getContinuous(j+1));
//					cout << ind << ": " << *p << endl;
					(*this)[ind] = *p;
					integrals++;
					p++;
				}
			}


			return core;
			delete[] buf;
		}
		break;

	default:
	Fort31FirstRecord	f31_1(_f31.name, _f31.format);
	Fort31SecondRecord	f31(_f31.name, _f31.format);
	TimeTicks	ticks;

		for ( INT checking=0 ; checking<=1 ; checking++ )
		{	if ( !checking )
				cout << "reading one electron integrals...";
			else
			{	cout << "checking one electron integrals...";
	//			check();
	//			exit(0);
			}
			cout.flush();
			ticks.start();

			switch ( _f31.format ) {
			case Fort31RecordFormatNew:
			{
				fort31.end();
			INT	recSize;
				fort31.previousRecord();
				fort31.previousRecord();
				fort31.previousRecord();
				fort31.previousRecord();
				recSize = f31_1.getN1el();
			OneElectronIntegralType* buf = new OneElectronIntegralType[recSize];

				fort31.read(buf, recSize*sizeof(OneElectronIntegralType));



		OneElectronIntegralType	*p = buf;
				for ( INT i=0 ; i<mrmos->getMaxMO() ; i++ )
				{	ind.setI(moMapping.getContinuous(i+1));
					for ( INT j=0 ; j<=i ; j++ )
					{	if ( mrmos->getProd(mrmos->getIrRep(i+1), mrmos->getIrRep(j+1)) )
							continue;
						ind.setJ(moMapping.getContinuous(j+1));
		//				cout << ind << ": " << *p << endl;
						if ( checking )
						{	if ( (*this)[ind] != *p )
							{	cout << "Fehler: " << ind << endl;
								cout << (*this)[ind] << " " << *p << endl;
							}
						}
						else
							(*this)[ind] = *p;
						integrals++;
						p++;
					}
				}
				delete[] buf;
			}
				break;

			case Fort31RecordFormatTRADPT:
			{
			INT	recSize = 2000;
			OneElectronIntegralType* buf = new OneElectronIntegralType[recSize];
			INT	pos = recSize;

				fort31.end();
				// skip VNUC
				fort31.previousRecord();

			INT	l = f31_1.getNBOX()*(f31_1.getNBOX()+1)/2;
				for ( INT i=0 ; i<(l+recSize-1)/recSize ; i++ )
				{
					fort31.previousRecord();
					fort31.previousRecord();
				}


				for ( INT i=0 ; i<mrmos->getMaxMO() ; i++ )
				{	ind.setI(moMapping.getContinuous(i+1));
					for ( INT j=0 ; j<=i ; j++ )
					{	if ( mrmos->getProd(mrmos->getIrRep(i+1), mrmos->getIrRep(j+1)) )
							continue;
						ind.setJ(moMapping.getContinuous(j+1));
						if ( pos>=recSize )
						{	pos = 0;
							fort31.read(buf, recSize*sizeof(OneElectronIntegralType));
						}
	//					cout << ind << ": " << buf[pos] << endl;
						if ( checking )
						{	if ( (*this)[ind] != buf[pos++] )
							{	cout << "Fehler: " << ind << endl;
								cout << (*this)[ind] << " " << buf[pos-1] << endl;
							}
						}
						else
							(*this)[ind] = buf[pos++];
						integrals++;
					}
				}
				delete[] buf;
			}
				break;
			default:
				exit(1);
			}


			ticks.stop();
			cout << ticks << endl;
		}	
		if ( verbosity.isActive(Verbosity::Integrals) )
			cout << "number of read one electron integrals: " << integrals << endl;

		fort31.end();
		fort31.previousRecord();
		fort31.read(&core, sizeof(double));
		return core;
	}
	
}
