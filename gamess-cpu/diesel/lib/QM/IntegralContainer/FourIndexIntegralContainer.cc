//***********************************************************************
//
//	Name:			FourIndexIntegralContainer.cc
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




#include "FourIndexIntegralContainer.h"

#include "../IO/Fortran/FortranFileIO.h"
#include "../IO/Fortran/Fort31SecondRecord.h"

#include "../IO/TimeTicks.h"

#include "../MO/MOMapping.h"
#include "../IO/Verbosity.h"

using namespace std;

static const INT N = 8;

FourIndexIntegralContainer::FourIndexIntegralContainer(const MRMOs & _mrmos) :
	BaseContainer<SymmetryContainer>(N)
{
	SymmetryContainer::mrmos = _mrmos;
	mrmos = _mrmos;
	
	index = new IndexTranslation(_mrmos.getMaxMO());

	for ( INT i=0 ; i<N ; i++ )
		p[i] = new SymmetryContainer(i);
}


FourIndexIntegralContainer::~FourIndexIntegralContainer()
{
	delete index;
}


void	FourIndexIntegralContainer::loadIntegrals(Fort31File _f31)
{
/*

	for ( INT i=0 ; i<8 ; i++ )
		cout << mrmos.getStartMO(i) << endl;
		
	getchar();
	
*/	

/*TwoElectronIntegralTriadeIndex	indSort1(48, 46, 2, 1, 0);
	cout << indSort1 << endl;
	(*this)[indSort1] = 5.98;
	return;
*/
	
FortranFileIO	fort31(_f31.name);


INT	recSize = 2000;
IntegralType* buf = new IntegralType[recSize];
INT	integrals = 0, intsContained;
IntegralType	a;


TwoElectronIntegralIndex<MOType>	ind;
TwoElectronIntegralTriadeIndex	indSort;


	switch ( _f31.format ) {
	case GAMESSC1:
		{
			fort31.nextRecord();
			fort31.nextRecord();
			fort31.nextRecord();
		const INT recsize = 15000;
		struct Record {
			INT n; 
			unsigned char label[recsize][4];
			double	a[recsize];
			} data;
			
		INT	startMO = 0;
			{
				fort31.read(&data, sizeof(Record));
				startMO = data.label[0][0] -1;
				fort31.previousRecord();
			}

			while ( 1 ) 
			{
				fort31.read(&data, sizeof(Record));
				data.n = -data.n;
				if ( data.n<1 )
					break;

				for ( INT i=0 ; i<data.n ; ++i )
				{
					ind[0] = moMapping.getContinuous(data.label[i][3] - startMO);
					ind[1] = moMapping.getContinuous(data.label[i][2] - startMO);
					ind[2] = moMapping.getContinuous(data.label[i][1] - startMO);
					ind[3] = moMapping.getContinuous(data.label[i][0] - startMO);
					indSort.set(ind);
//					cout << "(" << ind[0] << " " << ind[1] << "|" << ind[2] << " " << ind[3] << "),  " 
//						<< indSort << ": " << data.a[i] << endl;
					intsContained += this->set(indSort, data.a[i]);
					integrals++;
				}
			}
		}
		break;
	default:
	{

	Fort31SecondRecord	f31(_f31.name, _f31.format);
	TimeTicks	ticks;

		intsContained = 0;
		for ( INT checking=0 ; checking<=0 ; checking++ )
	//	for ( INT checking=0 ; checking<=1 ; checking++ )
		{	if ( !checking )
				cout << "reading two electron integrals...";
			else
			{	cout << "checking two electron integrals...";
	//			check();
	//			exit(0);
			}
			cout.flush();
			ticks.start();
			fort31.rewind();
			if ( _f31.format==Fort31RecordFormatTRADPT )
				fort31.nextRecord();
			fort31.nextRecord();
			fort31.nextRecord();
			integrals = 0;

		INT	pos = recSize;
		INT	is, ic = 0;

		for ( INT i=0 ; i<mrmos.getNumberOfIrReps() ; i++ )
			for ( INT j=0 ; j<=i ; j++ )
				for ( INT k=0 ; k<=i ; k++ )
					for ( INT l=0 ; l<=((i==k) ? j : k) ; l++ )
					{//	cout << "(" << i << " " << j << "|" << k << " " << l << ")" << endl;

						is = 0;
						ic++;


						if ( mrmos.getProd(
								mrmos.getProd(i, j), mrmos.getProd(k, l)) )
							continue;
						if ( i==l )	
						{
	//					cout << "(AA|AA)" << endl;
		// (AA|AA) ---------------------------------------------------------
		for ( INT ii=0 ; ii<mrmos.getInIrRep(i) ; ii++ )
		{	ind[0] = moMapping.getContinuous(ii + mrmos.getStartMO(i));
			for ( INT jj=0 ; jj<=ii ; jj++ )
			{	ind[1] = moMapping.getContinuous(jj + mrmos.getStartMO(i));
				for ( INT kk=0 ; kk<=ii ; kk++ )
				{	ind[2] = moMapping.getContinuous(kk + mrmos.getStartMO(i));
					for ( INT ll=0 ; ll<=((ii==kk) ? jj : kk) ; ll++ )
					{	ind[3] = moMapping.getContinuous(ll + mrmos.getStartMO(i));
						indSort.set(ind);

						if ( pos>=recSize )
						{	pos = 0;
							fort31.read(buf, recSize*sizeof(IntegralType));
						}
	//					cout << indSort << endl;
						if ( checking )
						{
							if ( this->get(indSort, a) )
								if ( a != buf[pos] )
									cout << "Fehler: " << indSort << ": " << a << "!=" << buf[pos] << endl;
							pos++;
						}
						else
							intsContained += this->set(indSort, buf[pos++]);
						integrals++;
						is++;
					}
				}
			}
		}
	//					cout << is << " <--> " << f31.getNIT(ic)-f31.getNIT(ic-1) << endl;
						}
						else
						if ( i==j )
						{
	//					cout << "(AA|BB)" << endl;
		// (AA|BB) ---------------------------------------------------------
		for ( INT ii=0 ; ii<mrmos.getInIrRep(i) ; ii++ )
		{	ind[0] = moMapping.getContinuous(ii + mrmos.getStartMO(i));
			for ( INT jj=0 ; jj<=ii ; jj++ )
			{	ind[1] = moMapping.getContinuous(jj + mrmos.getStartMO(i));
				for ( INT kk=0 ; kk<mrmos.getInIrRep(k) ; kk++ )
				{	ind[2] = moMapping.getContinuous(kk + mrmos.getStartMO(k));
					for ( INT ll=0 ; ll<=kk ; ll++ )
					{	ind[3] = moMapping.getContinuous(ll + mrmos.getStartMO(k));
						indSort.set(ind);

						if ( pos>=recSize )
						{	pos = 0;
							fort31.read(buf, recSize*sizeof(IntegralType));
						}

	//					cout << indSort << endl;
						if ( checking )
						{
							if ( this->get(indSort, a) )
								if ( a != buf[pos] )
									cout << "Fehler: " << indSort << ": " << a << "!=" << buf[pos] << endl;
							pos++;
						}
						else
							intsContained += this->set(indSort, buf[pos++]);
						integrals++;
						is++;
					}
				}
			}
		}
	//					cout << is << " <--> " << f31.getNIT(ic)-f31.getNIT(ic-1) << endl;
						}
						else
						if ( i==k )
						{
	//					cout << "(AB|AB)" << endl;
		// (AB|AB) ---------------------------------------------------------
		for ( INT ii=0 ; ii<mrmos.getInIrRep(i) ; ii++ )
		{	ind[0] = moMapping.getContinuous(ii + mrmos.getStartMO(i));
			for ( INT jj=0 ; jj<mrmos.getInIrRep(j) ; jj++ )
			{	ind[1] = moMapping.getContinuous(jj + mrmos.getStartMO(j));
				for ( INT kk=0 ; kk<=ii ; kk++ )
				{	ind[2] = moMapping.getContinuous(kk + mrmos.getStartMO(i));
					for ( INT ll=0 ; ll<mrmos.getInIrRep(j) ; ll++ )
					{	if ( kk==ii && ll>jj )
							continue;

						ind[3] = moMapping.getContinuous(ll + mrmos.getStartMO(j));
						indSort.set(ind);
	//					cout << "indSort=" << indSort << endl;

						if ( pos>=recSize )
						{	pos = 0;
							fort31.read(buf, recSize*sizeof(IntegralType));
						}

						if ( checking )
						{
							if ( this->get(indSort, a) )
								if ( a != buf[pos] )
									cout << "Fehler: " << indSort << ": " << a << "!=" << buf[pos] << endl;
							pos++;
						}
						else
							intsContained += this->set(indSort, buf[pos++]);
	//					cout << "o.k." << endl;
						integrals++;
						is++;
					}
				}
			}
		}
	//					cout << is << " <--> " << f31.getNIT(ic)-f31.getNIT(ic-1) << endl;
						}
						else
						{
	//					cout << "(AB|CD)" << endl;
		// (AB|CD) ---------------------------------------------------------
		for ( INT ii=0 ; ii<mrmos.getInIrRep(i) ; ii++ )
		{	ind[0] = moMapping.getContinuous(ii + mrmos.getStartMO(i));
			for ( INT jj=0 ; jj<mrmos.getInIrRep(j) ; jj++ )
			{	ind[1] = moMapping.getContinuous(jj + mrmos.getStartMO(j));
				for ( INT kk=0 ; kk<mrmos.getInIrRep(k) ; kk++ )
				{	ind[2] =moMapping.getContinuous( kk + mrmos.getStartMO(k));
					for ( INT ll=0 ; ll<mrmos.getInIrRep(l) ; ll++ )
					{	ind[3] = moMapping.getContinuous(ll + mrmos.getStartMO(l));
						indSort.set(ind);
	//					cout << "indSort=" << indSort << endl;

						if ( pos>=recSize )
						{	pos = 0;
							fort31.read(buf, recSize*sizeof(IntegralType));
						}

						if ( checking )
						{
							if ( this->get(indSort, a) )
								if ( a != buf[pos] )
									cout << "Fehler: " << indSort << ": " << a << "!=" << buf[pos] << endl;
							pos++;
						}
						else
							intsContained += this->set(indSort, buf[pos++]);
						integrals++;
						is++;
					}
				}
			}
		}

	//					cout << is << " <--> " << f31.getNIT(ic)-f31.getNIT(ic-1) << endl;

						} // end if
					}	// end for l
			ticks.stop();
			cout << ticks << endl;
		}	
		if ( verbosity.isActive(Verbosity::Integrals) )
		{
			cout << "number of read two electron integrals: " << integrals << endl;
			cout << "number of stored two electron integrals: " << intsContained << endl;
		}

	/*	
		getchar();

	double  h = 0;
	double	cb, ex;

		ind[0] = 1;
		ind[1] = 3;
		ind[2] = 2;
		ind[3] = 4;
		indSort.set(ind);

	TwoElectronIntegralCbExIndex	indCbEx;
	TwoElectronIntegralIndex<MOType>	ind2(1, 2, 3, 4);

	NTupelContainer	*popt = getNTupelContainer(indSort);
		for ( INT i=0 ; i<1000000 ; i++ )
			for ( INT j=0 ; j<10 ; j++ )
			{	ind[3] = 4 + j + (i & 15);
	//			h += ind[3];

	//			indSort.set(ind);
				h += (*popt)[indSort];
	//
				indCbEx.set(ind, ind2);
				popt->get(indCbEx, cb, ex);
				h += cb;
	//			cout << indSort << "    " << (*popt)[indSort] << "    " << (*this)[indSort] << endl;
	//			cout << ind << "    " << indSort << endl;
	//			h += (*this)[indSort];
			}


		cout << "fertig" << h << endl;
		getchar();
	*/	

	/*	
		indSort = TwoElectronIntegralTriadeIndex(100, 100, 99, 99, 2);
		cout << indSort << " " << (*this)[indSort] << endl;
		getchar();
	*/

	/*	alt = -1;
		for ( INT i=0 ; i<n ; i++ )
		{	neu = list[i].request / recSize;
			while ( neu!=alt )
			{	fort31.read(buf);
	//			for ( INT j=0 ; j<10 ; j++ ) cout << buf[j] << " ";
				alt++;
			}
		}
	*/	
	break;
	}
	}

delete[] buf;
}
