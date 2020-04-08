//***********************************************************************
//
//	Name:			MOEquivalenceProjected.cc
//
//	Description:	set of MOs to be handled equivalently
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			10.04.1998
//
//
//
//
//
//***********************************************************************




#include "MOEquivalenceProjected.h"
#include "MOEquivalence.h"
#include "../MOType.h"
#include "../../Configuration/ConfigurationSet.h"
#include "../MRMOs.h"


MOEquivalenceProjected::MOEquivalenceProjected(
	const MOEquivalence &e, const Configuration<MOType> &_conf,
	const MRMOs *_mrmos)
{
	conf = _conf;
	mrmos = _mrmos;
	irrep = conf.calcIrRep(*mrmos);
	
INT* flags1 = new INT[e.n];
INT* flags2 = new INT[e.n];
	memset(flags1, 0, e.n*sizeof(INT));
	memset(flags2, 0, e.n*sizeof(INT));

	n = 0;
INT	m = conf.getNumberOfOpenShells() + conf.getNumberOfClosedShells();	
MOType* mos = new MOType[m];
	memcpy(mos, conf.getOpenShellP(), conf.getNumberOfOpenShells()*sizeof(MOType));
	memcpy(mos + conf.getNumberOfOpenShells(),
		 conf.getClosedShellP(), conf.getNumberOfClosedShells()*sizeof(MOType));

	for ( INT i=0 ; i<e.n ; i++ )
	{
INT	nn = 0;
		for ( INT j=0 ; j<e.equiv[i].n ; j++ )
		{
			for ( INT k=0 ; k<m ; k++ )
				if ( e.equiv[i].mo[j]==mos[k] )
				{
					nn = 1;
					if ( k<_conf.getNumberOfOpenShells() )
					{
						flags1[i]++;
						conf.deleteSingleMO(mos[k]);
					}
					else
					{
						flags2[i]++;
						conf.deleteDoubleMO(mos[k]);
					}
					break;
				}
		}
		n += nn;
	}

//	cout << "n=" << n << endl;
	
	equiv = new TEquiv[n];
	doubly = new INT[n];
	singly = new INT[n];

INT	j = 0;
	for ( INT i=0 ; i<e.n ; i++ )
		if ( flags1[i] || flags2[i] )
		{
			equiv[j].n = e.equiv[i].n;
			equiv[j].mo = new MOType[equiv[j].n];
			memcpy(equiv[j].mo, e.equiv[i].mo, equiv[j].n*sizeof(MOType));
			singly[j] = flags1[i];
			doubly[j] = flags2[i];
//			cout << "j=" << j << ", " << singly[j] << " " << doubly[j] << endl;
			j++;
		}
	delete[] mos;
	delete[] flags2;
	delete[] flags1;
}




MOEquivalenceProjected::~MOEquivalenceProjected()
{
	delete doubly;
	delete singly;
}



void	MOEquivalenceProjected::generate(ConfigurationSet &confSet)
{
	generate(confSet, conf, 0);
}


void	MOEquivalenceProjected::generate(
	ConfigurationSet &confSet,
	Configuration<MOType> conf,
	INT level)
{
	if ( level<n )
	{

//		cout << ":::::::::" << level << ":::::::::" << equiv[level].n << " " <<  doubly[level] << endl;
//		cout << "level1:" << level << endl;
	Iterator	iter1(equiv[level].n, doubly[level]);
//		cout << "level2:" << level << endl;
		while ( iter1.next() )
		{
//		cout << "level3:" << level << endl;
		Configuration<MOType> c1(conf);
//			cout << level << ":::::::::" << "1: ";
			for ( INT j=0 ; j<doubly[level] ; j++ )
			{
//				cout << " " << iter1.ii[j];
				c1.insertClosedMO(equiv[level].mo[iter1.ii[j]]);
			}
//			cout << endl;

		// filter out already inserted doubly occupied MOs
		MOType* mo = new MOType[equiv[level].n-doubly[level]];
			{
			INT j = 0;
			INT	k = 0;
				for ( INT i=0 ; i<equiv[level].n ; i++ )
					if ( j<doubly[level] )
					{
						if ( i==iter1.ii[j] )
							j++;
						else
							mo[k++] = equiv[level].mo[i];
					}
					else
						mo[k++] = equiv[level].mo[i];
			}
			
//			cout << "singly=" << singly[level] << endl;
//			cout << equiv[level].n << " " << doubly[level] << endl;
		Iterator	iter2(equiv[level].n-doubly[level], singly[level]);
			while ( iter2.next() )
			{
			Configuration<MOType> c2(c1);
//				cout<< level << ":::::::::"  << "2: ";
				for ( INT j=0 ; j<singly[level] ; j++ )
				{
//					cout << " " << iter2.ii[j];
					c2.insertOpenMO(mo[iter2.ii[j]]);
				}
//				cout << endl;
				generate(confSet, c2, level+1);
			}
		delete[] mo;
		}
	}
	else
	{
//		cout << "conf=" << conf << endl;
		if ( conf.calcIrRep(*mrmos)==irrep )
				confSet.add(conf);
	}
}



void	MOEquivalenceProjected::distribute(INT n, MOType *mo, INT m)
{
Iterator	iter(n, m);

	while ( iter.next() )
	{
		for ( INT j=0 ; j<m ; j++ )
			cout << " " << iter.ii[j];
		cout << endl;
	}
}










MOEquivalenceProjected::Iterator::Iterator(INT _n, INT _m)
{
	n = _n;
	m = _m;
	i = 0;
	if ( m>0 )
	{
		i = 0;
		ii = new INT[m];
		ii[i] = -1;
	}
	else
	{
		i = -1;
		ii = NULL;
	}
//	cout << "n=" << n << ", m=" << m << ", i=" << i << endl;
}

MOEquivalenceProjected::Iterator::~Iterator()
{
	if ( ii )
		delete ii;
}

//	distribute "n" objects into "m" containers
INT	MOEquivalenceProjected::Iterator::next()
{
//	cout << "   i=" << i << endl;
	if ( i==-2 )
		return 0;
	if ( i==-1 )
	{
		i--;
		return 1;
	}
	ii[i]++;
	for ( ; ; )
	{
		if ( ii[i]>n-(m-i) )
		{
			i--;
			if ( i<0 )
				return 0;
			ii[i]++;
			continue;
		}
		if ( i<m-1 )
		{
			ii[i+1] = ii[i]+1;
			i++;
			continue;
		}
		return 1;
	}
}

