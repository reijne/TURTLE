//***********************************************************************
//
//	Name:			MainRefSet.h
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			25.08.1999
//
//
//
//
//
//***********************************************************************

#include "MainRefSet.h"

#include <algorithm>
#include <stdlib.h>

using std::sort;

static INT root = 0;

INT	cmp(const MainRefSet::Data p1, const MainRefSet::Data p2)
{
//	cout << root << endl;
//	cout << p1.conf << endl;
//	cout << p2.conf << endl;
//	cout << endl;
	
	return p1.ci2[root] > p2.ci2[root];
}


void	MainRefSet::add(Configuration<MOType> conf, const double *ci2)
{
Data	d;
	d.conf = conf;
	for ( INT i=0 ; i<roots ; i++ )
		d.ci2.push_back(ci2[i]);
	data.push_back(d);
}
	
void	MainRefSet::calc()
{
	for ( INT i=0 ; i<roots ; i++ )
	{
//		for ( unsigned INT j=0 ; j<data.size() ; j++ )
//		{
//			cout << "AAA" << j << " " << data[j].conf << " " << data[j].ci2[i] << endl;
//		}
		root = i;
		sort(data.begin(), data.end(), cmp);
//		for ( unsigned INT j=0 ; j<data.size() ; j++ )
//		{
//			cout << "BBB" << j << " " << data[j].conf << " " << data[j].ci2[i] << endl;
//		}
		
	double sum = 0;	
		for ( unsigned INT j=0 ; j<data.size() ; j++ )
		{
			cout << "root " << i+1 << ": adding " << data[j].conf
				<< " with ci^2 of " << data[j].ci2[i] << endl;
			sum += data[j].ci2[i];
			confs.add(data[j].conf);
			if ( sum>=limit && 
				!(j<data.size()-1 && fabs(data[j+1].ci2[i] - data[j].ci2[i])<1e-4) )
				break;
		}
		cout << "root " << i+1 << ": sum= " << sum << endl;
	}

	cout << endl;
	for ( INT i=0 ; i<roots ; i++ )
	{
	double sum = 0;
		for ( unsigned INT j=0 ; j<data.size() ; j++ )
			if ( confs.contains(data[j].conf) )
				sum += data[j].ci2[i];

		cout << "root " << i+1 << ": sum= " << sum << endl;
	}
	
	
	
}
