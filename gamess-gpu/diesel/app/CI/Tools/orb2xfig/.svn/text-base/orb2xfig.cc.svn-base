#include "../../../../config.h"
#include "../../../../VersionDate.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string>
#include <math.h>
#include <iostream>
#include <fstream>
#include <iomanip>

#include "../../../../lib/Container/String.h"

#include "../../../../lib/Container/GrepAwk.h"

using namespace std;

void	laber()
{
	cerr << "purpose:" << endl;
	cerr << "generate plotable orbital energies in XFIG format" << endl;
	cerr << "usage:" << endl;
	cerr << "orb2xfig  [min eV] [max eV] " << endl;
	cerr << endl;
	cerr << "options:" << endl;
}


double	convE(double e, String unit)
{
String	u(upcase(unit));
	if ( u=="AU" )
		return e;
	if ( u=="EV" )
		return e*27.211;
	if ( u=="KCAL/MOL" )
		return e*627.5;
	cerr << "unknown unit " << unit << endl;
	exit(1);
}


struct	TData {
	INT	irrep;
	INT	degen;
	double	occ;
	double	E;
	};


static int 	cmpTData(const void *p1, const void *p2)
{
	if ( ((TData *) p1)->E > ((TData *) p2)->E )
		return 1;
	if ( ((TData *) p1)->E < ((TData *) p2)->E )
		return -1;
	return 0;
}

static INT 	cmpDouble(const void *p1, const void *p2)
{
	if ( *((double *) p1) > *((double *) p2) )
		return 1;
	if ( *((double *) p1) < *((double *) p2) )
		return -1;
	return 0;
}



INT	main(INT argc, char **argv)
{
	cerr << "prop2xfig (Part of DIESEL-MR-CI), " << VERSION << ", " << DATE << endl << endl;
	if ( argc>3 )
	{
		laber();
		exit(1);
	}
	argc--;

	
ifstream	in(0);
GrepAwk	ga(in);
	ga.setIgnoreCase(1);
	
	ga.grep("Number of basis functions");
INT	nOrbs = 0;
INT	irrep[8];
INT	nIrreps = 0;
	for ( INT i=5 ; i<=ga.getNumberOfWords() ; i++ )
		nOrbs += (irrep[nIrreps++] = atoi(ga.getWord(i)));
	cout << "nIrreps=" << nIrreps << ",  nOrbs=" << nOrbs << endl;	

TData	data[nOrbs];
	memset(data, 0 , nOrbs*sizeof(TData));
	
INT	n = 0;	
	for ( INT i=0 ; i<nIrreps ; i++ )
	{
	char	s[100];
		ga.head();
		sprintf(s, "Molecular orbitals for symmetry species %d", i+1);
		ga.grep(s);
		while ( !ga.illegal() )
		{
			ga.grep("Energy");
			if ( ga.illegal() )
				break;
			for ( INT k=2 ; k<=ga.getNumberOfWords() ; k++ )
			{
				data[n].E = convE(atof(ga.getWord(k)), "eV");
				data[n].degen = 1;
				data[n++].irrep = i;
			}
			ga++;
			if ( ga.getWord(1)=="Occ." )
			{
				n -= ga.getNumberOfWords()-1;

				for ( INT k=2 ; k<=ga.getNumberOfWords() ; k++ )
					data[n++].occ = atof(ga.getWord(k));
			}
				
			ga.skipBehindNextBlankLine();
			ga.skipBehindNextBlankLine();
			if ( !ga.illegal() && ga.getWord(1)=="Molecular" )
				break;
			ga++;
		}
	}

/*	cout << n << endl;
	for ( INT i=0 ; i<n ; i++ )
		cout << data[i].irrep << ": " << data[i].E << " " << data[i].occ << endl;
		
	cout << "-------------------------" << endl;
*/	
	
/*	for ( INT i=0 ; i<n ; i++ )
		cout << dataSort[i].irrep << " " << dataSort[i].occ << " " << dataSort[i].root << ": " << dataSort[i].E << endl;
*/
TData	dataSort[n];
	memcpy(dataSort, data, n*sizeof(TData));
	qsort(dataSort, n, sizeof(TData), cmpTData);

	
double	min = dataSort[0].E;
double	max = dataSort[n-1].E;
	if ( argc+1>1 )
		min = atof(argv[1]);
	if ( argc+1>2 )
		max = atof(argv[2]);
		
	cout << "min= " << min << " , max=" << max << endl;

	for ( INT i=0 ; i<n ; i++ )
		cout << dataSort[i].irrep << " " << dataSort[i].occ << " " << ": " << dataSort[i].E << endl;

INT	nlow = 0;
INT	nlowIrrep[8];
INT	nhighIrrep[8];
	memset(nlowIrrep, 0, 8*sizeof(INT));
	memset(nhighIrrep, 0, 8*sizeof(INT));
	for ( INT i=0 ; i<n ; i++ )
	{
		if ( dataSort[i].E<min )
			nlowIrrep[dataSort[i].irrep]++;
		if ( dataSort[i].E>max )
			nhighIrrep[dataSort[i].irrep]++;
	}


	for ( INT i=0 ; i<n ; i++ )
		if ( dataSort[i].E>=min )
		{
			nlow = i;
			break;
		}

INT	nhigh = 0;
	for ( INT i=n-1 ; i>=nlow ; i-- )
		if ( dataSort[i].E<=max )
		{
			nhigh = i;
			break;
		}
		
	n = nhigh-nlow+1;


	cout << "nlow=" << nlow << endl;
	cout << "nhigh=" << nhigh << endl;

	cout << endl;

	if ( nlow )
		for ( INT i=0 ; i<n ; i++ )
			memcpy(&dataSort[i], &dataSort[i+(nlow)], sizeof(TData));

	for ( INT i=0 ; i<n ; i++ )
		cout << dataSort[i].irrep << " " << dataSort[i].occ << " " << ": " << dataSort[i].E << endl;


double	w = dataSort[n-1].E-dataSort[0].E;
cout << "w=" << w << endl;


//	w = 11;
//	w = 15;


double	e10 = exp(log(10.0)*floor(log10(w)));
INT	nScales = (INT) (w / e10) + 1;
	if ( nScales < 2 )
	{
		nScales *= 10;
		e10 /= 10;
	}

	cout << nScales << endl;



	for ( INT i=0 ; i<n-1 ; i++ )
	{
/*		if ( fabs(dataSort[i].E-dataSort[i+1].E)<1e-2 )
		{
			cout << i << endl;
			n--;
			dataSort[i].degen++;
			for ( INT k=i ; k<n ; k++ )
				memcpy(&dataSort[k], &dataSort[k+1], sizeof(TData));
		}
*/
		
	}






	// generate states
	{
	ofstream	states(String("states.fig"));

	INT		xoffset = 1200;
	double	xscale = 800;
	INT	width = 300;
	INT	height = 300;
	INT		yoffset = 4000;
	INT		yrand = 500;
	INT		ticklength = 50;
	double	yscale = -(yoffset-yrand)/w;
	double	yadd = -dataSort[0].E*yscale;

		states << 
			"#FIG 3.2" << endl <<
			"Landscape" << endl <<
			"Center" << endl <<
			"Metric" << endl <<
			"A4      " << endl <<
			"100.00" << endl <<
			"Single" << endl <<
			"-2" << endl <<
			"1200 2" << endl;

		// draw ordinate
		{
			states << "2 1 0 1 0 7 100 0 -1 4.000 0 0 -1 1 0 2" << endl;
			states << "\t0 0 1.00 60.00 120.00" << endl;
		INT	x = (INT) (0*xscale + xoffset - width);
		INT	y = (INT) (w*yscale + yoffset) - height;
			states << "\t" << x << " " << yoffset << " " << x << " " << y << endl;
			states << "4 1 0 100 0 0 12 0.0000 6 165 570 " << x << " " << y << " $E$/eV\\001" << endl;
		}
		// draw ordinate text
		{
		INT	x = (INT) (0*xscale + xoffset - width - (ticklength+100));
			for ( INT i=0 ; i<=nScales ; i++ )
			{
			double	E = floor(((i*w/nScales+dataSort[0].E)/e10)+0.5)*e10;
			INT	y = (INT) (E*yscale + yoffset +yadd);
				states << "4 2 0 100 0 0 12 0.0000 6 165 570 " << x << " " << y+50 << " "
					<< setw(4) << E << "\\001" << endl;
				states << "2 1 0 1 0 7 100 0 -1 4.000 0 0 -1 0 0 2" << endl
					<< "\t" << xoffset- width << " " << y << " " << xoffset- width-ticklength << " " << y << endl;
			}
		}
		// draw baseline
		{
			states << "2 1 1 1 0 7 100 0 -1 4.000 0 0 -1 0 0 2" << endl;
		INT	x1 = (INT) (0*xscale + xoffset - width);
		INT	x2 = (INT) ((nIrreps-1)*xscale + xoffset + 2*width);
			states << "\t" << x1 << " " << yoffset << " " << x2 << " " << yoffset << endl;
		}
		// draw baseline text
		{
		INT	y = (INT) yoffset + height;
			for ( INT i=0 ; i<nIrreps ; i++ )
			{
				INT	x = (INT) (i*xscale + xoffset + width/2);
				states << "4 1 0 100 0 0 12 0.0000 6 165 570 " << x << " " << y
					<< " " << i << "\\001" << endl;
				states << "4 1 0 100 0 0 12 0.0000 6 165 570 " << x << " " << y+200
					<< " " << nlowIrrep[i] << "/" << nhighIrrep[i] << "\\001" << endl;
			}
		}
		for ( INT i=0 ; i<n ; i++ )
		{
			if ( dataSort[i].E!=0 )
			{
				states << "2 1 0 2 0 7 100 0 -1 0.000 0 0 -1 0 0 2" << endl;
			INT	x = (INT) (dataSort[i].irrep*xscale + xoffset);
			INT	y = (INT) ((dataSort[i].E)*yscale + yoffset)+yadd;
				states << "\t" << x << " " << y << " " << x+width << " " << y << endl;
			}
		}
	}
		
}
