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
	cerr << "generate plotable spectra in XFIG format" << endl;
	cerr << "usage:" << endl;
	cerr << "prop2xfig   " << endl;
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
	INT	root;
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

static int 	cmpDouble(const void *p1, const void *p2)
{
	if ( *((double *) p1) > *((double *) p2) )
		return 1;
	if ( *((double *) p1) < *((double *) p2) )
		return -1;
	return 0;
}



int	main(INT argc, char **argv)
{
	cerr << "prop2xfig (Part of DIESEL-MR-CI), " << VERSION << ", " << DATE << endl << endl;
	if ( argc!=2 )
	{
		laber();
		exit(1);
	}
	argc--;


//-------------------------------------------------------------------------

/*
ifstream	in(0);
GrepAwk	ga(in);
	
	ga.grep(argv[1]);
	ga.grep("root");
	ga += 2;
Pix	top = ga.getIndex();
INT	n = 0;
INT	nn = 0;
INT	flag = 0;
INT	irrep[8];
	while ( !ga.illegal() ) 
	{
		if ( !ga.getNumberOfWords() && flag )
			break;
		if ( ga.getNumberOfWords() )
		{
			if ( ga.getWordPos(1)<6 )
				irrep[nn++] = atoi(ga.getWord(1));
			
			n++;
			flag = 0;
		}
		else
			flag = 1;
		ga++;
	}
	
	ga.setIndex(top);
	
TData	data[n];

	n = 0;
	nn = 0;
	while ( !ga.illegal() ) 
	{
		if ( !ga.getNumberOfWords() && flag )
			break;
		if ( ga.getNumberOfWords() )
		{
			if ( ga.getWordPos(1)<6 )
			{
				data[n].E = atof(ga.getWord(3));
				data[n].root = atoi(ga.getWord(2));
			}
			else
			{
				data[n].E = atof(ga.getWord(2));
				data[n].root = atoi(ga.getWord(1));
			}
			data[n].irrep = irrep[nn];
			n++;
			flag = 0;
		}
		else
		{
			flag = 1;
			nn++;
		}
		ga++;
	}

	ga.setIndex(top);
	while ( ga.getWord(1)!="oscillator" )
		ga--;
	ga.grep("root");
	ga += 2;


double	osc[(n-1)*n/2];
INT	nosc = 0;
	while ( !ga.illegal() ) 
	{
		if ( !ga.getNumberOfWords() && flag )
			break;
		if ( ga.getNumberOfWords() )
		{
			if ( ga.getWordPos(1)<6 )
				for ( INT i=4 ; i<=ga.getNumberOfWords() ; i++ )
					osc[nosc++] = fabs(atof(ga.getWord(i)));
			else
				for ( INT i=3 ; i<=ga.getNumberOfWords() ; i++ )
					osc[nosc++] = fabs(atof(ga.getWord(i)));
			flag = 0;
		}
		else
			flag = 1;
		ga++;
	}


*/


//-----------------------------------------------------------------------

	
ifstream	in(0);
GrepAwk	ga(in);
	
	ga.grep(argv[1]);
	ga.grep("root");
	ga += 2;
Pix	top = ga.getIndex();
INT	n = 0;
INT	nn = 0;
INT	flag = 0;
INT	irrep[8];
	while ( !ga.illegal() ) 
	{
		if ( !ga.getNumberOfWords() && flag )
			break;
		if ( ga.getNumberOfWords() )
		{
			if ( ga.getWordPos(1)<6 )
				irrep[nn++] = atoi(ga.getWord(1));
			
			n++;
			flag = 0;
		}
		else
			flag = 1;
		ga++;
	}
	
	ga.setIndex(top);
	
TData	data[n];

	n = 0;
	nn = 0;
	while ( !ga.illegal() ) 
	{
		if ( !ga.getNumberOfWords() && flag )
			break;
		if ( ga.getNumberOfWords() )
		{
			if ( ga.getWordPos(1)<6 )
			{
				data[n].E = atof(ga.getWord(3));
				data[n].root = atoi(ga.getWord(2));
			}
			else
			{
				data[n].E = atof(ga.getWord(2));
				data[n].root = atoi(ga.getWord(1));
			}
			data[n].irrep = irrep[nn];
			n++;
			flag = 0;
		}
		else
		{
			flag = 1;
			nn++;
		}
		ga++;
	}

	ga.setIndex(top);
	while ( ga.getWord(1)!="oscillator" )
		ga--;
	ga.grep("root");
	ga += 2;


double	osc[(n-1)*n/2];
INT	nosc = 0;
	while ( !ga.illegal() ) 
	{
		if ( !ga.getNumberOfWords() && flag )
			break;
		if ( ga.getNumberOfWords() )
		{
			if ( ga.getWordPos(1)<6 )
				for ( INT i=4 ; i<=ga.getNumberOfWords() ; i++ )
					osc[nosc++] = fabs(atof(ga.getWord(i)));
			else
				for ( INT i=3 ; i<=ga.getNumberOfWords() ; i++ )
					osc[nosc++] = fabs(atof(ga.getWord(i)));
			flag = 0;
		}
		else
			flag = 1;
		ga++;
	}







//----------------------------------------------------------







double	oscSort[nosc];
	memcpy(oscSort, osc, nosc*sizeof(double));
	qsort(oscSort, nosc, sizeof(double), cmpDouble);
	

/*	cout << n << endl;
	for ( INT i=0 ; i<n ; i++ )
		cout << data[i].irrep << " " << data[i].root << ": " << data[i].E << endl;
		
	cout << "-------------------------" << endl;
*/	
TData	dataSort[n];
	memcpy(dataSort, data, n*sizeof(TData));
	qsort(dataSort, n, sizeof(TData), cmpTData);
	
/*	for ( INT i=0 ; i<n ; i++ )
		cout << dataSort[i].irrep << " " << dataSort[i].root << ": " << dataSort[i].E << endl;
*/

double	w = dataSort[n-1].E-dataSort[0].E;
double	wx = oscSort[nosc-1];
cout << "w=" << w << endl;
cout << "wx=" << wx << endl;


//	w = 11;
//	w = 15;
//	wx = 4;


double	e10 = exp(log(10.0)*floor(log10(w)));
INT	nScales = (INT) (w / e10) + 1;
	if ( nScales < 2 )
	{
		nScales *= 10;
		e10 /= 10;
	}

	cout << nScales << endl;


	// generate spectrum
	{
	ofstream	spec(String("spec.")+String(argv[1])+String(".fig"));

	INT		yrand = 500;
	INT		xrand = 800;
	INT		xoffset = 1200;
	INT	width = 4000;
	double	xscale = width /w;
	INT	height = 300;
	INT		yoffset = 3500;
	INT		ticklength = 50;
	double	yscale = -(yoffset-yrand-height)/wx;

		spec << 
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
			spec << "2 1 0 1 0 7 100 0 -1 4.000 0 0 -1 1 0 2" << endl;
			spec << "\t0 0 1.00 60.00 120.00" << endl;
		INT	x = (INT) ( xoffset );
			spec << "\t" << x << " " << yoffset << " " << x << " " << yrand << endl;
			spec << "4 1 0 100 0 0 12 0.0000 6 165 570 " << x << " " << yrand-50 << " $\\\\nu$\\001" << endl;
		}
			
		// draw absciss
		{
			spec << "2 1 0 1 0 7 100 0 -1 4.000 0 0 -1 1 0 2" << endl;
			spec << "\t0 0 1.00 60.00 120.00" << endl;
		INT	x = (INT) ( xoffset );
			spec << "\t" << xoffset << " " << yoffset << " " << xoffset+width << " " << yoffset << endl;


			spec << "4 0 0 100 0 0 12 0.0000 6 165 570 " << xoffset+width+100 << " " << yoffset << " $E$/eV\\001" << endl;
		}

		// draw absciss text
		{
		INT	y = (INT) yoffset;
			for ( INT i=0 ; i<=10 ; i++ )
			{
			double	E = i;
			INT	x = (INT) (E*xscale + xoffset);
				spec << "4 1 0 100 0 0 12 0.0000 6 165 570 " << x << " " << y+250 << " "
					<< E << "\\001" << endl;
				spec << "2 1 0 1 0 7 100 0 -1 4.000 0 0 -1 0 0 2" << endl
					<< "\t" << x << " " << y << " " << x << " " << y+ticklength << endl;
			}
		}

INT	iy = 0;
INT	ix = 0;
		for ( INT i=0 ; i<nosc ; i++ )
		{
			if ( ++ix>=n )
			{
				iy++;
				ix = iy+1;
			}
			if ( osc[i]>wx/100 )
			{
			double	E = fabs(data[iy].E-data[ix].E);
			INT	x = (INT) (E*xscale + xoffset);
			INT	y = (INT) (osc[i]*yscale + yoffset);
				spec << "2 1 0 1 0 7 100 0 -1 4.000 0 0 -1 0 0 2" << endl
					<< "\t" << x << " " << yoffset << " " << x << " " << y << endl;
				if ( osc[i]>wx/10 )
				{
//					cout << ix << " " << iy << " " << E << " " << osc[i] << endl;
					spec << "4 1 0 100 0 0 12 0.0000 6 165 570 " << x << " " << y-100 << " "
						<< "I" << data[iy].irrep << "R" << data[iy].root
						<< "$\\\\leftrightarrow $"
						<< "I" << data[ix].irrep << "R" << data[ix].root
						 << "\\001" << endl;
				}
			}
		}

	}


	// generate states
	{
	ofstream	states(String("states.")+String(argv[1])+String(".fig"));

	INT		xoffset = 1200;
	double	xscale = 800;
	INT	width = 300;
	INT	height = 300;
	INT		yoffset = 4000;
	INT		yrand = 500;
	INT		ticklength = 50;
	double	yscale = -(yoffset-yrand)/w;


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
		INT	x = (INT) (irrep[0]*xscale + xoffset - width);
		INT	y = (INT) (w*yscale + yoffset) - height;
			states << "\t" << x << " " << yoffset << " " << x << " " << y << endl;
			states << "4 1 0 100 0 0 12 0.0000 6 165 570 " << x << " " << y << " $E$/eV\\001" << endl;
		}
		// draw ordinate text
		{
		INT	x = (INT) (irrep[0]*xscale + xoffset - width - (ticklength+100));
			for ( INT i=0 ; i<=nScales ; i++ )
			{
			double	E = floor((i*w/nScales/e10)+0.5)*e10;
			INT	y = (INT) (E*yscale + yoffset);
				states << "4 2 0 100 0 0 12 0.0000 6 165 570 " << x << " " << y+50 << " "
					<< setw(4) << E << "\\001" << endl;
				states << "2 1 0 1 0 7 100 0 -1 4.000 0 0 -1 0 0 2" << endl
					<< "\t" << xoffset- width << " " << y << " " << xoffset- width-ticklength << " " << y << endl;
			}
		}
		// draw baseline
		{
			states << "2 1 1 1 0 7 100 0 -1 4.000 0 0 -1 0 0 2" << endl;
		INT	x1 = (INT) (irrep[0]*xscale + xoffset - width);
		INT	x2 = (INT) (irrep[nn-1]*xscale + xoffset + 2*width);
			states << "\t" << x1 << " " << yoffset << " " << x2 << " " << yoffset << endl;
		}
		// draw baseline text
		{
		INT	y = (INT) yoffset + height;
			for ( INT i=0 ; i<nn ; i++ )
			{
				INT	x = (INT) (irrep[i]*xscale + xoffset + width/2);
				states << "4 1 0 100 0 0 12 0.0000 6 165 570 " << x << " " << y
					<< " " << irrep[i] << "\\001" << endl;
			}
		}
		for ( INT i=0 ; i<n ; i++ )
		{
			if ( dataSort[i].E-dataSort[0].E <=w )
			{
				states << "2 1 0 2 0 7 100 0 -1 0.000 0 0 -1 0 0 2" << endl;
			INT	x = (INT) (dataSort[i].irrep*xscale + xoffset);
			INT	y = (INT) ((dataSort[i].E-dataSort[0].E)*yscale + yoffset);
				states << "\t" << x << " " << y << " " << x+width << " " << y << endl;
			}
		}
	}

		
}
