#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>

#include "../../../../lib/Container/GrepAwk.h"

#include <ctype.h>

#include "../../../../config.h"
#include "../../../../VersionDate.h"

#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdlib.h>

#include <vector>
#include <string>
#include <cmath>

using namespace std;

//static const enum { quadratic, down } mode = down;
static const enum { quadratic, down } mode = quadratic;

typedef	vector<vector<vector<double> > >	V3;

struct TLine;

struct TCoord {


	double & operator [] (char xyz)
	{
		switch ( xyz ) {
		case 'x': 
			return x;

		case 'y': 
			return y;

		case 'z': 
			return z;

		}
		return x;
	}

	void	add(char xyz, double a)
	{
		(*this) [xyz] += a;
	}

	INT	isZero(char xyz)
	{
		return (*this) [xyz] == 0;
	}


string	pre;
double	x, y, z;
string	post;
vector<TLine>::iterator	linkedTo;
};

struct	TLine {
	TLine(string s) 
	{
		plain = s; 
		isEntry = 0; 
	}
	
	TLine(string pre, double x, double y, double z, string post)
	{ 
		entry.pre = pre;
		entry.x = x;
		entry.y = y;
		entry.z = z;
		entry.post = post;
		entry.linkedTo = (vector<TLine>::iterator) 0;
		isEntry = 1;
	}

string	plain;
TCoord	entry;
INT	isEntry;
};


ostream	& operator << (ostream	&s, const TLine &tl)
{
	if ( tl.isEntry )
	{
//		s << "[" << tl.entry.linkedTo << "] ";
		s << tl.entry.pre << " ";
		s << tl.entry.x << " ";
		s << tl.entry.y << " ";
		s << tl.entry.z << " ";
		s << tl.entry.post;
	}
	else
		s << tl.plain;
	return s;
}


class Coords {
public:
	Coords(string sewIn);
	
	void	step(INT atom, char xyz, double delta);
	INT	isZero(INT atom, char xyz);
	
	INT	update(double Ebase, V3 Estep, double delta);
	
	
	INT	getTriples() const { return entries-linked; }
	
	friend ostream & operator << (ostream &s, const Coords &c);
	
private:
	vector<TLine>::iterator operator () (INT atom, char xyz);

vector<TLine>	lines;
INT	entries;
INT	linked;


};


Coords::Coords(string sewIn)
{
ifstream	f(sewIn.c_str());
GrepAwk	ga(f);


INT	nl = 0;
	while ( !ga.illegal() )
	{
		lines.push_back(TLine(ga.getLine().chars()));
		if ( lines[nl].plain=="End of basis" )
		{
		INT	m = 0;
			while ( 1 )
			{
				ga--;
				m++;
			double	h;
			
				if ( ga.getNumberOfWords()<4 || ga.getNumberOfWords()>5 ||
					sscanf(ga.getWord(1).chars(), "%lf", &h)==1 )
					break;

			double	x, y, z;
			char	dummy[1000];
				if ( sscanf(ga.getLine().chars(), "%s %lf %lf %lf", 
					dummy, &x, &y, &z) == 4 )
					lines[nl-m] = TLine(
						ga.getWord(1).chars(), x, y, z, ga.getWord(5).chars());
				else
					break;
			}
			ga += m;
		}
		nl++;
		ga++;
	}

const double	eps = 1e-3;
	entries = 0;
	linked = 0;
	for ( vector<TLine>::iterator i=lines.begin() ; i<lines.end() ; i++ )
	{
		entries += !!i->isEntry;
		if ( i->isEntry && ! ( i->entry.linkedTo!=(vector<TLine>::iterator)0) )
			for ( vector<TLine>::iterator j=i+1 ; j<lines.end() ; j++ )
				if ( j->isEntry )
					if ( 
						fabs(i->entry.x - j->entry.x)<eps &&
						fabs(i->entry.y - j->entry.y)<eps &&
						fabs(i->entry.z - j->entry.z)<eps )
					{
						j->entry.linkedTo = i;
						linked++;
					}
	}
}	
	
	
ostream	& operator << (ostream	&s, const Coords &c)
{
	s << "* " << c.entries << endl;
	s << "* " << c.linked << endl;
	for ( vector<TLine>::const_iterator i=c.lines.begin() ; i<c.lines.end() ; i++ )
		s << *i << endl;
	return s;
}


vector<TLine>::iterator Coords::operator () (INT atom, char xyz)
{
vector<TLine>::iterator i;
INT	nl = 0;
	for ( i=lines.begin() ; i<lines.end() ; i++ , nl++ )
	{
		if ( i->isEntry && !(i->entry.linkedTo==(vector<TLine>::iterator)0) )
			atom--;
		if ( atom<0 )
			break;
	}
	if ( atom>=0 )
		return i;
	return i;
}

void	Coords::step(INT atom, char xyz, double delta)
{
vector<TLine>::iterator i = (*this)(atom, xyz);
	
	i->entry.add(xyz, delta);
	for ( vector<TLine>::iterator j=i ; j<lines.end() ; j++ )
		if ( j->isEntry && j->entry.linkedTo==i )
			j->entry.add(xyz, delta);

}

INT	Coords::isZero(INT atom, char xyz)
{
	return (*this)(atom, xyz)->entry.isZero(xyz);
}


INT	Coords::update(double Ebase, V3 Estep, double delta)
{
double	max = 0;
	for ( INT atom=0 ; atom<entries-linked ; atom++ )
	{
		for ( char xyz='x' ; xyz<='z' ; xyz++ )
		{
		vector<TLine>::iterator i = (*this)(atom, xyz);
		double d = 0;
			cout << endl << endl << *i << ":" << endl;
		
			if ( mode==quadratic )
			{
			double	y0 = Ebase;
			double	ym = Estep[atom][xyz-'x'][0];
			double	yp = Estep[atom][xyz-'x'][1];
				if ( ym==y0 || yp==y0 )
						continue;
			double	x0 = i->entry[xyz];
			double	xm = x0-delta;
			double	xp = x0+delta;

			double	b = ((y0-yp)*xm*xm + (yp-ym)*x0*x0 + (ym-y0)*xp*xp) /
						(-2*delta*delta*delta);
			double	a = (ym-y0-b*(-delta))/(xm*xm-x0*x0);
			double	c = y0-(a*x0*x0 + b*x0);

				cout << "a=" << a << ", b=" << b << ", c=" << c << endl;


				cout << xm << " " << a*xm*xm + b*xm +c << endl;
				cout << x0 << " " <<  a*x0*x0 + b*x0 +c << endl;
				cout << xp << " " <<  a*xp*xp + b*xp +c << endl;

			double x = -b/(2*a);
				d = x - i->entry[xyz];
				cout << "-b/(2*a)=" << -b/(2*a) << ", old=" << i->entry[xyz] << ", d=" << d << endl;

				if ( a<0 )
				{
					cout << "wrong curvature" << endl;
					if ( yp<ym )
						d = 0.002;
					else
						d = -0.002;
				}


				if ( fabs(d)>0.005 )
					d = (d>0 ? 0.005 : -0.005);

				if ( fabs(d)>max )
					max = fabs(d);
			}
			else
			{
			double	y0 = Ebase;
			double	yp = Estep[atom][xyz-'x'][1];
			
				if ( yp==y0 )
					continue;
				if ( yp<y0 )
					d = delta;
				else
					d = -delta;
			
				if ( fabs(yp-y0)>max )
					max = fabs(yp-y0);
			}
			cout << "d=" << d << "\t";
			i->entry[xyz] += d;
		}
		cout << endl;
	}


	return (max<0.0001);
}



int	main(INT argc, char **argv)
{
	cerr << "geomopt (Part of DIESEL-MR-CI), " << VERSION << ", " << DATE 
		<< endl << endl;;
	
	if ( argc!=1 )
	{
		cerr << endl;
		cerr << "purpose:" << endl;
		cerr << endl;
		cerr << "usage: geomopt -{c|s|o|t} -{m|i|d} file1 file2 ..." << endl;
		cerr << endl;
		cerr << endl;
		exit(1);
	}
	
	cout.precision(10);
	
Coords	coords("sew.in");	
//	cout << coords << endl;


const double	delta = 0.01;
const INT maxiter = 1000;
	for ( INT iter=0 ; iter<maxiter ; iter++ )
	{
	char	buf[10];
		sprintf(buf, "%d", iter);
		mkdir(buf, 0750);
		chdir(buf);

	
	double	Ebase = 0;
		{
			mkdir("base", 0750);
			chdir("base");

		ifstream	r("result");
			if ( r )
				r >> Ebase;
			else
			{
				system("cp ../../* . 2>/dev/null");
				{
					ofstream	sewIn("sew.in");
					sewIn << coords;
				}
				system("./job 2>/dev/null");
			ifstream	r("result");
				r >> Ebase;
			}

			chdir("..");
		}
	
	INT	triples = coords.getTriples();
	V3	Estep(triples);
	
		for ( INT atom=0 ; atom<triples ; atom++ )
		{
		char	buf[10];
			sprintf(buf, "%d", atom);
			mkdir(buf, 0750);
			chdir(buf);

			Estep[atom].resize(3);
			
			for ( char xyz='x' ; xyz<='z' ; xyz++ )
			{
				Estep[atom][xyz-'x'].resize(2);
				if ( coords.isZero(atom, xyz) )
				{
					Estep[atom][xyz-'x'][0] = Ebase;
					Estep[atom][xyz-'x'][1] = Ebase;
					continue;
				}
				
				for ( INT f=( mode==down ? 1 : 0 ) ; f<=1 ; f++ )
				{
				char	buf[10];
					sprintf(buf, "%c%c", xyz, (f?'+':'-'));
					mkdir(buf, 0750);
					chdir(buf);
					{
					ifstream	r("result");
						if ( r )
						{
							r >> Estep[atom][xyz-'x'][f];
							chdir("..");
							continue;
						}
					}

					system("cp ../../../* . 2>/dev/null");

					{
					Coords	step(coords);
						step.step(atom, xyz, (2*f-1)*delta);
					ofstream	sewIn("sew.in");
						sewIn << step;
					}
					system("./job 2>/dev/null");

				ifstream	r("result");
					r >> Estep[atom][xyz-'x'][f];

				
					chdir("..");
				}
			}
			chdir("..");
		
		}
		chdir("..");
		
		
		cout << Ebase << endl;
		for ( INT atom=0 ; atom<triples ; atom++ )
		{
			for ( char xyz='x' ; xyz<='z' ; xyz++ )
			{
				for ( INT f=0 ; f<=1 ; f++ )
					cout << Estep[atom][xyz-'x'][f] << "\t";
				cout << "\t";
			}
			cout << endl;
		}


		if ( coords.update(Ebase, Estep, delta) )
			break;
	
	}
	

	return 0;
}
