#include "../../../../config.h"
#include "../../../../VersionDate.h"

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include <dirent.h>
#include <string>

#include <math.h>

#include <iostream>

#include <iomanip>
#include "../../../../lib/Container/String.h"

#include "../../../../lib/Container/AVLMap.h"
#include "../../../../lib/Container/SLList.h"

#include "../../../../lib/Container/GrepAwk.h"

using namespace std;

void	laber()
{
	cerr << "purpose:" << endl;
	cerr << "collect and print data (energies, properties) from a DIESEL-CI calculation" << endl;
	cerr << "grouped by irreducible representations and roots." << endl << endl;
	cerr << "usage:" << endl;
	cerr << "pretty  -global options  -printing modes " << endl;
	cerr << endl;
	cerr << "global options:" << endl;
	cerr << "  -unit {au, eV, kcal/mol}" << endl;
	cerr << "   default: au" << endl;
	cerr << endl;
	cerr << "  -sort" << endl;
	cerr << "   order output by transition energies" << endl;
	cerr << endl;
	cerr << "  -mode {E_CI, E_extpolCI, E_Dav1, E_Dav2}" << endl;
	cerr << "   default: E_Dav2" << endl;
	cerr << endl;
	cerr << "  -thresh value" << endl;
	cerr << "   no default, must be specified when using printing modes \"-e\" or \"-osc\"" << endl;
	cerr << endl;
	cerr << "  -intpath IntegralPath" << endl;
	cerr << "   default: ../ONEINT" << endl;
	cerr << endl;
	cerr << "  -orbpath OrbitalPath" << endl;
	cerr << "   default: ../SCFORB" << endl;
	cerr << endl;
	cerr << "printing modes:" << endl;
	cerr << "  -e" << endl;
	cerr << "   print excitation energies" << endl;
	cerr << endl;
	cerr << "  -prop operator component" << endl;
	cerr << "   calculate properties (depends on density matrices)" << endl;
	cerr << endl;
	cerr << "  -osc" << endl;
	cerr << "   calculate oscillator strengths (depends on density matrices) " << endl;
	cerr << endl;
}



String	centerTo(String s, unsigned INT w)
{
	if ( w<s.length() )
		w = s.length();
unsigned INT	l = (w-s.length())/2;
	for ( unsigned INT i=0 ; i<l ; i++ )
		s.prepend(' ');

	while ( s.length()<w )
		s += ' ';
	return s;
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


void	getPropCol(String p, INT &col, INT &block)
{
String	u(upcase(p));
	if ( u=="E_CI" )
	{
		col = 4;
		block = 1;
		return;
	}
	if ( u=="E_EXTPOLCI" )
	{
		col = 4;
		block = 2;
		return;
	}
	if ( u=="E_DAV1" )
	{
		col = 3;
		block = 3;
		return;
	}
	if ( u=="E_DAV2" )
	{
		col = 5;
		block = 3;
		return;
	}
	cerr << "unknown energy class " << p << endl;
	exit(1);
}



struct	TProp {
	double	value;
	String	name;
	};
	
	

class PropList : public SLList<TProp> {
public:

	double	getProp(String);
	void	append(String, double);
};

double	PropList::getProp(String name)
{
Pix	pix = first();
	while ( pix )
	{
		if ( (*this)(pix).name == name )
			return (*this)(pix).value;
		next(pix);
	}
	return 0;
}

	
void	PropList::append(String name, double value)
{
struct TProp	prop;
	prop.value = value;
	prop.name = name;
	SLList<TProp>::append(prop);
}

typedef AVLMap<String, PropList> TMap;





struct	TSortData {
	double	e;
	Pix	pix;
	};

static int cmpFnc(const void *p1, const void *p2)
{
	if ( ((TSortData *) p1)->e > ((TSortData *) p2)->e )
		return 1;
	else
	if ( ((TSortData *) p1)->e < ((TSortData *) p2)->e )
		return -1;
	return 0;
}

//----------------------------------------------------------------------------
struct	TIrrep {
	INT	n;
	SLList<INT>	roots;
	};

typedef	SLList<TIrrep> TIrrepList;

void	Map2Irreps(TMap &map, TIrrepList &irrepList)
{
Pix	pix = map.first();
INT	oldIrrep1 = -1;
INT	oldRoot1 = -1;
TIrrep	ti;
	while ( pix )
	{
	INT	irrep1, irrep2;
	INT root1, root2;

		sscanf(map.key(pix).chars(), "I%dR%d_I%dR%d\n", &irrep1, &root1, &irrep2, &root2);

//		cout << "AAA " << map.key(pix).chars() << endl;
//			cout << irrep1 << " " << root1 << " " << irrep2 << " " << root2 << endl;

		map.next(pix);
		if ( irrep1!=oldIrrep1 && oldIrrep1!=-1 || !pix )
		{
			if ( !pix && root1!=oldRoot1 ) 
				ti.roots.append(root1);
			ti.n = oldIrrep1;
			irrepList.append(ti);
			ti.roots.clear();
			oldRoot1 = -1;
		}

		if ( root1!=oldRoot1 )
			ti.roots.append(root1);

		oldRoot1 = root1;

		oldIrrep1 = irrep1;
	}
}

void	prettyPrintSpec(String title, TMap &map, String prop)
{
Pix	pix = map.first();
struct	TSortData	data[map.length()];
	{
	INT	i = 0;
		while ( pix )
		{
			data[i].e = map.contents(pix).getProp("E");
			data[i].pix = pix;
			i++;
			map.next(pix);
		}
	}
	qsort(data, map.length(), sizeof(TSortData), cmpFnc);

const INT	width = 14;		

	cout << setw(width) << "excitation";
	pix = map.contents(data[0].pix).first();
	while ( pix )
	{
		cout << setw(width) << map.contents(data[0].pix)(pix).name;
		map.contents(data[0].pix).next(pix);
	}
	cout << endl;


	for ( INT i=0 ; i<map.length() ; i++ )
	{
		cout << setw(width) << map.key(data[i].pix);
		pix = map.contents(data[i].pix).first();
		while ( pix )
		{
			cout << setw(width) << map.contents(data[i].pix)(pix).value;
			map.contents(data[i].pix).next(pix);
		}
		cout << endl;

	}
}

void	GnuPlotE(TMap &map)
{
double	Emin = 1e10;
Pix	minPix = NULL;

	{
	Pix	pix = map.first();
		while ( pix )
		{
		INT	i1, i2, r1, r2;
			sscanf(map.key(pix), "I%dR%d_I%dR%d", &i1, &r1, &i2, &r2);
			if ( i1==i2 && r1==r2 )
			{
				if ( map.contents(pix).getProp("E")<Emin )
				{
					Emin = map.contents(pix).getProp("E");
					minPix = pix;
				}
			}
			map.next(pix);
		}
	}

INT	I1, I2, R1, R2;
	sscanf(map.key(minPix), "I%dR%d_I%dR%d", &I1, &R1, &I2, &R2);

	{
	Pix	pix = map.first();
		while ( pix )
		{
		INT	i1, i2, r1, r2;
			sscanf(map.key(pix), "I%dR%d_I%dR%d", &i1, &r1, &i2, &r2);
			if ( I1==i1 && R1==r1 )
				cout << i2 << " " << r2 << ": " << map.contents(pix).getProp("E") << endl;
			if ( I1==i2 && R1==r2 )
				cout << i1 << " " << r1 << ": " << -map.contents(pix).getProp("E") << endl;
			map.next(pix);
		}
	}
}

void	GnuPlotSpec(TMap &map)
{
Pix	pix = map.first();
	while ( pix )
	{
	INT	i1, i2, r1, r2;
		sscanf(map.key(pix), "I%dR%d_I%dR%d", &i1, &r1, &i2, &r2);
		if ( i1!=i2 && r1!=r2 )
			cout << map.key(pix) << ": " << fabs(map.contents(pix).getProp("E")) << " " << map.contents(pix).getProp("osc") << endl;
		map.next(pix);
	}
}

void	prettyPrint(String title, TMap &map, String prop, INT flag = 1)
{
SLList<TIrrep>	irreps;


	Map2Irreps(map, irreps);
	
	
	if ( flag ) 
	{
	INT	w = 80;
	INT	fw = 9;
	
	Pix	iIrrep = irreps.first();
		while ( iIrrep )
		{
		Pix	jIrrep = irreps.first();
			while ( jIrrep )
			{
				if ( irreps(iIrrep).n<=irreps(jIrrep).n )
				{
				cout << centerTo(title, w) << endl;
				char	s[1000];
					sprintf(s, "irrep %d <-> irrep %d\n", irreps(iIrrep).n, irreps(jIrrep).n);
				cout << centerTo(s, w) << endl;


			Pix	iRoot = irreps(iIrrep).roots.first();
				cout << "           ";
				 	while (iRoot)
				{
				char	s[10];
					sprintf(s, "%d", irreps(iIrrep).roots(iRoot));
					printf("%s", centerTo(s, fw+2).chars());
					irreps(iIrrep).roots.next(iRoot);
				}
				cout << endl;

			iRoot = irreps(iIrrep).roots.first();
		INT	iiRoot = 0;
			while ( iRoot )
			{
				{
				char	s[10];
					sprintf(s, "%d", iiRoot+1);
					printf("%s", centerTo(s, fw+2).chars());
				}


				Pix	jRoot = irreps(jIrrep).roots.first();
					while ( jRoot )
					{
						if ( irreps(iIrrep).n<irreps(jIrrep).n || 
								irreps(iIrrep).n==irreps(jIrrep).n && 
								irreps(iIrrep).roots(iRoot)<=irreps(jIrrep).roots(jRoot) )
						{
						char	state[1000];
							sprintf(state, "I%dR%d_I%dR%d", 
								irreps(iIrrep).n, irreps(iIrrep).roots(iRoot),
								irreps(jIrrep).n, irreps(jIrrep).roots(jRoot));
							if ( map.contains(state) )
								printf("%9.5f  ", map[state].getProp(prop));
							else
								printf("    ---    ");
						}
						else
							printf("%9s  ", " ");

						irreps(jIrrep).roots.next(jRoot);
					}
					printf("\n");
				irreps(iIrrep).roots.next(iRoot);
				iiRoot++;
				}
				}
				irreps.next(jIrrep);
				printf("\n\n");
			}
			printf("\n---------------------------------------------------------------------------------\n\n");
			irreps.next(iIrrep);
		}
		printf("\n\n=================================================================================\n\n");
		printf("\n\n\n");
	}
	else
	{
	INT	fw = 9;
		{
		INT	w = 15;
		Pix	iIrrep = irreps.first();
			while ( iIrrep )
			{
				w += (fw+3)*irreps(iIrrep).roots.length();
				irreps.next(iIrrep);
			}
			cout << centerTo(title, w) << endl << endl;
		}
		printf("irrep");
		for ( INT i=0 ; i<=fw ; i++ )
			printf(" ");
		{
		Pix	iIrrep = irreps.first();
			while ( iIrrep )
			{
			char	s[10];
				sprintf(s, "%d", irreps(iIrrep).n);
				printf("  %s  ", centerTo(s, (fw+2)*irreps(iIrrep).roots.length()).chars());
				irreps.next(iIrrep);
			}
		}

		printf("\n       root");
		for ( INT i=0 ; i<fw-3 ; i++ )
			printf(" ");
		{
		Pix	iIrrep = irreps.first();
			while ( iIrrep )
			{
			Pix	iRoot = irreps(iIrrep).roots.first();
				while (iRoot)
				{
				char	s[10];
					sprintf(s, "%d", irreps(iIrrep).roots(iRoot));
					printf("%s", centerTo(s, fw+2).chars());
					irreps(iIrrep).roots.next(iRoot);
				}
				irreps.next(iIrrep);
				printf("    ");
			}
		}
		printf("\n\n");

	Pix	iIrrep = irreps.first();
		while ( iIrrep )
		{
		Pix	iRoot = irreps(iIrrep).roots.first();
		INT	iiRoot = 0;
			while ( iRoot )
			{
				if ( iiRoot==irreps(iIrrep).roots.length()/2 )
				{
				char	s[10];
					sprintf(s, "%d", irreps(iIrrep).n);
					printf("%5s ", s);
				}
				else
					printf("%5s ", " ");

				{
				char	s[10];
					sprintf(s, "%d", iiRoot+1);
					printf("%s", centerTo(s, fw+2).chars());
				}


			Pix	jIrrep = irreps.first();
				while ( jIrrep )
				{
				Pix	jRoot = irreps(jIrrep).roots.first();
					while ( jRoot )
					{
						if ( irreps(iIrrep).n<irreps(jIrrep).n || 
								irreps(iIrrep).n==irreps(jIrrep).n && 
								irreps(iIrrep).roots(iRoot)<=irreps(jIrrep).roots(jRoot) )
						{
						char	state[1000];
							sprintf(state, "I%dR%d_I%dR%d", 
								irreps(iIrrep).n, irreps(iIrrep).roots(iRoot),
								irreps(jIrrep).n, irreps(jIrrep).roots(jRoot));
							if ( map.contains(state) )
								printf("%9.5f  ", map[state].getProp(prop));
							else
								printf("    ---    ");
						}
						else
							printf("%9s  ", " ");

						irreps(jIrrep).roots.next(jRoot);
					}
					printf("    ");
					irreps.next(jIrrep);
				}
				printf("\n");
				irreps(iIrrep).roots.next(iRoot);
				iiRoot++;
			}
			printf("\n");
			irreps.next(iIrrep);
		}
		printf("\n\n\n");
	}
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

//	global options
String	IntPath("../ONEINT");
String	OrbPath("../SCFORB");
String	Mode("E_Dav2");
String	Unit("au");
String	Thresh;


//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

void	calcE(TMap & map, String prop)
{
INT	propCol, propBlock;
	getPropCol(Mode, propCol, propBlock);


SLList<TIrrep>	irreps;
	Map2Irreps(map, irreps);
	
	
INT	dirs = irreps.length();
const char	**dir = new const char *[dirs];

	{
	Pix i = irreps.first();
	INT	j = 0;
		while ( i )
		{
		char	buf[100];
			sprintf(buf, "%d", irreps(i).n);
			dir[j] = new char[strlen(buf)+1];
			strcpy((char *) dir[j], buf);
			irreps.next(i);
			j++;
		}
	}
	

INT	*nRootsInIrrep = new INT[dirs];
double	**e = new double * [dirs];

	for ( INT iIrrep=0 ; iIrrep<dirs ; iIrrep++ )
	{
	int	des[2];
		pipe(des);
	pid_t	pid = fork();
		if ( !pid ) //	child
		{
			close(des[0]);
			close(1);
			dup(des[1]);
			close(des[1]);
		String	command(getenv("DIESEL_EXE_DIR"));
			command += "/dr";
			chdir(dir[iIrrep]);
			execl("/bin/SH", "SH", "-c", command.chars(), NULL);

		}
		close(des[1]);
	FILE	*fin = fdopen(des[0], "r");
	GrepAwk	ga(fin);
		nRootsInIrrep[iIrrep] = 0;
		while ( ga.grep("root #") )
		{
			nRootsInIrrep[iIrrep]++;
			ga++;
		}
		ga.head();
		e[iIrrep] = new double[nRootsInIrrep[iIrrep]];
		
		for ( INT iRoot=0; iRoot<nRootsInIrrep[iIrrep] ; iRoot++ )
		{
			ga.grep("root #");
			ga += 3 ;
			for ( INT i=0 ; i<propBlock ; i++ )
			{
				ga.grep("threshold");
				ga++;
			}
			while ( !ga.illegal() )
			{
				if ( fabs((atof(ga.getWord(1))/1000-atof(Thresh))/atof(Thresh))<1e-4 )
				{
					e[iIrrep][iRoot] = atof(ga.getWord(propCol));
					break;
				}
				ga++;
			}
			if ( ga.illegal() )
			{
				cout << "Threshold \"" << Thresh << "\" not found." << endl;
				exit(1);
			}
		}
		fclose(fin);
		wait(0);
	}
	
	for ( INT iIrrep=0 ; iIrrep<dirs ; iIrrep++ )
		for ( INT iRoot=0 ; iRoot<nRootsInIrrep[iIrrep] ; iRoot++ )
			for ( INT jIrrep=iIrrep ; jIrrep<dirs ; jIrrep++ )
				for ( INT jRoot=(iIrrep==jIrrep) ? iRoot : 0 ; jRoot<nRootsInIrrep[jIrrep] ; jRoot++ )
				{
				char	buf[1000];
					sprintf(buf, "I%sR%d_I%sR%d",
						dir[iIrrep], iRoot+1, dir[jIrrep], jRoot+1);
					map[buf].append(prop, (iRoot==jRoot && iIrrep==jIrrep) ?
						convE(e[iIrrep][iRoot], Unit):
						convE(e[jIrrep][jRoot]-e[iIrrep][iRoot], Unit));
				}
	
	for ( INT i=0 ; i<dirs ; i++ )
		delete e[i];
	delete e;
	delete nRootsInIrrep;
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

void	calcProp(TMap & map, const char *op, const char *comp)
{
SLList<String>	ls;


{
	DIR	*dir = opendir("Densities");
	struct dirent	*entry;

		while ( (entry=readdir(dir)) )
		{
			if ( !strncmp("Density.dat.", entry->d_name, 12) &&
				Thresh==entry->d_name+strlen(entry->d_name)-Thresh.length() )

			{
			char	d[1000];
				sprintf(d, "Densities/%s", entry->d_name);
			String	S(d);
				ls.append(S);
			}
		}
		closedir(dir);
}

/*	n = scandir("Densities", &namelist, 0, alphasort);
String	args;
	if (n < 0)
		perror("scandir");
	else
		while(n--)
		{
			if ( !strncmp("Density.dat.", namelist[n]->d_name, 12) &&
				Thresh==namelist[n]->d_name+strlen(namelist[n]->d_name)-Thresh.length() )

			{
			char	d[1000];
				sprintf(d, "Densities/%s", namelist[n]->d_name);
			String	S(d);
				ls.append(S);
			}
		}
*/

//	cout << args << endl;


	int	des[2];
		pipe(des);
	pid_t	pid = fork();
		if ( !pid ) //	child
		{
			close(des[0]);
			close(1);
			dup(des[1]);
			close(des[1]);
		String	command(getenv("DIESEL_EXE_DIR"));
			command += "/prop";
		const char	*argsV[ls.length()+5+1];
			argsV[0] = "prop";
			argsV[1] = op;
			argsV[2] = comp;
			argsV[3] = IntPath;
			argsV[4] = OrbPath;
		Pix pix = ls.first();
		INT	i = 5;
			while ( pix )
			{
				argsV[i++] = ls(pix);
				ls.next(pix);
			}
			argsV[ls.length()+5] = NULL;

			execv(command.chars(), (char *const *) argsV);

		}
		close(des[1]);
	FILE	*fin = fdopen(des[0], "r");
	char	buf[1000];
		fgets(buf, 1000, fin);
		fgets(buf, 1000, fin);
	Pix pix = ls.first();
		for ( INT i=0 ; i<ls.length() ; i++ )
		{
			fgets(buf, 1000, fin);
		char	b[1000];
			strcpy(b, ls(pix));
			b[strlen(ls(pix))-Thresh.length()-1] = 0;
//				cout << (String(op) + String(" ") + String(comp)) << endl;
			map[b+13+strlen("Densities")].append(String(op) + String(" ") + String(comp), atof(buf));
			ls.next(pix);
		}
		fclose(fin);
		wait(0);
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

void	calcOsc(TMap &map)
{
	{
	Pix	pix = map.first();
		while ( pix )
		{
		char	c[100];
			for ( INT i=1 ; i<=3 ; i++ )
			{
				sprintf(c, "Mltpl1 %d", i);
			double	y = map.contents(pix).getProp(c);
				y = 2.0/3.0*y*y*map.contents(pix).getProp("E");
				if ( fabs(y)>1e-6 )
				{
					map.contents(pix).append("osc", y);
					break;
				}
			}
			map.next(pix);
		}
	}
}


//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------


int	main(INT argc, char **argv)
{
	cerr << "prettyProp (Part of DIESEL-MR-CI), " << VERSION << ", " << DATE << endl << endl;
	if ( argc<2 )
	{
		laber();
		exit(1);
	}
	
	argc--;
String*	arg = new String[argc];
INT	sortFlag = 0;
INT	oscFlag = 0;
	for ( INT i=0 ; i<argc ; i++ )
		arg[i] = argv[i+1];
		
	for ( INT i=0 ; i<argc ; i++ )
	{
		if ( arg[i]=="-unit" )
		{
			Unit = arg[i+1];
			arg[i] = arg[i+1] = "";
		}
		if ( arg[i]=="-mode" )
		{
			Mode = arg[i+1];
			arg[i] = arg[i+1] = "";
		}
		if ( arg[i]=="-thresh" )
		{
			Thresh = arg[i+1];
			arg[i] = arg[i+1] = "";
		}
		if ( arg[i]=="-intpath" )
		{
			IntPath = arg[i+1];
			arg[i] = arg[i+1] = "";
		}
		if ( arg[i]=="-orbpath" )
		{
			OrbPath = arg[i+1];
			arg[i] = arg[i+1] = "";
		}
		if ( arg[i]=="-sort" )
		{
			sortFlag = 1;
			arg[i] = "";
		}
		if ( arg[i]=="-osc" )
			oscFlag = 1;
	}

	if ( Thresh.length()==0 )
	{
		cerr << "error: threshold missing" << endl;
		exit(1);
	}

PropList	propList;
TMap	map(propList);

	calcProp(map, "Mltpl1", "1");
	calcProp(map, "Mltpl1", "2");
	calcProp(map, "Mltpl1", "3");

	for ( INT i=0 ; i<argc ; i++ )
	{
		if ( arg[i]=="-e" && !oscFlag )
		{
			calcE(map, "E");
			continue;
		}

		if ( arg[i]=="-prop" && !oscFlag )
		{
//			calcProp(map, arg[i+1], arg[i+2]);
			i += 2;
			continue;
		}	


		if ( arg[i]=="-osc" )
		{
			calcE(map, "E");



			calcOsc(map);
		}
	}
	for ( INT i=0 ; i<argc ; i++ )
	{
		if ( arg[i]=="-e" )
		{
			if ( !sortFlag )
				prettyPrint(Mode + String(" in ") + Unit + 
					String(", thresh=") + Thresh, map, "E");
			continue;
		}

		if ( arg[i]=="-prop" )
		{
			if ( !sortFlag )
				prettyPrint(String("operator: ") + arg[i+1] + 
					String(", component: ") + arg[i+2] + String(" in ") + Unit +
					String(", thresh=") + Thresh, map, 
						String(arg[i+1]) + String(" ") + String(arg[i+2]));
			i += 2;
			continue;
		}	


		if ( arg[i]=="-osc" )
		{
			if ( !sortFlag )
				prettyPrint(String("oscillator strength")+ String(" in ") + Unit +
					String(", thresh=") + Thresh, map, "osc");
					
//			GnuPlotE(map);
//			GnuPlotSpec(map);
		}
	}
	if ( sortFlag )
		prettyPrintSpec(Thresh, map, "");

delete[] arg;
}





// template instanciation

#include "../../../../lib/Container/AVLMap.cc"
template class AVLMap<String, double>;
template <> String*   AVLMap<String, double>::_target_item = NULL;     // add/del_item target
template <> AVLNode<String, double>* 
	AVLMap<String, double>::_found_node = NULL; // returned added/deleted node

template class AVLMap<String, PropList>;
template <> String*   AVLMap<String, PropList>::_target_item = NULL;     // add/del_item target
template <> AVLNode<String, PropList>* 
	AVLMap<String, PropList>::_found_node = NULL; // returned added/deleted node

#include "../../../../lib/Container/Map.cc"
template class Map<String, double>;
template class Map<String, PropList>;


#include "../../../../lib/Container/SLList.cc"
template class SLList<String>;
