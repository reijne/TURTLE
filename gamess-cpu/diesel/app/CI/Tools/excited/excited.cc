#include "../../../../config.h"
#include "../../../../VersionDate.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <math.h>

#include <iostream>
#include <fstream>
#include <iomanip>

#include "../../../../lib/Container/String.h"

#include "../../../../lib/Container/SLList.h"

#include "../../../../lib/Container/GrepAwk.h"

#include <string>
#include <list>
#include <map>
#include <algorithm>
#include "../../../../lib/QM/IntegralContainer/DensityMatrix.h"
#include "../../../../lib/QM/IO/Fortran/FortranFileIO.h"
#include <sys/wait.h>
#include <sys/types.h>
#include <cctype>

using namespace std;

String unit = "eV";

double	Erel;
bool	RI = false;	// flag for RI-Calculation

void	laber()
{
	cerr << "purpose:" << endl;
	cerr << "generate excitation spectrum" << endl;
	cerr << "usage:" << endl;
	cerr << "excite diesel.out mult thresh MRCIextpol lambda Davidson [root irrep]" << endl;
	cerr << endl;
	cerr << "MRCIextpol: 0: none, 1: Pey/Buenker, 2: MRMP2, 3: MRMP3" << endl;
	cerr << "lambda    : 0: lambda=1, 1: lambda=lambdalin" << endl;
	cerr << "Davidonson: 0: none, 1: Davidson 1, 2: Davidson 2" << endl;
}


double	convE(double e)
{
	if ( unit=="au" )
		return e;
	if ( unit=="eV" )
		return e*27.211;
	if ( unit=="kcal/mol" )
		return e*627.5;
	cerr << "unknown unit " << unit << endl;
	exit(1);
}


struct	TData 
{
//	TData() : mult(-1), irrep(-1), root(-1), E(0) {}
// 	TData(INT mult, INT irrep, INT root, double E = 0) :
// 	mult(mult), irrep(irrep), root(root), E(E) {}
	
	TData() : mult(-1), irrep(-1), root(-1), Mx(0), My(0), Mz(0), Osc(-1), E(0) {}
	TData(INT mult, INT irrep, INT root,
		double Mx, double My, double Mz, double Osc = 0, double E = 0) :
	mult(mult), irrep(irrep), root(root), Mx(Mx), My(My), Mz(Mz), Osc(Osc), E(E) {}
	
	INT mult;
	INT	irrep;
	INT	root;
	double	Mx, My, Mz;
	double	Osc;
	double	E;
//	string  label;
	SLList<String>	Wave;
	
	
	
	friend ostream & operator << (ostream &, const TData &);
};



	
inline bool operator < (const TData &a, const TData &b) 
{	return a.E<b.E;	}

inline bool operator == (const TData &a, const TData &b) 
{	return a.mult==b.mult && a.irrep==b.irrep && a.root==b.root;	}

// ostream & operator << (ostream &s, const TData &d)
// {
// 	s << setprecision(4) << setw(3) << d.mult << setw(3) << d.root << "|" << d.irrep 
// 		<< setw(15) << setprecision(9) << d.E << ", " << setw(8) << setprecision(2) << convE(d.E - Erel) << " eV   ";
// 	s << endl << endl;
// 	return s;
//}



ostream & operator << (ostream &s, const TData &d)
{
	s << setprecision(4) << setw(3) << d.mult << setw(3) << d.root << "|" << d.irrep 
		<< setw(10);
	if ( d.E!=0 )
		s << d.E << ", " << setw(8) << setprecision(2) << convE(d.E - Erel) << " eV   ";
	else
		s << "---" << ", " << setw(8) << "---" << " eV   ";
	if ( d.Osc!=-1 )
	{	
		s << setprecision(5) << setw(9) << d.Osc
		<< "  (" << setw(0) << d.Mx << setw(9) << d.My << setw(9) << d.Mz << ")";
	}
	s << endl << endl;
Pix j = d.Wave.first();
	while ( j )
	{
		s << d.Wave(j) << endl;
		d.Wave.next(j);
	}
	s << endl << endl << endl;
	return s;
}




template <class Key, class Value>
map<Value, Key> invert(map<Key, Value> const & m)
{
        map<Value, Key>	r;

	for ( typename map<Key, Value>::const_iterator i=m.begin() ; i!=m.end() ; ++i )
		r[i->second] = i->first;

	return r;
}


class	RIDensityMatrix
{
public:
//	RIDensityMatrix()
// 	:label(""), Energy(0), nmos(0), isdens(0), DensityFileName("")
// 	{}
// 	
// 	
// 	RIDensityMatrix(string label, double Energy, INT nmos, INT isdens, string DensityFileName)
// 	:label(label), Energy(Energy), nmos(nmos), isdens(isdens), DensityFileName(DensityFileName)
// 	{}
	
	RIDensityMatrix(INT multi, INT Lirrep, INT Lroot, INT Rirrep, INT Rroot, double LEnergy, double REnergy,string T)
	:Lirrep(Lirrep), Lroot(Lroot), Rirrep(Rirrep), Rroot(Rroot), multi(multi), Thresh(T.c_str()),isSwitched(false)
	{
		char buf[12];		
		sprintf(buf,"%da%d   %da%d   ", Lroot, Lirrep, Rroot, Rirrep);
                sprintf(label,"%s",buf);
//FD		label = buf;


		if (Lirrep > Rirrep)
		{
			INT t = Rirrep;
			Rirrep = Lirrep;
			Lirrep = t;
			t = Rroot;
			Rroot = Lroot;
			Lroot = t;
			isSwitched = true;
		}
		
		if ((Lirrep == Rirrep) && (Lroot > Rroot))
		{
			INT t = Rroot;
			Rroot = Lroot;
			Lroot = t;
			isSwitched = true;
		}

		if ( Lirrep == Rirrep && Lroot == Rroot )
		{
			isdens = 1;
			Energy = LEnergy;
		}
		else
		{
			isdens = 0;
			Energy = convE(LEnergy - REnergy);
		}
		
		//char Thresh[] = "1e-5";
		
		char dbuf[100];
		sprintf(dbuf,"%d/Densities/Density.dat.I%dR%d_I%dR%d.%s",multi,Lirrep,Lroot,Rirrep,Rroot,Thresh);
                sprintf(DensityFileName,"%s",dbuf);
//FD		DensityFileName = dbuf;
		//cout << "HHALLOO: suche densFile: " << DensityFileName << endl;
// 		DensityFileName = "" + multi + "Densities/" + ".I" + Lirrep + "R" + Lroot;
// 		DensityFileName += "_I" + Rirrep + "R" + Rroot + "." + Thresh;
		if (!DensityFileName)
		{
			cerr << "cannot find file " << DensityFileName << endl << endl;
			exit(-1);
		}
		nmos = 0;
	}
	
	
	
	// data for input
	INT		Lirrep, Lroot, Rirrep, Rroot;
	INT 		multi;
	
	// data for output
	char 		label[12];
	const char*		Thresh;
	double		Energy;
	INT		nmos;
	INT 		isdens;
	char		DensityFileName[100];
	void		writeToStream(FortranFileIO &);
	bool		isSwitched;
};


void    RIDensityMatrix::writeToStream (FortranFileIO &wdens)
{
	ifstream	fDens(DensityFileName);
//	cout << "reading DensityMatrix " << DensityFileName << endl;
	DensityMatrix	densM(fDens);
//	cout << "... success" << endl;
	nmos = densM.getMaxMO();
	
	//Energy = abs(Energy);
	
//	cout << " install DensityMatrix" << endl;
	wdens.write(label,12*sizeof(label[0]));
	wdens.write(&Energy,sizeof(Energy));
	wdens.write(&nmos,sizeof(nmos));
	wdens.write(&isdens,sizeof(isdens));
//	cout << "transforming DensityMatrix now" << endl;
	double*	dummy = new double[nmos*nmos];
//	cout << "start now" << endl;
	for (INT j = 1; j <= nmos; j++)
	{
		for (INT k = 1; k <= nmos; k++)
		{

			if (!(isSwitched))		
			{
				dummy[(j-1)*nmos+(k-1)] = densM.operator()(j,k);
				//dummy[(j-1)*nmos+(k-1)] = 0;
			}
			// Mit Doppeldreher
			// dummy[(k-1)*nmos+(j-1)] = densM.operator()(k,j);
			else
			{
			// Mit Einfachdreher
				//dummy[(k-1)*nmos+(j-1)] = densM.operator()(j,k);
				dummy[(j-1)*nmos+(k-1)] = densM.operator()(k,j);
				//dummy[(k-1)*nmos+(j-1)] = densM.operator()(k,j);
				//dummy[(k-1)*nmos+(j-1)] = 0;
			}
//			cout << j-1 << "|" << k << " ";
		}
	}
 	wdens.write(dummy,nmos*nmos*sizeof(dummy[0]));
	delete [] dummy;
	fDens.close();
//	cout << "Transformation finished." << endl;
}




INT	main(INT argc, char **argv)
{
	cerr << "excited (Part of DIESEL-MR-CI), " << VERSION << ", " << DATE << endl << endl;
	if ( argc!=7 && argc!=9 )
	{
		laber();
		exit(1);
	}
	argc--;

INT	multSel = atoi(argv[2]);

String	Thresh = argv[3];
INT MRCIExtPol = atoi(argv[4]); // 0=none, 1=Pey/Buenker, 2=MRMP2, 3=MRMP3
INT	lambdaOpt = atoi(argv[5]);
INT	Davidson = atoi(argv[6]);	// 0=none, 1=Davidson1, 2=Davidson2

INT	block = -1;
INT	col = -1;


	switch ( MRCIExtPol )
	{
	case 0:
		block = 1;
		col = 4;
		break;
		
	case 1:
		switch ( Davidson ) {
		case 0:
			block = 2;
			col = lambdaOpt ? 4 : 2;
			break;

		case 1:
			block = 3;
			col = lambdaOpt ? 3 : 2;
			break;

		case 2:
			block = 3;
			col = lambdaOpt ? 5 : 4;
			break;
		}
		break;
		
	case 2:
		switch ( Davidson ) {
		case 0:
			block = 2;
			col = lambdaOpt ? 7 : 5;
			break;

		case 1:
			block = 3;
			col = lambdaOpt ? 7 : 6;
			break;

		case 2:
			block = 3;
			col = lambdaOpt ? 9 : 8;
			break;
		}
		break;
		
	case 3:
		switch ( Davidson ) {
		case 0:
			block = 2;
			col = lambdaOpt ?10 : 8;
			break;

		case 1:
			block = 3;
			col = lambdaOpt ?11 :10;
			break;

		case 2:
			block = 3;
			col = lambdaOpt ?13 :12;
			break;
		}
		break;
		
	}




//-----------------------------------------------------------------------
	
ifstream	in(argv[1]);
GrepAwk	ga(in);
	

list<TData>	states;
map<INT, INT>	rootNr[8];

		// MOLCASRootDir is used to switch between the directories in RI-Calculation
		ga.grep("MOLCASRootDir");
		String	MOLCASRootDir(ga.getWord(3));
		
		
		
		// search FileFormats and look for RIFormat
		ga.grep("MOIntegralFileFormat");
		if (ga.getWord(3) == "RIFormat")
		{
			cerr << "Found Flag for RI-Type Calculation." << endl;
			cerr << "Properties will be calculated with proper by S.Grimme." << endl;
			RI = true;
			cerr << "Expecting the turbomol-files in MOLCASRootDir = " << MOLCASRootDir << endl;
		}
	while ( !ga.illegal() )
	{


		ga.grep("Multiplicity");
	INT	mult = atoi(ga.getWord(3));
		ga.grep("Irrep");
	INT	irrep = atoi(ga.getWord(3));		
	INT	root = 0;
//	INT	rootNr[8][1000];
		while ( !ga.illegal() )
		{
			if ( ga.getWord(1)=="reference" && ga.getWord(2)=="energy:" )
			{
				ga++;
				root++;
	
				ga.grep("/mH");
			
				while ( !ga.illegal() )
				{
					if ( fabs(atof(ga.getWord(1))-atof(Thresh)*1000)<1e-6 )
						break;
					ga++;
				}
				rootNr[irrep][root] = atoi(ga.getWord(5));
				
				for ( INT b=1 ; b<block ; b++ )
				{
					ga.grep("/mH");
					ga++;
				}
				
				while ( !ga.illegal() )
				{
					if ( fabs(atof(ga.getWord(1))-atof(Thresh)*1000)<1e-6 )
						break;
					ga++;
				}
				
			TData	data;
				data.mult = mult;
				data.irrep = irrep;
				data.root = root;
				data.E = atof(ga.getWord(col));
				
			char	diagOut[100];
				sprintf(diagOut, "%d/%d/diag.out.%s", mult, irrep, Thresh.chars());
			ifstream	fdiagOut(diagOut);
				if ( diagOut )
				{
				GrepAwk	gaDiag(fdiagOut);
					gaDiag.grep("wave function with ci^2/nSAFS>");
					while ( !gaDiag.illegal() )
					{
						gaDiag.grep("root #");
						if ( atoi(gaDiag.getWord(3))==rootNr[irrep][root] )
						{
							gaDiag += 3;
							while ( gaDiag.getNumberOfWords()>0 )
							{
								if ( atof(gaDiag.getWord(1))>=0.004 && gaDiag.getNumberOfWords()>1 )
									do
									{
									String	s(gaDiag.getLine());
										data.Wave.append(s);
									gaDiag++;
									} while ( gaDiag.getNumberOfWords()==1 );
								else
									gaDiag++;
							}
							
							break;
						}
						gaDiag++;
					}
				}


				if ( multSel==0 || mult==multSel )
				{
					//cout << "found new state" << endl;	
					states.push_back(data);
//					cout << "new state inlcuded in list." << endl;
				}
				
			}
			if ( ga.getWord(1)=="Multiplicity" )
				break;

			ga++;
		}
	
	}

//	cout << "sorting the states." << endl;
	states.sort();
	
INT	Irel;
INT	Rrel;
INT	ics = 1;

//	cout << "searching for ground-state" << endl;
	if ( argc==8 )
	{
		double	Mx = 0, My = 0, Mz = 0;
		Rrel = atoi(argv[7]);
		Irel = atoi(argv[8]);

	TData	data(multSel, Irel, Rrel, Mx, My, Mz);
		list<TData>::iterator i = find(states.begin(), states.end(), data);
		if ( i==states.end() )
		{
			cerr << "no such irrep/root" << endl;
			exit(1);
		}
		Erel = i->E;
//		cout << Erel << endl;

		for ( list<TData>::const_iterator j=states.begin() ; j!=i ; ++j )
			--ics;
	}
	else
	{
		Erel = states.begin()->E;
		Irel = states.begin()->irrep;
		Rrel = states.begin()->root;
	}



map<INT, INT>	invRootNr[8];
	for ( INT i=0 ; i<8 ; ++i )
		invRootNr[i] = invert(rootNr[i]);

//	cout << "reading of diesel.out etc finished" << endl;

if (!RI)
{	
	char	buf[100];
		sprintf(buf, "%d/prop.dat.%s", multSel, Thresh.chars());
	ifstream prop(buf);
		
		while ( prop )
		{
		INT	i1, r1, i2, r2;
		double	Mx, My, Mz;
		
			prop >> i1 >> r1 >> i2 >> r2 >> Mx >> My >> Mz;
			
			if ( i1==Irel && invRootNr[i1].find(r1)->second==Rrel )
			{
				TData	data(multSel, i2, invRootNr[i2].find(r2)->second, Mx, My, Mz);
			list<TData>::iterator i = find(states.begin(), states.end(), data);
				i->Mx = Mx;
				i->My = My;
				i->Mz = Mz;
				if ( i->E!=0 )
					i->Osc = 2.0/3.0 * (Mx*Mx + My*My + Mz*Mz) * (i->E - Erel);
			}

			if ( i2==Irel && invRootNr[i2].find(r2)->second==Rrel )
			{
				TData	data(multSel, i1, invRootNr[i1].find(r1)->second, Mx, My, Mz);
			list<TData>::iterator i = find(states.begin(), states.end(), data);
				i->Mx = Mx;
				i->My = My;
				i->Mz = Mz;
				if ( i->E!=0 )
					i->Osc = 2.0/3.0 * (Mx*Mx + My*My + Mz*Mz) * (i->E - Erel);
			}
		}
	//}

// 	char	buf[100];
// 		sprintf(buf, "%d/prop.dat.%s", multSel, Thresh.chars());
// 	ifstream prop(buf);
// 		
// 		while ( prop )
// 		{
// 		INT	i1, r1, i2, r2;
// 		double	Mx, My, Mz;
// 		
// 			prop >> i1 >> r1 >> i2 >> r2 >> Mx >> My >> Mz;
// 			cout << "found something" << endl;			
// 			if (( i1==Irel && rootNr[i1].find(r1)->second==Rrel ))
// 			{
// 				TData	data(multSel, i2, rootNr[i2].find(r2)->second, Mx, My, Mz);
// 			list<TData>::iterator i = find(states.begin(), states.end(), data);
// 				i->Mx = Mx;
// 				i->My = My;
// 				i->Mz = Mz;
// 				i->Osc = 2.0/3.0 * (Mx*Mx + My*My + Mz*Mz) * (i->E - Erel);
// 			}
// 
// 			if (( i2==Irel && rootNr[i2].find(r2)->second==Rrel ))
// 			{
// 				TData	data(multSel, i1, rootNr[i1].find(r1)->second, Mx, My, Mz);
// 			list<TData>::iterator i = find(states.begin(), states.end(), data);
// 				i->Mx = Mx;
// 				i->My = My;
// 				i->Mz = Mz;
// 				i->Osc = 2.0/3.0 * (Mx*Mx + My*My + Mz*Mz) * (i->E - Erel);
// 			}
// 		}

}



else
{
// 	{
// 		//cout << "deleting files from previous calculations" << endl;
// 		char	buf[100];
// 		sprintf(buf,"%s/mos.bin",MOLCASRootDir.chars());
// 		INT errn = unlink(buf);
// 		sprintf(buf,"%s/mrci.cidens",MOLCASRootDir.chars());
// 		INT errn2 = unlink(buf);
// 		unlink("mrci.cidens");
// 		unlink("mrci.prp");
// 	}
//	cout << "start to do ri-specific things" << endl;
	unit = "au";
	// Searching for roots, which are set to 0
//	cout << "deleting roots, which are set to zero" << endl;
	for ( list<TData>::iterator i = states.begin() ; i!=states.end() ; i++ )
	{
		INT RefIr=i->irrep, RefRo=i->root;
		INT ciRo=rootNr[RefIr].find(RefRo)->second;
		//cout << "RefRo= " << RefRo << " ciRo= " << ciRo << endl;
		if (ciRo == 0)
		{
			list<TData>::iterator j = i;
			i--;
			states.erase(j);
	//		cout << "deleting state with root set to 0" << endl;
		}
	}	
	
//	cout << "converting density-matrices to format for proper" << endl;
	INT ndens = (INT)(states.size());
	INT ic = 0; // for counting the densities
	{
		FortranFileIO wdens("mrci.cidens");
		wdens.write(&ndens,sizeof(ndens));
		ndens = 2 * ndens - 1;
//		cout << "start now" << endl;
		for ( list<TData>::const_iterator i = states.begin() ; i!=states.end() ; i++)
			{
//				cout << setw(3) << ic <<  ".  " << *i;
				INT RefIr=i->irrep, RefRo=i->root;
				INT ciRo=rootNr[RefIr].find(RefRo)->second;
				INT ciIr=RefIr;
//				cout << "initlizing density" << endl;
				RIDensityMatrix rdens(multSel, ciIr, ciRo, ciIr, ciRo, i->E, i->E, string(Thresh.chars()));
//				cout << "prepare writing of density-matrices" << endl;
				rdens.writeToStream(wdens);
//				cout << "...success" << endl;

//				cout << "icounter=" <<  ++ic << endl
//					 << "Erel=    " <<  Erel << endl
//					 << "i->E=    " <<  i->E << endl
//					 << "Thresh=  " << Thresh  << endl << endl;
	
			}
//		INT rn = 0;
		for ( list<TData>::const_iterator i = states.begin() ; i!=states.end() ; i++)
			{
				if (!(i->E == Erel && i->irrep == Irel && i->root ==  Rrel))
				{
					//cout << setw(3) << ic << ".  " << *i;
					INT RefIr=i->irrep;
					INT RefRo=i->root;
					INT ciRo=rootNr[RefIr].find(RefRo)->second;
					INT ciIr=RefIr;
	
					INT ciRrel=rootNr[Irel].find(Rrel)->second;
// 					INT irL=Irel,	//	Grundzustand
// 						roL=ciRrel,
// 						irR=ciIr,	//	angeregter Zustand
// 						roR=ciRo;

					RIDensityMatrix rdens(multSel, Irel, ciRrel, ciIr, ciRo, Erel, i->E, string(Thresh.chars()));
					rdens.writeToStream(wdens);

// 				cout << "ic=      " <<  ++ic << endl
// 					 << "rn=      " <<  ++rn << endl
// 					 << "Erel=    " << setprecision(8) <<  Erel << endl
// 					 << "i->E=    " << setprecision(8) <<  i->E << endl
// 					 << "Thresh=  " << Thresh  << endl << endl;
				}
			}
	}
	//wdens.close();


	cerr << "Happy Landing." << endl << endl;
	cerr << "All Density Infromation now transformend." << endl;
	cerr << "Starting turbomole programs to continue." << endl;
	cerr << "Using tmwfn. -now" << endl;
	
	{	// Erzeugen der mos.bin Datei mittels tmwfn
		int des[2];
		pipe(des);
		pid_t	pid = fork();
		if (!pid)	//child
		{
			close(des[0]);
			close(1);
			dup(des[1]);
			close(des[1]);

			
			// linking file mrci.cidens from dieselDir to MOLCASRootDir
			char	dieselDir[100];
			getcwd(dieselDir,100);
			char	densFileImDieselDir[120];
			sprintf(densFileImDieselDir,"%s/mrci.cidens",dieselDir);
			char	densFileImMOLCASRootDir[120];
			sprintf(densFileImMOLCASRootDir,"%s/mrci.cidens",MOLCASRootDir.chars());
			
			INT er1=symlink(densFileImDieselDir, densFileImMOLCASRootDir);

			chdir(MOLCASRootDir.chars());
			
			cout << "MOLCASRootDir: " << MOLCASRootDir <<  endl;			
			////string	command(getenv("TURBODIR"));
			////command += "/../Intel/tmwfn";
			string	command(getenv("DIESEL_EXE_DIR"));
			command += "/tmwfn";
	
			////execl(command.c_str() , "tmwfn", "mos", NULL);
			execl(command.c_str() , command.c_str(), NULL);
			//execl(command.c_str() , "tmwfn", NULL);
			cerr << "problems with Program tmwfn" << endl;
		}
		close(des[1]);
		wait(0);
	}
	


	
	
//	cout << endl << endl << "Starting props.x" <<endl;	
	{	// aufruf von props.x
//		INT des[2];
//		pipe(des);
		//pid_t	pid = fork();
		//if (!pid)
		{
// 			close(des[0]);
// 			close(1);
// 			dup(des[1]);
// 			close(des[1]);
			
			char	dieselDir[100];
			getcwd(dieselDir,100);
			chdir(MOLCASRootDir.chars());
		
			//string pfad("/home/jan/bin/diesel0322/diesel/lib/QM/app/CI/Tools/excited/props.x");
		
			string pfad(getenv("DIESEL_EXE_DIR"));
			pfad += "/proper";
			
			int des[2];
			pipe(des);
			pid_t	pida = fork();
			if (!pida)
			{
	 			close(des[0]);
	 			close(1);
	 			dup(des[1]);
	 			close(des[1]);
				//execl(pfad.c_str(), "props.x", "mos", "mrci", NULL);
				execl(pfad.c_str(), pfad.c_str(), "mos.bin", "mrci", NULL);
				//execl(pfad.c_str(), "proper", "mos.bin", "mrci", NULL);
				cerr << "problems with Program proper" << endl;
			}
			wait(0);
	 		close(des[1]);
			cerr << "Program proper finished." << endl;
			
			
			char	propFileImDieselDir[120],propFileImMOLCASRootDir[120];
			
			sprintf(propFileImDieselDir,"%s/mrci.prp",dieselDir);
			sprintf(propFileImMOLCASRootDir,"%s/mrci.prp",MOLCASRootDir.chars());
			
			INT er=link(propFileImMOLCASRootDir, propFileImDieselDir);
			chdir(dieselDir);
			//execl("/bin/cp/", "-f", propFileImMOLCASRootDir, propFileImDieselDir, NULL);
				
			
			// Der Aufruf "props.x mos mrci"
			// l"asst props den File mos.bin lesen und schreibt den File mrci.prp
			// Der File mrci.prp kann dann von den nachfolgenden Programmen 
			// mittels GrepAwk.h weiterverarbeitet werden
		}
		//cout << "es kann weitergehen." << endl;
		//wait(0);
// 		close(des[1]);

	//FILE	*fin = fdopen(des[0], "r");
	}
	
	
	// lesen von mrci.prp
	cerr << "reading informations from file mrci.prp: Transition-Densities" << endl;
	ifstream	prop("mrci.prp");
	while (prop)
	{	
		cerr << "file exists, start reading" << endl;
		GrepAwk	gp(prop);
		INT i2=0,r2=0;
		double	Mx,My,Mz,Osc;
		gp.grep("T - D E N S I T Y - O U T P U T");
		//cout << gp.getLine() << endl;
		while (!gp.illegal())
		{
			if (gp.getWord(1) == "transition")
			{
				String	Z = gp.getWord(4);
				string  Zustand(Z.chars());
				
// 				INT	apos = Zustand.find("a");
// 				string	ro(Zustand, 0,apos);
// 				r2 = atoi(ro.c_str());
// 				string  ir(Zustand, apos++, Zustand.length()+1);
// 				i2 = atoi(ir.c_str());
// 				cout << "irrep: " << i2 << " root: " << r2 << endl;
				
				unsigned INT apos = 0;
				for(INT j=0; j<2; ++j)
				{
					while(apos < Zustand.size() && !isdigit(Zustand.at(apos)))
						++apos;
					
					INT number=0;
					char c='0';
					while(apos < Zustand.size() && isdigit(Zustand.at(apos)))
					{
						number = 10*number + (INT)(Zustand.at(apos)) - (INT)(c);
						++apos;
					}
					switch(j)
					{
						case 0: r2 = number; break;
						case 1: i2 = number;
					}
				}
 				//cout << "irrep: " << i2 << " root: " << r2 << endl;				
				
				gp.grep("dipole (L)");
				Mx = atof((gp.getWord(3)));
				My = atof((gp.getWord(4)).chars());
				Mz = atof((gp.getWord(5)).chars());
				
				
				gp.grep("osc.str. (L)");
 				Osc = atof((gp.getWord(3)).chars());
// 				cout << " Mx= " << Mx; 
// 				cout << " My= " << My;
// 				cout << " Mz= " << Mz;
// 				cout << " Osc= " << Osc; 
// 				cout << endl;
				

				//if (!( i2==Irel && rootNr[i2].find(r2)->second==Rrel ))
				{
					TData	data(multSel, i2, invRootNr[i2].find(r2)->second, Mx, My, Mz);
				list<TData>::iterator i = find(states.begin(), states.end(), data);
					i->Mx = Mx;
					i->My = My;
					i->Mz = Mz;
					i->Osc = Osc;
					//i->Osc = 2.0/3.0 * (Mx*Mx + My*My + Mz*Mz) * (i->E - Erel);
					//cout << "Zustand wird mit Momenten versorgt." << endl << *i << endl;
				}
			}
			//else { cerr << "nothing found" << endl; }
			
			
			gp++;
		}
	}
	prop.close();
	ifstream	bprop("mrci.prp");
	while (bprop)
	{	
		GrepAwk	gp(bprop);
		INT i2=0,r2=0;
		double	Mx,My,Mz;
		cerr << "Searching for state-densities..." << endl;
		gp.grep("T - D E N S I T Y - O U T P U T");
		gp.grep("D E N S I T Y - O U T P U T");
		while (!gp.illegal())
		{
			if (gp.getWord(1) == "state")
			{
				String	Z = gp.getWord(2);
				string  Zustand(Z.chars());
				
// 				INT	apos = Zustand.find("a");
// 				string	ro(Zustand, 0,apos);
// 				r2 = atoi(ro.c_str());
// 				string  ir(Zustand, apos++, Zustand.length()+1);
// 				i2 = atoi(ir.c_str());
// 				cout << "irrep: " << i2 << " root: " << r2 << endl;
				
				unsigned INT apos = 0;
				for(INT j=0; j<2; ++j)
				{
					while(apos < Zustand.size() && !isdigit(Zustand.at(apos)))
						++apos;
					
					INT number=0;
					char c='0';
					while(apos < Zustand.size() && isdigit(Zustand.at(apos)))
					{
						number = 10*number + (INT)(Zustand.at(apos)) - (INT)(c);
						++apos;
					}
					switch(j)
					{
						case 0: r2 = number; break;
						case 1: i2 = number;
					}
				}
 				//cout << "irrep: " << i2 << " root: " << r2 << endl;				
				
				gp.grep("dipole");
				Mx = atof((gp.getWord(3)));
				My = atof((gp.getWord(4)).chars());
				Mz = atof((gp.getWord(5)).chars());
				
				
// 				gp.grep("osc.str. (L)");
// 				double		Osc = atof((gp.getWord(3)).chars());
// 				cout << " Mx= " << Mx; 
// 				cout << " My= " << My;
// 				cout << " Mz= " << Mz;
// 				cout << " Osc= " << Osc; 
// 				cout << endl;
				

				map<INT, INT>	invRootNr[8];
				for ( INT i=0 ; i<8 ; ++i )
					invRootNr[i] = invert(rootNr[i]);
				if (( i2==Irel && rootNr[i2].find(r2)->second==Rrel ))
				{
					TData	data(multSel, i2, invRootNr[i2].find(r2)->second, Mx, My, Mz);
				list<TData>::iterator i = find(states.begin(), states.end(), data);
					i->Mx = Mx;
					i->My = My;
					i->Mz = Mz;
					//i->Osc = 0;
					i->Osc = 2.0/3.0 * (Mx*Mx + My*My + Mz*Mz) * (i->E - Erel);
					//cout << "Zustand wird mit Momenten versorgt." << endl << *i << endl;
				}
			}
			
			
			gp++;
		}
	}
	
	unlink("mrci.cidens");
	
unit = "eV";	
}	// end of else


	//cout << "Ausgabe beginnt. " << endl;
	cout.setf(ios::fixed);
	
INT	ic = ics;

	cout << "   **********   Threshold=" << Thresh << "   **********" << endl << endl;
	for ( list<TData>::const_iterator i = states.begin() ; 
				i!=states.end() ; i++ , ic++ )
	{
		if ( ic==1 )
			cout << "GS--> ";
		else
			cout << "      ";
		cout << setw(3) << ic << ".  " << *i;
		if ( ic==1 )
			cout << endl;
	}
	
} 




#include "../../../../lib/Container/SLList.cc"


template class SLList<String>;
template class SLList<TData>;
