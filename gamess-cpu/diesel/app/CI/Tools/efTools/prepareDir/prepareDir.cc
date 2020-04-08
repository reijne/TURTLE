#include <math.h>
//#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <iomanip>
#include <iostream>
#include <fstream>
//#include <sstream>
#include <string>
#include <list>
#include <vector>
#include <unistd.h>

using namespace std;

#include "../../../../../lib/Container/Nawk.h"
#include <dirent.h>


//global


string	getNewDir()
{
	char	bufa[200];
	getcwd(bufa,200);
	string	CurrentDir(bufa);
	
	//double	t;
	DIR	*dir = opendir(CurrentDir.c_str());
	struct dirent	*entry;
	INT iter=0;
		while ( (entry=readdir(dir)) )
		{
			string	EntryName(entry->d_name);
			if ((EntryName.find("iter")!=string::npos) )
			{
				iter++;
			}
		}
	closedir(dir);
	
	char	newDirbuf[200];
	sprintf(newDirbuf,"%s/iter%d/",CurrentDir.c_str(),iter++);
	
	string	newDir(newDirbuf);
	
	return newDir;
}


INT		cpFile(string File, string oldDir, string newDir)
{
	string	oldFile;
	oldFile = oldDir + "/" + File;
	string	newFile;
	newFile = newDir + "/" + File;
	string	Kommando("cp -f");
	Kommando += " " + oldFile + " " + newFile;
	system(Kommando.c_str());
    cout << Kommando << endl;
	//link(oldFile.c_str(),newFile.c_str());
	return 0;
}


INT		copyFiles(string newDir)
{
	char 	buf[200];
	getcwd(buf,200);
	string	CurrentDir(buf);
	
	string	BaseDir;
	BaseDir = CurrentDir + "/iter0/";
	
	
	mkdir(newDir.c_str(),0750);
	//	copy from CurrentDir
	cpFile("coord",CurrentDir,newDir);
	cpFile("gradinfo",CurrentDir,newDir);
	//	copy from BaseDir
	cpFile("basis",BaseDir,newDir);
	cpFile("auxbasis",BaseDir,newDir);
	cpFile("control",BaseDir,newDir);
	cpFile("mos",BaseDir,newDir);
	//	copy files for diesel;
	
	string	diesel_BaseDir;
	diesel_BaseDir = BaseDir + "CI/";
	
	string	diesel_newDir;
	diesel_newDir = newDir + "CI/";
	
	mkdir(diesel_newDir.c_str(),0750);
	cpFile("diesel.job",diesel_BaseDir,diesel_newDir);
	cpFile("grad.job",diesel_BaseDir,diesel_newDir);
		
	return 0;
}

INT		copyFiles(string newDir, INT Multi, INT Irrep)
{
	char 	buf[200];
	getcwd(buf,200);
	string	CurrentDir(buf);
	
	string	BaseDir;
	BaseDir = CurrentDir + "/iter0/";
	
	
	mkdir(newDir.c_str(),0750);
	//	copy from CurrentDir
	cpFile("coord",CurrentDir,newDir);
	cpFile("gradinfo",CurrentDir,newDir);
	//	copy from BaseDir
	cpFile("basis",BaseDir,newDir);
	cpFile("auxbasis",BaseDir,newDir);
	cpFile("control",BaseDir,newDir);
	cpFile("mos",BaseDir,newDir);
	//	copy files for diesel;
	
	string	diesel_BaseDir;
	diesel_BaseDir = BaseDir + "CI/";
	
	string	diesel_newDir;
	diesel_newDir = newDir + "CI/";
	
	mkdir(diesel_newDir.c_str(),0750);
	cpFile("diesel.job",diesel_BaseDir,diesel_newDir);
	cpFile("grad.job",diesel_BaseDir,diesel_newDir);
		
	
	// linken der COnfTrees und sel.in.all sel.out.all
	
	
	
	char	bbuf[200];
	sprintf(bbuf,"%s%d",diesel_newDir.c_str(),Multi);
	mkdir(bbuf,0750);
	//chdir(buf);
	sprintf(bbuf,"%s%d/%d",diesel_newDir.c_str(),Multi,Irrep);
	mkdir(bbuf,0750);
//	chdir(bbuf);
	
	char	MultIrrep[200];
	sprintf(MultIrrep,"%d/%d/",Multi,Irrep);
	string	CalcDir(diesel_BaseDir);
	CalcDir += string(MultIrrep);
	DIR	*dir = opendir(CalcDir.c_str());
	struct dirent	*entry;
	
//	char	jetztbuf[200];
//	getcwd(jetztbuf,200);
		cout << "copy files " << endl;
		cout << "from dir " << CalcDir << endl;
		cout << "to dir   " << string(bbuf) << endl;
		while ( (entry=readdir(dir)) )
		{
			//cout << entry->d_name << endl;
			//if( entry->d_name == "sel.in.all" ) cout << "***********************hurra"<<endl;
			string	EntryName(entry->d_name);
			if ((EntryName.find("ConfTree")!=string::npos) ||
				( EntryName == "sel.in.all" ) ||
				( EntryName == "sel.out.all" ) )
			{
			    string	FileName(CalcDir);
   				FileName += "/" + EntryName;
	    		char	hereDirbuf[200];
                        sprintf(hereDirbuf,"%s",bbuf);
// FD		    	hereDirbuf = bbuf;
                if (EntryName.find("ConfTree")!=string::npos)
				{
				//cout << "*****" << EntryName << endl;
			    	EntryName = string(hereDirbuf)+ "/" + EntryName;
				
				//cout << "ln -sf " << FileName << " " << EntryName << endl;
				//symlink(FileName.c_str(),EntryName.c_str());
				//cout << "ln -sf " << FileName << " " << EntryName << endl;
                //if (EntryName.find("ConfTree")!=string::npos)
				    symlink(FileName.c_str(),EntryName.c_str());
                }
                else
                {    
                    cout << "copy " << EntryName << endl;
			    	EntryName = string(hereDirbuf)+ "/" + EntryName;
                    string  Komm;
                    Komm = "cp -f " + FileName + " " + EntryName;
                    system(Komm.c_str());
                    //link(FileName.c_str(),EntryName.c_str());
                }
                    
			}
		}
	closedir(dir);
	
	// ende des linken
	
	return 0;
}


INT	executeNumGrad(string newDir, string Parameter)
{
	char buf[200];
	getcwd(buf,200);
	string	hereDir(buf);
	//mkdir(newDir.c_str(),0750);
	chdir(newDir.c_str());
		string	Kommando(getenv("DIESEL_EXE_DIR"));
		Kommando += "/intgrad ";
		Kommando += Parameter;
		system(Kommando.c_str());
	chdir(hereDir.c_str());
	return 0;
}

INT	excited(string newDir, string Parameter)
{
	char buf[200];
	getcwd(buf,200);
	string	hereDir(buf);
	//mkdir(newDir.c_str(),0750);
	chdir(newDir.c_str());
        chdir("CI");
			string	Kommando(getenv("DIESEL_EXE_DIR"));
			Kommando += "/excited ";
			Kommando += Parameter;
            Kommando += " > excited.out";
            cout << Kommando << endl;
			system(Kommando.c_str());
        chdir("../");
	chdir(hereDir.c_str());
	return 0;
}


INT	copyResults(string newDir)
{
	char buf[200];
	getcwd(buf,200);
	string	hereDir(buf);
	//mkdir(newDir.c_str(),0750);
	chdir(newDir.c_str());
		//cpFile("gradient.diesel",newDir,hereDir);
		char	KommA[200];
		sprintf(KommA,"cp -f %s/gradient.diesel %s/gradient.diesel",newDir.c_str(),hereDir.c_str());
		system(KommA);
		//cpFile("energy.diesel",newDir,hereDir);
		char	KommB[200];
		sprintf(KommB,"cp -f %s/energy.diesel %s/energy.diesel",newDir.c_str(),hereDir.c_str());
		system(KommB);
	chdir(hereDir.c_str());
	return 0;
}


INT readExcitedFile(string exFileName, 
                    INT const IrrepWanted, 
                    INT const RootWanted, 
                    double const fValueWanted, 
                    double const fDelta,
                    double & fValueOfRoot)
{
    ifstream  exFile(exFileName.c_str());
    if ( exFile )
    {
        Nawk  na(exFile);
        INT Nstates(0);
        INT RootWith_f(0);
        while( !na.illegal() )
        {
            na.grep("eV");
            Nstates++;
            if( !na.illegal() )
            {
                INT shift(0);
                cout << na.getLine() << endl;
                if ( 1 == Nstates ) 
                    shift = 1;
                string  RootIrrep(na.getWord(3+shift));
                string  fvalueString(na.getWord(7+shift));
                //cout << "RI= " << RootIrrep << " f= " << fvalueString << endl;
                size_t pos = RootIrrep.find("|");
                string  Rs = RootIrrep.substr(0,pos);
                INT     RootActual(atoi(Rs.c_str()));
                pos++;
                string  Is = RootIrrep.substr(pos);
                INT     IrrepActual(atoi(Is.c_str()));
                
                double  fValueActual(atof(fvalueString.c_str()));
                
//                 cout << "Root= " << RootActual 
//                      << " Irrep= " << IrrepActual 
//                      << " f= " << fValueActual << endl;
                if (   ( IrrepActual == IrrepWanted )
//                    && ( RootActual == RootWanted )
                    && ( fValueActual >= fValueWanted - fDelta )
                    && ( fValueActual <= fValueWanted + fDelta ) )
                {
                    RootWith_f++;
                    cout << "found state with irrep " << IrrepActual
                         << " and Fvalue in the interval ["
                         << fValueWanted - fDelta << " , "
                         << fValueWanted + fDelta << " ]" << endl;
//                    cout << "RootWith_f: " << RootWith_f << endl;
                    if ( RootWith_f == RootWanted ) 
                    {
                        cout << "**** Root: " << RootActual 
                             << " Irrep: " << IrrepActual << endl;
                        fValueOfRoot = fValueActual;
                        return RootActual;
                    }
                }
                na++;
            }
        }
        cout << "cannot find required root" << endl;
    }
    else
    {
        cerr << "ERROR: Problems with " << exFileName << endl;
    }
    
    return 0;
}




void	laber()
{
	cerr << "purpose:" << endl;
	cerr << "interfaces numerical gradient programm to ef.x" << endl;
	cerr << "usage:" << endl;
	cerr << "prepareDir < inputFile" << endl;
	cerr << endl;
	cerr << "MRCIextpol: 0: none, 1: Pey/Buenker, 2: MRMP2, 3: MRMP3" << endl;
	cerr << "lambda    : 0: lambda=1, 1: lambda=lambdalin" << endl;
	cerr << "Davidson  : 0: none, 1: Davidson 1, 2: Davidson 2" << endl;
}


INT	main()
{
//	cerr << "excited (Part of DIESEL-MR-CI), " << VERSION << ", " << DATE << endl << endl;
	cerr << "prepareDir (Part of DIESEL-MR-CI) "  << endl << endl;
	
	INT	 Mult(0) ;
	INT	 Irrep(0) ;
	INT	 Root(0) ; 
	string Thresh =  "1000";
	string MRCIextpol = "0"; 
	string lambda     = "0";
	string Davidson   = "0";
	string Geo_T   = Thresh     ;
	string Geo_MR  = MRCIextpol ;
	string Geo_lam = lambda     ;
	string Geo_Dav = Davidson   ; 
	
	bool shouldLink(true);

    bool    fValuesearch(false);
    double  fValue(0.25);
    double  deltaf(1.00);
    INT     GS_Root(1);
    INT     GS_Irrep(0);
    
    
	string  ParameterFileName("numgrad.in");
	cout << "reading parameters for intgrad from file " 
         << ParameterFileName << endl;
    ifstream  ParameterFile(ParameterFileName.c_str());
    if ( ParameterFile )
    {
        Nawk  na(ParameterFile);
        INT   NumberOfParameters(na.getNumberOfWords());
        string  Arguments(na.getLine());
        cout << Arguments << endl;

        if  (na.grep("f"))
        {
            fValuesearch = true;
            cout << "Mode fValue-Following switched on" << endl;
            
            // check for Groundstate
            if (na.grep("GS"))
            {
                string  s(na.getWord(NumberOfParameters));
                size_t  pos = s.find("GS=");
                pos += 3;
                s = s.substr(pos);
                pos = s.find("|");
                string  Rs = s.substr(0,pos);
                GS_Root  = atoi(Rs.c_str());
                pos++;
                string  Is = s.substr(pos);
                GS_Irrep = atoi(Is.c_str());
                
                NumberOfParameters--;
 
            }
            else
            {
                cout << "default: choosing Groundstate: "
                     << GS_Root << ".Root in Irrep " << GS_Irrep << endl;
            }
            na.head();

            if (na.grep("deltaf"))
            {
                string  s(na.getWord(NumberOfParameters));
                INT pos = s.find("-deltaf=");
                pos += 8;
                deltaf = atof((s.substr(pos)).c_str());
                cout << "choosing: deltaf=" << deltaf << endl;
                NumberOfParameters--;
            }
            else
            {
                
                cout << "defaults to deltaf=" << deltaf << endl;
            }
            na.head();
            string  s(na.getWord(NumberOfParameters));
            INT pos = s.find("-f=");
            pos += 3;
            fValue = atof((s.substr(pos)).c_str());
            NumberOfParameters--;
            
            
            cout << "fValue=" << fValue << "+" << deltaf << endl;
            deltaf = deltaf/2;
            fValue += deltaf;
        }

        na.head();        
        cout << "reading in " << NumberOfParameters << " parameters " << endl;
        if ((  7 == NumberOfParameters ) 
         || (  4 == NumberOfParameters ) 
         || ( 11 == NumberOfParameters ) 
         && ( !na.illegal() ) )
		{
            Mult           = atoi((na.getWord(1)).c_str());
            cout << "Mult  " << Mult << endl;
            Irrep          = atoi((na.getWord(2)).c_str());
            cout << "Irrep " << Irrep << endl;
			Root           = atoi((na.getWord(3)).c_str());
            cout << "Root  " << Root << endl;
			Geo_T = Thresh = na.getWord(4);
            cout << "Thresh " << Thresh << endl;
            if ( NumberOfParameters >= 7)
            {
    			Geo_MR     = MRCIextpol = na.getWord(5);
	    		Geo_lam    = lambda     = na.getWord(6);
		    	Geo_Dav    = Davidson   = na.getWord(7);
            }
            if (11 == NumberOfParameters)
            {
                cout << "still reading" << endl;
			    Geo_T      = na.getWord(8);
    			Geo_MR     = na.getWord(9);
	    		Geo_lam    = na.getWord(10);
		    	Geo_Dav    = na.getWord(11);
            }
        }
        else
        {
            cerr << "ERROR in prepareDir: wrong Number of Parameters" << endl;
        } 
    }
    ParameterFile.close();



	string	NewDir;
	NewDir = getNewDir();
	if (shouldLink)	
    {
        copyFiles(NewDir,Mult,Irrep);
        if ( fValuesearch )
        {
            copyFiles(NewDir,Mult,GS_Irrep);
        }
    }
	else 
    {
        copyFiles(NewDir);
    }
    
	char	buf[200];
    cout << "starting Program" << endl;
	sprintf(buf,"  %d %d %d %s %s %s %s %s %s %s %s",Mult, Irrep, Root, Thresh.c_str(), 
        MRCIextpol.c_str(), lambda.c_str(), Davidson.c_str(), Geo_T.c_str(), Geo_MR.c_str(), 
        Geo_lam.c_str(), Geo_Dav.c_str());
	string	Parameter(buf);
    cout << "Parameter" << Parameter << endl;
	
    if ( fValuesearch )
    {
        string  fParameter(Parameter);
        fParameter += " -b";
                    
        //  start base-Calculation
        executeNumGrad(NewDir, fParameter);    
//        return 0;
        //  start excited with correct input
        char    bufe[200];
// 		sprintf(bufe,"  diesel.out %d %s %s %s %s",Mult, Geo_T.c_str(),
//             Geo_MR.c_str(), Geo_lam.c_str(), Geo_Dav.c_str());
		sprintf(bufe,"  diesel.out %d %s %s %s %s %d %d",Mult, Geo_T.c_str(),
            Geo_MR.c_str(), Geo_lam.c_str(), Geo_Dav.c_str(), GS_Root, GS_Irrep);
		string	eParameter(bufe);
        excited(NewDir, eParameter);
            
        //  look up for correct fValue
        double  fValueOfRoot(0);
        string    excitedFileName;
        excitedFileName = NewDir + "/CI/excited.out";
        INT RootSearched = readExcitedFile(excitedFileName, Irrep, Root, fValue,  deltaf, fValueOfRoot);
        cout << "**** found Root " << RootSearched << " with fValue " << fValueOfRoot << endl;

        //  set correct Parameters for gradient calculation
        Root = RootSearched;
		sprintf(buf,"  %d %d %d %s %s %s %s %s %s %s %s",Mult, Irrep, Root, Thresh.c_str(),
            MRCIextpol.c_str(), lambda.c_str(), Davidson.c_str(), Geo_T.c_str(), Geo_MR.c_str(),
            Geo_lam.c_str(), Geo_Dav.c_str());
        Parameter = string(buf);
    }

	executeNumGrad(NewDir, Parameter);		
    copyResults(NewDir);
	
	return 0;
}

