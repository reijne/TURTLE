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

#include "../../../../../lib/Container/Nawk.h"
#include <dirent.h>

using namespace std;

// global **************
INT		Block,Column;
INT		Geo_Block,Geo_Column;
    bool    selMinusR ;
    bool    Geo_selMinusR;


void	executeProgram(string Directory, string Programm);
void	turbocalc();
void	dieselcalc();
// ****************

class   intCoord
{
public:
    intCoord();
    intCoord(INT _Nu, string _Var, double _Stoerung, string _efInput);

    INT     Nummer;
    string  Name;
    double  Betrag;
    string  efInput;
    double  Energie;
    double  Gradient;
};


intCoord::intCoord()
:Nummer(0)
,Name("")
,Betrag(1)
,efInput("")
,Energie(0)
,Gradient(0)
{}

intCoord::intCoord(INT _Nu, string _Name, double _Betrag, string _efInput)
:Nummer(_Nu)
,Name(_Name)
,Betrag(_Betrag)
,efInput(_efInput)
,Energie(0)
,Gradient(0)
{}

ostream&    operator<< (ostream& s, intCoord iC)
{
    s << iC.efInput << endl;
    return  s;
}

enum	workMode {prepare, calc, collectResults};


class   icListe
{
public:
    icListe();
    icListe(string _TMBaseDir, string _DieselDir, string _CalcDir
                    ,INT _Irrep,INT _Multi,INT _Root,string _Thresh);
    icListe(string _TMBaseDir
                    ,INT _Irrep,INT _Multi,INT _Root,string _Thresh);


    INT     readInfoFile(string InfoFileName);
    INT		calcBase();
//    void    calcGradient();
    bool	getTM(string File);
    bool	makeTurboDir(const intCoord& iC);
    bool	makeDieselDir();
    INT     writeCoordFile(const intCoord& iC);
    INT     writeGradientFile();
    INT     calcGradient();
    bool	navigateFileSystem(workMode wMode);

    bool            doublePrecision;
    intCoord        base;
    list<intCoord>  VarListe;
    
    string  (InfoFileName);
    
    string  TMBaseDir;
    string  DieselDir;
    string  CalcDir;
    INT     Irrep;
    INT     Multi;
    INT     Root;
    string  Thresh;
    double  GradientNorm;
    double  BaseEnergy;
    vector<double>  Gradient;
};



icListe::icListe()
:doublePrecision(true)
,InfoFileName("gradinfo")
,TMBaseDir()
,DieselDir()
,CalcDir()
,Irrep(0),Multi(1),Root(1),Thresh("1000")
{	
//	delta_xyz = 0.010;
}

icListe::icListe(string _TMBaseDir, string _DieselDir, string _CalcDir
                    ,INT _Irrep,INT _Multi,INT _Root,string _Thresh)
:doublePrecision(true)
,InfoFileName("gradinfo")
,TMBaseDir(_TMBaseDir)
,DieselDir(_DieselDir)
,CalcDir(_CalcDir)
,Irrep(_Irrep),Multi(_Multi),Root(_Root),Thresh(_Thresh)
{
//	delta_xyz = 0.010;
	
	cout << "setting defaults:" << endl;
	cout << "TMBaseDir:   " << TMBaseDir<< endl;
	cout << "DieselDir:   " << DieselDir<< endl;
	cout << "CalcDir:     " << CalcDir<< endl<< endl;
}

icListe::icListe(string _TMBaseDir
                    ,INT _Irrep,INT _Multi,INT _Root,string _Thresh)
:doublePrecision(true)
,InfoFileName("gradinfo")
,TMBaseDir(_TMBaseDir)
,Irrep(_Irrep),Multi(_Multi),Root(_Root),Thresh(_Thresh)
{
	DieselDir = TMBaseDir + "/CI/";
	char	buff[200];
	sprintf(buff,"%s/%d/%d/",DieselDir.c_str(),Multi,Irrep);
	CalcDir = string(buff);
//	delta_xyz = 0.010;
	
	cout << "setting defaults:" << endl;
	cout << "TMBaseDir:   " << TMBaseDir<< endl;
	cout << "DieselDir:   " << DieselDir<< endl;
	cout << "CalcDir:     " << CalcDir<< endl<< endl;
}


INT icListe::readInfoFile(string InfoFileName)
{
    doublePrecision = true;

    ifstream    InfoFile(InfoFileName.c_str());
    Nawk    na(InfoFile);
    INT     Zaehler(0);
    while( !na.illegal() )
    {
        string  act;
        act = na.getWord(1);
        cout << "1.Wort = " << act << endl;
        if (act == "var" || act == "base")
        {
            string  _Name;
            _Name   = na.getWord(1) + na.getWord(2);
            double  _Betrag;
            _Betrag = double(atof((na.getWord(3)).c_str()));
            na++;
            act = na.getWord(1);
            INT Zeile(0);
            string  efText;
            while (act != "var" 
                  && (!na.illegal()) 
                  && (Zaehler < 1000) )
            {
                Zeile++;
                if ( 1 == Zeile)
                {
                    // entferne die keywords DIESEL und INTERN 
                    // und fuege statt dessen 0scf als keyword ein
                    INT nwords = na.getNumberOfWords();
                    string  ZeileEins;
                    ZeileEins = " ";
                    for (INT i=1; i<=nwords; i++)
                    {
                        string  st = na.getWord(i);
                        if ( (st != "DIESEL") 
                          && (st != "INTERN") )
                        {
                            ZeileEins += st + " ";
                        }
                    }
                    ZeileEins += "0scf";
                    cout    << "1: " << ZeileEins << endl;
                    efText = ZeileEins + "\n";
                }
                else
                {
                    cout    << Zeile << ": " << na.getLine() << endl;
                    efText += na.getLine() + "\n";
                }
                na++;
                if(!na.illegal()) 
                    act = na.getWord(1);
            }
//             icoord.efInput = efText;
//             icoord.Name = _Name;
//             icoord.Betrag = _Betrag;
//             icoord.Nummer = Zaehler;
            intCoord    icoord(Zaehler, _Name, _Betrag, efText);
            if (0 <= Zaehler)
                VarListe.push_back(icoord);
            else
                base = icoord;
            
        }
        Zaehler++;
        cout << "Habe neue Auslenkung gefunden." << endl;
        cout << "Infos werden in File-Nr."<< Zaehler << " geschrieben" << endl;
        //na++;
    }
    cout << "Dimension des Gradienten: " << VarListe.size() << endl;
    Gradient = vector<double>(VarListe.size(),0);    
    
    cout << "Einlesevorgang beendet." << endl;
    cout << "Es wurden insgesamt " << Zaehler 
         << " Auslenkungen registriert." << endl;
    
    return 0;
}

INT		icListe::calcBase()
{
    cout << "starting turbomole-calculation for base-geometry" << endl;
	turbocalc();
	chdir(DieselDir.c_str());
        // vielleicht sel -r Option starten.
        if ( selMinusR )
        {
     		char	buf[200];
	    	sprintf(buf,"%d/%d/",Multi,Irrep);
		    mkdir(buf,0750);
    		chdir(buf);
            {
                // pruefe ob sel.out.all doch schon berechnet wurde
//                 bool    sel_schon_da(false);
//                 {
//         		    ifstream sel_out_all("sel.out.all");
//         		    if (sel_out_all)
//                         sel_schon_da = true;
//                 }
//                 if ( !sel_schon_da )
//                 {
                    // linking fort.31 Directory 
                    string  fort31("../../../fort.31");
                    symlink (fort31.c_str(),"fort.31");
                    // recalculate EpsteinNesbet Perturbation Sums
                    cout << "recalculating EpsteinNesbet Perturbation Sums" << endl;
                    string  sel(getenv("DIESEL_EXE_DIR"));
                    sel += "sel -r < sel.in.all > sel.out.all";
                    system(sel.c_str());
//                 }
//                 else
//                 {
//                     cout << "recalcuation not necessary" << endl;
//                 }
                
            }    
            chdir(DieselDir.c_str());
        }
        else
        {
            cout << "Pertubation Sums are not needed - using old sel.out.all file" << endl;
        }
        
		dieselcalc();
		char	buf[200];
		sprintf(buf,"%d/%d/",Multi,Irrep);
		mkdir(buf,0750);
		chdir(buf);
		
		{
			char	Parameter[200];
			{
				char	b[200];
				getcwd(b,200);
				cout << "now in directory " << b << endl;
			}
			
			sprintf(Parameter,"dr -w -h -T=%s -R=%d -C=%d,%d",Thresh.c_str(),Root,Block,Column);
			string	dr(getenv("DIESEL_EXE_DIR"));
			dr += string(Parameter) + " > dr.out";
	
			system(dr.c_str());
			ifstream	cfile("dr.out");
			cfile >> BaseEnergy;
			cout << "reading BaseEnergy: " << BaseEnergy << endl;
		}
		chdir("../../");
	chdir(TMBaseDir.c_str());
	
	return 0;
}

INT     icListe::writeCoordFile(const intCoord& iC)
{
    //cout << "not yet complete" << endl;
    // iC muss noch in einen File ef.in geschrieben werden
    string  efin("ef.in");
    ofstream    efinFile(efin.c_str());
    if ( efinFile )
    {
        efinFile << iC << endl;
        efinFile.close();
        string  efx(getenv("DIESEL_EXE_DIR"));
        efx += "ef.x < " + efin;
        system(efx.c_str());
        return 0;
    }
    cerr << "ERROR: cannot write coord File" << endl;
    return 1;
}

INT     icListe::writeGradientFile()
{
    string      GradFileName("gradient.diesel");
    ofstream    GradFile(GradFileName.c_str());
    
    if ( GradFile )
    {
        cout << "writing Gradient into File" << endl;
        for (size_t i = 0; i < VarListe.size(); i++)
        {
            //LI->irgendwas;
            cout      << "var" << i << "= " 
                      << setprecision(8) << setw(12) << Gradient[i] << endl;
            GradFile  << setprecision(8) << setw(12) << Gradient[i] << endl;
        }
        //return 0;
    }
    else
    {
        cerr << "ERROR: cannot write File " << GradFileName << endl;
        return 1;
    }

	// Energie in seperaten File schreiben
	string		en("energy.diesel");
	//fn += "." + Thresh;
	ofstream	efile(en.c_str());
	if ( efile )
	{
		efile << setprecision(8) << setw(12) << BaseEnergy << endl;
		cout  << "**************************BaseEnergy: " 
              << setprecision(8) << setw(12) << BaseEnergy << endl;
	}
    else
    {
        cerr << "ERROR: cannot write File " << en << endl;
        return 1;
    }

    return 0;
}


INT     icListe::calcGradient()
{
    cout << "not yet ready" << endl;
//    list<intCoord>::iterator    LI;
    INT counter(0);
    for (list<intCoord>::iterator LI = VarListe.begin(); LI != VarListe.end(); LI++)
    {
        cout << "calculating gradient component for var" << ++counter << ": ";
        LI->Gradient = 1/(LI->Betrag) * ( LI->Energie - BaseEnergy );
        cout << setprecision(8) << setw(12) << LI->Gradient << "= 1/(h="
             << setprecision(8) << setw(12) << LI->Betrag << ") * ( " 
             << setprecision(8) << setw(12) << LI->Energie
             << " - " 
             << setprecision(8) << setw(12) << BaseEnergy << " )" << endl;
    }
    
    if ( doublePrecision )
    {
        cout << "calculation for doublePrecision" << endl;
        INT     i(0);
        bool    first(true);
        for (list<intCoord>::iterator LI = VarListe.begin(); LI != VarListe.end(); LI++)
        {
            //LI->irgendwas;
            if ( !first )
            {
                list<intCoord>::iterator LJ = LI;
                LJ--;
                Gradient[i] = 0.5 * ( LI->Gradient + LJ->Gradient );
                cout << "Gradient[" << i << "]= " 
                     << setprecision(8) << setw(12)  << Gradient[i] 
                     << "= 0.5 * (" << LI->Gradient << " + "
                     << LJ->Gradient << " )" << endl;
                i++;
                first = true;
            }
            else
                first = false;
        }
    }
    else
    {
        INT i(0);
        for (list<intCoord>::iterator LI = VarListe.begin(); LI != VarListe.end(); LI++ )
        {
            //LI->irgendwas;
            Gradient[i] = LI-> Gradient;
            cout << "Gradient[" << i << "]= " << Gradient[i] << endl;
            i++;
        }
        
    }
    
    return 0;
}


bool	icListe::getTM(string File)
{
	string	Original;
	Original = TMBaseDir + "/" + File;
	char	hereDirBuf[200];
	getcwd(hereDirBuf,200);
	string	Kopie(hereDirBuf);
	Kopie += "/" + File;
	string	Kommando("cp -f");
	Kommando += " " + Original + " " + Kopie;
	system(Kommando.c_str());
	//link(Original.c_str(),Kopie.c_str());
	return	true;
}


bool	icListe::makeTurboDir(const intCoord& iC)
{
	writeCoordFile(iC);
	getTM("basis");
	getTM("auxbasis");
	getTM("control");
	getTM("mos");
	//writeCoordFile();
	return true;
}

bool	icListe::makeDieselDir()
{
	char	bufa[200];
	getcwd(bufa,200);
	string	CurrentDir(bufa);
	
	mkdir("CI",0750);
	chdir("CI");
	
	string	gradJob = "/grad.job";
	getcwd(bufa,200);
	string	Kopie_gradJob = string(bufa) + gradJob;
	string	Original_gradJob = DieselDir + gradJob;
	// hier macht link eigentlich nix b"oses
	link(Original_gradJob.c_str(),Kopie_gradJob.c_str());
	
	char	buf[200];
	sprintf(buf,"%d",Multi);
	mkdir(buf,0750);
	chdir(buf);
	sprintf(buf,"%d",Irrep);
	mkdir(buf,0750);
	chdir(buf);
	
	
	//double	t;
	DIR	*dir = opendir(CalcDir.c_str());
	struct dirent	*entry;
	
	char	jetztbuf[200];
	getcwd(jetztbuf,200);
		cout << "copy files " << endl;
		cout << "from dir " << CalcDir << endl;
		cout << "to dir   " << string(jetztbuf) << endl;
		while ( (entry=readdir(dir)) )
		{
			//cout << entry->d_name << endl;
			//if( entry->d_name == "sel.in.all" ) cout << "***********************hurra"<<endl;
			string	EntryName(entry->d_name);
			if ((EntryName.find("ConfTree")!=string::npos) ||
				( EntryName == "sel.in.all" ) ||
				( EntryName == "sel.out.all" ) )
			{
				//cout << "*****" << EntryName << endl;
       			string	FileName(CalcDir);
	    		FileName += "/" + EntryName;
		    	char	hereDirbuf[200];
			    getcwd(hereDirbuf,200);
                if (!( EntryName == "sel.out.all" ) || !(Geo_selMinusR) )
				{
				    EntryName = string(hereDirbuf)+ "/" + EntryName;
				
                    symlink(FileName.c_str(),EntryName.c_str());
				    cout << "ln -sf " << FileName << " " << EntryName << endl;
                }
                else
                {
//                     cout << "copy " << EntryName << endl;
// 			    	EntryName = string(hereDirbuf)+ "/" + EntryName;
//                     string  Komm;
//                     Komm = "cp -f " + FileName + " " + EntryName;
//                     system(Komm.c_str());
                    cout << "do not need to copy " << EntryName << endl;
                    
                }
//                 else
//                     link(FileName.c_str(),EntryName.c_str());
			}
		}
	closedir(dir);
	
	chdir(CurrentDir.c_str());
	return true;
}

//////////////////////////////////////
bool	icListe::navigateFileSystem(workMode wMode)
{	
	INT counter=0;
    cout << "Looping over " << VarListe.size() << " Variables." << endl;
	for (list<intCoord>::iterator IL(VarListe.begin()); IL != VarListe.end(); IL++, counter++)
	{
        string  DirName;
        DirName = IL->Name;
        char    hereDirbuf[200];
	    getcwd(hereDirbuf,200);
        


		// do special stuff: mkdir; chdir; cp und ln files
		switch(wMode)
		{
		case(prepare):
			{
			cout << "Perturbation-Nr " << counter << ": " << DirName << endl; 
			cout << "make new directory " << DirName << endl;
			mkdir(DirName.c_str(),0750);
			chdir(DirName.c_str());
				makeTurboDir(*IL);
				makeDieselDir();
			chdir("../");
            chdir(hereDirbuf);
			}
			break;
		case(calc):
			{
			chdir(DirName.c_str());
				turbocalc();
				chdir("CI/");
                    // vielleicht sel -r Option starten.
//                     if ( Geo_selMinusR )
//                     {
//                         char    buf[200];
// 	                    sprintf(buf,"%d/%d/",Multi,Irrep);
// 	                    mkdir(buf,0750);
//                         chdir(buf);
//                         {
//                             cout << "recalculating EpsteinNesbet Perturbation Sums" << endl;
//                             string  sel(getenv("DIESEL_EXE_DIR"));
//                             sel += "sel -r < sel.in.all > sel.out.all";
//                             system(sel.c_str());
//                         }
//                         chdir("../../");
//                     }
                    if ( Geo_selMinusR )
                    {
     	     	    	char	buf[200];
	    	    		sprintf(buf,"%d/%d/",Multi,Irrep);
					    mkdir(buf,0750);
    	    	    	chdir(buf);
                        {
                            // pruefe ob sel.out.all doch schon berechnet wurde
                            bool    sel_schon_da(false);
                            {
                    		    ifstream sel_out_all("sel.out.all");
                    		    if (sel_out_all)
                                    sel_schon_da = true;
                            }
                            if ( !sel_schon_da )
                            {
                                // linking fort.31 Directory 
                                string  fort31("../../../fort.31");
                                symlink (fort31.c_str(),"fort.31");
                                
                                cout << "recalculating EpsteinNesbet Perturbation Sums" << endl;
                                string  sel(getenv("DIESEL_EXE_DIR"));
                                sel += "sel -r < sel.in.all > sel.out.all";
                                system(sel.c_str());
                            }
                            else
                            {
                                cout << "recalcuation not necessary" << endl;
                            }
 
                        }
                        chdir("../../");
                    }
                    else
                    {
                        cout << "Pertubation Sums are not needed - using old sel.out.all file" << endl;
                    }
                    
                    
					//dieselcalc();
					system("./grad.job");
					char	bufa[200];
					sprintf(bufa,"%d/%d",Multi,Irrep);
					chdir(bufa);
						string	dr(getenv("DIESEL_EXE_DIR"));
						dr += "/dr";

						// Parameter f"ur dr einstellen
						// kann eventuell fr"uher geschehen
						// siehe excited
						
						char	Parameter[200];
						sprintf(Parameter," -w -h -T=%s -R=%d -C=%d,%d",Thresh.c_str(),Root,Geo_Block,Geo_Column);
						dr += string(Parameter) + " > dr.out";

						cout << "getting diesel results: " << dr << endl;
                        system(dr.c_str());

					chdir("../../");
				chdir("../");
			chdir("../");
            chdir(hereDirbuf);
			}
			break;
		case(collectResults):
			{
			chdir(DirName.c_str());
				chdir("CI/");
					char	bufa[200];
					sprintf(bufa,"%d/%d",Multi,Irrep);
					chdir(bufa);
						cout << "reading Gradient Information" << endl;
						ifstream drFile;
						drFile.open("dr.out");
						if (!drFile)
						{
							cerr << "ERROR: cannot read File dr.out" << endl;
						}
						else
						{
							double	r;
							drFile	>> r;
							IL->Energie = r;
							//drFile	>>  Energien[updown][xyz];
							cout << "Energy for Perturbation " 
                                 << IL->Name << ": " << IL->Energie << endl;
						}
					chdir("../../");
				chdir("../");
			chdir("../");
            chdir(hereDirbuf);
			}
			break;
		default:
			cerr << "ERROR: wrong workMode" << endl;
		}
	}
	return true;
}

void	executeProgram(string Directory, string Programm)
{
	string	command;
	command = Directory + "/" + Programm;
	string	output;
	output = Programm + ".out";
//	INT	des[2];
//	pipe(des);
//	pid_t	pid = fork();
//	if (!pid)
//	{
//		close(des[0]);
//		close(1);
//		dup(des[1]);
//		close(des[1]);
	    cout << "starting " << Programm << endl;
    	
//		execl(command.c_str(), Programm.c_str(), (string(" > ")).c_str() , " ", output.c_str(), NULL);
//        output = " > " + output;
//		execl(command.c_str(), Programm.c_str(), output.c_str(), NULL);
		string Kommando;
		Kommando = command + " > " + output;
        cout << "system(" << Kommando << ");" << endl;
		system(Kommando.c_str());
//		cout << "ERROR: !!! problems with " << Programm << "!!!" << endl;
//        exit(1);
//    }
//    close(des[1]);
//    wait(0);

	
	return;
}
	
void	turbocalc()
{
	// turbomol-Rechnung
	{
		ifstream dscf_out("dscf.out");
		if (dscf_out)
		{
			cout << "turbomole-calculation allready done" << endl;
			cout << "directly starting with CI-Calculation" << endl;
			return;
		}
	}
	
	string	TURBODIR(getenv("TURBODIR"));
	TURBODIR += "/../Intel";
	cout << "TURBODIR: " << TURBODIR << endl;
	executeProgram(TURBODIR,"dscf");
	
	executeProgram("/u2/jan/bin/","ritraf");
	executeProgram("/u2/jan/bin/","oneint");
	mkdir("fort.31",0750);
	char	hereDirbuf[200];
	getcwd(hereDirbuf,200);
	string	hereDir(hereDirbuf);
	string	bkji;
	bkji = hereDir + "/" + "bkji";
	symlink(bkji.c_str(),"fort.31/bkji");
	string	oneint;
	oneint = hereDir + "/" + "oneint";
	symlink(oneint.c_str(),"fort.31/oneint");
}	
	
void	dieselcalc()
{
	// diesel-Rechnung
    cout << "starting diesel" << endl;
	system("./diesel.job");
}


void	laber()
{
	cerr << "purpose:" << endl;
	cerr << "calculates numerical gradient in internal coords" << endl;
	cerr << "usage:" << endl;
	cerr << "intgrad Mult Irrep Root Thresh [MRCIextpol lambda Davidson]" << endl;
//	cerr << "intgrad Mult Irrep Root Thresh [MRCIextpol lambda Davidson]  [Geo_T [Geo_MR Geo_lam Geo_Dav]]" << endl;
	cerr << endl;
	cerr << "MRCIextpol: 0: none, 1: Pey/Buenker, 2: MRMP2, 3: MRMP3" << endl;
	cerr << "lambda    : 0: lambda=1, 1: lambda=lambdalin" << endl;
	cerr << "Davidson  : 0: none, 1: Davidson 1, 2: Davidson 2" << endl;
// 	cerr << "additional options with Geo_..." << endl;
// 	cerr << "Geo_Mu    : Multiplicity for Geometry-Point (default: Multi)" << endl;
// 	cerr << "Geo_Ir    : Irrep for Geometry-Point (default: Irrep)" << endl;
// 	cerr << "Geo_Ro    : Root for Geometry-Point (default: Root)" << endl;
// 	cerr << "Geo_T     : Thresh for Geometry-Point" << endl;
// 	cerr << endl;
}

INT	main(INT argc, char **argv)
{
//	cerr << "excited (Part of DIESEL-MR-CI), " << VERSION << ", " << DATE << endl << endl;
	cerr << "intgrad (Part of DIESEL-MR-CI) "  << endl << endl;
	cout << "intgrad (Part of DIESEL-MR-CI) "  << endl << endl;
    //  look if only base calculation is wanted
    bool    calcBaseOnly = false;
    for (INT i=1; i < argc; i++)
    {
        //cout << "argv["<<i<<"]= "<< argv[i] << endl;
        if (string(argv[i]) == "-b")
        {
            calcBaseOnly = true;
            argc--;
            cout << "calculating only Base" << endl;
        }
    }
    
    
	if ( argc!=5 && argc!=8 && argc!=9 &&argc!=12)
	{
		laber();
		exit(1);
	}
// 	INT		_Irrep		 =0;
//  	INT		_Multi		 =1;
//  	INT		_Root		 =1;
//  	string	_Thresh  	 ="1e-4";
	INT		_Irrep	= atoi(argv[2]);
 	INT		_Multi	= atoi(argv[1]);
 	INT		_Root	= atoi(argv[3]);
 	vector<string>	_Thresh(2,argv[4]);
	vector<INT> MRCIExtPol(2,0); // 0=none, 1=Pey/Buenker, 2=MRMP2, 3=MRMP3
	vector<INT>	lambdaOpt(2,0);
	vector<INT>	Davidson (2,0);   // 0=none, 1=Davidson1, 2=Davidson2
	if (argc>=8)
	{
		MRCIExtPol[1] = MRCIExtPol[0] = atoi(argv[5]); // 0=none, 1=Pey/Buenker, 2=MRMP2, 3=MRMP3
		lambdaOpt[1]  = lambdaOpt[0]  = atoi(argv[6]);
		Davidson [1]  = Davidson [0]  = atoi(argv[7]);	// 0=none, 1=Davidson1, 2=Davidson2
	}
	if (argc>=9)
	{
		_Thresh[1] = atoi(argv[8]);
	}
	
	if (argc==12)
	{
		MRCIExtPol[1] = atoi(argv[9]); // 0=none, 1=Pey/Buenker, 2=MRMP2, 3=MRMP3
		lambdaOpt[1]  = atoi(argv[10]);
		Davidson [1]  = atoi(argv[11]);    // 0=none, 1=Davidson1, 2=Davidson2
	}
		
	
	vector<INT>	block(2,-1);
	vector<INT>	col(2,-1);

	for (INT i=0; i<=1; i++)
	{
		switch ( MRCIExtPol[i] )
		{
		case 0:
			block[i] = 1;
			col[i] = 4;
			break;
	
		case 1:
			switch ( Davidson[i] ) {
			case 0:
				block[i] = 2;
				col[i] = lambdaOpt[i] ? 4 : 2;
				break;

			case 1:
				block[i] = 3;
				col[i] = lambdaOpt[i] ? 3 : 2;
				break;

			case 2:
				block[i] = 3;
				col[i] = lambdaOpt[i] ? 5 : 4;
				break;
			}
			break;
	
		case 2:
			switch ( Davidson[i] ) {
			case 0:
				block[i] = 2;
				col[i] = lambdaOpt[i] ? 7 : 5;
				break;

			case 1:
				block[i] = 3;
				col[i] = lambdaOpt[i] ? 7 : 6;
				break;

			case 2:
				block[i] = 3;
				col[i] = lambdaOpt[i] ? 9 : 8;
				break;
			}
			break;
	
		case 3:
			switch ( Davidson[i] ) {
			case 0:
				block[i] = 2;
				col[i] = lambdaOpt[i] ?10 : 8;
				break;

			case 1:
				block[i] = 3;
				col[i] = lambdaOpt[i] ?11 :10;
				break;

			case 2:
				block[i] = 3;
				col[i] = lambdaOpt[i] ?13 :12;
				break;
			}
			break;
		}
	}

	Block = block[0];
	Column = col[0];
	
	Geo_Block = block[1];
	Geo_Column = col[1];
    
    
    if ( Block > 1 ) selMinusR = true;
    else    selMinusR =  false;

    if ( Geo_Block > 1 ) Geo_selMinusR = true;
    else    Geo_selMinusR = false;

    cout << "Block = " << Block  << endl;
    cout << "Column = " << Column << endl;
    cout << "selMinusR = " << selMinusR << endl;
    cout << "Geo_Block = " << Geo_Block << endl;
    cout << "Geo_Column = " << Geo_Column << endl;
    cout << "Geo_selMinusR = " << Geo_selMinusR << endl;
    
    
 	char 	dirBuf[200];
 	getcwd(dirBuf,200);
 	string	_TMBaseDir(dirBuf);
 	icListe	cf( _TMBaseDir, _Irrep, _Multi, _Root, _Thresh[0]);
 

//      Noch ein Input-Flag einbauen, so dass entweder nur Base 
//      oder alles ausser der Base gerechnet werden kann
//      Dann beachte aber, dass die Dinge wie BaseEnergie an ihren
//      Platz kommen.
//      Dazu ist aber zu sagen, dass die BaseEnergie bei
//      doppelter Genauigkeit nicht unbedingt gebraucht wird
//      , sondern dort nur, um sie in den Energie-File zu schreiben
//      Ausserdem macht eine Neuberechnung keine Arbeit,
//      weil ja nur das gerechnet wird, was noch nicht existiert    
    if ( calcBaseOnly )
    {
      	cf.calcBase();    
    }
    else
    {
  	  	cf.calcBase();
 	 	cf.readInfoFile("gradinfo");
 	 	cf.navigateFileSystem(prepare);
 	 	cf.navigateFileSystem(calc);
 	 	cf.navigateFileSystem(collectResults);
 	 	cf.calcGradient();
 	 	cf.writeGradientFile();
	}
    
	return 0;
}




