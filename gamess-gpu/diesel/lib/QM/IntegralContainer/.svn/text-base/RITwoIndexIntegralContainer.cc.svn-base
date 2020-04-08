//***********************************************************************
//
//	Name:			TwoIndexIntegralContainer.cc
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


#include "RITwoIndexIntegralContainer.h"
#include "TwoIndexIntegralContainer.h"

#include "../IO/Fortran/FortranFileIO.h"
#include "../IO/Fortran/Fort31FirstRecord.h"
#include "../IO/Fortran/Fort31SecondRecord.h"

#include "../IO/TimeTicks.h"

#include "../MO/MOMapping.h"
#include "../IO/Verbosity.h"

#include <fstream>

using namespace std;

// TwoIndexIntegralContainer::TwoIndexIntegralContainer(const MRMOs & _mrmos) 
// {
// 	mrmos = &_mrmos;
// 	n = mrmos->getMaxMO()*(mrmos->getMaxMO()+1)/2;
// 	p = new OneElectronIntegralType[n];
// 	tab = new INT[mrmos->getMaxMO()];
// 	for ( INT i=0 ; i<mrmos->getMaxMO() ; i++ )
// 		tab[i] = i*(i+1)/2;
// }

RITwoIndexIntegralContainer::RITwoIndexIntegralContainer(const MRMOs & _mrmos) 
:TwoIndexIntegralContainer(_mrmos)
,oneintFilename("oneint")
{
}


RITwoIndexIntegralContainer::~RITwoIndexIntegralContainer()
{
// 	delete p;     // wird eigentlich schon in ~TwoIndexIntegralContainer() gemacht;
// 	delete tab;   // wird eigentlich schon in ~TwoIndexIntegralContainer() gemacht;
//	~TwoIndexIntegralContainer();
}




double RITwoIndexIntegralContainer::loadIntegrals(Fort31File _f31)
{
// 	cout << "*******************************************" << endl;
// 	cout << "Modul roneint" << endl;
// 	cout << "*******************************************" << endl;
	cout << "reading informations about core and one-electron-integrals...";
	
	// Beachte Format des Types "Fort31File", siehe in Datei
	if(_f31.format != RIFormat)
	{
		cerr << "Error in FileFormat" << endl;
		exit(-1);
	}
	// construct Filenames
	
	
	oneintFilename = string(_f31.name) + "/" + oneintFilename;
	ifstream Quelle;
	Quelle.open(oneintFilename.c_str(),ios::binary|ios::in);
	
	if(!Quelle)    //muss existieren
	{
		cerr << "ERROR: cannot open file " << oneintFilename << "!\n";
		exit(-1);
	}
	
	// 1.Zeile Zahlenformat
	// in diesem File (4e22.14)
	// bedeutet:	4 Zahlen pro Zeile
	// 		22 Zeichen je Feld fuer eine Zahl
	//		14 Nachkommastellen
	
	string Zahlenformat;
	Quelle >> Zahlenformat;
	
	// 2.Zeile: 3 Zahlen
	// 	1.Zahl:	core-Energie
	//	2.Zahl:	Kern-Kern-Abstossung
	//	3.Zahl: ohne Beudeutung
	
	double coreEnergie,KernKernAbstossung,nixWichtiges;
	Quelle >> coreEnergie >> KernKernAbstossung >> nixWichtiges;
//-------------------------------
/// set core-Energy
	core = coreEnergie;
//-------------------------------

	// 	cout << "coreEnergie:           " << coreEnergie << endl;
	// 	cout << "Kern-Kern-Abstossung:  " << KernKernAbstossung << endl;
		
	// 3.Zeile vier Zahlen
	// 1.Zahl:  Gesamtzahl aller Orbitale
	// 2.Zahl:  Anzahl der irreps
	// 3.Zahl:  Anzahl aller Einzentrenintegrale <i|h|j>, also 0.5 * (nmos[irrep]) * (nmos[irrep]+1) Stueck 
	// 4.Zahl:  weissNicht
	
	INT nmos, nirreps, nEinzentrenIntegrale, weissNicht;
	Quelle >> nmos >> nirreps >> nEinzentrenIntegrale >> weissNicht;
	
//-------------------------------
/// set number of contained integrals
	n = nEinzentrenIntegrale;
//-------------------------------
	
	// 4.Zeile Informationen "uber die MOs in den einzelnen Irreps
	//         und die Anzahl der aktiven MOs
	// 1.) Gesamtzahl der MOs pro Irrep, nach Irreps geordnet
	// 2.) Anzahl aktiver MOs, nach Irreps geordnet
	
	INT* mosInIrrep_ges = new INT[nirreps];
	INT* mosInIrrep_act = new INT[nirreps];
	INT* MOSymmetrieTabelle = new INT[nirreps];
	INT    kontroll_nmos = 0,     // Anzahl der MOs kontrollieren
	    kontroll_amos = 0;     // Anzahl aktiver MOs kontrollieren
	    
	for (INT i=0; i<nirreps; i++)
	{
		Quelle >> mosInIrrep_ges[i]; 
		kontroll_nmos += mosInIrrep_ges[i];
	}
	
	for (INT i=0; i<nirreps; i++)
	{
		Quelle >> mosInIrrep_act[i];
		MOSymmetrieTabelle[i] = kontroll_amos;
		kontroll_amos += mosInIrrep_act[i];
		
	}
//	INT nmos_act = kontroll_amos;
	
	if (kontroll_nmos != nmos)
	{
		cerr << "Fehler beim Lesen der MO-Zahlen" << endl;
		cerr << "nmos = " << nmos << "  Kontrolle: " << kontroll_nmos << endl;
	}
	
	// Zahlenblock mit Codes f"ur die einzelnen MOs
	// Alle MOs sind aufgelistet, nach Irreps geordnet
	// 16 Zahlen pro Zeile, 3 Leerzeichen, eine Zahl, drei Leerzeichen ...
	// Bedeutung der Codes:
	// 0 = virtuelles MO, aktiv
	// 1 = virtuelles MO, gefroren
	// 2 = doppelt besetztes MO, aktiv
	// 3 = doppelt besetztes MO, gefroren
	// 4 = einfach besetztes MO (offene Schale)
	
// 	INT* mocode = new INT[nmos+1];  
// 	
// 	for (INT i=1; i<=nmos; i++)
// 	{
// 		Quelle >> mocode[i];
// 	}
// 
// 	delete [] mocode;


	INT	 moNumber(0);//,moCode(0);
	INT* inIrRep = new INT[nirreps];
	INT** moCode = new INT* [nirreps];
	
	for (INT Ir=0; Ir < nirreps; Ir++)
	{
		inIrRep[Ir] = 0;
		INT frozenVirtual(0);
		INT frozenCore(0);
		moCode[Ir] = new INT[mosInIrrep_ges[Ir]];
		for (INT mo=0; mo < mosInIrrep_ges[Ir]; mo++)
		{
			moNumber++;
			Quelle >> moCode[Ir][mo];
			switch(moCode[Ir][mo])
			{
			case 1:
				frozenVirtual++;
				break;
			case 3:
				frozenCore++;
				break;
			default:
				inIrRep[Ir]++;
			}
			//cout << "moCode[Ir="<<Ir<<"][mo="<<mo<<"]="<<moCode[Ir][mo]<<endl;
		}
		//  check number of mos
		INT checkMOs = frozenVirtual + inIrRep[Ir];
		if (checkMOs != mosInIrrep_act[Ir])
		{
			cerr << "ERROR in reading integral information in RITwoIndexIntegralContainer.cc" << endl;
			cerr << "number of MOs (active + frozen virtual) is wrong" << endl;
		}
		checkMOs += frozenCore;
		if (checkMOs != mosInIrrep_ges[Ir])
		{
			cerr << "ERROR in reading integral information in RITwoIndexIntegralContainer.cc" << endl;
			cerr << "number of MOs (active + frozen virtual + frozen core) is wrong" << endl;
		}
		
	}
	


	// letzter Zahlenblock, der alle Einzentrenintegrale enthaelt
	// Format ist durch Zahlenformat gegeben
	// Es sind die Integrale fuer alle MOs angegeben, ausser den gefrorenen besetzten (mocode == 3)
	// beachte, dass die Integrale fuer die virtuellen gefrorenen MOs (mocode == 1) vorhanden sind.
	INT Blockbeginn = 0,
	  Blockende = 0;
	
	cout << "start reading one-electron integrals" << endl;
	INT r(0),s(0);
	INT rsBeginn(0);
	
	for (INT IR=0; IR< nirreps; IR++)
	  {
	    Blockbeginn = MOSymmetrieTabelle[IR];
	    Blockende = Blockbeginn + mosInIrrep_act[IR];
		
	    cout << "Irrep= " << IR << endl;
	    for (INT i=Blockbeginn; i < Blockende; i++)
	    {
		
			for (INT j=Blockbeginn; j<=i; j++)
			  {
			  	OneElectronIntegralType	EinzentrenIntegral;
				Quelle >> EinzentrenIntegral;
				//cout << "moCode[" << IR << "]["
				//     <<i-Blockbeginn+mosInIrrep_ges[IR]-mosInIrrep_act[IR] << "]= "
				//     << moCode[IR][(i-Blockbeginn+mosInIrrep_ges[IR]-mosInIrrep_act[IR])] << endl;
				//cout << "reading integral  ("<<i<<"|"<<j<<")= "<< EinzentrenIntegral;
				if ((moCode[IR][(i-Blockbeginn+mosInIrrep_ges[IR]-mosInIrrep_act[IR])] != 1)
				 && (moCode[IR][(i-Blockbeginn+mosInIrrep_ges[IR]-mosInIrrep_act[IR])] != 3)
				 && (moCode[IR][(j-Blockbeginn+mosInIrrep_ges[IR]-mosInIrrep_act[IR])] != 1)
				 && (moCode[IR][(j-Blockbeginn+mosInIrrep_ges[IR]-mosInIrrep_act[IR])] != 3))
				{
					//cout << " storing integral in ";
 					//cout << "p[tab[" <<r<< "]+" <<s<< "]=p["<<tab[r]<<"+"<<s<<"]" << endl;
					p[tab[r]+s] = EinzentrenIntegral;
// 					cout << "p[tab[" <<r<< "]" <<s<< "]: ";
// 					cout << EinzentrenIntegral << endl;
					if (r==s) { r++; s=rsBeginn; }
					else { s++; }
					if ((r-rsBeginn)==inIrRep[IR]) { s=r; }
				}
				else { cout << endl; }
// 				cout << "EinzentrenIntegral[" <<i<< "][" <<j<< "]: ";
// 				cout << EinzentrenIntegral << endl;
			}
	    }

		rsBeginn += inIrRep[IR];

	  }
	
	delete [] mosInIrrep_ges;
	delete [] mosInIrrep_act;
	delete [] MOSymmetrieTabelle;
	delete [] inIrRep;
	for (INT Ir=0; Ir < nirreps; Ir++)
	{
		delete [] moCode[Ir];
	}
	delete [] moCode;
	// noch drei worte zum Schluss, aber auf die kann man auch direkt verzichten
	
// 	cout << "*******************************************" << endl;
// 	cout << "oneint erfolgreich gelesen" << endl;
// 	cout << "*******************************************" << endl;
	cout << "successful" << endl;
	return core;
}
