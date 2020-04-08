#include "RISym.h"

using namespace std;

RISym::RISym(string Pfad)
{
  /// Einlesen der Daten

  Pfadname = Pfad;
  Filename = "oneint";
  // Dateiname = sprintf("%s/%s",Pfadname.c_str(),Filename.c_str());
  Dateiname = Pfadname + "/" + Filename;
  ifstream Quelle;
  Quelle.open(Dateiname.c_str(), ios::binary|ios::in);
  
  if(!Quelle) // Datei muss existieren
    {
      cerr << "ERROR: file " << Dateiname << " does not exist." << endl;
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
	
	double coreEnergie, KernKernAbstossung, nixWichtiges;
	Quelle >> coreEnergie >> KernKernAbstossung >> nixWichtiges;
	//cout << "coreEnergie:           " << coreEnergie << endl;
	//cout << "Kern-Kern-Abstossung:  " << KernKernAbstossung << endl;
	
	// 3.Zeile vier Zahlen
	// 1.Zahl:  Gesamtzahl aller Orbitale
	// 2.Zahl:  Anzahl der irreps
	// 3.Zahl:  Anzahl aller Einzentrenintegrale <i|h|j>, also 0.5 * (nmos) * (nmos-1) Stueck
	// 4.Zahl:  weissNicht
	
	INT nmos, nirreps, nEinzentrenIntegrale, weissNicht;
	Quelle >> nmos >> nirreps >> nEinzentrenIntegrale >> weissNicht;
	IrReps = nirreps;
	
	// 4.Zeile Informationen "uber die MOs in den einzelnen Irreps
	//         und die Anzahl der aktiven MOs
	// 1.) Gesamtzahl der MOs pro Irrep, nach Irreps geordnet
	// 2.) Anzahl aktiver MOs, nach Irreps geordnet
	
	INT *mosInIrrep_ges = new INT[nirreps];
	INT *mosInIrrep_act = new INT[nirreps];
	inIrRep = new INT[nirreps]; // Anzahl der aktiven MOs je irrep
	INT    kontroll_nmos = 0;     // Anzahl der MOs kontrollieren
	    
	for (INT i=0; i<nirreps; i++)
	{
		Quelle >> mosInIrrep_ges[i]; 
		kontroll_nmos += mosInIrrep_ges[i];
	}
	
	for (INT i=0; i<nirreps; i++)
	{
		Quelle >> mosInIrrep_act[i];
	}
	
	if (kontroll_nmos != nmos)
	{
	  cerr << "ERROR: number of MOs in file " << Dateiname << " might be wrong." << endl;
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
	
// 	INT* mocode = new INT[nmos];  
// 	for (INT i=1; i<=nmos; i++)
// 	{
// 		Quelle >> mocode[i];
// 	}
// 	delete [] mocode;


	INT	 moNumber(0),moCode(0);
	
	for (INT Ir=0; Ir < nirreps; Ir++)
	{
		inIrRep[Ir] = 0;
		INT frozenVirtual(0);
		INT frozenCore(0);
		for (INT mo=0; mo < mosInIrrep_ges[Ir]; mo++)
		{
			moNumber++;
			Quelle >> moCode;
			switch(moCode)
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
		}
		//  check number of mos
		INT checkMOs = frozenVirtual + inIrRep[Ir];
// 		cout << "checking irrep " << Ir 
// 		     << " mos (active=" << inIrRep[Ir]
// 			 << " frozen(virt)=" << frozenVirtual
// 			 << " frozen(occ.)=" << frozenCore
// 			 << ")" << endl;
		
		if (checkMOs != mosInIrrep_act[Ir])
		{
			cerr << "ERROR in reading integral information in RISym.cc" << endl;
			cerr << "number of MOs (active + frozen virtual) is wrong" << endl;
		}
		checkMOs += frozenCore;
		if (checkMOs != mosInIrrep_ges[Ir])
		{
			cerr << "ERROR in reading integral information in RISym.cc" << endl;
			cerr << "number of MOs (active + frozen virtual + frozen core) is wrong" << endl;
		}
		
	}
	
	
	delete [] mosInIrrep_ges;
	delete [] mosInIrrep_act;
	
	
	
/// Produkttafel erzeugen;
	
//	ACHTUNG:
//	Produkttafel gilt nur, wenn die Irreps in turbomol
//	so geordnet sind, dass 
//  ir1 und ir2 eine Untergruppe bilden
//	ir1 bis ir4 die naechst groessere
//	ir1 bis ir8 dann d2h bilden
//	eventuell Erweiterung noetig, falls
//	turbomol von dieser Konvention abweicht.
//	? vielleicht mit map Vertauschung entsprechender
//	Zeilen und Reihen vornehmen
//	Standardreihenfolge fuer jede Symmetriegruppe fest vorgeben
//	eingelesene Reihenfolge damit vergleichen
//	Reihen und Zeilen entsprechend vertauschen
//	weitererer Nachteil: es muessten weitere Files eingelesen werden

	INT n=(INT)(double(log((double)IrReps))/double(log(2.0))+.1);
	prodTab = new IrRep[IrReps*IrReps];
	prodTab[0] = 0;
	//INT n=INT(log(IrReps)/log(2));
	INT A = 0, q = 0;
	for (INT z=1; z<=n; z++)
	{	
	  for (INT i=(1<<(z-1)); i<(1<<z); i++)
		{
		  for (INT j=0; j<(1<<z); j++)
			{
			  q = (1<<(z-1));
				if (j<q) {A = prodTab[(i-q)*IrReps+j] + q;}
				else A = prodTab[(i-q)*IrReps+j-q];
				prodTab[i*IrReps+j] = prodTab[j*IrReps+i]= A;
			}
		}
	}
}
	
	
RISym::~RISym()
{
  delete [] inIrRep;
  delete [] prodTab;
}
	
