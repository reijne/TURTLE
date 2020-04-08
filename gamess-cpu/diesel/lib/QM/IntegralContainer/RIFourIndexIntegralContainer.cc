#include "RIFourIndexIntegralContainer.h"
#include "../IO/TimeTicks.h"
#include "IndexTranslation.h"
#include "../MO/MOMapping.h"
#include "../IO/Verbosity.h"

using namespace std;

RIFourIndexIntegralContainer::RIFourIndexIntegralContainer(const MRMOs & mrmos)
:FourIndexIntegralContainer(mrmos)
,precalcMode(storeNone)
,rikram(neu)
,InfoFileName("cisdata")
,IntegralFileName("bkji")
,nhilf(0)
,nmos_act(0)
,bkji(NULL)
,isAllprecalced(false)
{
	
}

RIFourIndexIntegralContainer::RIFourIndexIntegralContainer(const MRMOs & mrmos, INT pcMode)
:FourIndexIntegralContainer(mrmos)
,precalcMode(storeAll)
,rikram(neu)
,InfoFileName("cisdata")
,IntegralFileName("bkji")
,nhilf(0)
,nmos_act(0)
,bkji(NULL)
,isAllprecalced(false)
{
	
}

//RIFourIndexIntegralContainer::RIFourIndexIntegralContainer()
//{}




RIFourIndexIntegralContainer::~RIFourIndexIntegralContainer()
{
	// in loadInterals gelesen
  	for (INT i=1; i<=nmos_act; i++)
    	{
      		for (INT j=1; j<=i; j++)
		{
	  		delete [] bkji[i][j];
		}
      		delete [] bkji[i];
    	}
  	delete [] bkji;
}



void   RIFourIndexIntegralContainer::loadInfoFile(string InfoFileName)
{
	// reading Discription of IntegralFile.
	// Only necessary Information: "Number of RI-Integrals"
// ------------------------------------ 
// altes Modul rcisdata 
// ------------------------------------
// 	cout << "*******************************************" << endl;
// 	cout << "read rcisdata" << endl;
// 	cout << "*******************************************" << endl;
	
	ifstream Quelle;
	Quelle.open(InfoFileName.c_str(),ios::binary|ios::in);
	
	if(!Quelle)    //muss existieren
	{
		cerr << "Error: cannot open File " << InfoFileName << " !\n";
		exit(-1);
	}
	
	// 1.Zeile  Zahlenformat
	// eigentlich immer (4e20.12):
	// bedeutet:  4 Zahlen pro Zeile
	//           20 Zeichen je Feld fuer eine Zahl
	//           12 Nachkommastellen 
	string Zahlenformat;

	Quelle >> Zahlenformat;
	

	// 2.Zeile drei bedeutungslose Zahlen
	double bedeutungslos[3];

	for (INT i=0; i<=2; i++)
	{
		Quelle >> bedeutungslos[i];
		cout << bedeutungslos[i] << endl;
	}
	

	// 3.Zeile vier Zahlen
	// 1.Zahl:  Gesamtzahl aller Orbitale
	// 2.Zahl:  Anzahl der irreps
	// 3. und 4.Zahl:  ohne Bedeutung
        INT kontroll_nmos,      // Anzahl der MOs
            kontroll_nirreps,  // Anzahl der Irreps
            dummy1, dummy2;  // Dummyvariablen

	
	
	Quelle >> kontroll_nmos >> kontroll_nirreps >> dummy1 >> dummy2;

// ********************************
//Hier erstmal festlegen der Anzahl der MOs
	INT nmos = kontroll_nmos; // ACHTUNG: -- HIER gesamte MOs
	INT nirreps = kontroll_nirreps;
	// Vergleich mit der Zahl aus MOIrReps
	INT IrReps = mrmos.getNumberOfIrReps();
	if (nirreps != IrReps)
	{
		cerr << "ERROR: Number of IrReps is NOT OK." << endl;
		cerr << "       file:   " << nirreps << endl;
		cerr << "       IrReps: " << IrReps << endl;
		exit(-1);
	}
	//else { cout << "Number of IrReps is " << nirreps << endl; }
// ********************************
	
	// 4.Zeile Informationen "uber die MOs in den einzelnen Irreps
	//         und die Anzahl der aktiven MOs
	// 1.) Gesamtzahl der MOs pro Irrep, nach Irreps geordnet
	// 2.) Anzahl aktiver MOs, nach Irreps geordnet
	
	INT kontroll_mosInIrrep_ges;
	INT kontroll_mosInIrrep_act;
	    
	for (INT i=0; i<nirreps; i++)
	{
		Quelle >> kontroll_mosInIrrep_ges; 
	}
	
	for (INT i=0; i<nirreps; i++)
	{
		Quelle >> kontroll_mosInIrrep_act;
	}
	
	
	
	
	
	// Zahlenblock mit Codes f"ur die einzelnen MOs
	// Alle MOs sind aufgelistet, nach Irreps geordnet
	// 16 Zahlen pro Zeile, 3 Leerzeichen, eine Zahl, drei Leerzeichen ...
	// Bedeutung der Codes:
	// 0 = virtuelles MO, aktiv
	// 1 = virtuelles MO, gefroren
	// 2 = doppelt besetztes MO, aktiv
	// 3 = doppelt besetztes MO, virtuell
	// 4 = einfach besetztes MO (offene Schale)
	
	INT kontroll_mocode;  
	for (INT i=1; i<=nmos; i++)
	{
		Quelle >> kontroll_mocode;
	}
	
	
	
	
	
	// neue Zeile: Antwort aller Fragen
	// Zahl 42 in e-Darstellung
	
        double zweiundvierzig;
  		
	Quelle >> zweiundvierzig;
	
	cout << zweiundvierzig << endl;
	
	// Ein Satz ohne Bedeutung
	// Im Beispiel vier Woerter, vielleicht nicht immer so
	// Besser ganzen Satz/Zeile mit getline einlesen
	// Aber wie funktioniert das vom File aus ???
	
	//INT maxz;
	//char dummeZeile[maxz];
	
	string dumm1,dumm2,dumm3,dumm4;
	
	//Quelle.getline(dummeZeile,maxz);
	
	Quelle >> dumm1 >> dumm2 >> dumm3 >> dumm3;
	
	
	// ein Sternchen
	
	char Stern;

	Quelle >> Stern;
	
	
	// zwei Zahlen
	// 1.Zahl:  Anzahl der besetzten, aktiven MOs
	// 2.Zahl:  Anzahl der virtuellen, aktiven MOs
	
	INT nmos_occ,     // Anzahl der besetzten, aktiven MOs
	    nmos_virt;    // Anzahl der virtuellen, aktiven MOs
	
	Quelle >> nmos_occ >> nmos_virt;
// ********************************
// DAS ist die Anzahl aller aktiven MOs;
	nmos_act = nmos_occ + nmos_virt;
	cout << "Total-Number of active Orbitals: " << nmos_act << endl;
	//Vergleich mit der in MRMOs bestimmten Anzahl der aktiven MOs
	INT maxMO = mrmos.getMaxMO();
	//cout << "MaxMO: " << maxMO << " (as defined in MOIrReps)" << endl;
	if (nmos_act != maxMO)
	{
		cerr << "ERROR: Number of active Orbitals is NOT OK." << endl;
		cerr << "       file:   " << nmos_act << endl;
		cerr << "       IrReps: " << maxMO << endl;
		exit(-1);
	}
	//else { cout << "Number of active Orbitals is OK" << endl; }
// ********************************
	
	
	
	
	// Zahlenblock mit den Eigenwerten der aktiven MOs
	// Zahlenformat ist in der ersten Zeile angegeben
	// "ublicherweise:  (4e20.12)
	
	
	double* Eigenwert = new double[nmos+1];
	
	for (INT i=1; i<= nmos_act; i++)
	{
		Quelle >> Eigenwert[i];
		//	cout << i << "-ter Eigenwert: " << Eigenwert[i] << endl;
	}
	
	delete [] Eigenwert;
	
	// letzte Zeile  Anzahl der Hilfsbasisfunktionen
	// ********************************
	// Endlich die einzige wichtige Zeile in der ganzen Funktion
	// ********************************
		
	
	Quelle >> nhilf;
	cout << "Number of auxiliary basis functions: " << nhilf << endl;
	
	// Ende des Files
	
// 	cout << "*******************************************" << endl;
// 	cout << "read of " << InfoFileName << " successful." << endl;
// 	cout << "*******************************************" << endl;
	
}


void   RIFourIndexIntegralContainer::precalculateIntegrals()
{
	cout << "precalculating integrals" << endl;
	cout << "There are "<< mrmos.getMaxMO() <<"^3= "<< mrmos.getMaxMO()*mrmos.getMaxMO()*mrmos.getMaxMO()  
	     << " Integrals to calculate."<<endl;
	TwoElectronIntegralIndex<MOType> index(1,1,1,1);
	for (MOType i = 1; i <= mrmos.getMaxMO(); i++)
	{
		for (MOType j = 1; j <= i; j++)
		{
			//for (MOType k = 1; k <= j; k++)
			//{
				IntegralType TEI = 0;
				cout << "precalc ("<<i<<","<<i<<"|"<<j<<","<<j<<")"<<endl;
				index.set(i,i,j,j);
				TEI = this->operator[](index);
				cout << "  mit dem Wert:  " << TEI << endl;
				//if(this->FourIndexIntegralContainer::set(index,TEI))
				if(FourIndexIntegralContainer::set(index,TEI))
					cout << "es hat funktioniert" << endl;
				else
					cout << "es hat nicht geklappt." << endl;
				cout << "precalc ("<<i<<","<<j<<"|"<<i<<","<<j<<")"<<endl;
				index.set(i,j,i,j);
				TEI = this->operator[](index);
				cout << "  mit dem Wert:  " << TEI << endl;
				//this->FourIndexIntegralContainer::set(index,TEI);
				FourIndexIntegralContainer::set(index,TEI);
			//}
		}
	}
	cout << "precalculation of integrals finished" << endl;
}



void   RIFourIndexIntegralContainer::precalculateAllIntegrals()
{

	if ( verbosity.isActive(Verbosity::Integrals) )
	{
		cout << "precalculating integrals" << endl;
	}



//INT	recSize = 2000;
//IntegralType	buf[recSize];
INT	integrals = 0, intsContained;
IntegralType	TEI;


TwoElectronIntegralIndex<MOType>	ind;
TwoElectronIntegralTriadeIndex	indSort;



//Fort31SecondRecord	f31(_f31.name, _f31.format);
TimeTicks	ticks;

	intsContained = 0;
	for ( INT checking=0 ; checking<=0 ; checking++ )
	{	if ( !checking )
			cout << "reading two electron integrals...";
		else
		{	cout << "checking two electron integrals...";
		}
		cout.flush();
		ticks.start();
//		fort31.rewind();
//		if ( _f31.format==Fort31RecordFormatTRADPT )
//			fort31.nextRecord();
//		fort31.nextRecord();
//		fort31.nextRecord();
		integrals = 0;

//	INT	pos = recSize;
	INT	is, ic = 0;

	for ( INT i=0 ; i<mrmos.getNumberOfIrReps() ; i++ )
		for ( INT j=0 ; j<=i ; j++ )
			for ( INT k=0 ; k<=i ; k++ )
				for ( INT l=0 ; l<=((i==k) ? j : k) ; l++ )
				{//	cout << "(" << i << " " << j << "|" << k << " " << l << ")" << endl;
					
					is = 0;
					ic++;
					
					
					if ( mrmos.getProd(
							mrmos.getProd(i, j), mrmos.getProd(k, l)) )
						continue;
					if ( i==l )	
					{
					//cout << "(AA|AA)" << endl;
	// (AA|AA) ---------------------------------------------------------
	for ( INT ii=0 ; ii<mrmos.getInIrRep(i) ; ii++ )
	{	ind[0] = moMapping.getContinuous(ii + mrmos.getStartMO(i));
		for ( INT jj=0 ; jj<=ii ; jj++ )
		{	ind[1] = moMapping.getContinuous(jj + mrmos.getStartMO(i));
			for ( INT kk=0 ; kk<=ii ; kk++ )
			{	ind[2] = moMapping.getContinuous(kk + mrmos.getStartMO(i));
				for ( INT ll=0 ; ll<=((ii==kk) ? jj : kk) ; ll++ )
				{	ind[3] = moMapping.getContinuous(ll + mrmos.getStartMO(i));
					indSort.set(ind);
					
					TEI = this->operator[](indSort);	
					intsContained += this->set(indSort, TEI);

					integrals++;
					is++;
				}
			}
		}
	}
//					cout << is << " <--> " << f31.getNIT(ic)-f31.getNIT(ic-1) << endl;
					}
					else
					if ( i==j )
					{
					//cout << "(AA|BB)" << endl;
	// (AA|BB) ---------------------------------------------------------
	for ( INT ii=0 ; ii<mrmos.getInIrRep(i) ; ii++ )
	{	ind[0] = moMapping.getContinuous(ii + mrmos.getStartMO(i));
		for ( INT jj=0 ; jj<=ii ; jj++ )
		{	ind[1] = moMapping.getContinuous(jj + mrmos.getStartMO(i));
			for ( INT kk=0 ; kk<mrmos.getInIrRep(k) ; kk++ )
			{	ind[2] = moMapping.getContinuous(kk + mrmos.getStartMO(k));
				for ( INT ll=0 ; ll<=kk ; ll++ )
				{	ind[3] = moMapping.getContinuous(ll + mrmos.getStartMO(k));
					indSort.set(ind);
					
					
					TEI = this->operator[](indSort);	
					intsContained += this->set(indSort, TEI);
					integrals++;
					is++;
				}
			}
		}
	}
//					cout << is << " <--> " << f31.getNIT(ic)-f31.getNIT(ic-1) << endl;
					}
					else
					if ( i==k )
					{
					//cout << "(AB|AB)" << endl;
	// (AB|AB) ---------------------------------------------------------
	for ( INT ii=0 ; ii<mrmos.getInIrRep(i) ; ii++ )
	{	ind[0] = moMapping.getContinuous(ii + mrmos.getStartMO(i));
		for ( INT jj=0 ; jj<mrmos.getInIrRep(j) ; jj++ )
		{	ind[1] = moMapping.getContinuous(jj + mrmos.getStartMO(j));
			for ( INT kk=0 ; kk<=ii ; kk++ )
			{	ind[2] = moMapping.getContinuous(kk + mrmos.getStartMO(i));
				for ( INT ll=0 ; ll<mrmos.getInIrRep(j) ; ll++ )
				{	if ( kk==ii && ll>jj )
						continue;

					ind[3] = moMapping.getContinuous(ll + mrmos.getStartMO(j));
					indSort.set(ind);
					
					TEI = this->operator[](indSort);	
					intsContained += this->set(indSort, TEI);
					integrals++;
					is++;
				}
			}
		}
	}
//					cout << is << " <--> " << f31.getNIT(ic)-f31.getNIT(ic-1) << endl;
					}
					else
					{
					//cout << "(AB|CD)" << endl;
	// (AB|CD) ---------------------------------------------------------
	for ( INT ii=0 ; ii<mrmos.getInIrRep(i) ; ii++ )
	{	ind[0] = moMapping.getContinuous(ii + mrmos.getStartMO(i));
		for ( INT jj=0 ; jj<mrmos.getInIrRep(j) ; jj++ )
		{	ind[1] = moMapping.getContinuous(jj + mrmos.getStartMO(j));
			for ( INT kk=0 ; kk<mrmos.getInIrRep(k) ; kk++ )
			{	ind[2] =moMapping.getContinuous( kk + mrmos.getStartMO(k));
				for ( INT ll=0 ; ll<mrmos.getInIrRep(l) ; ll++ )
				{	ind[3] = moMapping.getContinuous(ll + mrmos.getStartMO(l));
					indSort.set(ind);

					
					TEI = this->operator[](indSort);	
					intsContained += this->set(indSort, TEI);
					integrals++;
					is++;
				}
			}
		}
	}

//					cout << is << " <--> " << f31.getNIT(ic)-f31.getNIT(ic-1) << endl;

					} // end if
				}	// end for l
		ticks.stop();
		cout << ticks << endl;
	}	
	if ( verbosity.isActive(Verbosity::Integrals) )
	{
		cout << "precalculation of integrals finished" << endl;
		cout << "number of read two electron integrals: " << integrals << endl;
		cout << "number of stored two electron integrals: " << intsContained << endl;
	}
	isAllprecalced = true;
}



void   RIFourIndexIntegralContainer::loadIntegrals(Fort31File _f31)
{
	// Beachte Format des Types "Fort31File", siehe in Datei
	if(_f31.format != RIFormat)
	{
		cerr << "Error in FileFormat" << endl;
		exit(-1);
	}
	// construct Filenames
	InfoFileName = string(_f31.name) + "/" + InfoFileName;
	IntegralFileName = string(_f31.name) + "/" + IntegralFileName;
	// Es m"ussen die Dateien "coredata" und "bkji" gelesen werden
	// Erstmal lesen vor coredata, um zu erfahren, wieviele Integrale vorhanden sind.
	// Dies geschieht in einer seperaten Funktion, um sp"atere Ver"anderungen in der 
	// Formatierung des Files leichter zu verarbeiten.
	
	switch (rikram)
	{
	case alt:
		loadInfoFile(InfoFileName);
		break;
	case neu:
		nmos_act = mrmos.getMaxMO();
		break;
	}
	
	
// *******************************************
	
// ------------------------------------ 
// altes Modul rbkji
// ------------------------------------

// 	cout << "*******************************************" << endl;
// 	cout << "Modul rbkji" << endl;
// 	cout << "*******************************************" << endl;
	
	
	

	
	ifstream bkjiQuelle;
	bkjiQuelle.open(IntegralFileName.c_str(),ios::binary|ios::in);
	
	if(!bkjiQuelle)    //muss existieren
	{
		cerr << "Cannot open "<< IntegralFileName << "!\n";
		exit(-1);
	}
	
	double f;
	INT Recordgroesse,
	    Zaehler = 0,
	    nbuf = 0,
	    Anzahl;
	    
//	In der neuen RI-Version: lese zuerst die Anzahl der Hilfsbasis-Funktionen
	if (rikram==neu)
	{
		INT	nmos_ges = 0;
		bkjiQuelle.read(reinterpret_cast<char*>(&Recordgroesse), sizeof(Recordgroesse));
		bkjiQuelle.read(reinterpret_cast<char*>(&nmos_act), sizeof(nmos_ges));
		bkjiQuelle.read(reinterpret_cast<char*>(&nhilf), sizeof(nhilf));
		bkjiQuelle.read(reinterpret_cast<char*>(&Recordgroesse), sizeof(Recordgroesse));
		cout << "Total Number of mos (including frozen):  " << nmos_act << endl;
		cout << "Number of auxiliary basis functions:     " << nhilf << endl;
	}	

// *******************************************

 	bkji = new double**[nmos_act+1];
        for (INT i=0; i<=nmos_act; i++)
        {
	        bkji[i] = new double*[i+1];
	        for (INT j=0; j<=i; j++)
	        {
		        bkji[i][j] = new double[nhilf+1];
	        }
        }

// *******************************************

//	Lese ersten integer mit Recordgroesse
	
	bkjiQuelle.read(reinterpret_cast<char*>(&Recordgroesse), sizeof(Recordgroesse));
	Anzahl = Recordgroesse / sizeof(f);
	

//	AAAAAAACHTUUUUUUUUNG
//	Der Leseprozess funktioniert nur dann, wenn nhilf < Anzahl ist
//	In den bisher betrachteten F"allen ist Anzahl = 8192,
//	d.h. es darf nicht mehr als 8192 RI Hilfsbasisfunktionen geben
//	!!!!  bitte beachten.
//	f"ur den Fall, dass nhilf > 8192, m"usste eine Modulo Division eigebaut,
//	die das Programm dazu veranlasst, die entsprechende Bufferzahl zu lesen
//	Falls also (Bufferzahl = nhilf div Anzahl) > 0, werden Bufferzahl Buffer
//	gelesen. Falls Bufferzahl = 0, werden halt keine Buffer gelesen.


//	Eine andere M"oglichkeit ist es nat"urlich, einen Hilfsbuffer zu benutzen,
//	in den jeweils ein ganzer Buffer gelesen wird, und der dann Indexweise in die
//	Matrix kopiert wird.
//	Dies ist aber sicherlich die aufwendigere Prozedur

	
	for (INT i = 1; i<= nmos_act; i++) 
	{
		for (INT j = 1; j <= i; j++)
		{
			Zaehler += nhilf;
			//cout << " Zaehler = " << Zaehler+1 - nhilf << " bis " << Zaehler+1 << endl;
			if (Zaehler <= Anzahl)
			{
				bkjiQuelle.read(reinterpret_cast<char*>(&bkji[i][j][1]), nhilf*sizeof(f));
				//cout << "immer noch Buffer-Nummer" << nbuf << endl;
			}
			else    // Buffergrenze wird "uberschritten
			{
				INT zuviel = Zaehler - Anzahl,    // Elemente, die bereits im neuen Buffer sind
				    noch_drin = nhilf - zuviel;  // Elemente, die noch im alten Buffer sind
				
				//cout << zuviel << " " << noch_drin << endl;
				// lese erstmal den alten Buffer zu Ende    
				bkjiQuelle.read(reinterpret_cast<char*>(&bkji[i][j][1]), noch_drin*sizeof(f));
				nbuf++;
				//cout << nbuf << " Buffer gelesen" << endl;
				// bearbeite die Buffergrenze des alten Buffers
				bkjiQuelle.read(reinterpret_cast<char*>(&Recordgroesse), sizeof(Recordgroesse));
				
				
				// lese Recordgroesse des neuen Buffers, eigentlich sollten aber alle Buffer gleich gross sein
				bkjiQuelle.read(reinterpret_cast<char*>(&Recordgroesse), sizeof(Recordgroesse));
				Anzahl = Recordgroesse / sizeof(f);
				
				// lese die restlichen Elemente im neuen Buffer
				bkjiQuelle.read(reinterpret_cast<char*>(&bkji[i][j][1+noch_drin]), zuviel*sizeof(f));
				
				// korrigiere den Zaehlerstand
				Zaehler = zuviel;
			}
		}
	}
	
//	cout << "*******************************************" << endl;
	cout << "Read successful RI-Integrals " << endl;
// 	cout << "from file << " << IntegralFileName << endl;
// 	cout << "*******************************************" << endl;

	if (precalcMode == storeAll) precalculateAllIntegrals();

}

