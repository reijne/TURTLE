#include "../../../config.h"
#ifndef __RIINTEGRALCONTAINER_H
#define __RIINTEGRALCONTAINER_H __RIINTEGRALCONTAINER_H

#include "../IntegralIndex/TwoElectronIntegralIndex.h"
#include "../IntegralIndex/TwoElectronIntegralTriadeIndex.h"
#include "../IntegralIndex/TwoElectronIntegralCbExIndex.h"

#include "FourIndexIntegralContainer.h"
#include "IndexTranslation.h"
#include "BaseContainer.h"
#include "SymmetryContainer.h"

#include "IntegralType.h"

#include "../IO/Fortran/Fort31File.h"

#include <fstream>
#include <iostream>
using std::cerr;
using std::cout;
using std::endl;

#include <string>

class RIFourIndexIntegralContainer : public FourIndexIntegralContainer {
public:
//	RIFourIndexIntegralContainer();
	RIFourIndexIntegralContainer(const MRMOs & mrmos);
	RIFourIndexIntegralContainer(const MRMOs & mrmos, INT pcMode);
	~RIFourIndexIntegralContainer();

	//	the folowing operators do not perform any check if then integrals
	//	are actually contained in tree
	IntegralType	operator[] (const TwoElectronIntegralIndex<MOType> &) const;
	IntegralType	operator[] (const TwoElectronIntegralTriadeIndex &) const;

	//	all get/set methods return if the requested integral 
	//	is currently contained in tree
	//INT	set(const TwoElectronIntegralIndex<MOType> &, IntegralType);  // braucht man erstmal nicht
	//INT	set(const TwoElectronIntegralTriadeIndex &, IntegralType);    // braucht man erstmal nicht

	INT	get(const TwoElectronIntegralIndex<MOType> &, IntegralType &) const;
	INT	get(const TwoElectronIntegralTriadeIndex &, IntegralType &) const;

	INT	get(const TwoElectronIntegralCbExIndex &,
		IntegralType &, IntegralType &) const;
//--------------------------------------------------------------------	

	void	loadIntegrals(Fort31File f31);   //Integrale laden // und Speicher alloziieren
	void	precalculateIntegrals();	//precalculation of some integrals
	void	precalculateAllIntegrals();	//precalculation of some integrals
	enum 	IntegralStorageMode { storeNone, storeAll } precalcMode;
	enum	riVersion	{alt,neu};

//--------------------------------------------------------------------	
//  tauchen nirgendwo auf -- ueberflussig ?? --
// 	INT	getNumberOfTotalIntegrals() const;    
// 
// 	INT	getNumberOfContainedIntegrals() const;
// 

//--------------------------------------------------------------------	

	const MRMOs	*getMRMOs() const;
	
// ******************
// weitere wichtige Funktionen, die mehrfach benoetigt werden:
//	- Abfrage der Symmetrie nicht noetig, weil immer nur solche Integrale abgerufen werden, die von sich aus ungleich 0 sind
//	- Routine zum Bilden des Skalarproduktes
//	- Beachte hierbei den Zugriff ueber verkuerzte Zeiger *p[b=0..ende] statt bkji[i][j][b=0..ende]
//	- Eigene Funktion zum Umwandeln von TriadenIndices in normale Indices
//	  ??? Namenskonfikte mit Konstruktoren bzw. Typumwandlungsoperatoren aus den entsprechenden Klassen ???
// ******************
	
	

private:
// ------------------------------------
// Geerbte Daten
//IndexTranslation	*index;
//MRMOs	mrmos;
// ------------------------------------
riVersion   rikram;
double		sdot(double* pa,double* pb) const;
void 		loadInfoFile(string InfoFileName);

string 		InfoFileName;
string 		IntegralFileName;
INT 		nhilf;		// Anzahl der Hilfsbasisfunktionen
INT 		nmos_act;	// Anzahl aller aktiven MOs
double***	bkji;

void		tausche(MOType& i,MOType& j) const;

// Debugfunktion
// void		Aus() const;
// void		Ausgabe(const INT i,const INT j, const INT k, const INT l, IntegralType TEI) const;
//weitere Testfunktionen
//ofstream	intstat;
void		Ausgabe(const INT i,const INT j, const INT k, const INT l, IntegralType TEI) const;
bool		iscontained(const MOType i,const  MOType j,const  MOType k,const  MOType l) const;
bool		isAllprecalced;
//INT		counter;
};


// inline 
// void		RIFourIndexIntegralContainer::Aus() const
// {
// 	cout << "Kontroll-";
// }


// inline 
// void		RIFourIndexIntegralContainer::Ausgabe(const INT i,const INT j,const INT k,const INT l,IntegralType TEI) const
// {
// 	cout << "Integral (" << i << ","<<j<<"|"<<k<<","<<l<<")= "<<TEI<<endl;
// }

static INT	counter;
//static	ofstream	intstat("intstat");

inline 
void		RIFourIndexIntegralContainer::Ausgabe(const INT i,const INT j,const INT k,const INT l,IntegralType TEI) const
{
	//counter++;
//	static INT	counter;
//	static	ofstream intstat("intstat");
	cout << counter++ << endl;
//	intstat << i << " "<<j<<" "<<k<<" "<<l<<" :    "<< counter <<endl;
}


inline
bool		RIFourIndexIntegralContainer::iscontained(const MOType i,const  MOType j,const  MOType k,const  MOType l) const
{
	if (i==j || i==k || i==l || j==k || j==l || k==l)
		return true;
	return false;
}


// vielleicht global
inline
double   RIFourIndexIntegralContainer::sdot(double* pa,double* pb) const
{
        double   s = 0;
        for (INT i = 1; i <= nhilf; i++)
        {
                s += pa[i] * pb[i];
        }
        return s;
}

inline
void    RIFourIndexIntegralContainer::tausche(MOType& i,MOType& j) const
{
        MOType temp = i;
        i = j;
        j = temp;
//      cout << "\n --------- \n Tausch wirklich notwendig! \n ----------- " << endl;
}


inline	
IntegralType	RIFourIndexIntegralContainer::operator[] (const TwoElectronIntegralIndex<MOType> &ind) const
{
	if (isAllprecalced) 
	{
		return FourIndexIntegralContainer::operator[](ind);
	}
	MOType i = ind.getI();
	MOType j = ind.getJ();
	MOType k = ind.getK();
	MOType l = ind.getL();
	if (j > i) tausche(i,j);
	if (l > k) tausche(k,l);
	IntegralType TEI = 0;
	//if (iscontained(i,j,k,l)) Ausgabe(i,j,k,l,TEI);
	TEI = sdot(bkji[i][j],bkji[k][l]);
	return TEI;
}


inline 
IntegralType	RIFourIndexIntegralContainer::operator[] (const TwoElectronIntegralTriadeIndex &ind) const
{
	if (isAllprecalced) 
	{
		return FourIndexIntegralContainer::operator[](ind);
	}
	INT h = ind.getM();
	MOType i,j,k,l;
	switch(h)
	{
		case 0	:
			i = ind.getI();
			j = ind.getJ();
			k = ind.getK();
			l = ind.getL();
			break;
		case 1	:
			i = ind.getI();
			j = ind.getL();
			k = ind.getK();
			l = ind.getJ();
			break;
		case 2	:
			i = ind.getI();
			j = ind.getK();
			k = ind.getJ();
			l = ind.getL();
			break;
		default	: cerr << "ERROR: Cannot convert TwoElectronIntegralTriadeIndex" << endl;
			exit(-1);
	}
	if (j > i) tausche(i,j);
	if (l > k) tausche(k,l);
	IntegralType TEI = 0;
	//if (iscontained(i,j,k,l)) Ausgabe(i,j,k,l,TEI);
	TEI = sdot(bkji[i][j],bkji[k][l]);
// 	TwoElectronIntegralIndex tind(i,j,k,l);
// 	
// 	IntegralType TEI = operator[](tind);
	return TEI;
}


inline	
INT	RIFourIndexIntegralContainer::get(const TwoElectronIntegralIndex<MOType> &ind, IntegralType &TEI) const
{
	if (isAllprecalced) 
	{
		return FourIndexIntegralContainer::get(ind,TEI);
	}
	INT irgendwas = 0;
	MOType i = ind.getI();
	MOType j = ind.getJ();
	MOType k = ind.getK();
	MOType l = ind.getL();
	if (j > i) tausche(i,j);
	if (l > k) tausche(k,l);
	//if (iscontained(i,j,k,l)) Ausgabe(i,j,k,l,TEI);
	TEI = sdot(bkji[i][j],bkji[k][l]);
	return irgendwas;
}


inline	
INT	RIFourIndexIntegralContainer::get(const TwoElectronIntegralTriadeIndex &ind, IntegralType &TEI) const
{
	if (isAllprecalced) 
	{
		return FourIndexIntegralContainer::get(ind,TEI);
	}
	INT irgendwas = 0;
	
	INT h = ind.getM();
	MOType i,j,k,l;
	switch(h)
	{
		case 0	:
			i = ind.getI();
			j = ind.getJ();
			k = ind.getK();
			l = ind.getL();
			break;
		case 1	:
			i = ind.getI();
			j = ind.getL();
			k = ind.getK();
			l = ind.getJ();
			break;
		case 2	:
			i = ind.getI();
			j = ind.getK();
			k = ind.getJ();
			l = ind.getL();
			break;
		default	: cerr << "ERROR: Cannot convert TwoElectronIntegralTriadeIndex" << endl;
			exit(-1);
	}
	if (j > i) tausche(i,j);
	if (l > k) tausche(k,l);
	//if (iscontained(i,j,k,l)) Ausgabe(i,j,k,l,TEI);
	TEI = sdot(bkji[i][j],bkji[k][l]);
	
// 	TwoElectronIntegralIndex<MOType> tind(i,j,k,l);
// 	
// 	get(tind,TEI);
	return irgendwas;
}


inline	
INT	RIFourIndexIntegralContainer::get(const TwoElectronIntegralCbExIndex &ind,
		IntegralType &CbTEI, IntegralType &ExTEI) const
{
	if (isAllprecalced) 
	{
		return FourIndexIntegralContainer::get(ind,CbTEI,ExTEI);
	}
//	cout << "\n RICbEx aufgerufen." << endl;
	INT irgendwas = 0;
	MOType i,j,k,l;
	
	INT m = ind.getCoulombTriade();
	switch(m)
	{
		case 0	:
			i = ind.getI();
			j = ind.getJ();
			k = ind.getK();
			l = ind.getL();
			break;
		case 1	:
			i = ind.getI();
			j = ind.getL();
			k = ind.getK();
			l = ind.getJ();
			break;
		case 2	:
			i = ind.getI();
			j = ind.getK();
			k = ind.getJ();
			l = ind.getL();
			break;
		default	: cerr << "ERROR: Cannot convert TwoElectronIntegralTriadeIndex" << endl;
			exit(-1);
	}

	if (j > i) tausche(i,j);
	if (l > k) tausche(k,l);
	//if (iscontained(i,j,k,l)) Ausgabe(i,j,k,l,CbTEI);
	CbTEI = sdot(bkji[i][j],bkji[k][l]);

	INT n = ind.getExchangeTriade();
	switch(n)
	{
		case 0	:
			i = ind.getI();
			j = ind.getJ();
			k = ind.getK();
			l = ind.getL();
			break;
		case 1	:
			i = ind.getI();
			j = ind.getL();
			k = ind.getK();
			l = ind.getJ();
			break;
		case 2	:
			i = ind.getI();
			j = ind.getK();
			k = ind.getJ();
			l = ind.getL();
			break;
		default	: cerr << "ERROR: Cannot convert TwoElectronIntegralTriadeIndex" << endl;
			exit(-1);
	}
	
	if (j > i) tausche(i,j);
	if (l > k) tausche(k,l);
	//if (iscontained(i,j,k,l)) Ausgabe(i,j,k,l,ExTEI);
	ExTEI = sdot(bkji[i][j],bkji[k][l]);

	
//------------------
// Funktioniert leider nicht fehlerfrei
// --> leider wird nicht die richtige Reihenfolge der Indices eingehalten

// 	TwoElectronIntegralTriadeIndex Cbind(i,j,k,l,m);
// 	get(Cbind,CbTEI);
// 	
// 	TwoElectronIntegralTriadeIndex Exind(i,j,k,l,n);
// 	get(Exind,ExTEI);
//------------------	

	return irgendwas;
}


#endif
