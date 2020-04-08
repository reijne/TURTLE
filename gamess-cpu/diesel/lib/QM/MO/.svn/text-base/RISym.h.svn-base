// Kommentare im Definitionsteil beachten !!
#ifndef __RISYM_H
#define __RISYM_H __RISYM_H

#include <math.h>
#include <iostream>
#include <string>
using std::string;

#include <stdio.h>
#include <fstream>

#include "../../../config.h"

#include "../Symmetry/IrRep.h"

class RISym {
public:
  RISym();  // fehlt noch
  RISym(const string Pfad);  // ACHTUNG statt Type string: Type Fort31File !!
  ~RISym();
  // Funktionen zur Verwaltung der Daten // f"urs erste Ueberfl"ussig
  //  INT getinIrreps();
  //  void setinIrreps();
  //  INT getIrReps();
  //  void setIrReps();
  // Funktionen, die im diesel vorkommen
  INT getN() const;
  INT   getLJ(INT i) const;
  IrRep   getJAB(INT i) const;
private:
  INT IrReps;
  INT *inIrRep;
  string Pfadname;
  string Filename;
  string Dateiname;
  IrRep	*prodTab;		// product table
};
  

inline
INT RISym::getN() const
{
  return IrReps;
}

inline 
INT  RISym::getLJ(INT i) const
{
  return inIrRep[i];
}

inline
IrRep  RISym::getJAB(INT i) const
{
  //  index von quadratisch auf dreieckig umrechnen
/*   INT h1 = i,delta = 0; */
/*   while (h1 > 0) */
/*     { */
/*       delta++; */
/*       h1 -= delta; */
/*     } */
/*   h1 += delta; */
/*   delta--; */
/*   INT neuIndex = h1*IrReps + delta; */
  return (prodTab[i]+1);
}



#endif
