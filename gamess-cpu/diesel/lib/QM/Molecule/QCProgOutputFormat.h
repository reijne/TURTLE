//***********************************************************************
//
//	Name:			QCProgOutputFormat.h
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			26. Feb 1998
//
//***********************************************************************

#ifndef __QCProgOutputFormat
#define __QCProgOutputFormat

//FD class	istream;

#include "../../../config.h"

#include "../../Container/String.h"

#include <iostream>
using std::istream;

class	QCProgOutputFormat {
public:

enum	QCPackage { undefinedPackage, Hondo, MOLPRO, MOLCAS,
		TURBOMOLE, Gaussian, GaussianWFN, IntCoord };

enum	MOLCASPart { undefinedMOLCASPart, seward, scf, rasscf };

union TPart {
	MOLCASPart	part;
	};

	QCProgOutputFormat();
	QCProgOutputFormat(istream &);
	~QCProgOutputFormat();

	String	getIdentifier() const;

	QCPackage	getPackage() const;
	TPart	getPart() const;

private:
	INT	isInteger(const String &) const;
	INT	isFloat(const String &) const;

QCPackage	package;
String	version;
TPart	part;
};











#endif
