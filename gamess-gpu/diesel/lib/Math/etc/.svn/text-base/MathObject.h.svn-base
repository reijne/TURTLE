//***********************************************************************
//
//	Name:			MathObject.h
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			27. Jan 1998
//
//***********************************************************************

#ifndef __MATHOBJ_H
#define __MATHOBJ_H

#include "../../../config.h"

#include <iostream>
using std::ostream;


class	MathObject {
public:
enum OutputMode {TerminalASCII, TeX};



	void	setOutputMode(OutputMode o);
	OutputMode	getOutputMode() const;

	friend ostream& operator<<(ostream& s, const MathObject);

protected:

static OutputMode	outputMode;	
};



inline
void	MathObject::setOutputMode(OutputMode o)
{	outputMode = o;	}

inline
MathObject::OutputMode	MathObject::getOutputMode() const
{	return outputMode;	}



#endif

