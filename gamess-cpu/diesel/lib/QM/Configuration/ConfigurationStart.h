//***********************************************************************
//
//	Name:			ConfigurationStart.h
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			10.06.1997
//
//
//
//
//
//***********************************************************************

#ifndef __ConfigurationStart_h
#define __ConfigurationStart_h

#include "../../../config.h"

#include "Configuration.h"




template <class TMOType>
class ConfigurationStart : 
	public Configuration<TMOType> {
public:
	ConfigurationStart() :
		Configuration<TMOType>() {	SAFStart = 0; ConfStart = 0;	}

	ConfigurationStart(const Configuration<TMOType> &conf) :
		Configuration<TMOType>(conf) {	SAFStart = 0; ConfStart = 0;	}

	INT	getSAFStart() const {	return SAFStart;	}
	void	setSAFStart(INT n) {	SAFStart = n;	}

	INT	getConfStart() const {	return ConfStart;	}
	void	setConfStart(INT n) {	ConfStart = n;	}

private:
INT	SAFStart;
INT	ConfStart;


};



#endif
