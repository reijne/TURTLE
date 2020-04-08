//***********************************************************************
//
//	Name:			ConfigurationStart.cc
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




#include "ConfigurationStart.h"




#include "../../Container/AVLSet.cc"
template class AVLSet<ConfigurationStart<MOType> >;

#include "../../Container/Set.cc"
template class Set<ConfigurationStart<MOType> >;
