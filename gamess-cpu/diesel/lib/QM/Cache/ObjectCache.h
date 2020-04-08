//***********************************************************************
//
//	Name:	ReadOnlyCache.h
//
//	Description:	
//
//	Author:	Michael Hanrath
//
//	Version:	0.0
//
//	Date:	21.08.1996
//
//
//
//
//***********************************************************************


#ifndef __READONLYCACHE_H
#define __READONLYCACHE_H

#include "../../../config.h"



template <class KeyType, class ObjectType>
class ReadOnlyCache {
public:
	ReadOnlyCache(INT entries);
	~ReadOnlyCache();


	const ObjectType & 	operator[] (KeyType);
	
	
protected:
	const ObjectType &	cacheMiss(KeyType);



private:
INT	entries;


};














#endif

