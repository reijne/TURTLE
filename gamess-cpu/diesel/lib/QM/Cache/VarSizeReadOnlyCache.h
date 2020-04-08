//***********************************************************************
//
//	Name:	VarSizeReadOnlyCache.h
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

#ifndef __VARSIZEREADONLYCACHE_H
#define __VARSIZEREADONLYCACHE_H

#include "../../../config.h"

#include "../../Container/AVLMap.h"

#include "CacheableObject.h"

#include "TAccessList.h"

#include <iostream>
using std::ostream;
using std::endl;

template <class KeyType, class ObjectType> class VarSizeReadOnlyCache;
template <class KeyType, class ObjectType> 
         ostream& operator<< (ostream& s, const VarSizeReadOnlyCache<KeyType, ObjectType> & cache);




template <class KeyType, class ObjectType>
class VarSizeReadOnlyCache {
public:
	VarSizeReadOnlyCache(INT entries, INT maxDataMemory);
	~VarSizeReadOnlyCache();

//----------------------------------------------------------------------	

	const ObjectType * 	operator[] (KeyType);

//----------------------------------------------------------------------	

	unsigned LONG_INT getNumberOfCacheHits() const;
	unsigned LONG_INT getNumberOfCacheMisses() const;
	unsigned LONG_INT getNumberOfCacheAccesses() const;
	
	INT	getNumberOfEntries() const;
	INT	getMaxNumberOfEntries() const;
	INT	getOccupiedDataMemory() const;
	INT	getMaxDataMemory() const;
	INT	getTotalOccupiedMemory() const;
	
	void	clearStatistics();

	void	clearCache();
	
//----------------------------------------------------------------------	
	
	friend ostream& operator<< <KeyType, ObjectType> (ostream& s,
		const VarSizeReadOnlyCache<KeyType, ObjectType> & cache);
	
//----------------------------------------------------------------------	

private:
	void				adjustTimeCodes();
	TAccessList *		findMinUsed();
	void				deleteLeastRecentlyUsed();
	const ObjectType *	addEntry(KeyType key);
	
	
INT	maxEntries;				// maximum number of entries in cache
INT	maxDataMemory;			// maximum amount of memory in byte used by data


INT	availableMemory;		// available memory

unsigned INT	timeCode;	// for determination of least recently used entry

unsigned LONG_INT	cacheHits;	// counter for cache hits
unsigned LONG_INT	cacheMisses;// counter for cache misses
AVLMap<KeyType, TAccessList>
	*keyAVLMap;				// AVLMap implementation of GNU C++ Library

KeyType	lastKey;			// key of last accessed object
const ObjectType
	*lastObject;			// pointer to last accessed object
};


template <class KeyType, class ObjectType>
ostream& operator<<(ostream& s,
		const VarSizeReadOnlyCache<KeyType, ObjectType> & cache);
	




#include <limits.h>



template <class KeyType, class ObjectType>
inline
const ObjectType *	VarSizeReadOnlyCache<KeyType, ObjectType>::operator[]
	(KeyType key)
{
//	cout << "key, lastkey " << key << " " << lastKey << endl;
	if ( key == lastKey )
	{	cacheHits++;
//		cout << "operator[] return" << endl;
		return lastObject;
	}

	if ( ++timeCode == INT_MAX )
		adjustTimeCodes();
		
	lastKey = key;

Pix	pix = keyAVLMap->seek(key);
TAccessList	*ind;
	if ( pix )
	{	// cache hit
//		cout << "hit" << endl;
		cacheHits++;
		ind = &keyAVLMap->contents(pix);
		ind->used = timeCode;
		return (lastObject = ((ObjectType *) ind->p));
	}
	// cache miss
	cacheMisses++;
//	cout << "miss" << endl;
	if ( keyAVLMap->length()<maxEntries )
		return (lastObject = addEntry(key));
	
	// no more free entries
	deleteLeastRecentlyUsed();
	return (lastObject = addEntry(key));
}


#endif

