//***********************************************************************
//
//	Name:	VarSizeReadOnlyCache.cc
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



#include "VarSizeReadOnlyCache.h"


#include <stdlib.h>


template <class KeyType, class ObjectType>
VarSizeReadOnlyCache<KeyType, ObjectType>::VarSizeReadOnlyCache
	(INT _maxEntries, INT _maxDataMemory)
{
	maxEntries = _maxEntries;
	availableMemory = maxDataMemory = _maxDataMemory;
TAccessList	l;
	l.p = NULL;	
	l.key = NULL;	
	keyAVLMap = new AVLMap<KeyType, TAccessList>(l);
		
	timeCode = 0;
	cacheHits = cacheMisses = 0;
}


template <class KeyType, class ObjectType>
VarSizeReadOnlyCache<KeyType, ObjectType>::~VarSizeReadOnlyCache()
{
Pix	i = keyAVLMap->first();
	while ( i )
	{	delete ((ObjectType *) keyAVLMap->contents(i).p);
		delete ((KeyType *) keyAVLMap->contents(i).key);
		keyAVLMap->next(i);
	}
	delete keyAVLMap;
}


template <class KeyType, class ObjectType>
void	VarSizeReadOnlyCache<KeyType, ObjectType>::adjustTimeCodes()
{
unsigned INT	sub = findMinUsed()->used;
Pix	i = keyAVLMap->first();
	while ( i )
	{	keyAVLMap->contents(i).used -= sub;
		keyAVLMap->next(i);
	}
}

template <class KeyType, class ObjectType>
TAccessList *	VarSizeReadOnlyCache<KeyType, ObjectType>::findMinUsed()
{
unsigned INT min = INT_MAX;
TAccessList *	pmin = NULL;
Pix	i = keyAVLMap->first();
	while ( i )
	{	if ( keyAVLMap->contents(i).used<min )
		{	min = keyAVLMap->contents(i).used;
			pmin =  &keyAVLMap->contents(i);
		}
		keyAVLMap->next(i);
	}
	return pmin;
}



template <class KeyType, class ObjectType>
unsigned LONG_INT	VarSizeReadOnlyCache<KeyType, ObjectType>::
	getNumberOfCacheHits() const
{	return cacheHits;	}


template <class KeyType, class ObjectType>
unsigned LONG_INT	VarSizeReadOnlyCache<KeyType, ObjectType>::
	getNumberOfCacheMisses() const
{	return cacheMisses;	}


template <class KeyType, class ObjectType>
unsigned LONG_INT	VarSizeReadOnlyCache<KeyType, ObjectType>::
	getNumberOfCacheAccesses() const
{	return cacheHits + cacheMisses;	}


template <class KeyType, class ObjectType>
INT	VarSizeReadOnlyCache<KeyType, ObjectType>::getNumberOfEntries() const
{	return keyAVLMap->length();	}

template <class KeyType, class ObjectType>
INT	VarSizeReadOnlyCache<KeyType, ObjectType>::getMaxNumberOfEntries() const
{	return maxEntries;	}

template <class KeyType, class ObjectType>
INT	VarSizeReadOnlyCache<KeyType, ObjectType>::getOccupiedDataMemory() const
{	return	maxDataMemory-availableMemory;	}

template <class KeyType, class ObjectType>
INT	VarSizeReadOnlyCache<KeyType, ObjectType>::getMaxDataMemory() const
{	return	maxDataMemory;	}

template <class KeyType, class ObjectType>
INT	VarSizeReadOnlyCache<KeyType, ObjectType>::getTotalOccupiedMemory() const
{	return	maxDataMemory-availableMemory + 
		keyAVLMap->length()*sizeof(TAccessList);	}

template <class KeyType, class ObjectType>
void	VarSizeReadOnlyCache<KeyType, ObjectType>::clearStatistics()
{	cacheHits = cacheMisses = 0;	}


template <class KeyType, class ObjectType>
void	VarSizeReadOnlyCache<KeyType, ObjectType>::clearCache()
{	
Pix	i = keyAVLMap->first();
	while ( i )
	{	keyAVLMap->contents(i).used = keyAVLMap->contents(i).size = 0;
		delete ((ObjectType *) keyAVLMap->contents(i).p);
		delete ((KeyType *) keyAVLMap->contents(i).key);
		keyAVLMap->next(i);
	}
	timeCode = 0;
	availableMemory = maxDataMemory;
	cacheHits = cacheMisses = 0;
	keyAVLMap->clear();
}


//----------------------------------------------------------------------	


template <class KeyType, class ObjectType>
void	VarSizeReadOnlyCache<KeyType, ObjectType>::
	deleteLeastRecentlyUsed()
{
TAccessList	*p = findMinUsed();
	availableMemory += p->size;
//	cout << "avail: " << availableMemory << ", ";
//	cout << "addEntry: delete " << p << " " << p->p << " " << p->used << endl;
	delete ((ObjectType *) p->p);
KeyType *pk = (KeyType *) p->key;
//	cout << *pk << endl;
	keyAVLMap->del(*pk);
	delete pk;
}


template <class KeyType, class ObjectType>
const ObjectType *	VarSizeReadOnlyCache<KeyType, ObjectType>::
	addEntry(KeyType key)
{
TAccessList	*p;
	p = &(*keyAVLMap)[key];
	p->used = timeCode;
	p->key = new KeyType(key);
//	cout << "key: " << key << ", timeCode: " << timeCode << " " << p << endl;

	// allocate and generate new object
	p->p = (void *) new ObjectType(key);
	p->size = ((ObjectType *) p->p)->getOccupiedMemory();

//	cout << "O.K." << endl;
	availableMemory -= p->size;
	while ( availableMemory<0 )
		deleteLeastRecentlyUsed();
//	cout << "O.K.2" << endl;
	
	return ((ObjectType *) p->p);
}


//----------------------------------------------------------------------	


template <class KeyType, class ObjectType>
ostream& operator<<(ostream& s,
		const VarSizeReadOnlyCache<KeyType, ObjectType> & cache)
{
	s << "cache statistics" << endl;
	s << "================" << endl << endl;
	s << "hits                 : " << cache.getNumberOfCacheHits() << endl;
	s << "misses               : " << cache.getNumberOfCacheMisses() << endl;
	s << "entries              : " << cache.getNumberOfEntries() << endl;
	s << "max entries          : " << cache.getMaxNumberOfEntries() << endl;
	s << "occupied data memory : " << cache.getOccupiedDataMemory() <<	" Byte" << endl;
	s << "max data memory      : " << cache.getMaxDataMemory() << " Byte" << endl;
	s << "total occupied memory: " << cache.getTotalOccupiedMemory() << " Byte" << endl;
	return s;
}


#include "TableKey.h"

#include "../RepresentationMatrices/RepresentationMatrices.h"
#include "../RepresentationMatrices/HMatElements.h"
#include "../RepresentationMatrices/MRMPH0MatElements.h"

template class VarSizeReadOnlyCache<TableKey, RepresentationMatrices<float> >;
template class VarSizeReadOnlyCache<TableKey, RepresentationMatrices<double> >;

template ostream& operator << (ostream& s, const
	VarSizeReadOnlyCache<TableKey, RepresentationMatrices<float> > & v);
template ostream& operator << (ostream& s, const
	VarSizeReadOnlyCache<TableKey, RepresentationMatrices<double> > & v);


template class VarSizeReadOnlyCache<TableKey, HMatElements<float, float> >;
template class VarSizeReadOnlyCache<TableKey, HMatElements<double, float> >;
template class VarSizeReadOnlyCache<TableKey, HMatElements<float, double> >;
template class VarSizeReadOnlyCache<TableKey, HMatElements<double, double> >;

template ostream& operator << (ostream& s, const
	VarSizeReadOnlyCache<TableKey, HMatElements<float, float> > & v);
template ostream& operator << (ostream& s, const
	VarSizeReadOnlyCache<TableKey, HMatElements<double, float> > & v);
template ostream& operator << (ostream& s, const
	VarSizeReadOnlyCache<TableKey, HMatElements<float, double> > & v);
template ostream& operator << (ostream& s, const
	VarSizeReadOnlyCache<TableKey, HMatElements<double, double> > & v);


template class VarSizeReadOnlyCache<TableKey, MRMPH0MatElements<float, float> >;
template class VarSizeReadOnlyCache<TableKey, MRMPH0MatElements<double, float> >;
template class VarSizeReadOnlyCache<TableKey, MRMPH0MatElements<float, double> >;
template class VarSizeReadOnlyCache<TableKey, MRMPH0MatElements<double, double> >;

template ostream& operator << (ostream& s, const
	VarSizeReadOnlyCache<TableKey, MRMPH0MatElements<float, float> > & v);
template ostream& operator << (ostream& s, const
	VarSizeReadOnlyCache<TableKey, MRMPH0MatElements<double, float> > & v);
template ostream& operator << (ostream& s, const
	VarSizeReadOnlyCache<TableKey, MRMPH0MatElements<float, double> > & v);
template ostream& operator << (ostream& s, const
	VarSizeReadOnlyCache<TableKey, MRMPH0MatElements<double, double> > & v);


//**********************************************************************
//	explicit instantiation for g++ 2.7.2 
//	("template <class TMOType> void calcInteraction(...)" does not work)
/*
void	VarSizeReadOnlyCacheInstanciateTemplates()
{
	VarSizeReadOnlyCache<TableKey, RepresentationMatrices<float> >	cache2(0, 0);
		cout << cache2;
}
*/
	
#undef TEMPLATE_INSTANTIATE

//--------------------------------------------------------------------------
#include "../../Container/AVLMap.cc"


template class AVLMap<TableKey, TAccessList>;
template <> TableKey*   AVLMap<TableKey, TAccessList>::_target_item = NULL;     // add/del_item target
template <> AVLNode<TableKey, TAccessList>* 
	AVLMap<TableKey, TAccessList>::_found_node = NULL; // returned added/deleted node


class EnergyEntry;
template class AVLMap<Configuration<MOType>, EnergyEntry *>;
template <> Configuration<MOType>*   AVLMap<Configuration<MOType>, EnergyEntry *>::_target_item = NULL ;     // add/del_item target
template <> AVLNode<Configuration<MOType>, EnergyEntry *>* 
	AVLMap<Configuration<MOType>, EnergyEntry *>::_found_node = NULL ; // returned added/deleted node

//--------------------------------------------------------------------------
#include "../../Container/Map.cc"

template class Map<TableKey, TAccessList>;


class EnergyEntry;
template class Map<Configuration<MOType>, EnergyEntry *>;

