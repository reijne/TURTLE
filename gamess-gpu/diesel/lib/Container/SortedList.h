#ifndef __SortedList_h
#define __SortedList_h

#include "../../config.h"

#include <algorithm>
#include <vector>
using std::vector;

template <class T>
class SortedList : public vector<T> 
{
public:
	SortedList() :
		maxEntries(0) {};
	SortedList(unsigned INT maxEntries) :
		maxEntries(maxEntries) {};

	void insert(const T & v)
	{
		if ( this->size()==maxEntries && max<v )
			return;
		typename vector<T>::iterator	p = lower_bound(this->begin(), this->end(), v);
//		cout << "v.conf= " << v.conf << ", v.e= " << v.e << "  p=" << p << endl;
//		cout << "end()" << this->end() << endl;
		if ( p!=this->end() )
		{
//			cout << "::::::::   " << p->conf << " " << p->e << endl;
		typename vector<T>::iterator	pu = upper_bound(this->begin(), this->end(), v);
			for ( ; p!=pu ; ++p )
			{
//				cout << "!!!!!!!!" << p->conf << " " << p->e << endl;
				if ( p->conf==v.conf )
					return;
			}
			
		}
//		cout << "insert!!" << endl;
		vector<T>::insert(p, v);
		if ( this->size()>maxEntries )
			erase(this->end()-1);
		max = *(this->end()-1);
	}
	
	void setMaxEntries(INT n)
	{
		maxEntries = n;
	}

private:
unsigned INT	maxEntries;
T	max;	
};

#endif
	
