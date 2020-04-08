//***********************************************************************
//
//	Name:	TableKey.h
//
//	Description:	key for addressing representation matrices
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

#ifndef __TABLEKEY_H
#define __TABLEKEY_H

#include "../../../config.h"

#include "../Configuration/TableCase.h"

#include <iostream>

class TableKey {
public:
typedef unsigned LONG_INT TableKeyType;
	TableKey() : key(0)  {}
	~TableKey() {}
	
	// type conversion
	TableKey(const TableCase<MOType> &);
	TableKey(const TableCase<GeneralizedMO> &);
	
	

	INT	operator == (const TableKey &) const;
	INT	operator != (const TableKey &) const;
	INT	operator <(const TableKey &) const;
	INT	operator <= (const TableKey &) const;
	INT	operator > (const TableKey &) const;
	INT	operator >= (const TableKey &) const;

//----------------------------------------------------------------------	

	TableKeyType	getKey() const;
	TableCase<MOType> getTableCase() const;

//----------------------------------------------------------------------	
	
	friend ostream& operator<<(ostream& s, const TableKey & key);
	
private:
//	coding of key:
//
//		parameter		possible values		next 2^n	n	start
//-----------------------------------------------------------------
//		qr				[1, ..., 999]		1024		10	   0
//		ql				[1, ..., 999]		1024		10	  10
//		inopen			[0, ..., 15]		16			4	  20
//		pfa				[1, ..., 4]			4			2	  24
//      rfa				[1, ..., 6]			8			3	  26
//      dsk				[-2, -1, 0, 1, 2]	8			3	  29
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!!!! Coding of q-cases in ten bits is sufficient to 10 open shells !!!!
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

TableKeyType	key;
};









inline
TableKey::TableKey(const TableCase<MOType> & tc)
{
//	cout << "TableCase<MOType>: " << tc << endl;
	if ( tc.getP()!=5 )
		key = 	 ((TableKeyType) tc.getqR()) |
				(((TableKeyType) tc.getqL()) << 10) |
				(((TableKeyType) tc.getNumberOfMoreOpenShells()) << 20) |
				(((TableKeyType) (tc.getP()-1)) << 24) |
				(((TableKeyType) tc.getR()) << 26) |
				(((TableKeyType) tc.getdK()) << 29);
	else
		key = 	(((TableKeyType) tc.getNumberOfMoreOpenShells()) << 20) |
				(((TableKeyType) 7) << 26);
	
//	cout << *this << endl;
}


inline
TableKey::TableKey(const TableCase<GeneralizedMO> & tc)
{
//	cout << "TableCase<GeneralizedMO>: " << tc << endl;
	if ( tc.getP()!=5 )
		key = 	 ((TableKeyType) tc.getqR()) |
				(((TableKeyType) tc.getqL()) << 10) |
				(((TableKeyType) tc.getNumberOfMoreOpenShells()) << 20) |
				(((TableKeyType) (tc.getP()-1)) << 24) |
				(((TableKeyType) tc.getR()) << 26) |
				(((TableKeyType) tc.getdK()) << 29);
	else
		key = 	(((TableKeyType) tc.getNumberOfMoreOpenShells()) << 20) |
				(((TableKeyType) 7) << 26);
//	cout << *this << endl;
}




inline
TableCase<MOType>	TableKey::getTableCase() const
{
/*	cout << "B: " << key << " " << ((INT) key) << " " <<
		(((INT) key) & ((INT) 0xe0000000)) << " " <<
		((((INT) key) & ((INT) 0xe0000000)) >> 29) << endl;
*/	
	if ( ((key >> 26) & 7) < 7 )
		return TableCase<MOType>(
			(key >> 20) & 15,
			(((int) key) & ((int) 0xe0000000)) >> 29,
			((key >> 24) & 3) + 1,
			(key >> 26) & 7,
			key & 1023,
			(key >> 10) & 1023);
	else
		return TableCase<MOType>(
			(key >> 20) & 15,
			0,
			5,
			1,
			0,
			0);
}



inline
TableKey::TableKeyType	TableKey::getKey() const
{	return	key;	}


inline
INT	TableKey::operator == (const TableKey & k) const
{	return key==k.getKey();	}

inline
INT	TableKey::operator != (const TableKey & k) const
{	return key!=k.getKey();	}

inline
INT	TableKey::operator < (const TableKey & k) const
{	return key<k.getKey();	}

inline
INT	TableKey::operator <= (const TableKey & k) const
{	return key<=k.getKey();	}

inline
INT	TableKey::operator > (const TableKey & k) const
{	return key>k.getKey();	}

inline
INT	TableKey::operator >= (const TableKey & k) const
{	return key>=k.getKey();	}



#endif
