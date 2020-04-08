//***********************************************************************
//
//	Name:			IntegerStream.h
//
//	Description:	integer stream IO
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			29.04.1998
//
//
//
//
//
//***********************************************************************

#ifndef __IntegerStream_h
#define __IntegerStream_h

#include "../../../config.h"


#include "../../Container/SLList.h"
#include "../../Container/String.h"

using std::ostream;
using std::istream;

class IntegerStream : public SLList<INT> {
public:
	IntegerStream(char separator = ' ');
	IntegerStream(String s, char separator = ' ');
	~IntegerStream();


	friend ostream & operator << (ostream &, const IntegerStream &);
	friend istream & operator >> (istream &, IntegerStream &);


private:
char	separator;
INT	compaction;			// flag if consecutive number should be grouped by "1-10"
};



#endif


