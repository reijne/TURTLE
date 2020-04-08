//***********************************************************************
//
//	Name:			IntegerStream.cc
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



#include "IntegerStream.h"

#include <stdlib.h>
#include <stdio.h>


IntegerStream::IntegerStream(char _separator)
{
	separator = _separator;
	compaction = 1;
}



IntegerStream::IntegerStream(String s, char _separator)
{
	separator = _separator;
	compaction = 1;

const INT	maxn = 1000;
String	res[maxn];
INT	n = split(s, res, maxn, Regex(String(separator)+"+"));
	for ( INT i=0 ; i<n ; i++ )
	{
		if ( res[i].length()==0 )
			continue;
	String	res2[maxn];
	INT	n2 = split(res[i], res2, maxn, String("-"));
	INT	i1;
        int     i1_2;
		if ( sscanf(res2[0], "%d", &i1_2)==1 )
		{
                        i1=i1_2;
		INT	i2 = i1;
			if ( n2>1 )
				i2 = atoi(res2[1]);
			for ( INT j=i1 ; j<=i2 ; j++ )
				append(j);
		}
	}
}


IntegerStream::~IntegerStream()
{
}

ostream & operator << (ostream & s, const IntegerStream & is)
{
	if ( is.compaction )
	{
	Pix	pix = is.first();
	INT last = -1000;
	INT	n = 0;
		while ( pix )
		{
		INT	i = is(pix);
			if ( i-1>last || i<last )
			{
				if ( n>1 )
					s << "-" << last;
				s << is.separator;
				s << i;
				n = 0;
			}
			n++;
			last = i;
			
			is.next(pix);
		}
		if ( n>1 )
			s << "-" << last;
	}
	else
	{
	Pix	pix = is.first();
		while ( pix )
		{
			s << is(pix) << " ";
			is.next(pix);
		}
	}
	return s;
}


istream & operator >> (istream & s, IntegerStream & is)
{
char	buf[10000];
	s.getline(buf, 10000);
	is = IntegerStream(String(buf),' ');
/*INT	i;

	while ( s >> i )
		is.append(i);
*/	return s;
}



