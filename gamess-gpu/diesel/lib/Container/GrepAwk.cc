//***********************************************************************
//
//	Name:			GrepAwk.cc
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			28. Jul 1998
//
//***********************************************************************




#include "GrepAwk.h"

#include <iostream>
#include <string>



GrepAwk::GrepAwk(std::istream &i, INT maxLines)
{
const INT bufsize = 10000;
char	buf[bufsize];
INT	iline = 0;
	while ( i.getline(buf, bufsize) && (maxLines==-1 || iline<maxLines) )
	{
		for ( INT j=0 ; j<((INT) strlen(buf)) ; j++ )
			if ( buf[j]=='\t' )
				buf[j] = ' ';
	String	*S = new String(buf);
		append(S);
		iline++;
	}
	line = first();
	ignoreCase = 0;
}


GrepAwk::GrepAwk(FILE *fin, INT maxLines)
{
const INT bufsize = 10000;
char	buf[bufsize];
INT	iline = 0;
	while ( !feof(fin) )
	{
		if( !fgets(buf, bufsize, fin) )
			break;
	String	*S = new String(buf);
		append(S);
		iline++;
	}

	line = first();
	ignoreCase = 0;
}



GrepAwk::~GrepAwk()
{
Pix	pix;
	pix = first();
	while ( pix )
	{
		delete (*this)(pix);
		next(pix);
	}
}

INT	GrepAwk::illegal() const
{
	return	line==NULL;
}



Pix	GrepAwk::getIndex() const
{
	return	line;
}


void	GrepAwk::setIndex(Pix _line)
{
	line = _line;
}



String	GrepAwk::getLine() const
{
	return *(*((GrepAwk *) this))(line);
}

INT	GrepAwk::getLineLength() const
{
	return (*((GrepAwk *) this))(line)->length();
}

INT	GrepAwk::getNumberOfWords() const
{
INT	p = 0;
INT	n = 0;
String	s = getLine();
	while ( p<((INT) s.length()) && s[p]==' ' )
		p++;
	while ( p<((INT) s.length()) )
	{
		n++;
		while ( p<((INT) s.length()) && s[p]!=' ' )
			p++;
		while ( p<((INT) s.length()) && s[p]==' ' )
			p++;
	}
	return n;
}

INT	GrepAwk::getWordPos(INT i) const
{
INT	pos = 0;
	getWords(i, pos, 1);
	return pos;
}

String	GrepAwk::getWord(INT i) const
{
INT	pos;
	return getWords(i, pos, 1);
}

String	GrepAwk::getWord(INT i, INT &pos) const
{
	return getWords(i, pos, 1);
}

String	GrepAwk::getWords(INT i, INT n) const
{
INT	pos;
	return getWords(i, pos, n);
}

String	GrepAwk::getWords(INT i, INT &pos, INT n) const
{
String	s(*(*((GrepAwk *) this))(line));
String	r;
INT	p = 0;

	while ( p<((INT) s.length()) && s[p]==' ' )
		p++;
	for ( INT j=0 ; j<i-1 ; j++ )
	{
		p = s.index(' ', p);
		if ( p<0 )
			return r;
		else
			while ( p<((INT) s.length()) && s[p]==' ' )
				p++;
	}
	
INT	pold = p;
	pos = p;
	for ( INT j=0 ; j<n ; j++ )
	{
		p = s.index(' ', pold);
		if ( p<0 )
		{
			r += s.at(pold, s.length()-pold);
			return r;
		}
		r += s.at(pold, p-pold);
		while ( p<((INT) s.length()) && s[p]==' ' )
			p++;
		pold = p;
	}
	
	return r;
}

INT	GrepAwk::grep(String subString, GrepAwk::Direction dir)
{
	if ( ignoreCase )
		subString.capitalize();
	while ( line )
	{
		if (
			(*((GrepAwk *) this))(line)->index(subString)>=0 && !ignoreCase
			||
			capitalize(*(*((GrepAwk *) this))(line)).index(subString)>=0 && ignoreCase
			)
			return 1;
			
		if ( dir==Forward )
			next(line);
		else
			prev(line);
	}
	return 0;
}


void	GrepAwk::setIgnoreCase(INT i)
{
	ignoreCase = i;
}

void	GrepAwk::head()
{
	line = first();
}


void	GrepAwk::tail()
{
	line = last();
}
	
void	GrepAwk::skipBehindNextBlankLine()
{
	while ( !illegal() && getLine().length()>0 )
		(*this)++;
		
	while ( !illegal() && getLine().length()==0 )
		(*this)++;
}

GrepAwk &	GrepAwk::operator++(int)
{
	next(line);
	return *this;
}


GrepAwk &	GrepAwk::operator--(int)
{
	prev(line);
	return *this;
}


GrepAwk & GrepAwk::operator+=(INT i)
{
	for ( ; i>0 && line ; i-- )
		next(line);
	return *this;
}
	 

GrepAwk & GrepAwk::operator-=(INT i)
{
	for ( ; i>0 && line ; i-- )
		prev(line);
	return *this;
}


#include "DLList.cc"
template class DLList<String *>;
