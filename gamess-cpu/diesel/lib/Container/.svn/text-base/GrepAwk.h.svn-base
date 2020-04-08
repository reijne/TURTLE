//***********************************************************************
//
//	Name:			GrepAwk.h
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




#ifndef __GrepAwk_h
#define __GrepAwk_h

#include "../../config.h"

#include "DLList.h"
#include "String.h"
#include <stdio.h>

//class	istream;

class GrepAwk :
	private DLList<String *> {
public:
friend class String;
friend class SubString;
enum	Direction { Forward, Backward };

	GrepAwk(std::istream & i, INT maxLines = -1);
	GrepAwk(FILE *, INT maxLines = -1);
	~GrepAwk();


	INT	illegal() const;
	
	void	setIgnoreCase(INT i);
	
	Pix	getIndex() const;
	void	setIndex(Pix);
	
	String getLine() const;

	INT	getNumberOfWords() const;
	INT	getLineLength() const;

	INT getWordPos(INT i) const;
	String getWord(INT i) const;
	String getWord(INT i, INT &pos) const;
	String getWords(INT i, INT n=1) const;
	String getWords(INT i, INT &pos, INT n) const;

	INT	grep(String, Direction dir = Forward);
	
	void	head();
	void	tail();
	void	skipBehindNextBlankLine();
	
	GrepAwk & operator++(int);
	GrepAwk & operator--(int);
	GrepAwk & operator+=(INT i);
	GrepAwk & operator-=(INT i);




private:
	GrepAwk(const GrepAwk &);
Pix	line;
INT	ignoreCase;
};




#endif
