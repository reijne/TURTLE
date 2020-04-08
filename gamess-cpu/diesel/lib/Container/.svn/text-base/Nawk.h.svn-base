#ifndef __Nawk_h
#define __Nawk_h


#include <list>
#include <string>
#include <sstream>
#include <iostream>
#include <stdio.h>

#include "../../config.h"

using std::string;
using std::list;
using std::istream;

class Nawk : private list<string> 
{
public:
	enum	Direction { Forward, Backward };

	Nawk(istream & i, INT maxLines = -1);
	Nawk(FILE *, INT maxLines = -1);
	Nawk(string s, INT maxLines = -1);
	~Nawk();


//	INT		illegal() const;
	bool	illegal();
	
	void	setIgnoreCase(INT i);
	
	list<string>::iterator	getIndex() const;
	void	setIndex(list<string>::iterator);
	
	string getLine() const;

	INT	getNumberOfWords() const;
	INT	getLineLength() const;

	INT getWordPos(INT i) const;
	string getWord(INT i) const;
	string getWord(INT i, INT &pos) const;
	string getWords(INT i, INT n=1) const;
	string getWords(INT i, INT &pos, INT n) const;

	INT	grep(string, Direction dir = Forward);
	
	void	head();
	void	tail();
	void	skipBehindNextBlankLine();
	
	Nawk & operator++(int);
	Nawk & operator--(int);
	Nawk & operator+=(INT i);
	Nawk & operator-=(INT i);

	void writeToStream();



private:
	Nawk(const Nawk &);
	list<string>::iterator	line;
	INT	ignoreCase;
};




#endif
