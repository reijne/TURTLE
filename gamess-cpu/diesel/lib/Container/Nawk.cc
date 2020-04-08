#include "Nawk.h"

using namespace std;

string&	capitalize(string s)
{
	cerr << "Function capitalize() not yet installed!!!" << endl;
	return s;
}

Nawk::Nawk(istream &i, INT maxLines)
{
	const INT 	bufsize = 10000;
	char		buf[bufsize];
	INT			iline = 0;
	while ( i.getline(buf, bufsize) && (maxLines==-1 || iline<maxLines) )
	{
		for ( INT j=0 ; j<((INT)( strlen(buf))) ; j++ )
			if ( buf[j]=='\t' )
				buf[j] = ' ';
		string	S(buf);
		push_back(S);
		iline++;
	}
	//cout << "ilines= " << iline << endl;
	line = begin();
	//cout << "*line= " << *line << endl;
	ignoreCase = 0;
}


Nawk::Nawk(string i, INT maxLines)
{
	istringstream is(i.c_str(), ostringstream::in | ostringstream::out);
	Nawk(is,maxLines);
}


Nawk::Nawk(FILE *fin, INT maxLines)
{
	const INT 	bufsize = 10000;
	char		buf[bufsize];
	INT			iline = 0;
	while ( !feof(fin) )
	{
		if( !fgets(buf, bufsize, fin) )
			break;
		string	S(buf);
		push_back(S);
		iline++;
	}

	line = begin();
	ignoreCase = 0;
}

void Nawk::writeToStream()
{
	for (list<string>::iterator	pix = begin(); pix != end(); pix++)
		cout << *pix << endl;
}


Nawk::~Nawk()
{
	list<string>::iterator	pix;
	pix = begin();
	while ( pix != end())
	{
		pop_front();
		pix = begin();
	}
}

// INT	Nawk::illegal() const
// {
// 	return	INT(!(line!=end()));
// }

bool	Nawk::illegal() 
{
	return	(!(line!=end()));
}



list<string>::iterator	Nawk::getIndex() const
{
	return	line;
}


void	Nawk::setIndex(list<string>::iterator _line)
{
	line = _line;
}



string	Nawk::getLine() const
{
	//cout <<"*line= " << *line << endl;
	return *line;
}

INT	Nawk::getLineLength() const
{
	return line->length();
}

INT	Nawk::getNumberOfWords() const
{
INT	p = 0;
INT	n = 0;
// 	cout << "gok1a"<<endl;
// 	cout << "getLine(): "<<endl;
// 	cout <<getLine()<<endl;
// 	cout << "gok1a2"<<endl;
// 	cout << *begin() << endl;
// 	cout << "gok1b"<<endl;
// 	cout << *line << endl;
// 	cout << "gok2a"<<endl;
string	s = getLine();
//	cout << "gok2b"<<endl;
	while ( (p<(((INT)( s.length())))) && s[p]==' ' )
		p++;
//	cout << "gok3"<<endl;
	while ( p<((INT) (s.length())))
	{
		n++;
		while ( (p<((INT)( s.length()))) && s[p]!=' ' )
			p++;
		while ( (p<((INT)(s.length()))) && s[p]==' ' )
			p++;
	}
	return n;
}

INT	Nawk::getWordPos(INT i) const
{
INT	pos = 0;
	getWords(i, pos, 1);
	return pos;
}

string	Nawk::getWord(INT i) const
{
INT	pos;
	return getWords(i, pos, 1);
}

string	Nawk::getWord(INT i, INT &pos) const
{
	return getWords(i, pos, 1);
}

string	Nawk::getWords(INT i, INT n) const
{
INT	pos;
	return getWords(i, pos, n);
}

string	Nawk::getWords(INT firstWord, INT &pos, INT numberOfWords) const
{
	string	s(*line);
	string	r;
	INT	p = 0;

	while ((p<(s.length()) && s[p]==' '))
		p++;
	for ( INT j=0 ; j<firstWord-1 ; j++ )  //gehe zum richtigen Wort
	{
		while ((p<(s.length()) && s[p]!=' '))
			p++;
		if ( (p>=s.length()) || p < 0 )//Ende des strings fr"uher erreicht als gedacht
			return r;
		else  //"uberspringe die blanks
			while ((p<(s.length()) && s[p]!=' '))
				p++;
	}
	
	unsigned INT	pold = p;
	pos = p;
	for ( INT j=0 ; j<numberOfWords ; j++ )  // nehme die verlangte Anzahl an W"ortern
	{
		p = s.find(" ", pold);  
		if ( (p>=s.length()) || p < 0)  //Ende des strings fr"uher erreicht als gedacht
		{
			r += s.substr(pold, s.length()-pold);  //gesamter verbleibender string wird zur"uckgegeben
			return r;
		}
		r += s.substr(pold, p-pold);   // Der Teil zwischen pold und p wird angeh"angt
		while ((p<(s.length()) && s[p]==' '))
			p++;
		pold = p;
	}
	
	return r;
}

INT	Nawk::grep(string substring, Nawk::Direction dir)
{
	//if ( ignoreCase )
	//	substring=capitalize(substring);
	while ( line != end())
	{
		//cout << **line << endl;
		if (
			line->find(substring)<(line->length())//&& !ignoreCase
// 			||
// 			capitalize(**line)).find(substring)>=0 && ignoreCase
			)
			return 1;
			
		if ( dir==Forward )
			line++;
		else
			line--;
	}
	return 0;
}


void	Nawk::setIgnoreCase(INT i)
{
	ignoreCase = i;
}

void	Nawk::head()
{
	line = begin();
}


void	Nawk::tail()
{
	line = end();
}
	
void	Nawk::skipBehindNextBlankLine()
{
	while ( !illegal() && getLine().length()>0 )
		(*this)++;
		
	while ( !illegal() && getLine().length()==0 )
		(*this)++;
}

Nawk &	Nawk::operator++(int)
{
	line++;
	return *this;
}


Nawk &	Nawk::operator--(int)
{
	line--;
	return *this;
}


Nawk & Nawk::operator+=(INT i)
{
	for ( ; i>0 && line != end(); i-- )
		line++;
	//line+=i;
	return *this;
}
	 

Nawk & Nawk::operator-=(INT i)
{
	for ( ; i>0 && line != begin(); i-- )
		line--;
	//line+=i;
	return *this;
}


template class list<string>;
