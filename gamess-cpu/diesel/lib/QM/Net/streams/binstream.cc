using namespace std;

#include <iostream>


class	BinStream {

public:
	BinStream() { n=0; }
	
	unsigned char	operator [] (INT i) const {	return buf[i];	} 
	
	INT	getN() const	{	return n;	} 
	
	BinStream & operator <<(const INT);
	BinStream & operator <<(const float);
	BinStream & operator <<(const double);

	BinStream & operator >>(const INT);
	BinStream & operator >>(const float);
	BinStream & operator >>(const double);

	friend ostream& operator<<(ostream & s, const BinStream &);

private:
unsigned char	buf[100];
INT	n;
};




BinStream & BinStream::operator << (const INT a)
{
	for ( INT i=0 ; i<sizeof(INT) ; i++ )
		buf[n++] = ((unsigned char *) &a)[i];

	return *this;
}

BinStream & BinStream::operator << (const float a)
{
	for ( INT i=0 ; i<sizeof(float) ; i++ )
		buf[n++] = ((unsigned char *) &a)[i];

	return *this;
}

BinStream & BinStream::operator << (const double a)
{
	for ( INT i=0 ; i<sizeof(double) ; i++ )
		buf[n++] = ((unsigned char *) &a)[i];

	return *this;
}


ostream& operator<<(ostream & s, const BinStream & b)
{
	for ( INT i=0 ; i<b.getN() ; i++ )
		s << ((INT) b[i]) << " "; 
	return s;
}

#include <fstream>

main()
{
BinStream	b;

INT	i = 8765;
float	a = 3.1415926;
double	aa = a;

	b << i << a << aa;
	
	
	cout << b << endl;

	cout << "Ende" << endl;

ofstream	of("test.of");
	of << b;


	


}
