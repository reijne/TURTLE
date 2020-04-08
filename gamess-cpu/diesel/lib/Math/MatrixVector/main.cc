#include "../../../config.h"

#include "Matrix.h"
#include "BufferedVector.h"


int main()
{
//	BufferedVector<double>	hhh1(40000);

const INT	dim = 10000000;

BufferedVector<double>	v(dim);

	v[0] = 10;
	v[40] = 15;

//	cout << v << endl;

	cout << "START" << endl;
	cout << v*v << endl;


BufferedVector<double>	v2(dim);

	v2 = v;

//	cout << v2 << endl;

/*
double	p[dim];

	v2.get(p);
	for ( INT i=0 ; i<dim ; i++ )
		cout << p[i] << " ";
	cout << endl;
*/
}

/*
main()
{
const INT	dim = 3;
Matrix<double>	A(dim, dim);


	A(1, 1) = 3;
	A(1, 2) = 2;
	A(1, 3) = 5;
	A(2, 1) = 8;
	A(2, 2) = 4;
	A(2, 3) = 2;
	A(3, 1) = 9;
	A(3, 2) = 5;
	A(3, 3) = 7;
	cout << A << endl;

Matrix<double>	B = A.invert();
	
	cout << B << endl;
	
	cout << A*B << endl;
	


}
*/
