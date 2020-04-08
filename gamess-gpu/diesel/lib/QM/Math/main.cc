#include "../../../config.h"

#include "SpinEigenFunctionDegeneration.h"



int	main()
{
SpinEigenFunctionDegeneration	*s;

	for ( INT i=0 ; i<10 ; i++ )
	{	s = new SpinEigenFunctionDegeneration(2, 20);
//		delete s;
	}
	
	cout << *s << endl;
}
