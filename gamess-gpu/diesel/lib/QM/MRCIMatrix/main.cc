#include "../../../config.h"

#include <iostream>
#include <stdlib.h>

#include "../../Container/PackedIntStorage.h"

using namespace std;

int	main()
{
const INT n = 102402;
INT	ldBits = 3;


PackedIntStorage	pack(n, ldBits);

unsigned INT buf[n];


	for ( INT i=0 ; i<n ; i++ )
	{
		buf[i] = rand() & ((1<<(ldBits+1))-1);
//		cout << buf[i] << endl;
		pack.set(i, buf[i]);
	}

	for ( INT i=0 ; i<n ; i++ )
	{
		if ( buf[i] != pack.get(i) )
		{
			cout << n << ", " << buf[i] << " " << pack.get(i) << endl;
		}
	}
	

}

