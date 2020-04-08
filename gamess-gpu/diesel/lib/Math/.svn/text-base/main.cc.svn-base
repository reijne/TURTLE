using namespace std;

#include "../../config.h"
#include "Histogram.h"

#include <stdlib.h>


int	main()
{
Histogram<double>	hist(1e-5, 1.0, Histogram<double>::Logarithmic);

	for ( INT i=0 ; i<100000 ; i++ )
		hist.count(1.0*rand()/RAND_MAX)-0.5;
		
	hist.printPrettyASCII();
	


}
