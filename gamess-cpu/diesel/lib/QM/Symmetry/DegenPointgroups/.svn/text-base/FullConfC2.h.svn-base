#ifndef __FullConfC2_h
#define __FullConfC2_h

#include <vector>

#include "../../../../config.h"
#include "FullConf.h"


class FullConfC2 : public FullConf {
public:
	FullConfC2(const FullConf &fc, INT irreps, INT roots) :
		FullConf(fc)
	{
		c2.resize(irreps);
		for ( INT i=0 ; i<irreps ; i++ )
			c2[i].resize(roots);
	}
	
	FullConfC2(const FullConf &fc) :
		FullConf(fc) {}

vector<vector<double> >	c2;			// c2[irrep][root]
private:

};

ostream & operator << (ostream &s, const FullConfC2 &f);



#endif
