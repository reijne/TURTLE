#include "DieselInput.h"

#include "../../../common/StartUp.h"

#include "../../../../lib/Container/String.h"
#include <sstream>
#include <fstream>
#include <stdlib.h>

#include <iomanip>

#include "../../../../lib/QM/IO/TimeTicks.h"

#include "../../../../lib/QM/IO/Verbosity.h"

#include "../../../../Configured.h"
#include "Compiled.h"
#include "../../../common/Banner.h"

#include "../../../../config.h"
#include "../../../../VersionDate.h"

using namespace std;

void	MakeBanner()
{
const INT w = 80;

	MakeBannerTop(w, "D I E S E L   M R - C I   D r i v e r");
	center(VERSION, w);
	center(DATE, w);
	MakeBannerBottom(w);
}




TimeTicks	globalTime;
//==========================================================================




int main(INT argc, char **argv)
{
	globalTime.start();

	MakeBanner();


	StartUp();
/*
*/

			
DieselInput	mrinp;


	cin >> mrinp;
	cout << mrinp << endl;

	mrinp.start();

}


