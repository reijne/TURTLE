#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <ctype.h>

#include "../../../../lib/QM/Configuration/Configuration.h"

#include "../../../../config.h"

using namespace std;

void	del(Configuration<MOType> & conf, INT k, INT n)
{
INT	mo;
	for ( INT i=0 ; i<n ; i++ )
	{
		conf.annihilate(1+i+k);
		conf.annihilate(1+i+k);
	}

	for ( INT i=0 ; i<conf.getNumberOfOpenShells() ; i++ )
		conf.setOpenShell(i, (mo=conf.getOpenShell(i))>=k ? mo-n : mo);
	for ( INT i=0 ; i<conf.getNumberOfClosedShells() ; i++ )
		conf.setClosedShell(i, (mo=conf.getClosedShell(i))>=k ? mo-n : mo);

}

void	ins(Configuration<MOType> & conf, INT k, INT n, INT f)
{
INT	mo;
	for ( INT i=0 ; i<conf.getNumberOfOpenShells() ; i++ )
		conf.setOpenShell(i, (mo=conf.getOpenShell(i))>=k ? mo+n : mo);
	for ( INT i=0 ; i<conf.getNumberOfClosedShells() ; i++ )
		conf.setClosedShell(i, (mo=conf.getClosedShell(i))>=k ? mo+n : mo);
	if ( f )
		for ( INT i=0 ; i<n ; i++ )
		{
			conf.create(1+i+k);
			conf.create(1+i+k);
		}
}


int	main(INT argc, char **argv)
{

	if ( !strcmp(argv[1], "-f") )
	{
	Configuration<MOType>	conf;


	INT	map[1000];
		for ( INT i=0 ; i<1000 ; i++ )
			map[i] = i;

	INT	irreps;

		cin >> irreps;

	INT	inIrrep[irreps];

		for ( INT i=0 ; i<irreps ; i++ )
			cin >> inIrrep[i];

	INT	frozen1[irreps];	
	INT	deleted1[irreps];	
	INT	frozen2[irreps];	
	INT	deleted2[irreps];	

		for ( INT i=0 ; i<irreps ; i++ )
			cin >> frozen1[i];
		for ( INT i=0 ; i<irreps ; i++ )
			cin >> deleted1[i];
		for ( INT i=0 ; i<irreps ; i++ )
			cin >> frozen2[i];
		for ( INT i=0 ; i<irreps ; i++ )
			cin >> deleted2[i];


		while ( cin )
		{
			cin >> conf;
		INT	s = 0;
			for ( INT i=0 ; i<irreps ; i++ )
			{
				ins(conf, s, frozen1[i], 1);
				s += inIrrep[i];
				ins(conf, s-deleted1[i], deleted1[i], 0);
			}

			for ( INT i=irreps-1 ; i>=0 ; i-- )
			{
				del(conf, s-deleted2[i], deleted2[i]);
				s -= inIrrep[i];
				del(conf, s, frozen2[i]);

			}
			cout << conf << endl;
		}

		return 0;
	}


	if ( !strcmp(argv[1], "-m") )
	{

	INT n;
	cin >> n;

	INT	map[1000];
		for ( INT i=0 ; i<1000 ; i++ )
			map[i] = i;
	for ( INT i=0 ; i<n ; i++ )
	{
	INT	j, k;
		cin >> j >> k;
		map[j] = k;
	}
	Configuration<MOType>	conf;

		while ( !cin.eof() )
		{
		char c = cin.get();
			if ( !c || c==EOF )
				break;
			if ( !isdigit(c) )
			{
			char	buf[10000];
				cin.getline(buf, 10000);
				cout << "#" << buf << endl;
				continue;
			}
			cin.putback(c);

			cin >> conf;

			for ( INT i=0 ; i<conf.getNumberOfOpenShells() ; i++ )
				conf.setOpenShell(i, map[conf.getOpenShell(i)]);

			for ( INT i=0 ; i<conf.getNumberOfClosedShells() ; i++ )
				conf.setClosedShell(i, map[conf.getClosedShell(i)]);

			conf.sort();

			cout << conf << endl;
		}

		return 0;
	}

}
