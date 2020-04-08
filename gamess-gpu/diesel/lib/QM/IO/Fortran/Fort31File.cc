//***********************************************************************
//
//	Name:			Fort31File.cc
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			10. Jun 1998
//
//***********************************************************************

#include "Fort31File.h"
#include "Fort31RecordFormat.h"
#include <stdio.h>
#include <fcntl.h>
#include <unistd.h>
#include <iostream>

using namespace std;

Fort31File::Fort31File(String _name)
{
	name = _name;
cout<<"remco name "<<name<<endl;
	format = Fort31RecordFormatUndefined;
	
INT	fd = open(name, O_RDONLY);
int	l;
	if ( fd==-1 )
	{
		cout << "no file " << name << "." << endl;
		exit(1);
	}
	if ( ::read(fd, &l, sizeof(int))==-1 )
	{
		cout << "file " << name << " is empty." << endl;
		cout << "but it might be a dummy file for RI-Calculation." << endl;
		//exit(1);
	}

	if ( l>10000 )
	{
		cout << "endian notation of file " << name << " is not compatible to this machine." << endl;
		cout << "please use \"f31endian\" to convert." << endl;
		exit(1);
	}
	if ( l==3634 )
	{
		format = Fort31RecordFormatNew;
		return;
	}
	if ( l==23434 )
	{
		format = Fort31RecordFormatOld;
		return;
	}
	lseek(fd, l + sizeof(int), SEEK_CUR);
	::read(fd, &l, sizeof(int));
	if ( l==230 || l==366 )
	{
		format = Fort31RecordFormatTRADPT;
		return;
	}
}
