//***********************************************************************
//
//	Name:			FortranFileIO.cc
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			27. Jan 1998
//
//***********************************************************************



#include "FortranFileIO.h"

#include <stdio.h>
#include <fcntl.h>
#include <unistd.h>
#include <string>

#include <iostream>

FortranFileIO::FortranFileIO(INT n, INT _bigEndian)
{
char	s[100];
	sprintf(s, "fort.%d", n);
printf("remco in FortranFileIO %s\n",s);
	bigEndian = _bigEndian;
	recPos = 0;
	init(s);
}

FortranFileIO::FortranFileIO(const char *s, INT _bigEndian)
{
	bigEndian = _bigEndian;
	recPos = 0;
	init(s);
}

FortranFileIO::~FortranFileIO()
{
	close(fd);
}


void	FortranFileIO::init(const char *s)
{
	strcpy(filename, s);
	fd = open(filename, O_RDONLY);
	if ( fd==-1 )
		fd = creat(filename, 0x180 | 0x20 | 0x4 );
}

void	FortranFileIO::swap(int & i)
{
char	h, *p = (char *) &i;
	h = ((char *) p)[0];
	((char *) p)[0] = ((char *) p)[3];
	((char *) p)[3] = h;

	h = ((char *) p)[1];
	((char *) p)[1] = ((char *) p)[2];
	((char *) p)[2] = h;
}

void	FortranFileIO::rewind()
{
	lseek(fd, 0, SEEK_SET);
}

void	FortranFileIO::end()
{
	lseek(fd, 0, SEEK_END);
}

void	FortranFileIO::nextRecord()
{
int	l;
	if ( ::read(fd, &l, sizeof(int))==-1 )
		l = -1;
	if ( bigEndian )
		swap(l);
	if ( l<0 )
		return;
	lseek(fd, l + sizeof(int), SEEK_CUR);
}

void    FortranFileIO::previousRecord()
{
int     l;
INT     pos = lseek(fd, 0, SEEK_CUR);
        lseek(fd, pos-sizeof(int), SEEK_SET);
        if ( ::read(fd, &l, sizeof(int))==-1 )
                l = -1;
        if ( bigEndian )
                swap(l);
        if ( l<0 )
                return;
        pos = lseek(fd, 0, SEEK_CUR);
        lseek(fd, pos-(l + 2*sizeof(int)), SEEK_SET);
}

/*
void	FortranFileIO::previousRecord()
{
int	l;
	lseek(fd, -sizeof(int), SEEK_CUR);
	if ( ::read(fd, &l, sizeof(int))==-1 )
		l = -1;
	if ( bigEndian )
		swap(l);
	if ( l<0 )
		return;
	lseek(fd, -(l + 2*sizeof(int)), SEEK_CUR);
}
*/

INT	FortranFileIO::read(void *p, INT max)
{
int	l;
	recPos = 0;
	if ( ::read(fd, &l, sizeof(int))==-1 )
		l = -1;
	if ( bigEndian )
		swap(l);
	if ( l<0 )
		return 0;
	if ( max>l )
		max = l;
INT	count = ::read(fd, p, max);
	lseek(fd, sizeof(int)+(l-max), SEEK_CUR);
	return count;	
}

INT	FortranFileIO::read(void *p)
{
int	l;
	recPos = 0;
	if ( ::read(fd, &l, sizeof(int))==-1 )
		l = -1;
	if ( bigEndian )
		swap(l);
	if ( l<0 )
		return 0;
INT	count = ::read(fd, p, l);
	lseek(fd, sizeof(int), SEEK_CUR);
	return count;	
}

INT     FortranFileIO::readInRec(void *p, INT max)
{
int     l = 0;
INT     pos = lseek(fd, 0, SEEK_CUR);
        if ( ::read(fd, &l, sizeof(int))==-1 )
                l = -1;
        if ( bigEndian )
                swap(l);
        l -= recPos;
        if ( l<0 )
                return 0;
        if ( max>l )
                max = l;

        lseek(fd, recPos, SEEK_CUR);
INT     count = ::read(fd, p, max);
        recPos += count;
        lseek(fd, pos, SEEK_SET);
        return count;
}

/*INT	FortranFileIO::readInRec(void *p, INT max)
{
int	l;
	if ( ::read(fd, &l, sizeof(int))==-1 )
		l = -1;
	if ( bigEndian )
		swap(l);
	l -= recPos;
	if ( l<0 )
		return 0;
	if ( max>l )
		max = l;

	lseek(fd, recPos, SEEK_CUR);		
INT	count = ::read(fd, p, max);
	recPos += count;
	lseek(fd, -(sizeof(int)+recPos), SEEK_CUR);
	return count;	
}
*/


INT	FortranFileIO::write(void *p, INT n)
{
int	l = n;
	if ( bigEndian )
		swap(l);
	::write(fd, &l, sizeof(int));
INT	count = ::write(fd, p, n);
	::write(fd, &l, sizeof(int));
	return count;	
}

