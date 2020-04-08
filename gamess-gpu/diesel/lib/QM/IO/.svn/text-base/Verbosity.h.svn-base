//***********************************************************************
//
//	Name:			Verbosity.h
//
//	Description:	levels of verbosity for output
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			24.02.1998
//
//
//
//
//
//***********************************************************************

#ifndef __VERBOSITY_H
#define __VERBOSITY_H


#include "../../../config.h"

#include <iostream>
using std::ostream;
using std::istream;

//FD class	istream;
//FD class	ostream;


class Verbosity {
public:

enum	Section {
			Input,						//     1
			Integrals,					//     2
			MOs,						//     4
			RefGuess,					//     8
			SGA, 						//    16
			RefMat,						//    32
			RefMatEigenValues,			//    64
			RefMatEigenVectors,			//   128
			IterationBlocks,			//   256
			WaveFunction,				//   512
			SelectionPerRoot,			//  1024
			CacheStatistics,			//  2048
			DegenGuess,					//  4096
			DiagHist					//  8192
		};

enum	Level { Quiet, Default, Verbose };

	Verbosity(Level level = Default);
	Verbosity(INT flags);
	Verbosity(istream &);


	INT	isActive(Section) const;
	void	setActive(Section);

	friend ostream & operator << (ostream &, const Verbosity &);

private:
static char	*SectionName[];
INT	flags;
static const INT nFlags = 14;
};


extern Verbosity verbosity;


inline
INT	Verbosity::isActive(Section section) const
{	return (flags & (1<<section)) ? 1 : 0;	}

inline
void	Verbosity::setActive(Section section)
{	flags |= 1<<section;	}


#endif
