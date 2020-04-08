//***********************************************************************
//
//	Name:			FortranLinkage.h
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




#ifdef __UNDERBAR
#define FORTRAN_LINKAGE(x) x ## _	/*	Linux, ...	*/
#else
#define FORTRAN_LINKAGE(x) x		/*	RS6000/AIX, PA/HP-UX, ... */
#endif

