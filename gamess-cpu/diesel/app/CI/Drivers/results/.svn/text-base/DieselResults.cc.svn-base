//***********************************************************************
//
//	Name:			DieselResults.cc
//
//	Description:	stores input for diesel driver 
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			09.11.1998
//
//
//
//
//
//***********************************************************************

#include "DieselResults.h"

#include "../../../../lib/Container/GrepAwk.h"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>
#include <string>

#include <dirent.h>

#include <math.h>

using namespace std;

const double DieselResults::imp = 0.01;

static INT feq(Number<double> a, Number<double> b)
{
	if ( !a.isNumber() || !b.isNumber() )
		return 0;
	return fabs(a.value()-b.value())<fabs(a.value())/1000.0;
}

DieselResults::DieselResults(ostream &s, Number<double> Thresh, INT Root, 
	INT noWave, INT noHeaders, String MRPTThresh) :
	Thresh(Thresh), Root(Root), 
	noWave(noWave), noHeaders(noHeaders)
{
ifstream	sel("sel.out.all");

	if ( !sel )
	{
		cerr << "no file \"sel.out.all\"" << endl;
		exit(1);
	}
	
	{
	GrepAwk	ga(sel);
		ga.grep("R e s u l t s");
		ga.grep("threshold");
		ga.grep("======");
		ga++;
		nRoots = ga.getNumberOfWords()-3;
	Pix	pos = ga.getIndex();
		nThresholds = 0;
		while ( ga.getNumberOfWords() )
		{
			ga++;
			nThresholds++;
		}
		ga.setIndex(pos);

		threshData = new TThreshData[nThresholds];
		for ( INT i=0 ; i<nThresholds ; i++ )
		{
			threshData[i].threshold = atof(ga.getWord(1));
			threshData[i].nConfs = atoi(ga.getWord(2));
			threshData[i].nSAFs = atoi(ga.getWord(3));
			threshData[i].root = new TThreshRootData[nRoots];
			for ( INT j=0 ; j<nRoots ; j++ )
			{
				threshData[i].root[j].rootNr_sel2ref = -1;
				threshData[i].root[j].rootNr_ref2sel = -1;
				
				threshData[i].root[j].EPTsum_EN = atof(ga.getWord(4+j));
			}
			threshData[i].overlap = Matrix<double>(nRoots, nRoots);
			threshData[i].overlap.setOne();
			ga++;
		}
		ga.head();
		if ( !ga.grep("selected reference configuration") )
			ga.head();
		ga.grep("eigenvalues of reference matrix:");
		ga++;
		rootData = new TRootData[nRoots];
		memset(rootData, 0, nRoots*sizeof(TRootData));
		for ( INT j=0 ; j<nRoots ; j++ )
		{
			rootData[j].RefE = atof(ga.getWord(1));
			ga++;
		}
	}
		
	for ( INT i=0 ; i<nThresholds ; i++ )
	{
	double	t;
	char	buf[1000];

	DIR	*dir = opendir(".");
	struct dirent	*entry;

		while ( (entry=readdir(dir)) )
			if ( sscanf(entry->d_name, "diag.out.%lf%s", &t, buf) == 1 &&
				feq(t, threshData[i].threshold/1000.0) )
			{
			ifstream	f(entry->d_name);

			GrepAwk	ga(f);
				if ( ga.grep("R e s u l t s") )
				{
					if ( ga.grep("eigenvector overlap with reference space:") )
					{
						ga += 2;
						threshData[i].overlap = Matrix<double>(nRoots, nRoots);
						for ( INT j=0 ; j<nRoots ; j++ )
						{
							for ( INT k=0 ; k<nRoots ; k++ )
								threshData[i].overlap(j, k) = atof(ga.getWord(k+1));
							ga++;
						}

						if ( ga.grep("renormalized eigenvector overlap with reference space:") )
						{
							ga += 2;
							threshData[i].overlapRenorm = Matrix<double>(nRoots, nRoots);
							threshData[i].overlapMax = Matrix<double>(nRoots, nRoots);
							threshData[i].overlapMax.setFill(0);
							for ( INT j=0 ; j<nRoots ; j++ )
							{
								for ( INT k=0 ; k<nRoots ; k++ )
									threshData[i].overlapRenorm(j, k) = atof(ga.getWord(k+1));
								ga++;
							INT	maxk = -1;
							double	max = 0;
								for ( INT k=0 ; k<nRoots ; k++ )
									if ( fabs(threshData[i].overlap(j, k))>max )
									{
										max = fabs(threshData[i].overlap(j, k));
										maxk = k;
									}
								threshData[i].overlapMax(j, maxk) = max;
							}
							for ( INT j=0 ; j<nRoots ; j++ )
							{
							double	max = 0;
							INT	maxk = -1;
								for ( INT k=0 ; k<nRoots ; k++ )
								{
									if ( fabs(threshData[i].overlapMax(k, j))>max )
									{
										max = fabs(threshData[i].overlapMax(k, j));
										maxk = k;
									}
								}
											
								if ( maxk!=-1 )
								{
									threshData[i].root[maxk].rootNr_sel2ref = j;
									threshData[i].root[j].rootNr_ref2sel = maxk;
								}
							
							}
							
							s << endl << endl << endl << "threshold " << threshData[i].threshold << " mH" << endl;
							s << "eigenvector overlap with reference space "
							 << "(column: reference root, row: MR-CI root):" << endl;

							s.precision(5);
							s.setf(ios::fixed);
							s.width(12);
							for ( INT j=0 ; j<nRoots ; j++ )
							{
							double	s2 = 0;
								for ( INT k=0 ; k<nRoots ; k++ )
								{
									s << (threshData[i].overlapMax(j, k)>0 ? '[' : ' ') <<
										setw(7) << setprecision(4) << threshData[i].overlap(j, k)
										<< (threshData[i].overlapMax(j, k)>0 ? ']' : ' ') << "  ";
									s2 += threshData[i].overlap(j, k)*threshData[i].overlap(j, k);
								}
								s << "   |   "  << setw(7) << setprecision(4) << s2;
								s << endl;
							}

							for ( INT j=0 ; j<nRoots ; j++ )
							{
								threshData[i].root[j].EPTsum_ENweighted = 0;
							INT	na = 0;
							INT	maxk = -1;
							double	max = 0;
								for ( INT k=0 ; k<nRoots ; k++ )
									if ( threshData[i].overlapMax(k, j)>0 )
									{
										na ++;
										if ( fabs(threshData[i].overlap(k, j))>max )
										{
											max = fabs(threshData[i].overlap(k, j));
											maxk = k;
										}
										threshData[i].root[maxk].overlap = max;	
									}

								if ( na>1 )
									s << "*** warning: reference root " << j+1 << " is multiply assigned "
										<< "(assigning to MR-CI root " << maxk+1 << ") ***" << endl
										<< "*** please check your references! ***" << endl;
								if ( na==0 )
									s << "*** warning: reference root " << j+1 << " is never assigned "
										<< "(probably new root in MR-CI) ***" << endl
										<< "*** please check your references! ***" << endl;
								
							}


							for ( INT j=0 ; j<nRoots ; j++ )
							{
								if ( threshData[i].root[j].rootNr_sel2ref != -1 )
								{
									for ( INT k=0 ; k<nRoots ; k++ )
											threshData[i].root[threshData[i].root[j].rootNr_sel2ref].EPTsum_ENweighted += 
												threshData[i].overlapRenorm(j, k)*
												threshData[i].overlapRenorm(j, k)*
												threshData[i].root[k].EPTsum_EN;
								}
							}
						}
					}
					else
					{
						cerr << "no \"eigenvector overlap\" in \"" << entry->d_name << "\"" << endl;
						cerr << "using diagonal" << endl;
						
						for ( INT j=0 ; j<nRoots ; j++ )
						{
							for ( INT k=0 ; k<nRoots ; k++ )
							{
								threshData[i].overlapRenorm(j, k) = (j==k);
								threshData[i].overlapMax(j, k) = (j==k);
								threshData[i].root[j].overlap = 1;	
								threshData[i].root[j].rootNr_sel2ref = j;
								threshData[i].root[j].rootNr_ref2sel = j;
							}

						}
					}

//----------------------------------------------------------

					ga.head();
					ga.grep("R e s u l t s");

					ga += 7;
					for ( INT j=0 ; j<nRoots ; j++ )
					{
//						cout << "!!!!!: " << threshData[i].root[j].rootNr_sel2ref << " " << atof(ga.getWord(4)) << endl;
						if ( threshData[i].root[j].rootNr_sel2ref!=-1 )
							threshData[i].root[threshData[i].root[j].rootNr_sel2ref].ECIvar = atof(ga.getWord(4));
						ga++;
					}
					for ( INT j=0 ; j<nRoots ; j++ )
					{
						if ( ga.grep("root #") )
						{
							ga += 3;
							if ( i==nThresholds-1 )
							{
							INT	flag = 0;
								while ( !ga.illegal() && ga.getNumberOfWords() )
								{
								String	s(ga.getLine());
									if ( ga.getNumberOfWords()>1 )
									{
									double	v;

										sscanf(s.chars(), "%lf", &v);
										flag = v>imp;

									}

									if ( flag && threshData[i].root[j].rootNr_sel2ref!=-1 )
										rootData[threshData[i].root[j].rootNr_sel2ref].waveFunction.append(s);
									ga++;
								}
							char	buf[100];
								sprintf(buf, "%.2e mH", threshData[i].threshold.value());
								waveFunctionThreshold = String(buf);
							}
						}
						else
							cerr << "no \"root #\" in \"" << entry->d_name << "\"" << endl;

						if ( ga.grep("ci^2 of references:") )
						{
							if ( threshData[i].root[j].rootNr_sel2ref!=-1 )
								threshData[i].root[threshData[i].root[j].rootNr_sel2ref].RefCISqr = atof(ga.getWord(4));
						}
						else
							cerr << "no \"ci^2 of references\" in \"" << entry->d_name << "\"" << endl;
						ga++;
					}
				}
			}
		closedir(dir);
	}


	{
	String	fn("mrpt.out");
//		if ( MRPTThresh.length()>0 )
//			fn += String(".") + MRPTThresh;
	ifstream	f(fn);
		MRMP = (f!=0);
		if ( f )
		{
		GrepAwk	ga(f);
			if ( ga.grep("MR-MP2 Results:") )
			{
				ga += 3;

				while ( ga.getNumberOfWords() )
				{
					for ( INT i=0 ; i<nThresholds ; i++ )
						if ( feq(threshData[i].threshold, atof(ga.getWord(1))) )
							for ( INT j=0 ; j<nRoots ; j++ )
								if ( threshData[i].root[j].rootNr_sel2ref!=-1 )
									threshData[i].root[threshData[i].root[j].rootNr_sel2ref].EPTsum_MRMP2 = atof(ga.getWord(4+j));
					ga++;
				}
				if ( ga.grep("MR-MP3 Results:") )
				{
					ga += 5;

					while ( ga.getNumberOfWords() )
					{
						for ( INT i=0 ; i<nThresholds ; i++ )
							if ( feq(threshData[i].threshold, atof(ga.getWord(1))) )
								for ( INT j=0 ; j<nRoots ; j++ )
									if ( threshData[i].root[j].rootNr_sel2ref!=-1 )
										threshData[i].root[threshData[i].root[j].rootNr_sel2ref].EPTsum_MRMP3 = atof(ga.getWord(4+j));
						ga++;
					}
				}
			}
			
		}
	}	

	for ( INT j=0 ; j<nRoots ; j++ )
	{
		for ( INT i=0 ; i<nThresholds ; i++ )
		{
			if ( threshData[i].root[threshData[i].root[j].rootNr_ref2sel].overlap.value()<0.5 )
				threshData[i].root[j].remark += "  *** small overlap";
			else
			{
				threshData[i].root[j].MRCIextPol[0] = 
					threshData[i].root[j].ECIvar + threshData[i].root[j].EPTsum_ENweighted/1000.0;
				threshData[i].root[j].MRCIextPol[3] = 
					threshData[i].root[j].ECIvar + threshData[i].root[j].EPTsum_MRMP2/1000.0;
				threshData[i].root[j].MRCIextPol[6] = 
					threshData[i].root[j].MRCIextPol[3] + threshData[i].root[j].EPTsum_MRMP3/1000.0;

				if ( i )
				{
				Number<double>	l;
					l = -1000.0*(threshData[i].root[j].ECIvar-threshData[i-1].root[j].ECIvar)/
						(threshData[i].root[j].EPTsum_ENweighted-threshData[i-1].root[j].EPTsum_ENweighted);
					threshData[i].root[j].MRCIextPol[1] = l;
					if ( l.value()<0.5 || l.value()>1.0 )
						threshData[i].root[j].remark += "  *** bad lambda EN";
					threshData[i].root[j].MRCIextPol[2] = threshData[i].root[j].ECIvar + 
						l*threshData[i].root[j].EPTsum_ENweighted/1000.0;

					if ( MRMP )
					{
						l = -1000.0*(threshData[i].root[j].ECIvar-threshData[i-1].root[j].ECIvar)/
							(threshData[i].root[j].EPTsum_MRMP2-threshData[i-1].root[j].EPTsum_MRMP2);
						threshData[i].root[j].MRCIextPol[4] = l;
						if ( l.value()<0.5 || l.value()>1.0 )
							threshData[i].root[j].remark += "  *** bad lambda MP2";

						threshData[i].root[j].MRCIextPol[5] = threshData[i].root[j].ECIvar + 
							l*threshData[i].root[j].EPTsum_MRMP2/1000.0;

						l = -1000.0*(threshData[i].root[j].ECIvar-threshData[i-1].root[j].ECIvar)/
							(threshData[i].root[j].EPTsum_MRMP2-threshData[i-1].root[j].EPTsum_MRMP2 +
							threshData[i].root[j].EPTsum_MRMP3-threshData[i-1].root[j].EPTsum_MRMP3);
						threshData[i].root[j].MRCIextPol[7] = l;
						if ( l.value()<0.5 || l.value()>1.0 )
							threshData[i].root[j].remark += "  *** bad lambda MP3";

						threshData[i].root[j].MRCIextPol[8] = threshData[i].root[j].ECIvar + 
							l*(threshData[i].root[j].EPTsum_MRMP2+threshData[i].root[j].EPTsum_MRMP3)/1000.0;
					}

				}

				for ( INT k=0 ; k<(MRMP?3:1) ; k++ )
				{
					threshData[i].root[j].fullCIextPol[k*4+0] = threshData[i].root[j].MRCIextPol[k*3+0];
					threshData[i].root[j].fullCIextPol[k*4+1] = threshData[i].root[j].MRCIextPol[k*3+2];
					threshData[i].root[j].fullCIextPol[k*4+2] = threshData[i].root[j].MRCIextPol[k*3+0];
					threshData[i].root[j].fullCIextPol[k*4+3] = threshData[i].root[j].MRCIextPol[k*3+2];

	//			s << "j=" << j << ", i=" << i << ", k=" << k << endl;
	//			INT l = threshData[i].root[j].rootNr_sel2ref;
				INT l = j;
	//				for ( INT l=0 ; l<nRoots ; l++ )
					if ( l!=-1 )
					{
					Number<double>	cisqr = threshData[i].root[j].RefCISqr;
					Number<double>	RefE = rootData[l].RefE;
					Number<double>	w = 1; //threshData[i].overlap(j, l)*threshData[i].overlap(j, l);
	//					s << w << endl;
	//					cout << "j=" << j << ", l=" << l << endl;
	//					cout << "RefE=" << RefE << endl;


						threshData[i].root[j].fullCIextPol[k*4+0] += 
							-w*(1.0-cisqr)*(RefE-threshData[i].root[j].MRCIextPol[k*3+0]);
						threshData[i].root[j].fullCIextPol[k*4+1] += 
							-w*(1.0-cisqr)*(RefE-threshData[i].root[j].MRCIextPol[k*3+2]);
						threshData[i].root[j].fullCIextPol[k*4+2] += 
							-w*(1.0-cisqr)*(RefE-threshData[i].root[j].MRCIextPol[k*3+0])
								/cisqr;
						threshData[i].root[j].fullCIextPol[k*4+3] += 
							-w*(1.0-cisqr)*(RefE-threshData[i].root[j].MRCIextPol[k*3+2])
								/cisqr;
					}
				}
				if ( threshData[i].root[j].RefCISqr.value()<0.7 || threshData[i].root[j].RefCISqr.value()>1.0 )
					threshData[i].root[j].remark += "  *** bad ci^2";
			}
		}
	}
}


DieselResults::~DieselResults()
{
	for ( INT i=0 ; i<nThresholds ; i++ )
		delete[] threshData[i].root;
	delete[] threshData;
	delete[] rootData;

}



ostream & operator << (ostream &s, const DieselResults &dr)
{
	s.precision(5);
	s.setf(ios::fixed);
	s << endl << endl;
	for ( INT j=0 ; j<dr.nRoots ; j++ )
	{
		if ( dr.Root!=-1 && j!=dr.Root )
			continue;
		if ( !dr.noWave )
		{
			s << "++++++++" << endl;
			s << "root #" << j+1 << ":" << endl;
			s << "reference energy: " << dr.rootData[j].RefE << endl;
			s << "character (at threshold " << dr.waveFunctionThreshold << ", ci^2>" << dr.imp << "): " << endl;
		Pix	pix = dr.rootData[j].waveFunction.first();
			while ( pix )
			{
				s << dr.rootData[j].waveFunction(pix) << endl;
				dr.rootData[j].waveFunction.next(pix);
			}
		}
		const INT w = 12;
		if ( !dr.noHeaders )
		{
			s << endl;
			s 
				<< setw(w) << "threshold"
				<< setw(w) << "sel. confs"
				<< setw(w) << "sel. CSFs"
				<< setw(w) << "CI energy"
				<< setw(w) << "CI root"
				<< setw(w) << "ref ci^2"
				<< setw(w) << "overlap"
				<< setw(w) << "PT(EN)"
				<< setw(w) << "PT(EN)";
			if ( dr.MRMP )
				s << setw(w) << "PT(MRMP2)"
					<< setw(w) << "PT(MRMP3)";
			s << endl;
			s 
				<< setw(w) << "/mH"
				<< setw(w) << ""
				<< setw(w) << ""
				<< setw(w) << "/H"
				<< setw(w) << ""
				<< setw(w) << ""
				<< setw(w) << "max"
				<< setw(w) << "/mH"
				<< setw(w) << "weighted/mH";
			if ( dr.MRMP )
				s << setw(w) << "/mH"
					<< setw(w) << "/mH";
			s << endl;
			s << "--------------------------------------------------------------------------------------------------------------";
			if ( dr.MRMP )
				s << "------------------------";
			s << endl;
		}
		for ( INT i=0 ; i<dr.nThresholds ; i++ )
		{
			if ( dr.Thresh.isNumber() && !feq(dr.Thresh, dr.threshData[i].threshold) )
				continue;
		iostream::fmtflags	f = s.flags();
			s.setf(iostream::scientific, iostream::floatfield);
			s << setw(w) << setprecision(3) << dr.threshData[i].threshold;
			s.flags(f);
			s << setw(w) << setprecision(5) << dr.threshData[i].nConfs
			<< setw(w) << dr.threshData[i].nSAFs
			<< setw(w) << dr.threshData[i].root[j].ECIvar
			<< setw(w) << dr.threshData[i].root[j].rootNr_ref2sel+1
			<< setw(w) << dr.threshData[i].root[j].RefCISqr
			<< setw(w) << dr.threshData[i].root[dr.threshData[i].root[j].rootNr_ref2sel].overlap
			<< setw(w) << setprecision(3) << dr.threshData[i].root[j].EPTsum_EN
			<< setw(w) << dr.threshData[i].root[j].EPTsum_ENweighted;
			s.precision(5);
			if ( dr.MRMP )
				s << setw(w) << dr.threshData[i].root[j].EPTsum_MRMP2
					<< setw(w) << dr.threshData[i].root[j].EPTsum_MRMP3;
			s << endl;
		}
		if ( !dr.noHeaders )
		{
			s << endl;
			s << endl;
			s << "extrapolation to full MRCI (all Energies in H):" << endl;
			s << endl;
			s 
				<< setw(w) << "threshold"
				<< setw(3*w) << "_______________ EN ________________";
			if ( dr.MRMP )
				s << setw(3*w) << "______________ MRMP2 ______________"
					<< setw(3*w) << "______________ MRMP3 ______________";
			s << endl;
			s 
				<< setw(w) << "/mH"
				<< setw(w) << "E(l=1)"
				<< setw(w) << "l"
				<< setw(w) << "E(l)";
			if ( dr.MRMP )
				s << setw(w) << "E(l=1)"
					<< setw(w) << "l"
					<< setw(w) << "E(l)"
					<< setw(w) << "E(l=1)"
					<< setw(w) << "l"
					<< setw(w) << "E(l)";
			s << endl;
			s << "------------------------------------------------";
			if ( dr.MRMP )
				s << "------------------------------------------------------------------------";
			s << endl;
		}

		for ( INT i=0 ; i<dr.nThresholds ; i++ )
		{
			if ( dr.Thresh.isNumber() && !feq(dr.Thresh, dr.threshData[i].threshold) )
				continue;
		iostream::fmtflags	f = s.flags();
			s.setf(iostream::scientific, iostream::floatfield);
			s << setw(w) << setprecision(3) << dr.threshData[i].threshold;
			s.flags(f);
			s.precision(5);
			for ( INT k=0 ; k<(dr.MRMP ? 9 : 3) ; k++ )
				s << setw(w) << dr.threshData[i].root[j].MRCIextPol[k];
			s << "    " << dr.threshData[i].root[j].remark << endl;
		}
		if ( !dr.noHeaders )
		{
			s << endl;
			s << endl;
			s << "extrapolation to full CI (all Energies in H):" << endl;
			s << endl;
			s << setw(w) << "threshold"
				<< setw(4*w) << "____________________ EN _____________________";
			if ( dr.MRMP )
				s << setw(4*w) << "___________________ MRMP2 ___________________"
				<< setw(4*w) << "___________________ MRMP3 ___________________";
			s << endl;




			s	<< setw(w) << "/mH"
				<< setw(2*w) << "____ Davidson 1 _____"
				<< setw(2*w) << "____ Davidson 2 _____";
			if ( dr.MRMP )
				s << setw(2*w) << "____ Davidson 1 _____"
					<< setw(2*w) << "____ Davidson 2 _____"
					<< setw(2*w) << "____ Davidson 1 _____"
					<< setw(2*w) << "____ Davidson 2 _____";
			s << endl;

			s 
				<< setw(w) << ""
				<< setw(w) << "E(l=1)"
				<< setw(w) << "E(l)"
				<< setw(w) << "E(l=1)"
				<< setw(w) << "E(l)";
			if ( dr.MRMP )
				s << setw(w) << "E(l=1)"
					<< setw(w) << "E(l)"
					<< setw(w) << "E(l=1)"
					<< setw(w) << "E(l)"
					<< setw(w) << "E(l=1)"
					<< setw(w) << "E(l)"
					<< setw(w) << "E(l=1)"
					<< setw(w) << "E(l)";
			s << endl;
			s << "------------------------------------------------------------";
			if ( dr.MRMP )
				s << "-----------------------------------------------------------------------------------------------";
			s << endl;
		}
		for ( INT i=0 ; i<dr.nThresholds ; i++ )
		{
			if ( dr.Thresh.isNumber() && !feq(dr.Thresh, dr.threshData[i].threshold) )
				continue;
		iostream::fmtflags	f = s.flags();
			s.setf(iostream::scientific, iostream::floatfield);
			s << setw(w) << setprecision(3) << dr.threshData[i].threshold;
			s.flags(f);
			s.precision(5);
			for ( INT k=0 ; k<(dr.MRMP ? 12 : 4) ; k++ )
				s << setw(w) << dr.threshData[i].root[j].fullCIextPol[k];
			s << "    " << dr.threshData[i].root[j].remark << endl;
		}
		if ( !dr.noWave )
		{
			s << "++++++++" << endl;
			s << endl;
			s << endl;
			s << endl;
		}
	}		
	return s;
}




#include "../../../../lib/Container/SLList.cc"

template class SLList<String>;
