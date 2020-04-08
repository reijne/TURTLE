using namespace std;

#include "Permutator.h"
#include "Spin.h"

#include <iostream>
#include <sstream>


extern "C" {
	void Jakobi(INT dim, double *hamilton, double *eigenvectors);

	void dsytrd_(
		char	*UPLO,
		INT		*N, 
		double	**A,
		INT		*LDA,
		double	**D,
		double	**E,
		double	**TAU,
		double	**WORK,
		INT		*LWORK,
		INT		*INFO
		);

	void dsteqr_(
		char	*COMPZ,
		INT		*N, 
		double	**D,
		double	**E,
		double	**Z,
		INT		*LDZ,
		double	**WORK,
		INT		*INFO
		);
}

/*
void Diag(INT dim, double *hamilton, double *eigenvectors)
{
char	UPLO = 'U';
char	COMPZ[2] = "V";
INT		LWORK = dim*32;
INT		INFO;
double	*A = new double[(dim+1)*dim];
double	*D = new double[dim];
double	*E = new double[dim-1];
double	*TAU = new double[dim-1];
double	*WORK = new double[LWORK];


double	*p = A;
	for ( INT i=0 ; i<dim ; i++ )
		for ( INT j=i ; j<dim ; j++ )
			*p++ = hamilton[i*dim + j];

INT	n = dim;
	dsytrd_(&UPLO,
		&dim, 
		&A,
		&n,
		&D,
		&E,
		&TAU,
		&WORK,
		&LWORK,
		&INFO
	);

	cout << "INFO=" << INFO << endl;

	n = dim;
	dsteqr_(COMPZ,
		&dim, 
		&D,
		&E,
		&eigenvectors,
		&n,
		&WORK,
		&INFO
	);

	cout << "INFO=" << INFO << endl;


	delete WORK;
	delete TAU;
	delete E;
	delete D;
	delete A;
}
*/

INT	main()
{

Spin	spin(1), spinSS;

INT	n;
INT	m;

	cout << "n= ";
	cin >> n;
INT	n2 = (1 << n);

double	pSS[n2*n2];
double	pv[n2];
double	pV[n2*n2];


	memset(pSS, 0, n2*n2*sizeof(double));

	for ( INT i=0 ; i<n2 ; i++ )
	{
		m = 0;
		for ( INT j=0 ; j<n ; j++ )
			m |= (( i & (1 << j) ) ? 2 : 1) << (2 * j);
		spin.setSpin(m);
		spinSS = spin.SS().simplify();

	INT k;
		for ( INT j=0 ; j<spinSS.getNumberOfTerms() ; j++ )
			if ( (k=spinSS.getSpinCompact(j))>=i )
			{	pSS[k*n2 + i] = real(spinSS.getC(j));
				pSS[i*n2 + k] = real(spinSS.getC(j));
			}
	}

	cout << "Gen fertig!\n";

double	pA[n2*n2];
	memcpy(pA, pSS, n2*n2*sizeof(double));

	Jakobi(n2, pA, pV);
//	Diag(n2, pA, pV);

	for ( INT i=0 ; i<n2 ; i++ )
	{	pv[i] = pA[i*n2 + i];
//		for ( INT j=0 ; j<n2 ; j++ )
//			pV[i*n2+j] = V(j, i);
	}
	for ( INT i=0 ; i<n2 ; i++ )
	{	for ( INT j=0 ; j<n2 ; j++ )
			cout << pSS[i*n2 + j] << " ";
		cout << endl;
	}
	for ( INT i=0 ; i<n2 ; i++ )
	{	
		cout << pv[i];
		cout << endl;
	}
	for ( INT i=0 ; i<n2 ; i++ )
	{	for ( INT j=0 ; j<n2 ; j++ )
			cout << pV[i*n2 + j] << " ";
		cout << endl;
	}
	

INT	NumberOfEigenFunctions = n2;
Spin	**EigenFunction = new Spin * [NumberOfEigenFunctions];
	for ( INT i=0 ; i<NumberOfEigenFunctions ; i++ )
	{
	INT	nn0 = 0;
		for ( INT j=0 ; j<n2 ; j++ )
			if ( fabs(pV[j*n2 + i])>1e-10 )
				nn0++;
				
		EigenFunction[i] = new Spin(nn0);
		nn0 = 0;
		for ( INT j=0 ; j<n2 ; j++ )
		{//	cout << "i, j= " << i << ", " << j << ", V(j,i) = " << pV[i*n2 + j] << endl;
			if ( fabs(pV[j*n2 + i])>1e-10 )
			{	m = 0;
				for ( INT k=0 ; k<n ; k++ )
					m |= (( j & (1 << k) ) ? 2 : 1) << (2 * k);
//		cout << *EigenFunction[i] << endl;
				EigenFunction[i]->setSpin(nn0, m);
				EigenFunction[i]->setC(nn0++, pV[j*n2 + i]);
			}
		}
		cout << *EigenFunction[i] << endl;
		cout << (*EigenFunction[i]).SS().simplify() << endl;
		cout << "=================" << endl;
	}
	
	NumberOfEigenFunctions = 1 << n;

INT	openShells = n;

Permutator	p(openShells);
Permutator	p1(openShells), p2(openShells);

//	p1.setTransposition(2, 0);
//	p2.setTransposition(3, 0);

Spin	spinPspin;


	for ( INT j=0 ; j<openShells-1 ; j++ )
	{	for ( INT i=j+1 ; i<openShells ; i++ )
		{	cout << i << " " << j << endl;
			p.setTransposition(i, j);
//			p = p2*p1;
			cout << p << endl;
			for ( INT ii=0 ; ii<NumberOfEigenFunctions ; ii++ )
			{	for ( INT jj=0 ; jj<NumberOfEigenFunctions ; jj++ )
				{	spin = *EigenFunction[ii];
					spin.permute(p);
					spinPspin = *EigenFunction[jj] * spin;
					spinPspin.simplify();
				//	cout << spin << endl;
				//	cout << spinPspin << endl;
				//	cout << "<..|P|..> = " << spinPspin.Integrate() 
				//		<< endl << endl << endl;
					cout << p.getParity()*real(spinPspin.Integrate()) << " ";
				}
				cout << endl;
			}
			cout << endl;
		}
	}

/*
Permutator	p(4);

Spin	spinPspin;

Complex	Sum = 0;

	p.resetIndex();
	for ( INT i=0 ; i<24 ; i++ )
	{
		spin = *EigenFunction[0];
		spin.permute(p);
		cout << p << endl;
		spinPspin = *EigenFunction[0] * spin;
		spinPspin.simplify();
		cout << spin << endl;
		cout << spinPspin << endl;
		cout << "<..|P|..> = " << spinPspin.Integrate() << endl << endl << endl;
		Sum += p.getParity()*spinPspin.Integrate();
		p.nextIndex();
	}
	
	cout << "Sum = " << Sum << endl;
*/
//		cout << V(i, 0) << endl;


/*
Spin	s[4], Ps[4];
Spin	e1, e2, e3;

	cout << "============================================" << endl;

	s[0] = alpha(1)*alpha(2);
	s[1] = beta(1)*beta(2);
	s[2] = 1/sqrt(2) * (alpha(1)*beta(2) + beta(1)*alpha(2));
	s[3] = 1/sqrt(2) * (alpha(1)*beta(2) - beta(1)*alpha(2));

	Ps[0] = alpha(2)*alpha(1);
	Ps[1] = beta(2)*beta(1);
	Ps[2] = 1/sqrt(2) * (alpha(2)*beta(1) + beta(2)*alpha(1));
	Ps[3] = 1/sqrt(2) * (alpha(2)*beta(1) - beta(2)*alpha(1));


	for ( INT i=0 ; i<4; i++ )
	{	for ( INT j=0 ; j<4 ; j++ )
			cout << (s[i]*Ps[j]).Integrate() << "\t";
		cout << endl;
	}

*/
	cout << endl;
	NumberOfEigenFunctions = n2;
	
	for ( INT i=0 ; i<NumberOfEigenFunctions ; i++ )
	{
	Complex	SS, Sz;
	INT	iSS, iSz;
		iSS = EigenFunction[i]->isMultipleOf(EigenFunction[i]->SS().simplify(), SS);
		iSz = EigenFunction[i]->isMultipleOf(EigenFunction[i]->Sz().simplify(), Sz);
		if ( iSS )
			cout << real(SS);
		else
			cout << "X";
		cout << " ";
		if ( iSz )
			cout << real(Sz);
		else
			cout << "X";

		cout << endl;

	}

	exit(0);
	cout << endl;

Spin	sp;

	for ( INT i=0 ; i<NumberOfEigenFunctions ; i++ )
	{	for ( INT j=0 ; j<NumberOfEigenFunctions ; j++ )
		{
			sp = *EigenFunction[i] + *EigenFunction[j];
//			cout << sp.Sz().simplify();
	Complex	SS, Sz;
	INT	iSS, iSz;
	
			iSS = sp.isMultipleOf(sp.SS().simplify(), SS);
			iSz = sp.isMultipleOf(sp.Sz().simplify(), Sz);

			cout << "[";
			if ( iSS )
				cout << real(SS);
			else
				cout << "X";
			cout << " ";
			if ( iSz )
				cout << real(Sz);
			else
				cout << "X";
			cout << "]\t";
		}
		cout << endl;
	}


/*

Spin	s = alpha(1)*beta(2)*alpha(3);
Spin	s2;
Permutator	p(3);

//	p[0] = 0;
//	p[1] = 1;
//	p[2] = 2;
	
	cout << s << endl;
	s2 = s;
	s2.symmetrize();
	s2.simplify();
	cout << s2 << endl;
	
	cout << s << endl;
	s2 = s;
	s2.antisymmetrize();
	s2.simplify();
	cout << s2 << endl;
	
//	cout << p << endl;
*/
/*	p.resetIndex();
	for ( INT i=0 ; i<6 ; i++ ) 
	{	cout << p << " " << p.getIndex() << "  ";
		cout << p.nextIndex() << endl;
	}	
*/
/*	p.setIndex(0);
	for ( INT i=0 ; i<6 ; i++ ) 
	{	cout << p << endl << p.getParity() << endl << endl;
		p++;
	}	
*/	
//	s2 = s.permute(p);
//	cout << s2 << endl; 
}
