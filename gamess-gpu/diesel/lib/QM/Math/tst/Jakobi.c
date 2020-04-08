#include <math.h>

#include <stdlib.h>



static int	comp(const void *p1, const void * p2)
{
	if ( *((double *) p1) > *((double *) p2))
		return 1;
	else
	if ( *((double *) p1) < *((double *) p2))
		return -1;
	else
		return 0;
}


void Jakobi(INT dim, double *hamilton, double *eigenvectors)
{
 double *h = hamilton;
 double *vv = eigenvectors;
 double *h1;
 double *hhp, *hhq, *php, *phq, *vvp, *vvq;
  INT a, b, c, p_, q, pp, qq, hp, hq;
  double m, br, max, si, co, ha, hb, hs;

	memset(vv, 0, sizeof(double) * dim * dim);
	h1 = vv;
	for ( a = 0 ; a < dim ; a++)
	{	*h1++ = 1.0;
		h1 += dim;
	}

_Lloop:
  max = 0.0;
  b = 0;
  for (a = 1; a <= dim; a++) {
	for (c = a + 1; c <= dim; c++) {
	  if (fabs(h[b + c - 1]) > fabs(max)) {
	max = h[b + c - 1];
	p_ = a;
	q = c;
	  }
	}
	b += dim;
  }
	if (fabs(max) <= 1.0e-10)
	{
/*	normalize eigenvectors	*/
/*	for ( a=0 ; a<dim ; a++ )
	{	
		m = 0;
		for ( b=0 ; b<dim ; b++ )
			m += vv[b*dim + a]*vv[b*dim + a];
		printf("%f\n", m);
		m = 1/sqrt(m);
		for ( b=0 ; b<dim ; b++ )
			vv[b*dim + a] *= m;
	}
*/

/*	sort by eigenvalues	*/
		h1 = (double *) malloc(dim*(dim+1)*sizeof(double));
		hhp = h1;
		for ( a=0 ; a<dim ; a++ )
		{	*hhp++ = h[a*dim + a];
			for ( b=0 ; b<dim ; b++ )
				*hhp++ = vv[b*dim + a];
		}

		qsort(h1, dim, (dim+1)*sizeof(double), comp);

		hhp = h1;
		for ( a=0 ; a<dim ; a++ )
		{	h[a*dim + a] = *hhp++;
			for ( b=0 ; b<dim ; b++ )
				vv[b*dim + a] = *hhp++;
		}

		free(h1);
		return;
	}
  hp = (p_ - 1) * dim;
  hq = (q - 1) * dim;
  pp = hp + p_;
  qq = hq + q;
  m = 0.5 * (h[pp - 1] - h[qq - 1]);
  br = sqrt(m * m + max * max);
  si = sqrt(fabs(1.0 - m / br)) * 7.071067811865e-1;
  if (max > 0)
	si = -si;
  co = sqrt(1.0 - si * si);
  hs = m + h[qq - 1];
  h[pp - 1] = hs + br;
  h[qq - 1] = hs - br;
  h[hp + q - 1] = 0.0;
  h[hq + p_ - 1] = 0.0;
  for (a = 1; a <= dim; a++) {
	hp++;
	hq++;
	if (hp != pp && hq != qq) {
	  h[hp - 1] = h[p_ - 1] * co - h[q - 1] * si;
	  h[hq - 1] = h[p_ - 1] * si + h[q - 1] * co;
	  h[p_ - 1] = h[hp - 1];
	  h[q - 1] = h[hq - 1];
    }
	ha = vv[p_ - 1];
    hb = vv[q - 1];
    vv[p_ - 1] = ha * co - hb * si;
	vv[q - 1] = ha * si + hb * co;
    p_ += dim;
    q += dim;
  }
  goto _Lloop;



/*_Lloop:
  printf("!");
  max = 0.0;

	h1 = h+1;
	for (a = 0; a < dim-1; a++)
	{	for ( c = a+1 ; c < dim ; c++)
		{	if (fabs(*h1) > fabs(max))
			{	max = *h1;
				p_ = a;
				q = c;
			}
			h1++;
		}
		h1 += 2+a;
	}
  if (fabs(max) <= 1.0e-5)
	return;
  hp = p_ * dim;
  hq = q * dim;
  pp = hp + p_;
  qq = hq + q;
  m = 0.5 * (h[pp] - h[qq]);
  br = sqrt(m * m + max * max);
  si = sqrt(fabs(1.0 - m / br)) * 7.071067811865e-1;
  if (max > 0)
	si = -si;
  co = sqrt(1.0 - si * si);
  hs = m + h[qq];
  h[pp] = hs + br;
  h[qq] = hs - br;
  h[hp + q] = 0.0;
  h[hq + p_] = 0.0;
  hhp = h + hp;
  hhq = h + hq;
  php = h + p_;
  phq = h + q;
  vvp = vv + p_;
  vvq = vv + q;
  for (a = 0; a < dim; a++) {
	if (hp != pp && hq != qq) {
	  *hhp = *php * co - *phq * si;
	  *hhq = *php * si + *phq * co;
	  *php = *hhp;
	  *phq = *hhq;
	}
	hp++;
	hq++;
	ha = *vvp;
	hb = *vvq;
	*vvp = ha * co - hb * si;
	*vvq = ha * si + hb * co;
	hhp++;
	hhq++;
	vvp += dim;
	vvq += dim;
	php += dim;
	phq += dim;
  }
  goto _Lloop;
*/
}
