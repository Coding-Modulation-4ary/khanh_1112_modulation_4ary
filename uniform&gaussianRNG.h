#include <math.h>

#define	IA		16807			
#define	IM		2147483647		
#define	AM		(1.0/IM)		
#define	IQ		127773			
#define	IR		2836			
#define	NTAB	32				
#define	NDIV	(1+(IM-1)/NTAB)	
#define	EPS		1.2e-7			
#define	RNMX	(1.0-EPS)		

float rng(long *idum);			
float gaussian(float m, float s, long *idum);


float rng(long *idum)
{
	int j;
	long k;
	static long iy = 0;
	static long iv[NTAB];
	float temp;

	if ( *idum <= 0 || !iy )
	{
		if (-(*idum) < 1 ) *idum = 1;
		else *idum = -(*idum);

		for (j = NTAB+7 ; j >= 0 ; j--)
		{
			k = (*idum)/IQ;
			*idum = IA*(*idum-k*IQ)-IR*k;

			if (*idum < 0)	*idum += IM;

			if (j < NTAB)	iv[j] = *idum;
		}

		iy = iv[0];
	}

	k = (*idum)/IQ;
	*idum = IA*(*idum-k*IQ)-IR*k;

	if (*idum < 0)	*idum += IM;

	j = iy/NDIV;
	iy = iv[j];
	iv[j] = *idum;

	if ((temp=(float)AM*iy) > RNMX)
		return (float)RNMX;
	else
  		return temp;
}

/* returns a Gaussian random number with a normal distribution,
   mean = m, standard deviation = s  */
float gaussian(float m, float s, long *idum)
{
	float rng(long *idum);
	float rsq=2, v1, v2, fac, gset=0, iset;

	if (*idum < 0)
	{
		iset = 0; 
		
		return (m + s * gset);
	}
	else 
	{
		while (rsq > 1)
		{
			v1 = 2 * rng(idum) - 1;
			v2 = 2 * rng(idum) - 1;
			rsq = v1 * v1 + v2 * v2;
		}
		fac = (float)sqrt((-2*log(rsq)) / rsq);
		iset = 1;
		gset = v2 * fac;
    
		return ( m + s * (v1 * fac));
	}
}