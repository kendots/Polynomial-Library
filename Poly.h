#ifndef POLY_H
#define POLY_H


#include<stdio.h>
#include<stdlib.h>
#include<stdarg.h>
#include"Vector.h"


typedef struct
{
	int n;
	double *t;
} poly;

poly pnull={-5,NULL};

void PolyPrint (poly);
void PolyInit (poly *);
void PolyFree (poly *);

void PolyInit (poly * P)
{
	if (P->t)
		free(P->t);
	P->t=calloc(P->n+1,sizeof(double));
}


void PolyFree(poly * P)
{
	if (P->t)
	{
		free(P->t);
		*P=pnull;
	}
}


//input polynomial
poly Poly (int n, ...)
{
	va_list valist;
	int i;
	poly P={n,NULL};
	PolyInit(&P);
	va_start(valist, n);
	for (i=0; i<=n; i++)
		P.t[i]=(double) va_arg(valist, double);
	va_end(valist);
	return P;
}


double PolySub (poly P, double x)
{
	int i,j;
	double s=0,r,z=1;
	for (i=0; i<=P.n; i++)
	{
		r=P.t[i];
		r*=z;
		s+=r;
		z*=x;
	}
	return s;
}


/*void fill (poly * P)
{
	double func (double x) 
	{
		return PolySub(*P,x);
	}

	P->f=func;
}*/



int max (int a, int b)
{
	if (a<=b) return b;
	return a;
}


int Degree (poly * P)
{
	int i;
	for (i=P->n; i>=0; i--)
	{
		if (ab(P->t[i])>eps)
			break;
	}
	P->n=i;
	return i;
}


void PolyAdd (poly P, poly Q, double x, double y, poly * R)
{
	int i,l,k;
	double z;
	l=max(P.n, Q.n);
	k=P.n+Q.n-l;

	if (R->t==NULL)
	{
		R->n=l;
		R->t=calloc(l+1,sizeof(double));
	}

	if (P.n==l)
	{
		for (i=k+1; i<=l; i++)
			R->t[i]=x*P.t[i];
	}
	else
	{
		for (i=k+1; i<=l; i++)
			R->t[i]=y*Q.t[i];
	}
	for (i=0; i<=k; i++)
		R->t[i]=x*P.t[i]+y*Q.t[i];
	Degree(R);
}


void PolyProd (poly P, poly Q, double x, poly * R)
{
	int i,j;
	if (R->t==NULL)
	{
		R->n=P.n+Q.n;
		R->t=calloc(1+R->n,sizeof(double));
	}

	for (i=0; i<=P.n; i++)
	{
		for (j=0; j<=Q.n; j++)
			R->t[i+j]+=P.t[i]*Q.t[j];
	}

	for (i=0; i<=R->n; i++)
		R->t[i]=x*R->t[i];
}


void PolyScal (poly P, double x, poly * R)
{
	if (R->t==NULL)
	{
		R->n=P.n;
		R->t=calloc(1+P.n,sizeof(double));
	}

	int i;
	R->n=P.n;
	for (i=0; i<=P.n; i++)
		R->t[i]=x*P.t[i];
}


int PolyDiv (poly P, poly S, poly * Q, poly * R )
{
	int i;
	double x,y,z;

	Degree(&S);
	if (S.n<0)
	{
		fputs("Error",stdout);
		return -1;
	}

	if (R->t==NULL)
	{
		R->n=P.n;
		R->t=calloc(1+P.n,sizeof(double));
	}

	if (Q->t==NULL)
	{
		Q->n=P.n-S.n;
		Q->t=calloc(1+Q->n,sizeof(double));
	}

	for (i=0; i<=P.n; i++)
		R->t[i]=P.t[i];

	while (R->n>=S.n)
	{
		z=R->t[R->n]/S.t[S.n];
		printf("z = %g\n",z);
		Q->t[R->n-S.n]=z;
		R->t[R->n]=0;
		for (i=1; i<=S.n; i++)
			R->t[R->n-i]=R->t[R->n-i]-z*S.t[S.n-i];
		Degree(R);
	}

	if (R->n<0)
		return 0;
	else
		return 1;
}


void PolyMake (vector v, poly * P)
{
	int i,j,n=v.n;
	if (P->t==NULL)
	{
		P->n=n;
		P->t=calloc(1+P->n,sizeof(double));	
	}

	P->t[n]=1;
	for (i=0; i<n; i++)
		P->t[i]=0;
	for (i=n; i>=0; i--)
	{
		for (j=i; j<n; j++)
			P->t[j]=P->t[j]-v.t[i]*P->t[j+1];
	}
}


void PolyGet (poly P)
{
	int i;
	for (i=P.n-1; i>=0; i--)
		scanf("%lf",&P.t[i]);
}


void PolyPrint (poly P)
{
	Degree(&P);
	if (P.n<0)
	{	
		puts("0");
		return;
	}

	int i;
	double z=P.t[0],x;
	if (ab(z)>eps) 
		printf("%g",z);

	for (i=1; i<=P.n; i++)
	{
		z=P.t[i];
		if (z>eps) fputs("+",stdout);
		else if (z>-eps) continue;
		x=ab(z-1);
		if (i>1)
		{
			if (x<eps)
				printf("X^%d",i);
			else
				printf("%.15gX^%d",z,i);
		}
		else
		{
			if (x<eps)	
				fputs("X",stdout);
			else
				printf("%gX",z);
		}
	}

	puts("");
}


void PolyDer (poly P, poly * R)
{
	int i,n=P.n;
	if (R->t==NULL)
	{
		R->n=P.n;
		PolyInit(R);
	}

	for (i=1; i<=n; i++)
		R->t[i-1]=i*P.t[i];
}

void PolyAntiDer (poly P, poly * R)
{
	int i,n=P.n+1;
	if (R->t==NULL)
	{
		R->n=P.n;
		R->t=calloc(1+P.n,sizeof(double));
	}

	R->t[0]=0;
	for (i=1; i<=n; i++)
		R->t[i]=P.t[i-1]/i;
}


double PolyInt (poly P, double a, double b)
{
	int i,n=P.n+1;
	double x=a,y=b,r,s=0;
	for (i=1; i<=n; i++)
	{
		r=P.t[i-1]/i;
		s+=r*(y-x);
		x=x*a;
		y=y*b;
	}
	return s;
}

#endif
