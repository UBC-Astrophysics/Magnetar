#include <math.h>
#include <stdio.h>

#define MAXSTP 10000
#define TINY 1.0e-30
#ifdef RHL
#define SMALL 1e-6			/* don't ask for tiny step */
#endif
#define RETURN				/* cleanup and exit */		\
   for (i=1;i<=nvar;i++) ystart[i]=y[i];				\
   if (kmax) {								\
      xp[++kount]=x;							\
      for (i=1;i<=nvar;i++) yp[i][kount]=y[i];				\
   }									\
   free_vector(dydx,1,nvar);						\
   free_vector(y,1,nvar);						\
   free_vector(yscal,1,nvar);						\
   return;


int kmax=0,kount=0;  /* defining declaration */
double *xp=0,**yp=0,dxsav=0;  /* defining declaration */

void odeint(ystart,nvar,x1,x2,eps,h1,hmin,nok,nbad,derivs,rkqc)
double ystart[],x1,x2,eps,h1,hmin;
int nvar,*nok,*nbad;
void (*derivs)();	/* ANSI: void (*derivs)(double,double *,double *); */
void (*rkqc)(); 	/* ANSI: void (*rkqc)(double *,double *,int,double *,double,
				double,double *,double *,double *,void (*)()); */
{
	int nstp,i;
	double xsav,x,hnext,hdid,h;
	double *yscal,*y,*dydx,*vector();
	void nrerror(),free_vector();

	yscal=vector(1,nvar);
	y=vector(1,nvar);
	dydx=vector(1,nvar);
	x=x1;
	h=(x2 > x1) ? fabs(h1) : -fabs(h1);
	*nok = (*nbad) = kount = 0;
	for (i=1;i<=nvar;i++) y[i]=ystart[i];
	if (kmax > 0) xsav=x-dxsav*2.0;
	for (nstp=1;nstp<=MAXSTP;nstp++) {
		(*derivs)(x,y,dydx);
		for (i=1;i<=nvar;i++)
			yscal[i]=fabs(y[i])+fabs(dydx[i]*h)+TINY;
		if (kmax > 0) {
			if (fabs(x-xsav) > fabs(dxsav)) {
				if (kount < kmax-1) {
					xp[++kount]=x;
					for (i=1;i<=nvar;i++) yp[i][kount]=y[i];
					xsav=x;
				}
			}
		}
		if ((x+h-x2)*(x+h-x1) > 0.0) {
			h=x2-x;
#ifdef RHL
		   if(h > SMALL*(x2 - x1)) {
		      (*rkqc)(y,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext,derivs);
		   } else {
		      RETURN;
		   }
		} else {
		   (*rkqc)(y,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext,derivs);
		}
#else
		}
		(*rkqc)(y,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext,derivs);
#endif
		if (hdid == h) ++(*nok); else ++(*nbad);
		if ((x-x2)*(x2-x1) >= 0.0) {
			RETURN;
		}
		if (fabs(hnext) <= hmin) {
		  fprintf(stderr,"Step size too small in ODEINT\n");
		  RETURN;
		}
		h=hnext;
	}
        fprintf(stderr,"Too many steps in routine ODEINT\n");
}

#undef MAXSTP
#undef TINY
