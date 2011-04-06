#include <stdio.h>
#include <math.h>

#include "illinois.h"

double illinoisInitialize(ILLINOIS *ctx,double r,double fr,double s,double fs) {
    ctx->r = r;
    ctx->s = s;
    ctx->fr = fr;
    ctx->fs = fs;
    ctx->t = (s*fr - r*fs)/(fr - fs);
    return ctx->t;
    }

int illinoisDone(ILLINOIS *ctx,double ft,double xacc,double yacc) {
    if (fabs(ctx->t-ctx->s)<=xacc) return 1;
    if (fabs(ft) <= yacc) return 1;
    return 0;
    }

double illinoisIterate(ILLINOIS *ctx,double ft) {
    double phis,phir,gamma;
    if (ft==0.0) return ctx->t;
    if (ft*ctx->fs < 0) {
	/*
	** Unmodified step.
	*/
	ctx->r = ctx->s;
	ctx->s = ctx->t;
	ctx->fr = ctx->fs;
	ctx->fs = ft;
	}
    else {
	/*
	** Modified step to make sure we do not retain the 
	** endpoint r indefinitely.
	*/
#if 1
	phis = ft/ctx->fs;
	phir = ft/ctx->fr;
	gamma = 1 - (phis/(1-phir));  /* method 3 */
	if (gamma < 0) gamma = 0.5;
#else
	gamma = 0.5;    /* illinois */
#endif
	ctx->fr *= gamma;
	ctx->s = ctx->t;
	ctx->fs = ft;
	}
    ctx->t = (ctx->s*ctx->fr - ctx->r*ctx->fs)/(ctx->fr - ctx->fs);
    return ctx->t;
    }


double illinois(double (*func)(double,void *),void *ctx,double r,double s,double xacc,double yacc,int *pnIter) {
    const int maxIter = 100;
    double v;
    ILLINOIS ictx;
    int nIter = 0;

    v = illinoisInitialize(&ictx,r,(*func)(r,ctx),s,(*func)(s,ctx));
    while(!illinoisDone(&ictx,v,xacc,yacc) && ++nIter <= maxIter) {
	v = illinoisIterate(&ictx,(*func)(v,ctx));
	}
    if (pnIter) *pnIter = nIter;
    return v;
    }

#ifdef TEST_ILLINOIS
double f(double v) {
    return v*v - 4.0;
    }

int main() {
    ILLINOIS ctx;
    double v;

    v = illinoisInitialize(&ctx,0.0,f(0.0),11125.0,f(11125.0));

    while(!illinoisDone(&ctx,v,1e-7,1e-7)) {
	printf("%.20g\n", v);
	v = illinoisIterate(&ctx,f(v));
	}

    printf("%g !!\n", v);





    return 0;
    }
#endif
