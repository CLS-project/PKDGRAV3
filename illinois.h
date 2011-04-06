#ifndef ILLINOIS_H
#define ILLINOIS_H

typedef struct {
    double r, s;
    double fr, fs;
    double t;
    } ILLINOIS;

double illinoisInitialize(ILLINOIS *ctx,double r,double fr,double s,double fs);
double illinoisIterate(ILLINOIS *ctx,double ft);
double illinois(double (*func)(double,void *),void *ctx,double r,double s,double xacc,double yacc,int *pnIter);

#endif
