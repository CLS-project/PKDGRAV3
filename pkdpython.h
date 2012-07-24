#ifndef PYTHON_H
#define PYTHON_H
typedef void *PPY;
void ppyInitialize(PPY *ppy, MSR msr,double dTime);
void ppyFinish(PPY ppy);
void ppyRunScript(PPY ppy,int argc,char *argv[]);
#endif
