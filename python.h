#ifndef PYTHON_H
#define PYTHON_H
#define PYTHON_H_MODULE_ID "$Id$"
typedef void *PPY;
void ppyInitialize(PPY *ppy, MSR msr,double dTime);
void ppyFinish(PPY ppy);
void ppyRunScript(PPY ppy,const char *achFilename);
#endif
