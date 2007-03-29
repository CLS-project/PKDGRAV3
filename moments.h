#ifndef MOMENTS_INCLUDED
#define MOMENTS_INCLUDED

#ifdef QUAD
typedef long double momFloat;
#define sqrt(x)	sqrtl(x)
#else
typedef double momFloat;
#endif

/*
 ** moment tensor components for reduced multipoles.
 */
typedef struct momReduced {
	momFloat m;
	momFloat xx,yy,xy,xz,yz;
	momFloat xxx,xyy,xxy,yyy,xxz,yyz,xyz;
	momFloat xxxx,xyyy,xxxy,yyyy,xxxz,yyyz,xxyy,xxyz,xyyz;
	} MOMR;

/*
 ** moment tensor components for complete multipoles.
 */
typedef struct momComplete {
	momFloat m;
	momFloat xx,yy,xy,xz,yz;
	momFloat xxx,xyy,xxy,yyy,xxz,yyz,xyz;
	momFloat xxxx,xyyy,xxxy,yyyy,xxxz,yyyz,xxyy,xxyz,xyyz;
	momFloat zz;
	momFloat xzz,yzz,zzz;
	momFloat xxzz,xyzz,xzzz,yyzz,yzzz,zzzz;
	} MOMC;


/*
 ** moment tensor components for reduced local expansion.
 */
typedef struct locReduced {
	momFloat m;
	momFloat x,y,z;
	momFloat xx,yy,xy,xz,yz;
	momFloat xxx,xyy,xxy,yyy,xxz,yyz,xyz;
	momFloat xxxx,xyyy,xxxy,yyyy,xxxz,yyyz,xxyy,xxyz,xyyz;
	} LOCR;

/*
 ** moment tensor components for complete local expansion.
 */
typedef struct locComplete {
	momFloat m;
	momFloat x,y,z;
	momFloat xx,yy,xy,xz,yz;
	momFloat xxx,xyy,xxy,yyy,xxz,yyz,xyz;
	momFloat xxxx,xyyy,xxxy,yyyy,xxxz,yyyz,xxyy,xxyz,xyyz;
	momFloat zz;
	momFloat xzz,yzz,zzz;
	momFloat xxzz,xyzz,xzzz,yyzz,yzzz,zzzz;
	} LOCC;


void momAddMomc(MOMC *,MOMC *);
void momAddMomr(MOMR *,MOMR *);
void momMulAddMomc(MOMC *,momFloat,MOMC *);
void momMulAddMomr(MOMR *,momFloat,MOMR *);
void momSubMomc(MOMC *,MOMC *);
void momSubMomr(MOMR *,MOMR *);
void momMakeMomc(MOMC *,momFloat,momFloat,momFloat,momFloat);
momFloat momMakeMomr(MOMR *,momFloat,momFloat,momFloat,momFloat);
void momOldMakeMomr(MOMR *,momFloat,momFloat,momFloat,momFloat);
void momShiftMomc(MOMC *,momFloat,momFloat,momFloat);
void momShiftMomr(MOMR *,momFloat,momFloat,momFloat);
void momShiftLocr(LOCR *,momFloat,momFloat,momFloat);
void momShiftLocr2(LOCR *,momFloat,momFloat,momFloat);
void momReduceMomc(MOMC *,MOMR *);
void momReduceLocc(LOCC *,LOCR *);
void momEvalMomr(MOMR *,momFloat,momFloat,momFloat,momFloat,
				 momFloat *,momFloat *,momFloat *,momFloat *);
void momClearLocc(LOCC *);
void momClearLocr(LOCR *);
void momClearMomr(MOMR *);
void momMomr2Momc(MOMR *,MOMC *);
void momLoccAddMomr(LOCC *,MOMC *,momFloat,momFloat,momFloat,momFloat);
void momLocrAddMomr(LOCR *,MOMR *,momFloat,momFloat,momFloat,momFloat);
void momEvalLocc(LOCC *,momFloat,momFloat,momFloat,
				 momFloat *,momFloat *,momFloat *,momFloat *); 
void momEvalLocr(LOCR *,momFloat,momFloat,momFloat,
				 momFloat *,momFloat *,momFloat *,momFloat *); 
void momSubLocc(LOCC *,LOCC *);
void momPrintMomc(MOMC *);
void momPrintMomr(MOMR *);
void momPrintLocc(LOCC *);
void momPrintLocr(LOCR *);

/*
** from newloc.c
*/
void momSymLocrAddMomr(LOCR *l1,LOCR *l2,MOMR *q1,MOMR *q2,momFloat dir,momFloat x,momFloat y,momFloat z);
void momNewLocrAddMomr(LOCR *l,MOMR *q,momFloat dir,momFloat x,momFloat y,momFloat z);

#endif
