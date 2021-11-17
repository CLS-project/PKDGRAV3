
void msrSetFeedbackParam(MSR msr);
void smSNFeedback(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf) ;
void initSNFeedback(void *vpkd, void *vp);
void combSNFeedback(void *vpkd, void *p1,void *p2);
void pkdAddFBEnergy(PKD pkd, PARTICLE *p, SPHFIELDS *psph);
float SNFeedbackEfficiency(PKD pkd, float Z, float rho);
