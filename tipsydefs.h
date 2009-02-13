#ifndef TIPSYDEFS_HINCLUDED
#define TIPSYDEFS_HINCLUDED
#define TIPSYDEFS_H_MODULE_ID "$Id$"

struct gas_particle {
    float mass;
    float pos[3];
    float vel[3];
    float rho;
    float temp;
    float hsmooth;
    float metals ;
    float phi ;
    } ;

struct dark_particle {
    float mass;
    float pos[3];
    float vel[3];
    float eps;
    float phi ;
    } ;

struct star_particle {
    float mass;
    float pos[3];
    float vel[3];
    float metals ;
    float tform ;
    float eps;
    float phi ;
    } ;

struct dump {
    double time ;
    unsigned nbodies ;
    unsigned ndim ;
    unsigned nsph ;
    unsigned ndark ;
    unsigned nstar ;
    unsigned pad ;
    } ;

#endif





