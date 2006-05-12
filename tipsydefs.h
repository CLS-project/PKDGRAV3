#ifndef TIPSYDEFS_HINCLUDED
#define TIPSYDEFS_HINCLUDED

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
    int nbodies ;
    int ndim ;
    int nsph ;
    int ndark ;
    int nstar ;
} ;

#endif





