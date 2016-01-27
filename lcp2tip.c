#include <stdlib.h>
#include <stdio.h>
#include <errno.h>

typedef struct partLightCone {
    float pos[3];
    float vel[3];
    } LIGHTCONEP;

typedef struct {
    double dTime;
    unsigned nBodies;
    unsigned nDim;
    unsigned nSph;
    unsigned nDark;
    unsigned nStar;
    unsigned nPad;
    } tipsyHdr;

typedef struct {
    float mass;
    float pos[3];
    float vel[3];
    float eps;
    float phi;
    } tipsyDark;


int main(int argc, char *argv[]) {
    LIGHTCONEP p;
    FILE *fp;
    int j;

    tipsyHdr h;
    tipsyDark d;

    if (argc!=3) {
	fprintf(stderr,"Usage: cat files | %s outfile mass\n", argv[0]);
	return EINVAL;
	}

    h.dTime = 1.0;
    h.nBodies = 0;
    h.nDim = 3;
    h.nSph = h.nDark = h.nStar = 0;
    h.nPad = 0;

    d.mass = atof(argv[2]);
    d.eps = 0.0;
    d.phi = 0.0;

    fp = fopen(argv[1],"wb");
    if (fp==NULL) {
	perror(argv[1]);
	abort();
	}

    fseek(fp,sizeof(tipsyHdr),SEEK_SET);
    while(fread(&p,sizeof(p),1,stdin)) {
	for(j=0; j<3; ++j) {
	    d.pos[j] = p.pos[j];
	    d.vel[j] = p.vel[j];
	    }
	fwrite(&d,sizeof(d),1,fp);
	++h.nBodies;
	}
    h.nDark = h.nBodies;

    fseek(fp,0,SEEK_SET);
    fwrite(&h,sizeof(h),1,fp);

    fclose(fp);

    }

