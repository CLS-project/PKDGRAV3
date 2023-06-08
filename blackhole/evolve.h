#include "smooth/smooth.h"
#include "pkd.h"

#define NOT_ACCRETED -1

struct bhEvolvePack {
    blitz::TinyVector<double,3> position;
    blitz::TinyVector<double,3> velocity;
    float c;
    float fMass;
#ifdef ENTROPY_SWITCH
    float fDensity;
#endif
    float fPotential;
    uint8_t uRung;
    uint8_t iClass;
};

struct bhEvolveFlush {
#ifdef OLD_FB_SCHEME
    double E;
    double Uint;
#ifdef ENTROPY_SWITCH
    double S;
#endif
#else // OLD_FB_SCHEME
    float fAccFBEnergy;
#endif
};

void smBHevolve(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf);
void packBHevolve(void *vpkd,void *dst,const void *src);
void unpackBHevolve(void *vpkd,void *dst,const void *src);
void initBHevolve(void *vpkd,void *dst);
void flushBHevolve(void *vpkd,void *dst,const void *src);
void combBHevolve(void *vpkd,void *dst,const void *src);

void pkdBHIntegrate(PKD pkd, particleStore::ParticleReference &p, double dTime, double dDelta, double dBHRadiativeEff);
void pkdBHAccretion(PKD pkd, double dScaleFactor);

class BHAccretionCache : public mdl::CACHEhelper {
protected:
    virtual void init(void *dst) override {
        auto p = static_cast<remoteID *>(dst);
        p->iPid = NOT_ACCRETED;
    }
    virtual void combine(void *dst,const void *src,const void *key) override {
        auto p1 = static_cast<remoteID *>(dst);
        auto p2 = static_cast<const remoteID *>(src);

        // Is this the first accretion attempt for this particle?
        if (p1->iPid == NOT_ACCRETED && p2->iPid != NOT_ACCRETED) {
            p1->iPid = p2->iPid;
            p1->iIndex = p2->iIndex;
        }
        // If not, just keep the previous attempt
    }
public:
    explicit BHAccretionCache() : CACHEhelper(sizeof(remoteID),true) {}
};

class BHAccretor {
public:
    BHAccretor() = delete;

    explicit BHAccretor(PKD pkd,int32_t iIndex)
        : pkd(pkd),rID(fetch(iIndex)),bMDL(false),bAcquire(false) {}

    explicit BHAccretor(PKD pkd,int cid,int32_t iIndex,int32_t iPid,bool bAcquire)
        : pkd(pkd),cid(cid),bMDL(pkd->Self() != iPid),bAcquire(bAcquire) {
        if (bMDL) {
            if (bAcquire) rID = acquire(iIndex,iPid);
            else rID = fetch(iIndex,iPid);
        }
        else rID = fetch(iIndex);
    }

    BHAccretor(const BHAccretor &o) = delete;
    BHAccretor &operator=(const BHAccretor &o) = delete;
    BHAccretor(BHAccretor &&o) = default;
    BHAccretor &operator=(BHAccretor &&o) = delete;
    ~BHAccretor() { if (bAcquire) mdlRelease(pkd->mdl,cid,rID); }

    auto get_pid() const { return rID->iPid; }
    auto get_index() const { return rID->iIndex; }
    void set_pid(int32_t iPid) { rID->iPid = iPid; }
    void set_index(int32_t iIndex) { rID->iIndex = iIndex; }

    bool is_remote() { return rID->iPid != pkd->Self(); }
    bool has_accreted() { return rID->iPid != NOT_ACCRETED; }

private:
    remoteID *fetch(int32_t iIndex) {
        return &static_cast<remoteID *>(pkd->pLite)[iIndex];
    }
    remoteID *fetch(int32_t iIndex,int32_t iPid) {
        return static_cast<remoteID *>(mdlFetch(pkd->mdl,cid,iIndex,iPid));
    }
    remoteID *acquire(int32_t iIndex,int32_t iPid) {
        return static_cast<remoteID *>(mdlAcquire(pkd->mdl,cid,iIndex,iPid));
    }

    PKD pkd;
    remoteID *rID;
    int cid;
    bool bMDL, bAcquire;
};

