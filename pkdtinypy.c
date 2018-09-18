#include "master.h"
#include "pkdtinypy.h"

static tp_obj missing(const char *func, const char *parm) {
    fprintf(stderr, "%s: missing parameter %s\n", func, parm);
    return tp_None;
    }


static tp_obj msr_tp_OutputOrbits(TP) {
    tp_obj v = TP_OBJ();
    int type = v.type;
    uint64_t ID[GET_PARTICLES_MAX];
    int nP=0;

    if (type == TP_NUMBER) { ID[nP++] = v.number.val; }
    else if (type == TP_LIST) {
	tp_obj i = tp_number(0);
	while(i.number.val < tp_len(tp,v).number.val) {
	    tp_obj id = tp_iter(tp,v,i);
	    if (id.type == TP_NUMBER) ID[nP++] = id.number.val;
	    i.number.val += 1;
	    }
	}
    else {}
    return tp_None;
    }

static tp_obj msr_tp_GenerateIC(TP) {
    tp_obj o, v = TP_DEFAULT(tp_None);
    int nGrid = 0;


    tp_obj msr_module = tp_get(tp, tp->modules, tp_string("msr"));
    tp_obj msr_data = tp_get(tp,msr_module,tp_string("__MSR__"));
    MSR msr = msr_data.data.val;

    /* We have named parameters: fetch them */
    if (v.type == TP_DICT) {
	if (tp_iget(tp,&o,v,tp_string("nGrid")) && o.type == TP_NUMBER)
	    nGrid = o.number.val;
	else return missing("GenerateIC","nGrid");
	}
    else {
    	fprintf(stderr,"GenerateIC: missing named parameters\n");
	}

    msrGenerateIC(msr);


    printf("========================================\n");
    printf("GenerateIC\n");
    printf("========================================\n");

    return tp_None;
    }
static tp_obj msr_tp_DomainDecompose(TP) {
    return tp_None;
    }
static tp_obj msr_tp_TreeBuild(TP) {
    return tp_None;
    }
static tp_obj msr_tp_Gravity(TP) {
    return tp_None;
    }
static tp_obj msr_tp_Fof(TP) {
    return tp_None;
    }


#define SERVICE(n) {#n, msr_tp_ ## n}
static struct {
    char name[32];
    tp_obj (*fcn)(TP);
    } msr_services[] = {
	SERVICE(GenerateIC),
	SERVICE(DomainDecompose),
	SERVICE(TreeBuild),
	SERVICE(Gravity),
	SERVICE(OutputOrbits),
	SERVICE(Fof),
    };

void tpyInitialize(TP, MSR msr) {
    tp_obj msr_module = tp_dict(tp);
    int i;

    for(i=0; i<sizeof(msr_services)/sizeof(msr_services[0]); ++i)
	tp_set(tp,msr_module, tp_string(msr_services[i].name), tp_fnc(tp,msr_services[i].fcn));
    tp_set(tp,msr_module,tp_string("__MSR__"),tp_data(tp,1,msr));
    tp_set(tp,tp->modules, tp_string("msr"), msr_module);  // Register the module

    }
