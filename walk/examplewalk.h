#include "../TraversePST.h"
#include "basewalk.h"
#include "core/simd.h"
#include "pkd.h"

// The structure here is akin to the example interaction, but
// without an evaluation block

class ExampleWalkTraits {
public:
    struct name {
        static constexpr std::string_view value = "ExampleWalk";
    };

    struct Options {
        float time;
        float delta;

        int nReps;

        uint8_t uRungHi;
        uint8_t uRungLo;

        int iLocalRoot1;
        int iLocalRoot2;
    };

    struct output {
        double flops;
    };
};

template <typename Interaction>
class ExampleWalk : public BaseTreeWalk<ExampleWalkTraits, Interaction> {

public:
    explicit ExampleWalk(PST pst)
        : BaseTreeWalk<ExampleWalkTraits, Interaction>(pst) {}

private:
    // This would be a good moment to eliminate the magic numbers
    // and use an enum instead..
    enum outcome {
        OPEN_CELL,
    };

    // Mock opening criteria, only for compiling tests
    void openOutcomeBlock(treeStore::NodePointer k, int n, clBlock &blk) {
        n = (n + fvec::mask()) / fvec::width();
        for (auto i = 0; i < n; ++i) {
            blk.iOpen.v[i] = OPEN_CELL;
        }
    }

    void selectOutcome(clBlock blk, int outcome) {
        switch (outcome) {
        case OPEN_CELL:
            return;
        }
    }

    bool testAddNodeCL(treeStore::NodePointer k) {
        return (k->min_rung() <= this->options.uRungHi &&
                k->max_rung() >= this->options.uRungLo);
    }
};
