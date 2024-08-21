#include "test.h"

int main() {
    testEnergyWCA();
    testEnergyLJ();
    // testEnergyYukawa();
    testBoundaryEnergy();
    testPBC();
    testMonteCarloEnergyReduction();
    testBoundaryForce();
    testForceWCA();
    return 0;
}
