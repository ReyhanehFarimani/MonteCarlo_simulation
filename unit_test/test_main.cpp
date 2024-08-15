#include "test.h"

int main() {
    testEnergyWCA();
    testEnergyLJ();
    // testEnergyYukawa();
    testBoundaryEnergy();
    testPBC();
    testMonteCarloEnergyReduction();

    return 0;
}
