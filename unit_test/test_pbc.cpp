// test_PBC.cpp
#include "../src/initial.h"
#include "../src/input.h"
#include "../src/logging.h"
#include "../src/potential.h"
#include "../src/simulation.h"
#include <iostream>
#include "test.h"

void testPBC() {
    SimulationBox box(10.0, 10.0);
    Particle p(5.0, 5.0);

    for (int i = 0; i < 100000; ++i) {
        double dx = (rand() / double(RAND_MAX) - 0.5) * 0.1; // Random small step
        double dy = (rand() / double(RAND_MAX) - 0.5) * 0.1;

        p.updatePosition(dx, dy);
        box.applyPBC(p);

        if (p.x < 0 || p.x >= box.getLx() || p.y < 0 || p.y >= box.getLy()) {
            std::cout << "PBC Test Failed." << std::endl;
            return;
        }
    }

    std::cout << "PBC Test Passed." << std::endl;
}