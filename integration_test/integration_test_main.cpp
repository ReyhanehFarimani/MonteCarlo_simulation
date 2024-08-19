#include "integration_test.h"

int main() {
    // Run the seed reproducibility test
    seed_test();
    cell_test_NVT();
    return 0;
}
