#include <random>
#include <vector>
#include "../serial_src/initial.h"
#include "../serial_src/logging.h"

// int main(){
//     const int N = 100;
//     const int nsteps = 10;
//     SimulationBox box(10.0, 10.0);
//     std::mt19937 gen(42);
//     std::uniform_real_distribution<> disx(0.0, box.getLx());
//     std::uniform_real_distribution<> disy(0.0, box.getLy());

//     std::vector<Particle> particles;
//     particles.reserve(N);
//     // initial random placement
//     for(int i = 0; i < N; ++i){
//         particles.emplace_back(disx(gen), disy(gen));
//     }

//     Logging logger("dump.log", "data.log");
//     for(int t = 0; t < nsteps; ++t){
//         // random walk and PBC
//         for(auto &p : particles){
//             double dx = disx(gen) - box.getLx()/2.0;
//             double dy = disy(gen) - box.getLy()/2.0;
//             p.updatePosition(dx, dy);
//             box.applyPBC(p);
//         }
//         logger.logPositions_dump(particles, box, t);
//     }
//     logger.close();
//     return 0;
// }