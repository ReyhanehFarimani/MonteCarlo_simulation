// MC_parallel.cpp  (optimized lightweight MPI profiler + coarser hot-path scopes)

#include "MC_parallel.h"
#include <array>
#include <algorithm>
#include <random>
#include <cmath>
#include <cassert>
#include <iostream>
#include <unordered_map>
#include <vector>
#include <string>
#include <cstdio>
#include <cstring>

// --------------------------- built-in headerless MPI profiler ------------------
// Disable entirely by compiling with -DNPROFILE.
#ifndef NPROFILE
namespace {
    struct __ProfRec { double sec{0.0}; unsigned long long calls{0}; };

    // Hash/eq for literal-pointer keys (labels must be stable string literals)
    struct __PtrHash { size_t operator()(const char* s) const noexcept {
        return std::hash<const void*>()((const void*)s);
    }};
    struct __PtrEq { bool operator()(const char* a, const char* b) const noexcept { return a==b; }};

    class __ProfilerMPI {
    public:
        explicit __ProfilerMPI(MPI_Comm comm=MPI_COMM_WORLD): comm_(comm) {
            MPI_Comm_rank(comm_, &rank_); MPI_Comm_size(comm_, &size_);
        }
        static __ProfilerMPI& get(){ static __ProfilerMPI inst; return inst; }

        struct Scope {
            const char* name; double t0;
            Scope(const char* n): name(n), t0(MPI_Wtime()) {}
            ~Scope(){ auto& P=__ProfilerMPI::get(); auto& r=P.data_[name]; r.sec+=MPI_Wtime()-t0; r.calls+=1; }
        };

        void begin(const char* n){ begins_[n]=MPI_Wtime(); }
        void end  (const char* n){
            auto it=begins_.find(n); if(it==begins_.end()) return;
            auto& r=data_[n]; r.sec += MPI_Wtime()-it->second; r.calls += 1; begins_.erase(it);
        }

        void dump(const char* txt="prof_functions_ranked.txt",
                  const char* csv="prof_functions_ranked.csv"){
            // collect labels (as literal pointers) and sort by content for stable cross-rank order
            std::vector<const char*> labels; labels.reserve(data_.size());
            for(auto& kv: data_) labels.push_back(kv.first);
            std::sort(labels.begin(), labels.end(),
                      [](const char* a, const char* b){ return std::strcmp(a,b) < 0; });

            // pack local
            std::vector<double> lsec(labels.size(),0.0), gsec(labels.size(),0.0);
            std::vector<unsigned long long> lc(labels.size(),0), gc(labels.size(),0);
            for(size_t i=0;i<labels.size();++i){
                auto it=data_.find(labels[i]);
                if(it!=data_.end()){ lsec[i]=it->second.sec; lc[i]=it->second.calls; }
            }
            // reduce
            MPI_Reduce(lsec.data(), gsec.data(), (int)labels.size(), MPI_DOUBLE, MPI_SUM, 0, comm_);
            MPI_Reduce(lc.data(),   gc.data(),   (int)labels.size(), MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, comm_);

            if(rank_==0){
                double total=0.0; for(double s: gsec) total+=s;
                struct Row{ const char* n; double s; unsigned long long c; double p; };
                std::vector<Row> rows; rows.reserve(labels.size());
                for(size_t i=0;i<labels.size();++i)
                    rows.push_back({labels[i], gsec[i], gc[i], total>0?100.0*gsec[i]/total:0.0});
                std::sort(rows.begin(), rows.end(), [](const Row&a,const Row&b){return a.s>b.s;});

                if(FILE* f=std::fopen(txt,"w")){
                    std::fprintf(f,"# Total aggregated wall-time (all ranks): %.6f s\n", total);
                    std::fprintf(f,"%-36s %12s %12s %8s\n","Function/Block","Time[s]","Calls","%");
                    for(auto& r: rows) std::fprintf(f,"%-36s %12.6f %12llu %7.2f\n", r.n, r.s, r.c, r.p);
                    std::fclose(f);
                }
                if(FILE* c=std::fopen(csv,"w")){
                    std::fprintf(c,"name,time_seconds,calls,percent\n");
                    for(auto& r: rows) std::fprintf(c,"%s,%.6f,%llu,%.2f\n", r.n, r.s, r.c, r.p);
                    std::fclose(c);
                }
            }
            MPI_Barrier(comm_);
        }
    private:
        MPI_Comm comm_; int rank_{0}, size_{1};
        std::unordered_map<const char*,__ProfRec,__PtrHash,__PtrEq> data_;
        std::unordered_map<const char*,double,__PtrHash,__PtrEq>    begins_;
    };
} // anonymous

  #define PROF_INIT()          do{ (void)__ProfilerMPI::get(); }while(0)
  #define PROF_SCOPE(name)     __ProfilerMPI::Scope __prof_scope_##__LINE__(name)
  #define PROF_BEGIN(name)     __ProfilerMPI::get().begin(name)
  #define PROF_END(name)       __ProfilerMPI::get().end(name)
  #define PROF_DUMP(txt,csv)   __ProfilerMPI::get().dump(txt,csv)
#else
  #define PROF_INIT()        do{}while(0)
  #define PROF_SCOPE(x)      do{}while(0)
  #define PROF_BEGIN(x)      do{}while(0)
  #define PROF_END(x)        do{}while(0)
  #define PROF_DUMP(a,b)     do{}while(0)
#endif
// ------------------------------------------------------------------------------

// ---- lightweight counters (aggregated at end) --------------------------------
namespace {
    struct Counters {
        unsigned long long parity_matches{0};
        unsigned long long parity_misses{0};
        unsigned long long map_failures{0};
        unsigned long long k_exhausted{0};
        unsigned long long rng_index_draws{0};
        unsigned long long rng_displacements{0};
        unsigned long long pot_evals_owned{0};
        unsigned long long pot_evals_point{0};
    };

    static Counters __cnt;

    inline void dump_counters(MPI_Comm comm, int rank) {
        Counters g{};
        MPI_Reduce(&__cnt.parity_matches, &g.parity_matches, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, comm);
        MPI_Reduce(&__cnt.parity_misses,  &g.parity_misses,  1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, comm);
        MPI_Reduce(&__cnt.map_failures,   &g.map_failures,   1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, comm);
        MPI_Reduce(&__cnt.k_exhausted,    &g.k_exhausted,    1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, comm);
        MPI_Reduce(&__cnt.rng_index_draws,&g.rng_index_draws,1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, comm);
        MPI_Reduce(&__cnt.rng_displacements,&g.rng_displacements,1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, comm);
        MPI_Reduce(&__cnt.pot_evals_owned,&g.pot_evals_owned,1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, comm);
        MPI_Reduce(&__cnt.pot_evals_point,&g.pot_evals_point,1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, comm);

        if (rank==0) {
            if (FILE* f=std::fopen("prof_counters.txt","w")) {
                std::fprintf(f, "# aggregated counters across ranks\n");
                std::fprintf(f, "parity_matches     %llu\n", g.parity_matches);
                std::fprintf(f, "parity_misses      %llu\n", g.parity_misses);
                std::fprintf(f, "map_failures       %llu\n", g.map_failures);
                std::fprintf(f, "k_exhausted        %llu\n", g.k_exhausted);
                std::fprintf(f, "rng_index_draws    %llu\n", g.rng_index_draws);
                std::fprintf(f, "rng_displacements  %llu\n", g.rng_displacements);
                std::fprintf(f, "pot_evals_owned    %llu\n", g.pot_evals_owned);
                std::fprintf(f, "pot_evals_point    %llu\n", g.pot_evals_point);
                std::fclose(f);
            }
            if (FILE* c=std::fopen("prof_counters.csv","w")) {
                std::fprintf(c,"name,value\n");
                std::fprintf(c,"parity_matches,%llu\n", g.parity_matches);
                std::fprintf(c,"parity_misses,%llu\n", g.parity_misses);
                std::fprintf(c,"map_failures,%llu\n", g.map_failures);
                std::fprintf(c,"k_exhausted,%llu\n", g.k_exhausted);
                std::fprintf(c,"rng_index_draws,%llu\n", g.rng_index_draws);
                std::fprintf(c,"rng_displacements,%llu\n", g.rng_displacements);
                std::fprintf(c,"pot_evals_owned,%llu\n", g.pot_evals_owned);
                std::fprintf(c,"pot_evals_point,%llu\n", g.pot_evals_point);
                std::fclose(c);
            }
        }
        MPI_Barrier(comm);
    }
}
// ------------------------------------------------------------------------------

// --- full rebuild + migrate + refresh ---
static inline void rebuild_migrate_refresh_(CellListParallel& cl,
                                            ParticleExchange&  pex,
                                            std::vector<Particle>& owned,
                                            SimulationBox& /*box*/)
{
    // Rebuild bins for current owned state
    {
        PROF_SCOPE("CL::buildInterior(full-rebuild)");
        cl.buildInterior(owned);
    }
    // Migrate ownership across ranks (based on current positions)
    {
        PROF_SCOPE("PX::migrate");
        pex.migrate(owned, cl);
    }
    // Rebuild bins again to reflect the post-migration owned set
    {
        PROF_SCOPE("CL::buildInterior(post-migrate)");
        cl.buildInterior(owned);
    }
    // Refresh ghost layer (and its bins)
    {
        PROF_SCOPE("PX::refreshGhosts");
        pex.refreshGhosts(owned, cl);
    }
}


MonteCarloNVT_MPI::MonteCarloNVT_MPI(MPI_Comm comm,
                                     SimulationBox& box,
                                     CellListParallel& cl,
                                     ParticleExchange& pex,
                                     ThermodynamicCalculatorParallel& thermo,
                                     std::vector<Particle>& owned,
                                     RNG_parallel& rng,
                                     LoggingTrajMPI* traj,
                                     LoggingDataMPI* data,
                                     const Params& p)
: comm_(comm)
, box_(box)
, cl_(cl)
, pex_(pex)
, thermo_(thermo)
, owned_(owned)
, rng_(rng)
, p_(p)
, traj_(traj)
, data_(data)
{
    PROF_INIT(); // ensure profiler constructed

    MPI_Comm_rank(comm_, &rank_);
    MPI_Comm_size(comm_, &size_);

    beta_ = 1.0 / thermo_.getTemperature();

    // --- Init order: migrate first, then build, then refresh ghosts ---
    // (cl_ already knows geometry; bins can be empty for migration)
    {
        PROF_SCOPE("PX::migrate(init)");
        pex_.migrate(owned_, cl_);       // drop particles leaving this rank, receive arrivals
    }
    {
        PROF_SCOPE("CL::buildInterior(init)");
        cl_.buildInterior(owned_);       // build bins for current owned
    }
    {
        PROF_SCOPE("PX::refreshGhosts(init)");
        pex_.refreshGhosts(owned_, cl_); // construct halo layer (ghost bins built inside)
    }
}

double MonteCarloNVT_MPI::run(std::size_t nsweeps, int start_timestep)
{
    PROF_SCOPE("MC::run");

    long long acc_loc = 0, att_loc = 0;
    int timestep = start_timestep;

    const std::array<Parity,4> base = {
        Parity::EvenEven, Parity::EvenOdd, Parity::OddEven, Parity::OddOdd
    };

    const int nloc_target = static_cast<int>(owned_.size());
    const int micro_per_attempt = 4;

    for (std::size_t sweep = 0; sweep < nsweeps; ++sweep) {
        PROF_BEGIN("MC::sweep");

        int attempts_done = 0;

        while (attempts_done < nloc_target) {

            std::array<Parity,4> order = base;
            // Keeping deterministic parity order; shuffling is commented out for now.
            // If you re-enable, do it once per sweep, not per attempt.

            for (Parity par : order) {
                PROF_SCOPE("MC::parity_step");
                (void)do_one_parity_step_(par, acc_loc, att_loc);
            }

            // Optional: occasional lightweight ghost refresh (without migration)
            if (p_.halo_every > 0 && (((attempts_done / micro_per_attempt) + 1) % p_.halo_every) == 0) {
                PROF_SCOPE("PX::refreshGhosts(periodic)");
                pex_.refreshGhosts(owned_, cl_);
            }

            // Periodic full rebuild + migrate + refresh
            if (p_.rebuild_every_attempts > 0 &&
                ((attempts_done + 1) % p_.rebuild_every_attempts) == 0)
            {
                PROF_SCOPE("MC::rebuild_migrate_refresh");
                rebuild_migrate_refresh_(cl_, pex_, owned_, box_);
            }

            attempts_done += 1; // one "attempt" = 4 parity micro-steps
        }

        if (p_.out_every && ((sweep + 1) % p_.out_every == 0)) {
            PROF_SCOPE("MC::logging");
            if (traj_) traj_->log_dump(owned_, box_, timestep);
            if (data_) data_->log_step(owned_, box_, cl_, pex_, thermo_, timestep);
        }

        ++timestep;
        // if (rank_ == 0) std::cout<<"timestep:\t"<<timestep<<std::endl;

        PROF_END("MC::sweep");
    }

    long long acc_glob = 0, att_glob = 0;
    {
        PROF_SCOPE("MC::allreduce_accept");
        MPI_Allreduce(&acc_loc, &acc_glob, 1, MPI_LONG_LONG_INT, MPI_SUM, comm_);
    }
    {
        PROF_SCOPE("MC::allreduce_attempt");
        MPI_Allreduce(&att_loc, &att_glob, 1, MPI_LONG_LONG_INT, MPI_SUM, comm_);
    }

    // Dump aggregated timing once per run (rank 0 writes files)
    PROF_DUMP("prof_functions_ranked.txt", "prof_functions_ranked.csv");

    // Dump counters (rank 0 writes files)
    dump_counters(comm_, rank_);

    return (att_glob > 0) ? static_cast<double>(acc_glob) / static_cast<double>(att_glob) : 0.0;
}


bool MonteCarloNVT_MPI::do_one_parity_step_(Parity par, long long& acc, long long& att)
{
    PROF_SCOPE("MC::do_one_parity_step");

    if (owned_.empty()) { ++att; return false; }

    // Try up to K times to find a particle whose cell has the requested parity
    constexpr int K = 4;
    int pidx = -1, ix = 0, iy = 0;

    // Coarsened candidate-selection timing (remove tiny per-op scopes)
    {
        PROF_SCOPE("MC::candidate_selection");
        for (int t = 0; t < K; ++t) {
            __cnt.rng_index_draws++;
            const int cand = (int)rng_.randint(0, (int)owned_.size()-1);

            bool mapped = cl_.mapToLocalCell(owned_[cand].x, owned_[cand].y, ix, iy);
            if (!mapped) { __cnt.map_failures++; continue; }
            if (ix < 1 || ix > cl_.nxInterior() || iy < 1 || iy > cl_.nyInterior()) {
                __cnt.map_failures++;
                continue;
            }

            const bool ex = ((ix-1) % 2) == 0;
            const bool ey = ((iy-1) % 2) == 0;
            const Parity cur = ex ? (ey ? Parity::EvenEven : Parity::EvenOdd)
                                  : (ey ? Parity::OddEven  : Parity::OddOdd);
            const bool match = (cur == par);

            if (match) {
                __cnt.parity_matches++;
                pidx = cand; break;
            } else {
                __cnt.parity_misses++;
            }
        }
    }

    if (pidx < 0) __cnt.k_exhausted++;

    bool ok = false;
    if (pidx >= 0) {
        PROF_SCOPE("MC::try_displacement_call");
        ok = try_displacement_(pidx);
    }

    ++att; if (ok) ++acc;
    return ok;
}


bool MonteCarloNVT_MPI::try_displacement_(int i)
{
    PROF_SCOPE("MC::try_displacement");

    // current local contribution for particle i
    PROF_BEGIN("MC::local_energy_of");
    const double U_old = local_energy_of_(i);
    PROF_END("MC::local_energy_of");

    // propose trial
    double dx = 0.0, dy = 0.0;
    if (p_.delta > 0.0) {
        PROF_SCOPE("RNG::displacement");
        __cnt.rng_displacements += 2; // dx and dy
        dx = rng_.uniform(-p_.delta, p_.delta);
        dy = rng_.uniform(-p_.delta, p_.delta);
    }

    Particle old = owned_[i];
    Particle trial = old;
    trial.x += dx; trial.y += dy;

    {
        PROF_SCOPE("Box::applyPBC");
        box_.applyPBC(trial);
    }

    // new local contribution at trial position
    PROF_BEGIN("MC::local_energy_of_point");
    const double U_new = local_energy_of_point_(i, trial.x, trial.y);
    PROF_END("MC::local_energy_of_point");

    const double dU    = U_new - U_old;

    // Metropolis
    bool accept = false;
    {
        PROF_SCOPE("MC::metropolis");
        accept = (dU <= 0.0) || (rng_.uniform01() < std::exp(-beta_ * dU));
    }
    if (accept) {
        owned_[i] = trial;
        {
            PROF_SCOPE("CL::onAcceptedMove");
            cl_.onAcceptedMove(i, old, trial);
        }
        // If you later add halo queuing, keep it coarse-grained.
        // pex_.queueGhostUpdateCandidate(trial);
    }
    return accept;
}

double MonteCarloNVT_MPI::local_energy_of_(int i) const
{
    PROF_SCOPE("MC::local_energy_of");

    double U = 0.0;

    const auto type   = thermo_.getPotentialType();
    const float fp    = static_cast<float>(thermo_.getFPrime());
    const float fpd   = static_cast<float>(thermo_.getFPrimeAttraction());
    const float kappa = static_cast<float>(thermo_.getKappa());
    const float alpha = static_cast<float>(thermo_.getAlpha());

    PROF_BEGIN("CL::neighborsOfOwned");
    const auto nbs = cl_.neighborsOfOwned(i, owned_);
    PROF_END("CL::neighborsOfOwned");

    {
        PROF_SCOPE("MC::pair_loop_owned");
        for (const auto& pr : nbs) {
            const double r2 = pr.second;
            if (r2 <= 1e-24) continue; // paranoid guard
            __cnt.pot_evals_owned++;
            U += computePairPotential(r2, type, fp, fpd, kappa, alpha);
        }
    }
    return U;
}

double MonteCarloNVT_MPI::local_energy_of_point_(int i_moved, double x, double y) const
{
    PROF_SCOPE("MC::local_energy_of_point");

    double U = 0.0;

    const auto type   = thermo_.getPotentialType();
    const float fp    = static_cast<float>(thermo_.getFPrime());
    const float fpd   = static_cast<float>(thermo_.getFPrimeAttraction());
    const float kappa = static_cast<float>(thermo_.getKappa());
    const float alpha = static_cast<float>(thermo_.getAlpha());

    PROF_BEGIN("CL::neighborsOfPoint");
    const auto nbs = cl_.neighborsOfPoint(x, y, owned_);
    PROF_END("CL::neighborsOfPoint");

    {
        PROF_SCOPE("MC::pair_loop_point");
        for (const auto& pr : nbs) {
            const int    j  = pr.first;   // >=0 owned, <0 ghost
            const double r2 = pr.second;

            // exclude self if it appears at zero distance
            if (j >= 0 && j == i_moved) continue;
            if (r2 <= 1e-24) continue;

            __cnt.pot_evals_point++;
            U += computePairPotential(r2, type, fp, fpd, kappa, alpha);
        }
    }
    return U;
}

std::uint64_t MonteCarloNVT_MPI::bcast_seed_()
{
    std::uint64_t seed = 0;
    if (rank_ == 0) {
        // compose 64 bits from two 31-bit draws (simple, reproducible)
        const std::uint64_t a = static_cast<std::uint64_t>(rng_.randint(0, 0x7fffffff));
        const std::uint64_t b = static_cast<std::uint64_t>(rng_.randint(0, 0x7fffffff));
        seed = (a << 33) ^ (b << 1) ^ 0x9E3779B97F4A7C15ull;
        if (seed == 0) seed = 0xD1B54A32D192ED03ull;
    }
    MPI_Bcast(&seed, 1, MPI_UINT64_T, 0, comm_);
    return seed;
}
