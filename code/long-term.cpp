#include "SpaceHub/src/spaceHub.hpp"
#include "SpaceHub/src/taskflow/taskflow.hpp"
using namespace hub;
using namespace unit;
using namespace callback;
using namespace force;
using namespace orbit;
using f = Interactions<NewtonianGrav, Tidal, PN1>;                   // add static tidal force
using Solver = methods::DefaultMethod<f, particles::TideParticles>;  // use the corresponding tidal particles
using Particle = Solver::Particle;

// Q_max: max closest approach
double calc_max_impact_parameter(double Q_max, double v_inf, double M_tot) {
    return Q_max * sqrt(1 + 2 * M_tot / v_inf / v_inf / Q_max);
}

void job(std::string job_name, size_t thread_id, size_t scattering_num) {
    double v_inf = 0.1 * 3.66_kms;
    double r_start = 1500_AU;
    double a_j1 = 5_AU;
    double a_j2 = 15_AU;
    double a_n = 45_AU;

    std::fstream post_flyby_file5("full-5-" + job_name + "-" + std::to_string(thread_id) + ".txt", std::ios::out);
    std::fstream post_flyby_file6("full-6-" + job_name + "-" + std::to_string(thread_id) + ".txt", std::ios::out);
    std::fstream post_flyby_file7("full-7-" + job_name + "-" + std::to_string(thread_id) + ".txt", std::ios::out);
    std::fstream post_flyby_file8("full-8-" + job_name + "-" + std::to_string(thread_id) + ".txt", std::ios::out);

    std::fstream* files[4] = {&post_flyby_file5, &post_flyby_file6, &post_flyby_file7, &post_flyby_file8};

    for (auto f : files) {
        print(*f, std::setprecision(16));
    }

    print(std::cout, "starting job on thread ", thread_id, " \n");

    for (size_t i = 0; i < scattering_num; ++i) {
        Particle star{1_Ms, 1_Rs, 0, 0};
        Particle jupiter1{1_Mj, 7.1492e4_km, 0.25, 0.66_sec};
        Particle jupiter2{1_Mj, 7.1492e4_km, 0, 0};
        Particle neptune{17.147_Me, 24622_km, 0, 0};
        Particle intruder{1_Ms, 1_Rs, 0, 0};

        double coll_R = star.radius + jupiter1.radius;

        // create planetary system with two giant planets
        double planet_inc = random::Uniform(0, consts::pi);

        auto jupiter1_orb = Elliptic(star.mass, jupiter1.mass, a_j1, 0.0, planet_inc, isotherm, isotherm, isotherm);

        auto jupiter2_orb =
            Elliptic(star.mass + jupiter1.mass, jupiter2.mass, a_j2, 0.0, planet_inc, isotherm, isotherm, isotherm);

        auto neptune_orb =
            Elliptic(M_tot(star, jupiter1, jupiter2), neptune.mass, a_n, 0.0, planet_inc, isotherm, isotherm, isotherm);

        move_particles(jupiter1_orb, jupiter1);

        // move_to_COM_frame(star, jupiter1);

        move_particles(jupiter2_orb, jupiter2);

        // move_to_COM_frame(star, jupiter1, jupiter2);

        move_particles(neptune_orb, neptune);

        move_to_COM_frame(star, jupiter1, jupiter2, neptune);

        // create binary star
        // create scattering hyperbolic orbit
        double b_max = calc_max_impact_parameter(a_n * 5, v_inf, M_tot(star, jupiter1, jupiter2, neptune, intruder));

        auto incident_orb =
            scattering::incident_orbit(M_tot(star, jupiter1, jupiter2, neptune), intruder.mass, v_inf, b_max, r_start);

        move_particles(incident_orb, intruder);

        move_to_COM_frame(star, jupiter1, jupiter2, neptune, intruder);

        double scattering_t_end = 2 * time_to_periapsis(group(star, jupiter1, jupiter2), intruder);

        Solver sim{0, star, jupiter1, jupiter2, neptune, intruder};

        bool remain_three = false;

        Solver::RunArgs args;

        args.add_stop_condition(scattering_t_end);

        args.add_stop_point_operation([&](auto& ptc, auto h) {
            auto aj1 =
                calc_semi_major_axis(ptc.mass(0) + ptc.mass(1), ptc.pos(0) - ptc.pos(1), ptc.vel(0) - ptc.vel(1));

            auto aj2 =
                calc_semi_major_axis(ptc.mass(0) + ptc.mass(2), ptc.pos(0) - ptc.pos(2), ptc.vel(0) - ptc.vel(2));

            auto an = calc_semi_major_axis(ptc.mass(0) + ptc.mass(3), ptc.pos(0) - ptc.pos(3), ptc.vel(0) - ptc.vel(3));

            if (aj1 > 0 && aj2 > 0 && an > 0) {
                remain_three = true;
            }
        });

        sim.run(args);

        if (remain_three == true) {  // find proper initial condition;
            auto new_init_condition = sim.particles().to_AoS();

            new_init_condition.resize(4);  // delete the intruder

            move_to_COM_frame(new_init_condition);

            Solver long_sim{0, new_init_condition};

            Solver::RunArgs long_args;

            long_args.add_stop_condition(1e8_year);

            int event_tag = 0;

            auto init_ptc = long_sim.particles();

            auto stop_check = [&](auto& ptc, auto h) -> bool {
                auto [a, e] = calc_a_e(ptc.mass(0) + ptc.mass(1), ptc.pos(0) - ptc.pos(1), ptc.vel(0) - ptc.vel(1));

                auto [aj2, ej2] = calc_a_e(ptc.mass(0) + ptc.mass(2), ptc.pos(0) - ptc.pos(2), ptc.vel(0) - ptc.vel(2));

                auto [an, en] = calc_a_e(ptc.mass(0) + ptc.mass(3), ptc.pos(0) - ptc.pos(3), ptc.vel(0) - ptc.vel(3));

                if (a < 0 || aj2 < 0 || an < 0) {  // planet ejection
                    event_tag = 3;
                    return true;
                }

                if (0 < a && a * (1 - e) <= coll_R) {  // collision
                    event_tag = 2;
                    return true;
                }

                if (0 < a && a * (1 - e) <= 0.05_AU && e < 0.1) {  // hot jupiter candidate
                    event_tag = 1;
                    return true;
                }
                return false;
            };

            long_args.add_stop_condition(StepSlice(stop_check, 1000));

            double time_stamp[4] = {1e5_year, 1e6_year, 1e7_year, 1e10_year};
            auto writing_file = files;
            double* out_time = time_stamp;

            auto write_callback = [&](auto& ptc, auto h) {
                if (ptc.time() > *out_time) {
                    print(**writing_file, jupiter1_orb, ',', jupiter2_orb, ',', neptune_orb, ',', incident_orb, '\n');
                    print(**writing_file, init_ptc);
                    print(**writing_file, event_tag, ',', ptc);
                    **writing_file << std::endl;
                    writing_file++;
                    out_time++;
                }
            };

            long_args.add_operation(write_callback);

            long_args.add_stop_point_operation([&](auto& ptc, auto h) {
                print(**writing_file, jupiter1_orb, ',', jupiter2_orb, ',', neptune_orb, ',', incident_orb, '\n');
                print(**writing_file, init_ptc);
                print(**writing_file, event_tag, ',', ptc);
                **writing_file << std::endl;
            });

            long_sim.run(long_args);
        } else {
            continue;
        }
    }
}

int main(int argc, char** argv) {
    size_t n = 1000000;  // total scattering number
    size_t job_num = 40;

    tf::Executor executor;

    std::string job_name(argv[1]);

    for (size_t i = 0; i < job_num; ++i) {
        executor.silent_async(job, job_name, i, n / job_num);
    }
    executor.wait_for_all();  // wait all jobs to be finished
    return 0;
}
