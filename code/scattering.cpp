#include "/home/yihanw/repositories/SpaceHub/src/spaceHub.hpp"
#include "/home/yihanw/repositories/SpaceHub/src/taskflow/taskflow.hpp"
using namespace hub;
using namespace unit;
using namespace callback;
using namespace force;
using namespace orbit;
using f = Interactions<NewtonianGrav>;
using Solver = methods::DefaultMethod<f>;
using Particle = Solver::Particle;

// Q_max: max closest approach
double calc_max_impact_parameter(double Q_max, double v_inf, double M_tot) {
    return Q_max * sqrt(1 + 2 * M_tot / v_inf / v_inf / Q_max);
}

void job(size_t thread_id, size_t scattering_num) {
    double v_inf = 3.66_kms;
    double r_start = 500_AU;
    double a_j1 = 5_AU;
    double a_j2 = 15_AU;
    double a_n = 45_AU;

    std::fstream post_flyby_file("post-flyby-" + std::to_string(thread_id) + ".txt", std::ios::out);

    print(post_flyby_file, std::setprecision(16));

    print(std::cout, "starting job on thread ", thread_id, " \n");

    for (size_t i = 0; i < scattering_num; ++i) {
        Particle star{1_Ms};
        Particle jupiter1{1_Mj};
        Particle jupiter2{1_Mj};
        Particle neptune{17.147_Me};
        Particle intruder{1_Ms};

        // create planetary system with two giant planets
        double planet_inc = random::Uniform(0, consts::pi);

        auto jupiter1_orb = Elliptic(star.mass, jupiter1.mass, a_j1, 0.0, planet_inc, isotherm, isotherm, isotherm);

        auto jupiter2_orb =
            Elliptic(star.mass + jupiter1.mass, jupiter2.mass, a_j2, 0.0, planet_inc, isotherm, isotherm, isotherm);

        auto neptune_orb =
            Elliptic(M_tot(star, jupiter1, jupiter2), neptune.mass, a_n, 0.0, planet_inc, isotherm, isotherm, isotherm);

        move_particles(jupiter1_orb, jupiter1);

        move_to_COM_frame(star, jupiter1);

        move_particles(jupiter2_orb, jupiter2);

        move_to_COM_frame(star, jupiter1, jupiter2);

        move_particles(neptune_orb, neptune);

        // create binary star
        // create scattering hyperbolic orbit

        double b_max = calc_max_impact_parameter(a_n * 2, v_inf, M_tot(star, jupiter1, jupiter2, neptune, intruder));

        auto incident_orb =
            scattering::incident_orbit(M_tot(star, jupiter1, jupiter2, neptune), intruder.mass, v_inf, b_max, r_start);

        move_particles(incident_orb, intruder);

        move_to_COM_frame(star, jupiter1, jupiter2, neptune, intruder);

        double scattering_t_end = 2 * time_to_periapsis(group(star, jupiter1, jupiter2), intruder);

        Solver sim{0, star, jupiter1, jupiter2, neptune, intruder};

        Solver::RunArgs args;

        args.add_stop_condition(scattering_t_end);

        args.add_stop_point_operation([&](auto& ptc, auto h) {
            auto [aj1, ej1, Lj1] =
                calc_a_e_L(ptc.mass(0), ptc.mass(1), ptc.pos(0) - ptc.pos(1), ptc.vel(0) - ptc.vel(1));

            auto [aj2, ej2, Lj2] =
                calc_a_e_L(ptc.mass(0), ptc.mass(2), ptc.pos(0) - ptc.pos(2), ptc.vel(0) - ptc.vel(2));

            auto [an, en, Ln] = calc_a_e_L(ptc.mass(0), ptc.mass(3), ptc.pos(0) - ptc.pos(3), ptc.vel(0) - ptc.vel(3));

            auto L_p = calc::calc_isolated_angular_momentum(ptc, std::array{0, 1, 2, 3});
            print(post_flyby_file, jupiter1_orb, ',', jupiter2_orb, ',', neptune_orb, ',', incident_orb, ',', aj1, ',',
                  ej1, ',', Lj1, ',', aj2, ',', ej2, ',', Lj2, ',', an, ',', en, ',', Ln, ',', L_p, '\n');
        });

        sim.run(args);
    }
}

int main() {
    size_t n = 5000000;  // total scattering number
    size_t job_num = 40;

    tf::Executor executor;

    for (size_t i = 0; i < job_num; ++i) {
        executor.silent_async(job, i, n / job_num);
    }
    executor.wait_for_all();  // wait all jobs to be finished
    return 0;
}
