#include "rebound.h"
#include "reboundx.h"
#include "math.h"

int main(int argc, char* argv[]){
	struct reb_simulation* sim = reb_simulation_create();
	struct reb_particle star = {0};

	// Sun
	star.m = 1.;
	reb_simulation_add(sim, star);

	// Spacecraft
	struct reb_particle planet = {0};  // add a planet on a circular orbit (with default units where G=1)
    planet.x = 1.;
    planet.vy = 1.1;
    reb_simulation_add(sim, planet);

	// Add external force
	struct rebx_extras* rebx = rebx_attach(sim);  // first initialize rebx
    struct rebx_force* sphere = rebx_load_force(rebx, "solid_sphere"); // add our new force
    rebx_add_force(rebx, sphere);

	// Set parameters and enable it
	rebx_set_param_double(rebx, &sphere->ap, "rho", 10.);
	rebx_set_param_double(rebx, &sphere->ap, "rad", 0.1);
	rebx_set_param_double(rebx, &sphere->ap, "x_cen", 0.1);
	rebx_set_param_double(rebx, &sphere->ap, "y_cen", 0.);
	rebx_set_param_double(rebx, &sphere->ap, "z_cen", 0.);
	rebx_set_param_int(rebx, &sim->particles[1].ap, "ext_enable", 1);

	// Integrate
	double tmax = 2.*M_PI*10000;
	reb_simulation_integrate(sim, tmax);
}