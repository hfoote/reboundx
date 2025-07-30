/** * @file solid_sphere.c
 * @brief   An external force from a constant-density sphere.
 * @author  Hayden R. Foote <haydenfoote@arizona.edu>
 * 
 * @section     LICENSE
 * Copyright (c) 2015 Dan Tamayo, Hanno Rein
 *
 * This file is part of reboundx.
 *
 * reboundx is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * reboundx is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rebound.  If not, see <http://www.gnu.org/licenses/>.
 *
 * The section after the dollar signs gets built into the documentation by a script.  All lines must start with space * space like below.
 * Tables always must be preceded and followed by a blank line.  See http://docutils.sourceforge.net/docs/user/rst/quickstart.html for a primer on rst.
 * $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
 *
 * $External Point Mass$       // Effect category (must be the first non-blank line after dollar signs and between dollar signs to be detected by script).
 *
 * ======================= ===============================================
 * Authors                 H. R. Foote
 * Implementation Paper    None
 * Based on                Central Force
 * C Example               
 * Python Example          
 * ======================= ===============================================
 * 
 * Adds an external acceleration due to a constant-density sphere. 
 *
 * **Effect Parameters**
 * 
 * ============== ========================
 * rad            Radius of sphere
 * rho            Density of source
 * b_max          Max distance from center of simulation volume to integrate the force
 * x_cen          x-position of source
 * y_cen          y-position of source
 * x_cen          z-position of source
 * ============== ========================
 * **Particle Parameters**
 *
 * None
 * 
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "rebound.h"
#include "reboundx.h"

static void rebx_calculate_sphere_force(struct reb_simulation* const sim, struct reb_particle* const particles, const int N, 
                                        const double rad, const double rho, const double x_cen, 
                                        const double y_cen, const double z_cen){
    const double G = sim->G; 
    const double R3 = rad*rad*rad;
    const double M = rho * 4.0/3.0 * M_PI * R3;

    for (int i=0; i<N; i++){
        const struct reb_particle p = particles[i];
        const int* const ext_enable = rebx_get_param(sim->extras, p.ap, "ext_enable");
        
        if (ext_enable != NULL){
            const double dx = p.x - x_cen;
            const double dy = p.y - y_cen;
            const double dz = p.z - z_cen;
            const double r2 = dx*dx + dy*dy + dz*dz;
            double prefac;

            if (r2 < (rad*rad)){
                prefac = -G*M/R3;
            } 
            else {
                prefac = -G*M*pow(r2, (-1.5));
            }

            particles[i].ax += prefac*dx;
            particles[i].ay += prefac*dy;
            particles[i].az += prefac*dz;
        }
    }
}

void rebx_sphere_force(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* const particles, const int N){

    const double* const rad = rebx_get_param(sim->extras, force->ap, "rad");
    const double* const rho = rebx_get_param(sim->extras, force->ap, "rho");
    const double* const x_cen = rebx_get_param(sim->extras, force->ap, "x_cen");
    const double* const y_cen = rebx_get_param(sim->extras, force->ap, "y_cen");
    const double* const z_cen = rebx_get_param(sim->extras, force->ap, "z_cen");

    if ((rad != NULL) & (rho != NULL) & (x_cen != NULL) & (y_cen != NULL) & (z_cen != NULL)){
        rebx_calculate_sphere_force(sim, particles, N, *rad, *rho, *x_cen, *y_cen, *z_cen);  
    }
}

