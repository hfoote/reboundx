/** * @file ext_point_mass.c
 * @brief   An external force from a stationary point mass that can be placed at an arbitrary location.
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
 * Adds an external acceleration due to a stationary point mass, which is not attached to a particle. User specifies the location
 * and mass of the "ghost" source particle at creation. 
 *
 * **Effect Parameters**
 * 
 * ============== ========================
 * M_ext          Mass of source 
 * x_ext          x-position of source
 * y_ext          y-position of source
 * x_ext          z-position of source
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

static void rebx_calculate_point_mass_force(struct reb_simulation* const sim, struct reb_particle* const particles, const int N, const double ext_M, const double ext_x, const double ext_y, const double ext_z){
    const double G = sim->G; 

    for (int i=0; i<N; i++){
        const struct reb_particle p = particles[i];
        const int* const ext_enable = rebx_get_param(sim->extras, p.ap, "ext_enable");

        if (ext_enable != NULL){
            const double dx = p.x - ext_x;
            const double dy = p.y - ext_y;
            const double dz = p.z - ext_z;
            const double r2 = dx*dx + dy*dy + dz*dz;
            const double prefac = -G*ext_M*pow(r2, (-1.5));

            particles[i].ax += prefac*dx;
            particles[i].ay += prefac*dy;
            particles[i].az += prefac*dz;
        }
    }
}

void rebx_external_point_mass_force(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* const particles, const int N){

    const double* const ext_M = rebx_get_param(sim->extras, force->ap, "ext_M");
    const double* const ext_x = rebx_get_param(sim->extras, force->ap, "ext_x");
    const double* const ext_y = rebx_get_param(sim->extras, force->ap, "ext_y");
    const double* const ext_z = rebx_get_param(sim->extras, force->ap, "ext_z");

    if ((ext_M != NULL) & (ext_x != NULL) & (ext_y != NULL) & (ext_z != NULL)){
        rebx_calculate_point_mass_force(sim, particles, N, *ext_M, *ext_x, *ext_y, *ext_z);  
    }
}

