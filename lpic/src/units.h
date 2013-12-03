/*
   This file is part of LPIC++, a particle-in-cell code for
   simulating the interaction of laser light with plasma.

   Copyright (C) 1994-1997 Roland Lichters

   LPIC++ is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
*/

//////////////////////////////////////////////////////////////////////////////////////////
//
// dimensionless quantities used in lpi
//
//////////////////////////////////////////////////////////////////////////////////////////
//
//  coordinate                    x / lambda_0  ---> x
//
//  velocity                             v / c  ---> v
//                                                   u = gamma v
//
//  field                  e E / ( m omega c )  ---> E
//                           e B / ( m omega )  ---> B
//
//
//  potential                   e Phi / m c^2   ---> Phi
//
//  vector potential                 e A / m c  ---> A
//
//  density                            n / n_c  ---> n
//             n_c = m epsilon_0 omega^2 / e^2
//
//  charge density                  rho / e n_c ---> rho = Z_e n_e/n_c + Z_i n_i/n_c
//
//  current density             j / ( e c n_c ) ---> j
//
//
//
//
//  in these units, the 'longitudinal' 1D Maxwell Equations read:
//
//       div E = rho/epsilon_0      --->       div E = 2pi rho
//
//  grad^2 Phi = - rho/epsilon_0    ---> grad^2 Phi = - (2pi)^2 rho
//
//           E = - grad Phi         --->          E = - 1/2pi grad Phi
//
//
//
//
//  The total energy is given in the following dimensionless units:
//
//  total energy/(m_e c^2 n_c A dx) = 1/2 SUM (E^2+B^2) + SUM (m/m_e * n/n_c * (gamma-1))
//                                        all cells       all macro particles
//                                    = field energy      = kinetic energy
//
//
//
//
//  lambda_0 is the laser wavelength in the laboratory frame L
//  omega and fields are both given in L or the moving frame M
//  the scaled field amplitudes defined above do not depend on the frame
//
//
/////////////////////////////////////////////////////////////////////////////////////////

