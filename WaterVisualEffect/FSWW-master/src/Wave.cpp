/* 
 * File: Wave.cpp
 *
 * Copyright (C) 2019  Camille Schreck
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "Wave.hpp"
#include "definitions.hpp"

Wave::Wave() {
  phase_exp = 1;
  phase = 0;
}

bool Wave::isActive() const {
  return is_active;
}

void Wave::setActive(int t) const {
  is_active = true;
}

void Wave::setPhase(FLOAT p) {
  phase_exp = exp(definitions::i_*p);
  phase = p;
}
COMPLEX Wave::getPhaseExp() const {
  return phase_exp;
}

FLOAT Wave::getPhase() const {
  return phase;
}
