/**
 * A converter from noncommutative polynomial optimization problems
 * to sparse SDPA input format
 *
 * Copyright (C) 2013 Peter Wittek
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
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

/*
 * An example that exports to sparse SDPA format for scalable computation. The
 * description of the example is in the following paper:
 * Pironio, S.; Navascués, M. & Acín, A. Convergent relaxations of
 * polynomial optimization problems with noncommuting variables SIAM Journal
 * on Optimization, SIAM, 2010, 20, 2157-2180.
 *
 */

#include "SdpRelaxation.h"

int main(void) {
	short int nVars = 2;
	short int order = 2;
    char filename[] = "examplenc.dat-s";

    // Declaring noncommutative variables
	Symbolic X("X", nVars);
	X = ~X;

	// Setting objective function
	Symbolic objective = X(0) * X(1) + X(1) * X(0);

	// Defining inequalities
	vector<Symbolic> inequalities;
	inequalities.push_back(-X(1)*X(1) + X(1) + 0.5);

	// Defining equalities
	vector<Symbolic> equalities;

	// Defining monomial substitutions
	unordered_map<Symbolic, Symbolic, hashMonomial> substitutions;
	substitutions[X(0)*X(0)] = X(0);

	// Obtaining relaxation and writing file
    SdpRelaxation *sdpRelaxation = new SdpRelaxation(substitutions);
	sdpRelaxation->getRelaxation(X, objective, inequalities, equalities, order);
    sdpRelaxation->writeToSdpa(filename);
    delete sdpRelaxation;
  
	return 0;
}
