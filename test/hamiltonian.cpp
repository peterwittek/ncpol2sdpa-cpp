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
 * Exporting a Hamiltonian ground state problem to SDPA. The Hamiltonian
 * is described in the following paper:
 * Corboz, P.; Evenbly, G.; Verstraete, F. & Vidal, G. (2009),
 * Simulation of interacting fermions with entanglement renormalization.
 * arXiv:0904.4151
 *
 */

#include <sys/time.h>
#include "SdpRelaxation.h"

Index linearToLattice(short int r, short int latticeDimensions) {
	Index index;
	index.first = r % latticeDimensions;
	index.second = r / latticeDimensions;
	return index;
}

vector<Symbolic> getNeighbors(Symbolic X, short int r,
		short int latticeDimensions) {
	vector<Symbolic> neighbors;
	Index index = linearToLattice(r, latticeDimensions);
	int j = index.first;
	int i = index.second;
	if (i > 0) {
		neighbors.push_back(X(r - latticeDimensions));
	}
	if (i < latticeDimensions - 1) {
		neighbors.push_back(X(r + latticeDimensions));
	}
	if (j > 0) {
		neighbors.push_back(X(r - 1));
	}
	if (j < latticeDimensions - 1) {
		neighbors.push_back(X(r + 1));
	}

	return neighbors;
}

int main(void) {
	short int latticeDimensions = 4;
	short int nVars = latticeDimensions * latticeDimensions;
	short int order = 2;

	double gamma = 1.0;
	double lambda = 2.0;

	char filename[] = "hamiltonian.dat-s";

	Symbolic C("C", nVars);

	C = ~C;
	struct timeval start, end;

	Symbolic hamiltonian = 0;
	for (short int r = 0; r < nVars; ++r) {
		hamiltonian -= 2 * lambda * C(r) * C(r);
		vector<Symbolic> neighbors = getNeighbors(C, r, latticeDimensions);
		for (vector<Symbolic>::const_iterator Cs = neighbors.begin();
				Cs != neighbors.end(); ++Cs) {
			hamiltonian += C(r) * (*Cs) + (*Cs) * C(r);
			hamiltonian -= gamma * (C(r) * (*Cs) + (*Cs) * C(r));
		}
	}
	unordered_map<Symbolic, Symbolic, hashMonomial> substitutions;
	vector<Symbolic> inequalities;
	vector<Symbolic> equalities;

	for (short int r = 0; r < nVars; ++r) {
		for (short int s = r; s < nVars; ++s) {
			if (r != s) {
//				equalities.push_back(C(r) * C(s) + C(s) * C(r));
				// Monomial substitutions result in much sparser
				// SDPs, but they are slower to generate
				 substitutions[C(r)*C(s)]=-C(s)*C(r);
			} else {
				equalities.push_back(C(r) * C(s) + C(s) * C(r) - 1);
			}
		}
	}

	SdpRelaxation *sdpRelaxation = new SdpRelaxation(substitutions);

	gettimeofday(&start, NULL);
	sdpRelaxation->getRelaxation(C, hamiltonian, inequalities, equalities,
			order);
	sdpRelaxation->writeToSdpa(filename);
	gettimeofday(&end, NULL);
	cout << latticeDimensions << " " << end.tv_sec - start.tv_sec << " s"
			<< endl;
	delete sdpRelaxation;

	return 0;
}
