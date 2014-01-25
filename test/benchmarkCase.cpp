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

#include <sys/time.h>
#include "SdpRelaxation.h"

int main(void) {
	short int nVars = 10;
	short int order = 1;
    char filename[] = "benchmark.dat-s";

    // Declaring noncommutative variables
	Symbolic X("X", nVars);
	X = ~X;

	// Setting objective function
	Symbolic objective = 0;
  for (int i=0; i < nVars; ++i) {
    for (int j=0; j < nVars; ++j) {
      objective += X(i)*X(j);
    }
  }

	// Defining inequalities
	vector<Symbolic> inequalities;
  for (int i=1; i < nVars; ++i) {
    inequalities.push_back( X(i)*X(i-1) - 0.5);
  }

	// Defining equalities
	vector<Symbolic> equalities;

	// Defining monomial substitutions
	unordered_map<Symbolic, Symbolic, hashMonomial> substitutions;
  for (int i=0; i < nVars; ++i) {
    substitutions[X(i)*X(i)] = X(i);
  }

	// Obtaining relaxation and writing file
  struct timeval start, end;
  gettimeofday(&start, NULL);
  SdpRelaxation *sdpRelaxation = new SdpRelaxation(substitutions);
  sdpRelaxation->getRelaxation(X, objective, inequalities, equalities, order);
  sdpRelaxation->writeToSdpa(filename);
  gettimeofday(&end, NULL);
  cout << nVars << " " << end.tv_sec - start.tv_sec << " s"
    << endl;
  delete sdpRelaxation;
  
	return 0;
}
