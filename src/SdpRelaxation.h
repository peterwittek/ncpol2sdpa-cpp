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

#include <iostream>
#include <unordered_map>
#include "symbolicc++.h"
#include "ncUtils.h"

#ifndef SDP_RELAXATION
#define SDP_RELAXATION

typedef pair<int, int> Index;

struct Entry {

	int blockIndex;
	int row;
	int column;
	double value;

};

class SdpRelaxation {

private:
	const unordered_map<Symbolic, Symbolic, hashMonomial> substitutions;
	unordered_map<Symbolic, Index, hashMonomial> monomialDictionary;
	int nMonomials;
	int nElements;
	vector<int> blockStruct;
	double *objFacVar;
	list<Entry> *F;

	Symbolic applySubstitution(Symbolic monomial);
	vector<Symbolic> getNcMonomials(const Symbolic variables, short int degree);
	double *getFacVar(const Symbolic polynomial);
	void generateMomentMatrix(const vector<Symbolic> monomials, int *blockIndex,
			int *nEq);
	void processInequalities(const vector<Symbolic> inequalities,
			const vector<Symbolic> monomials, const int blockIndex, const int order);
	void pushFacVarSparse(const Symbolic polynomial, const int blockIndex,
			const int i, const int j);

public:
	SdpRelaxation(
			const unordered_map<Symbolic, Symbolic, hashMonomial> substitutions);
	~SdpRelaxation();
	void getRelaxation(const Symbolic variables, const Symbolic objective,
			vector<Symbolic> inequalities, const vector<Symbolic> equalities,
			const short int order);
	void writeToSdpa(const char *filename);
};

#endif
