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

#include "symbolicc++.h"

#ifndef NC_UTILS
#define NC_UTILS

struct hashMonomial {
	size_t operator()(const Symbolic &monomial) const {
		stringstream ss;
		ss << monomial;
		return hash<string>()(ss.str());
	}
};

Symbolic conjugate(const Symbolic monomial);
int countNcMonomials(const vector<Symbolic> monomials, const short int degree);
Symbolic fastSubstitute(Symbolic monomial, Symbolic oldSub, Symbolic newSub);
double getCoefficient(const Symbolic monomial);
int index2linear(const int i, const int j, const int nMonomials);
int ncDegree(const Symbolic monomial);
vector<Symbolic> unique(const list<Symbolic> l);

#endif
