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

#include <unordered_map>
#include "ncUtils.h"

Symbolic conjugate(const Symbolic monomial) {
	Symbolic result = Symbolic(1);
	if (monomial.type() == typeid(Product)) {
		list<Symbolic> factors = CastPtr<const Product>(monomial)->factors;
		list<Symbolic>::const_iterator i = factors.end();
		for (i--; i != factors.end(); --i) {
			result *= *i;
		}
	} else if ((monomial.type() == typeid(Power))
			|| (monomial.type() == typeid(Symbol))
			|| (monomial.type() == typeid(Numeric))) {
		result = monomial;
	} else {
		cerr << "Not a monomial: " << monomial << endl;
	}
	return result;
}

int ncDegree(const Symbolic monomial) {
	int degree = 0;
	if (monomial.type() == typeid(Product)) {
		list<Symbolic> factors = CastPtr<const Product>(monomial)->factors;
		for (list<Symbolic>::const_iterator i = factors.begin();
				i != factors.end(); ++i) {
			if ((*i).type() == typeid(Power)) {
				CastPtr<const Power> p = *i;
				degree += (int) p->parameters.back();
			} else if ((*i).type() == typeid(Symbol)) {
				++degree;
			}
		}
	} else if (monomial.type() == typeid(Power)) {
		CastPtr<const Power> p = monomial;
		degree += (int) p->parameters.back();
	} else if (monomial.type() == typeid(Symbol)) {
		degree = 1;

	} else if (monomial.type() == typeid(Numeric)) {

	} else {
		cerr << "Not a monomial: " << monomial << endl;
	}
	return degree;
}

int countNcMonomials(const vector<Symbolic> monomials, const short int degree) {
	int count = 0;
	for (vector<Symbolic>::const_iterator i = monomials.begin();
			i != monomials.end(); ++i) {
		if (ncDegree(*i) <= degree) {
			++count;
		} else {
			break;
		}
	}
	return count;
}

vector<Symbolic> unique(const list<Symbolic> l) {
	vector<Symbolic> result;
	unordered_map<Symbolic, int, hashMonomial> seen;
	for (list<Symbolic>::const_iterator i = l.begin(); i != l.end(); ++i) {
		if (seen.count(*i) > 0) {
			continue;
		}
		seen[*i] = 1;
		result.push_back(*i);
	}
	return result;
}

int index2linear(const int i, const int j, const int nMonomials) {
	if (i == 0) {
		return j + 1;
	}
	int skew = i * (i + 1) / 2;
	return i * nMonomials - skew + j + 1;
}

double getCoefficient(const Symbolic monomial) {
	double result = 1.0;
	if (monomial.type() == typeid(Product)) {
		Symbolic front = CastPtr<const Product>(monomial)->factors.front();
		if (front.type() == typeid(Numeric)) {
			result = front;
		}
	} else if (monomial.type() == typeid(Numeric)) {
		result = monomial;
	}
	return result;
}
