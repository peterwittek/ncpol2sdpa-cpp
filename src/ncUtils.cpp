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

/**
 * A simple routine of conjugating a monomial of Hermitian variables
 */
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

/**
 * Experimental fast substitution routine that considers only restricted 
 * cases of noncommutative algebras. In rare cases, it fails to find a
 * substitution. Use it with proper testing.
 *
 * Arguments:
 * @param monomial - the monomial with parts need to be substituted
 * @param oldSub - the part to be replaced
 * @param newSub - the replacement
 */
Symbolic fastSubstitute(Symbolic monomial, Symbolic oldSub, Symbolic newSub) {
	bool isOldSubProduct = false;
	list<Symbolic> oldSubFactors;
	if (oldSub.type() == typeid(Product)) {
		isOldSubProduct = true;
		oldSubFactors = CastPtr<const Product>(oldSub)->factors;
	}
	if (monomial.type() == typeid(Product)) {
		bool changed = false;
		list<Symbolic> result;
		Symbolic remainder;
		list<Symbolic> factors = CastPtr<const Product>(monomial)->factors;
		list<Symbolic>::const_iterator factor = factors.begin();
		while ( factor != factors.end()) {
			if (!isOldSubProduct){
				if (*factor == oldSub){
					result.push_back(newSub);
					++factor;
					changed = true;
					break;
				} else {
					result.push_back(*factor);
				}
			} else {
				remainder = 1;
				list<Symbolic>::const_iterator it1 = factor;
				list<Symbolic>::const_iterator it2 = oldSubFactors.begin();
				bool match = false;
			    while (true) {
			    	if (it1!=factors.end() && it2!=oldSubFactors.end()) {
			    		if ((it1->type() == typeid(Symbol) && it2->type() == typeid(Symbol))){
			    			if (*it1!=*it2) {
			    				break;
			    			}
			    		} else if (it1->type() == typeid(Power)) {
							CastPtr<const Power> p = *it1;
							Symbolic x = (Symbolic) p->parameters.front();
							int degree1 = (int) p->parameters.back();
							Symbolic y; int degree2=1;
							if (it2->type() == typeid(Symbol)) {
								y = *it2;
							} else {
								CastPtr<const Power> p2 = *it2;
								y = (Symbolic) p2->parameters.front();
								degree2 = (int) p2->parameters.back();
							}
							if ( x != y) {
								break;
							}
							if (degree2 > degree1) {
								break;
							} else if (degree2 < degree1) {
								++it2;
								if (it2!=oldSubFactors.end()) {
									break;
									--it2;
								} else {
									remainder = x^(degree1-1);
									++it1;
									match = true;
									break;
								}
							}
			    		} else {
			    			break;
			    		}
			    	}
			        if (it2==oldSubFactors.end()) {
			        	match = true;
			        	break;
			        }
			        if (it1==factors.end()) {
			        	break;
			        }
			    	++it1; ++it2;
			    }
			    if (match) {
			    	changed = true;
			    	result.push_back(newSub);
			    	if (remainder != 1) {
			    		result.push_back(remainder);
			    	}
					factor=it1;
			    	break;
			    } else {
			    	result.push_back(*factor);
			    }
			}
			++factor;
		}
		if (changed) {
			Product newMonomial;
			for(list<Symbolic>::const_iterator i=result.begin();i!=result.end();++i)
			  newMonomial.factors.push_back(*i);
			while ( factor != factors.end()) {
				newMonomial.factors.push_back(*factor);
				++factor;
			}
			return newMonomial;
		} else {
			return monomial;
		}
	} else if (monomial==oldSub){
		return newSub;
	} else {
		return monomial;
	}
}

/**
 * Returns the degree of a noncommutative polynomial.
 */ 
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

/**
 * Given a list of monomials, it counts those that have a certain degree,
 * or less. The function is useful when certain monomials were eliminated
 * from the basis.
 * 
 * Arguments:
 * @param variables - the noncommutative variables making up the monomials
 * @param degree -- maximum degree to count

 * Returns the count of appropriate monomials.
 */
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

/**
 * Helper function to include only unique monomials in a basis.
 */ 
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
  return i * nMonomials + j + 1;
}

/**
 * Helper function to get the coefficient of a monomial.
 */ 
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
