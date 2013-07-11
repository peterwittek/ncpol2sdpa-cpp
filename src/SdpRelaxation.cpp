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

#include <fstream>
#include "SdpRelaxation.h"

using namespace std;

SdpRelaxation::SdpRelaxation(
		const unordered_map<Symbolic, Symbolic, hashMonomial> substitutions) :
		substitutions(substitutions) {
}

SdpRelaxation::~SdpRelaxation() {
	delete[] F;
}

Symbolic SdpRelaxation::applySubstitution(Symbolic monomial) {
	Symbolic originalMonomial;
	bool changed = true;
	while (changed) {
		originalMonomial = monomial;
		for (auto ii = substitutions.begin(); ii != substitutions.end(); ii++) {
			monomial = monomial.subst_all(ii->first, ii->second);
		}
		if (originalMonomial == monomial) {
			changed = false;
		}
	}
	return monomial;
}

vector<Symbolic> SdpRelaxation::getNcMonomials(const Symbolic variables,
		short int degree) {
	list<Symbolic> ncMonomials;
	short int nVars = variables.rows();
	if (degree > 0) {
		for (short int i = 0; i < nVars; ++i) {
			ncMonomials.push_back(variables(i));
		}
		--degree;
	}
	while (degree > 0) {
		list<Symbolic> temp;
		for (short int i = 0; i < nVars; ++i) {
			for (list<Symbolic>::const_iterator newVar = ncMonomials.begin();
					newVar != ncMonomials.end(); ++newVar) {
				Symbolic monomial = variables(i) * (*newVar);
				if (substitutions.count(monomial) == 0) {
					temp.push_back(monomial);
				}
			}
		}
		ncMonomials.splice(ncMonomials.end(), temp);
		--degree;
	}
	ncMonomials.push_front(Symbolic(1));
	return unique(ncMonomials);
}

void SdpRelaxation::generateMomentMatrix(const vector<Symbolic> monomials,
		int *blockIndex, int *nEq) {
	Entry entry;
	*blockIndex = 1;
	*nEq = 1;
	// Defining top left corner of momentum matrix
	entry.blockIndex = *blockIndex;
	entry.row = *nEq;
	entry.column = *nEq;
	entry.value = 1;
	F[0].push_back(entry);
	F[index2linear(0, 0, nMonomials)].push_back(entry);
	++(*nEq);
	entry.row = *nEq;
	entry.column = *nEq;
	entry.value = -1;
	F[0].push_back(entry);
	F[index2linear(0, 0, nMonomials)].push_back(entry);
	++(*nEq);
	// Generating the rest of the matrix
	Index index;
	Symbolic monomial;
	for (int i = 0; i < nMonomials; ++i) {
		for (int j = i; j < nMonomials; ++j) {
			monomial = conjugate(monomials[i]) * monomials[j];
			monomial = applySubstitution(monomial);
			if (getCoefficient(monomial) < 0) {
				monomial = -monomial;
			}
			index = monomialDictionary[monomial];
			if (index.first == 0 && index.second == 0) {
				index.first = i;
				index.second = j;
				monomialDictionary[monomial] = index;
			} else {
				entry.row = *nEq;
				entry.column = *nEq;
				entry.value = 1;
				F[index2linear(index.first, index.second, nMonomials)].push_back(
						entry);
				entry.value = -1;
				F[index2linear(i, j, nMonomials)].push_back(entry);
				++(*nEq);
				entry.row = *nEq;
				entry.column = *nEq;
				entry.value = -1;
				F[index2linear(index.first, index.second, nMonomials)].push_back(
						entry);
				entry.value = 1;
				F[index2linear(i, j, nMonomials)].push_back(entry);
				++(*nEq);
			}
		}
	}
	--(*nEq);
}

double *SdpRelaxation::getFacVar(const Symbolic polynomial) {
	double *facVar = new double[nElements];
	for (int i = 0; i < nElements; ++i) {
		facVar[i] = 0;
	}
	list<Symbolic> factors;
	if (polynomial.type() == typeid(Sum)) {
		factors = CastPtr<const Sum>(polynomial)->summands;
	} else {
		factors.push_back(polynomial);
	}
	for (list<Symbolic>::const_iterator monomial = factors.begin();
			monomial != factors.end(); ++monomial) {
		double coeff = getCoefficient(*monomial);
		Symbolic newMonomial = applySubstitution(*monomial / coeff);
		if (getCoefficient(newMonomial) < 0) {
			newMonomial = -1 * newMonomial;
			coeff = -1 * coeff;
		}
		Index index = monomialDictionary[newMonomial];
		facVar[index2linear(index.first, index.second, nMonomials) - 1] +=
				coeff;
	}
	return facVar;
}

void SdpRelaxation::pushFacVarSparse(const Symbolic polynomial,
		const int blockIndex, const int i, const int j) {
	list<Symbolic> factors;
	if (polynomial.type() == typeid(Sum)) {
		factors = CastPtr<const Sum>(polynomial)->summands;
	} else {
		factors.push_back(polynomial);
	}
	Entry entry;
	for (list<Symbolic>::const_iterator monomial = factors.begin();
			monomial != factors.end(); ++monomial) {
		double coeff = getCoefficient(*monomial);
		Symbolic newMonomial = applySubstitution(*monomial / coeff);
		if (getCoefficient(newMonomial) < 0) {
			newMonomial = -1 * newMonomial;
			coeff = -1 * coeff;
		}
		Index index = monomialDictionary[newMonomial];
		entry.blockIndex = blockIndex;
		entry.row = i + 1;
		entry.column = j + 1; //index.second+1;
		entry.value = coeff;
		int k = index2linear(index.first, index.second, nMonomials);
		F[k].push_back(entry);
	}
}

void SdpRelaxation::processInequalities(const vector<Symbolic> inequalities,
		const vector<Symbolic> monomials, int *blockIndex, const int order) {
	int nIneqMonomials = countNcMonomials(monomials, order - 1);
	for (auto ineq = inequalities.begin(); ineq != inequalities.end(); ++ineq) {
		blockStruct.push_back(nIneqMonomials);
		++(*blockIndex);
		for (int i = 0; i < nIneqMonomials; ++i) {
			for (int j = i; j < nIneqMonomials; ++j) {
				Symbolic polynomial = conjugate(monomials[i]) * (*ineq)
						* monomials[j];
				pushFacVarSparse(polynomial, *blockIndex, i, j);
			}
		}
	}
}

/** Obtain SDP relaxation
 * @param variables - the noncommutative variables
 * @param objective - the objective function to minimize
 * @param inequalities - the list of inequality constraints
 * @param equalities - the list of equality constraints
 * @param order - the order of the relaxation
 */
void SdpRelaxation::getRelaxation(const Symbolic variables,
		const Symbolic objective, vector<Symbolic> inequalities,
		const vector<Symbolic> equalities, const short int order) {

	vector<Symbolic> monomials = getNcMonomials(variables, order);
	nMonomials = monomials.size();
	nElements = nMonomials * (nMonomials + 1) / 2;

	F = new list<Entry> [nElements + 1];
	cout << "Generating moments..." << endl;
	int blockIndex, nEq;
	generateMomentMatrix(monomials, &blockIndex, &nEq);
	blockStruct.push_back(-nEq);

	// Objective function needs dense representation
	objFacVar = getFacVar(objective);

	cout << "Transforming " << equalities.size() << " equalities to "
			<< 2 * equalities.size() << " inequalities..." << endl;
	for (vector<Symbolic>::const_iterator eq = equalities.begin();
			eq != equalities.end(); ++eq) {
		inequalities.push_back(*eq);
		inequalities.push_back(-*eq);
	}

	cout << "Processing " << inequalities.size() << " inequalitites..." << endl;
	processInequalities(inequalities, monomials, &blockIndex, order);
	// Finally add the SDP constraint on the moment matrix
	++blockIndex;
	blockStruct.push_back(nMonomials);
	Entry entry;
	for (int i = 0; i < nMonomials; ++i) {
		for (int j = i; j < nMonomials; ++j) {
			int k = index2linear(i, j, nMonomials);
			entry.blockIndex = blockIndex;
			entry.row = i + 1;
			entry.column = j + 1;
			entry.value = 1;
			F[k].push_back(entry);
		}
	}
}

/** Write an SDP relaxation to SDPA format
 * @param filename - the name of the file
 */
void SdpRelaxation::writeToSdpa(const char *filename) {
	ofstream outfile(filename);

	// Writing header
	outfile << "\"file " << filename << " generated by ncpol2sdpa\"\n";
	cout << "writing problem in " << filename << endl;
	outfile << nElements << " = number of vars\n";
	outfile << blockStruct.size() << " = number of blocs\n";
	outfile << "(";
	for (unsigned int i = 0; i < blockStruct.size(); ++i) {
		outfile << blockStruct[i];
		if (i != blockStruct.size() - 1) {
			outfile << ", ";
		} else {
			outfile << ") = BlocStructure\n";
		}
	}
	// Objective function
	outfile << "{";
	for (int i = 0; i < nElements; ++i) {
		outfile << objFacVar[i];
		if (i != nElements - 1) {
			outfile << ", ";
		} else {
			outfile << "}\n";
		}

	}
	// Writing entries
	for (int k = 0; k < nElements + 1; ++k) {
		for (list<Entry>::const_iterator e = F[k].begin(); e != F[k].end();
				++e) {
			outfile << k << "\t" << e->blockIndex << "\t" << e->row << "\t"
					<< e->column << "\t" << e->value << "\n";
		}
	}
	outfile.close();
}
