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

/** 
 * Helper function to remove monomials from the basis.
 */
Symbolic SdpRelaxation::applySubstitution(Symbolic monomial) {
	Symbolic originalMonomial;
	bool changed = true;
	while (changed) {
		originalMonomial = monomial;
		for (auto ii = substitutions.begin(); ii != substitutions.end(); ii++) {
      // The fast substitution routine still fails on some rare 
      // conditions. In production environments, it is safer to use
      // the default substitution routine that comes with SymPy.      

			monomial = fastSubstitute(monomial, ii->first, ii->second);
      // monomial = monomial.subst_all(ii->first, ii->second);
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
				temp.push_back(monomial);
			}
		}
		ncMonomials.splice(ncMonomials.end(), temp);
		--degree;
	}
	ncMonomials.push_front(Symbolic(1));
	return unique(ncMonomials);
}

/**
 * Generate the moment matrix of monomials
 * 
 * Arguments:
 * @param monomials - |W_d| set of words of length up to the relaxation order d
 * @param block_index - current block index in the constraints matrices of the
 *                      SDP relaxation
 */
void SdpRelaxation::generateMomentMatrix(const vector<Symbolic> monomials,
		int *blockIndex) {
	Entry entry;
	*blockIndex = 1;
	int nEq = 1;
	// Defining top left corner of momentum matrix
	entry.blockIndex = *blockIndex;
	entry.row = nEq;
	entry.column = nEq;
	entry.value = 1;
	F[0].push_back(entry);
	F[index2linear(0, 0, nMonomials)].push_back(entry);
	++nEq;
	entry.row = nEq;
	entry.column = nEq;
	entry.value = -1;
	F[0].push_back(entry);
	F[index2linear(0, 0, nMonomials)].push_back(entry);
	blockStruct.push_back(-nEq);
  ++(*blockIndex);
  entry.blockIndex = *blockIndex;
	Index index;
	Symbolic monomial, monomial_dagger;
  int k, k_dagger;
  float value;
  // Generating the rest of the matrix.
  // This is potentially done in parallel, albeit the symbolic library used is 
  // not thread-safe.
	#pragma omp parallel default(shared) private(index, entry, k, k_dagger, monomial, monomial_dagger, value)
	{
	#pragma omp for schedule(runtime)
  // We process (u,w) elements of the matrix
	for (int row = 0; row < nMonomials; ++row) {
		for (int column = row; column < nMonomials; ++column) {
      // Calculate the monomial u*v
			monomial = conjugate(monomials[row]) * monomials[column];
      // Apply substitutions if any
			monomial = applySubstitution(monomial);
			if (getCoefficient(monomial) < 0) {
				monomial = -monomial;
			}
      // Look up the index of the monomial in the dictionary built so far
			index = monomialDictionary[monomial];
			#pragma omp critical(update)
			{
        // Have we seen this monomial before?
        if (index.first == 0 && index.second == 0) {
          // If not, we add a new entry to the monomial dictionary.
          index.first = row;
          index.second = column;
          monomialDictionary[monomial] = index;
          k = index2linear(row, column, nMonomials);
        } else {
          // If yes, then we improve sparsity by reusing the 
          // previous variable to denote this entry in the matrix
          k = index2linear(index.first, index.second, nMonomials);
        }
        if (row == column) {
          value = 1;
        } else {
          // Special care must be taken so that the resulting
          // constraint matrices are symmetric, not just 
          // Hermitian. The procedure is essentially the same as 
          // above.          
          value = 0.5;
          monomial_dagger = conjugate(monomials[column]) * monomials[row];
			    monomial_dagger = applySubstitution(monomial_dagger);
          if (getCoefficient(monomial_dagger) < 0) {
            monomial_dagger = -monomial_dagger;
          }
          index = monomialDictionary[monomial_dagger];
          if (index.first == 0 && index.second == 0) {
            index.first = column;
            index.second = row;
            monomialDictionary[monomial] = index;
            k_dagger = index2linear(column, row, nMonomials);
          } else {
            k_dagger = index2linear(index.first, index.second, nMonomials);
          }
          if (k_dagger == k) {
            value = 1;
          } else {
            entry.row = row + 1;
            entry.column = column + 1;
            entry.value = value;
            F[k_dagger].push_back(entry);
          }
        }
        entry.row = row + 1;
        entry.column = column + 1;
        entry.value = value;
        F[k].push_back(entry);
			}
		}
	}
	}
	#pragma omp barrier
  blockStruct.push_back(nMonomials);
  ++(*blockIndex);
}

/* 
 * Calculate the sparse vector representation of a polynomial
 * and pushes it to the F structure.
 */
void SdpRelaxation::pushFacVarSparse(const Symbolic polynomial,
		const int blockIndex, const int i, const int j) {
	list<Symbolic> factors;
  // Preprocess the polynomial for uniform handling later
	if (polynomial.type() == typeid(Sum)) {
		factors = CastPtr<const Sum>(polynomial)->summands;
	} else {
		factors.push_back(polynomial);
	}
	Entry entry;
  // Identify its constituent monomials    
	for (list<Symbolic>::const_iterator monomial = factors.begin();
			monomial != factors.end(); ++monomial) {
		double coeff = getCoefficient(*monomial);
		Symbolic newMonomial = applySubstitution(*monomial / coeff);
		if (getCoefficient(newMonomial) < 0) {
			newMonomial = -1 * newMonomial;
			coeff = -1.0 * coeff;
		}
    // Given the monomial, we need its mapping L_y(w) to push it into
    // a corresponding constraint matrix
		Index index = monomialDictionary[newMonomial];
		entry.blockIndex = blockIndex;
		entry.row = i + 1;
		entry.column = j + 1; 
		entry.value = coeff;
    // k identifies the mapped value of a word (monomial) w
		int k = index2linear(index.first, index.second, nMonomials);
		#pragma omp critical(pushFacVarSparse)
		{
		F[k].push_back(entry);
		}
	}
}

/*
 * Return dense vector representation of a polynomial. This function is
 * nearly identical to push_facvar_sparse, but instead of pushing
 * sparse entries to the constraint matrices, it returns a dense 
 * vector.
 */
double *SdpRelaxation::getFacVar(const Symbolic polynomial) {
	double *facVar = new double[nElements];
	for (int i = 0; i < nElements; ++i) {
		facVar[i] = 0;
	}
	list<Symbolic> factors;
  // Preprocess the polynomial for uniform handling later
	if (polynomial.type() == typeid(Sum)) {
		factors = CastPtr<const Sum>(polynomial)->summands;
	} else {
		factors.push_back(polynomial);
	}
  // Identify its constituent monomials    
	for (list<Symbolic>::const_iterator monomial = factors.begin();
			monomial != factors.end(); ++monomial) {
		double coeff = getCoefficient(*monomial);
		Symbolic newMonomial = applySubstitution(*monomial / coeff);
		if (getCoefficient(newMonomial) < 0) {
			newMonomial = -1 * newMonomial;
			coeff = -1 * coeff;
		}
    // Given the monomial, we need its mapping L_y(w) to find its 
    // location in the dense vector needed by the objective function.
		Index index = monomialDictionary[newMonomial];
		facVar[index2linear(index.first, index.second, nMonomials) - 1] +=
				coeff;
	}
	return facVar;
}

/** 
 * Generate localizing matrices
 *
 * Arguments:
 * @param inequalities - inequality constraints
 * @param monomials - monomials in the set |W_d| with d being the relaxation order
 * @param block_index - the current block index in constraint matrices of the 
 *                      SDP relaxation
 * @param - the order of the relaxation        
 */
void SdpRelaxation::processInequalities(const vector<Symbolic> inequalities,
		const vector<Symbolic> monomials, const int blockIndex, const int order) {
  // Identify the correct set of monomials
	int nIneqMonomials = countNcMonomials(monomials, order - 1);
  // Mark length of block in the constraint matrices
	for (int k = 0; k<inequalities.size(); ++k) {
		blockStruct.push_back(nIneqMonomials);
	}
  // Process M_y(gy)(u,w) entries. Technically this can be done in parallel.
	#pragma omp parallel default(shared)
	{
    #pragma omp for schedule(runtime)
    for (int k = 0; k<inequalities.size(); ++k) {
      int localBlockIndex = blockIndex + k;
      for (int row = 0; row < nIneqMonomials; ++row) {
        for (int column = row; column < nIneqMonomials; ++column) {
          // Calculate the moments of polynomial entries
          Symbolic polynomial = conjugate(monomials[row]) * inequalities[k]
              * monomials[column];
          if (row == column) {
              pushFacVarSparse(polynomial, localBlockIndex, row, column);
          } else {
              // Special care must be taken so that the resulting
              // constraint matrices are symmetric, not just 
              // Hermitian. The procedure is essentially the same as 
              // above.            
              Symbolic polynomial_dagger = conjugate(monomials[column]) * inequalities[k] * monomials[row];
              Symbolic poly = 0.5*polynomial_dagger + 0.5*polynomial;
              pushFacVarSparse(poly, localBlockIndex, row, column);          
          }
        }
      }
    }
	} // End pragma
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

  // Generate the set W_d containing words (monomials) of length up to d,
  // where d is the relaxation order
	vector<Symbolic> monomials = getNcMonomials(variables, order);
  
  // Initialize some helper variables, including the offsets of monomial
  // blocks if there is more than one.
	nMonomials = monomials.size();
	nElements = nMonomials * nMonomials; 
  int blockIndex;
  // Initialize sparse constant matrices in the target SDP
	F = new list<Entry> [nElements + 1];

  // Generate moment matrices for each blocks of variables 
	cout << "Generating moments..." << endl;
  generateMomentMatrix(monomials, &blockIndex);
	
  
  // Objective function needs dense representation
	objFacVar = getFacVar(objective);

  // Equalities are converted to pairs of inequalities
	cout << "Transforming " << equalities.size() << " equalities to "
			<< 2 * equalities.size() << " inequalities..." << endl;
	for (vector<Symbolic>::const_iterator eq = equalities.begin();
			eq != equalities.end(); ++eq) {
		inequalities.push_back(*eq);
		inequalities.push_back(-*eq);
	}
  
  // Process inequalities
	cout << "Processing " << inequalities.size() << " inequalitites..." << endl;
	processInequalities(inequalities, monomials, blockIndex, order);

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
