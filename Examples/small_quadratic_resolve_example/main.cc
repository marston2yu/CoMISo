/*===========================================================================*\
 *                                                                           *
 *                               CoMISo                                      *
 *      Copyright (C) 2008-2009 by Computer Graphics Group, RWTH Aachen      *
 *                           www.rwth-graphics.de                            *
 *                                                                           *
 *---------------------------------------------------------------------------*
 *  This file is part of CoMISo.                                             *
 *                                                                           *
 *  CoMISo is free software: you can redistribute it and/or modify           *
 *  it under the terms of the GNU General Public License as published by     *
 *  the Free Software Foundation, either version 3 of the License, or        *
 *  (at your option) any later version.                                      *
 *                                                                           *
 *  CoMISo is distributed in the hope that it will be useful,                *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of           *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            *
 *  GNU General Public License for more details.                             *
 *                                                                           *
 *  You should have received a copy of the GNU General Public License        *
 *  along with CoMISo.  If not, see <http://www.gnu.org/licenses/>.          *
 *                                                                           *
\*===========================================================================*/



#include <CoMISo/Utils/StopWatch.hh>
#include <gmm/gmm.h>
#include <vector>
#include <CoMISo/Solver/ConstrainedSolver.hh>
#include <CoMISo/Solver/MISolver.hh>
#include <CoMISo/Solver/GMM_Tools.hh>
#include <CoMISo/Utils/Tools.hh>


/// function to initialize a simple system of linear equations
template<class MatrixT>
void init_les( MatrixT& _A, std::vector< double >& _b)
{
  _A(0,0) = 25  ; _A(0,1) = 0 ; _A(0,2) = -15; _A(0,3) = 0 ;
  _A(1,0) = 0   ; _A(1,1) = 20; _A(1,2) = -8 ; _A(1,3) = -4;
  _A(2,0) = -15 ; _A(2,1) = -8; _A(2,2) = 17 ; _A(2,3) = 0 ;
  _A(3,0) = 0   ; _A(3,1) = -4; _A(3,2) = 0  ; _A(3,3) = 4 ;

  _b[0] = 0; _b[1] = 4; _b[2] = -2; _b[3] = 0;
}

/// function to print the equations corresponding to the matrices of an equation system
template<class MatrixT>
void print_equations( const MatrixT& _B)
{
  size_t m = gmm::mat_nrows( _B);
  size_t n = gmm::mat_ncols( _B);
  for( size_t i = 0; i < m; ++i)
  {
    for( size_t j = 0; j < n-1; ++j)
    {
      if( _B(i,j) != 0.0)
        std::cout << _B(i,j) << "*x" << j;
      else
        std::cout << "   0 ";
      if( j < n-2 ) std::cout << " + ";
    }
    std::cout << " = " << -_B(i, n-1) << std::endl;
  }
}


// Example main
int main(void)
{
  using COMISO_GMM::operator<<;

  std::cout << "---------- 1) Setup small (symmetric) test equation system Ax=b..." << std::endl;
  int n = 4;
  gmm::col_matrix< gmm::wsvector< double > > A(n,n);
  std::vector< double > x(n);
  std::vector< double > b(n);

  std::vector<double> x_bak;

  std::cout << "---------- 1) Set up problem..." << std::endl;

  init_les( A, b);

  // create an empty constraint matrix (will be used later)
  gmm::row_matrix< gmm::wsvector< double > > constraints(0,n+1); //n+1 because of right hand side
  // create an empty vector of variable indices to be rounded (will be used later)
  std::vector< int > ids_to_round;

  std::cout << A << std::endl << b << std::endl;

  // setup constraints
  gmm::resize( constraints, 3, n+1);
  constraints( 0, 0 ) = 1.0;
  constraints( 0, 1 ) = -1.0;
  constraints( 0, n ) = 2.0;
  constraints( 1, 3 ) = 1.0;
  constraints( 1, n ) = -1.0;
  // add one redundant constraint (this will be filtered out during Gaussian elimination)
  constraints( 2, 0 ) = 1.0;
  constraints( 2, 1 ) = -1.0;
  constraints( 2, n ) = 2.0;
  std::cout << "           the constraint equations looks like this:" << std::endl;
  print_equations( constraints);

  std::cout << "---------- 2) Solve full ..." << std::endl;
  COMISO::ConstrainedSolver cs;
  cs.solve_const( constraints, A, x, b, ids_to_round, 0.0, false);
  x_bak = x;

  // first test: resolve with identical rhs's
  std::vector<double> constraint_rhs(3);
  std::vector<double> b_new = b;
  constraint_rhs[0] = -2.0;
  constraint_rhs[1] =  1.0;
  constraint_rhs[2] = -2.0;

  std::cout << "---------- 2) Solve same rhs pre-factorized ..." << std::endl;
  cs.resolve(x, &constraint_rhs, &b_new);
  std::cout << "orig    result:    " << x_bak << std::endl;
  std::cout << "resolve result:    " << x     << std::endl;

  // second test: resolve with changed rhs
  constraint_rhs[0] =  4.0;
  constraint_rhs[1] = -2.0;
  constraint_rhs[2] =  4.0;
  b_new[0] =  1.0;
  b_new[1] = -2.0;
  b_new[2] =  3.0;
  b_new[3] = -5.0;

  std::cout << "---------- 3) Solve different rhs pre-factorized ..." << std::endl;
  cs.resolve(x, &constraint_rhs, &b_new);


  // solve with new factorization
  constraints( 0, n ) = -4.0;
  constraints( 1, n ) =  2.0;
  constraints( 2, n ) = -4.0;
  std::cout << "---------- 4) Solve different rhs full ..." << std::endl;
  cs.solve_const( constraints, A, x_bak, b_new, ids_to_round, 0.0, false);

  std::cout << "orig     result (with different rhs's):    " << x_bak << std::endl;
  std::cout << "resolve  result (with different rhs's):    " << x     << std::endl;

  return 0;
}

