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


//=============================================================================
//
//  CLASS CholmodSolver
//
//=============================================================================

#ifndef COMISO_CHOLMOD_SOLVER_HH
#define COMISO_CHOLMOD_SOLVER_HH

//== COMPILE-TIME PACKAGE REQUIREMENTS ========================================
#include <CoMISo/Config/config.hh>
#if COMISO_SUITESPARSE_AVAILABLE

//== INCLUDES =================================================================


#include <CoMISo/Config/CoMISoDefines.hh>
#include <Base/Utils/StopWatch.hh>

#include <iostream>
#include <vector>

#include "cholmod.h"
//#include "cholmod_matrixops.h"

// typedef struct cholmod_common_struct cholmod_common;
// typedef struct cholmod_factor_struct cholmod_factor;


//== NAMESPACES ===============================================================

namespace COMISO {

//== CLASS DEFINITION =========================================================
class COMISODLLEXPORT CholmodSolver
{
public:

    // _size is maximal size this instance can handle (smaller problems are possible!!!)
    CholmodSolver();
    ~CholmodSolver();

    bool calc_system( const std::vector<int>&    _colptr,
		      const std::vector<int>&    _rowind,
		      const std::vector<double>& _values );

    bool calc_system_prepare_pattern( const std::vector<int>&    _colptr,
                                      const std::vector<int>&    _rowind,
                                      const std::vector<double>& _values,
                                      const std::vector<int>&    _colptr2,
                                      const std::vector<int>&    _rowind2,
                                      const std::vector<double>& _values2 );

#if COMISO_GMM_AVAILABLE
    template< class GMM_MatrixT>
    bool calc_system_gmm( const GMM_MatrixT& _mat);
#endif

    template< class Eigen_MatrixT>
    bool calc_system_eigen( const Eigen_MatrixT& _mat);

    // factorize matrix _mat by optimizing the permutation for the pattern of _mat_pattern
    template< class Eigen_MatrixT>
    bool calc_system_eigen_prepare_pattern( const Eigen_MatrixT& _mat, const Eigen_MatrixT& _mat_pattern);


    bool update_system( const std::vector<int>&    _colptr,
 			const std::vector<int>&    _rowind,
 			const std::vector<double>& _values );

#if COMISO_GMM_AVAILABLE
    template< class GMM_MatrixT>
    bool update_system_gmm(const GMM_MatrixT& _mat);
#endif

    template< class Eigen_MatrixT>
    bool update_system_eigen( const Eigen_MatrixT& _mat);


    bool update_downdate_factor( const std::vector<int>&    _colptr,
                                 const std::vector<int>&    _rowind,
                                 const std::vector<double>& _values,
                                 const bool                 _upd);

    template< class Eigen_MatrixT>
    bool update_downdate_factor_eigen( const Eigen_MatrixT& _mat, const bool _upd = true);

    bool solve ( double *             _x0, double *             _b);

    bool solve ( std::vector<double>& _x0, std::vector<double>& _b);

    template <class Eigen_VectorT>
    bool solve(Eigen_VectorT& _x, const Eigen_VectorT& _b);

    bool& show_timings();

    int dimension();

private:

    cholmod_common * mp_cholmodCommon;

    cholmod_factor * mp_L;

    std::vector<double> values_;
    std::vector<int>    colptr_;
    std::vector<int>    rowind_;

    bool show_timings_;
    Base::StopWatch sw_;
};

//=============================================================================
} // namespace COMISO
//=============================================================================
#if defined(INCLUDE_TEMPLATES) && !defined(COMISO_CHOLMOD_SOLVER_TEMPLATES_C)
#define COMISO_CHOLMOD_SOLVER_TEMPLATES
#include "CholmodSolverT_impl.hh"
#endif
//=============================================================================
#else  // COMISO_SUITESPARSE_AVAILABLE
#include<CoMISo/Solver/EigenLDLTSolver.hh>

//#warning "CholmodSolver not available, fallback to EigenLDLTSolver..."

namespace COMISO
{
  typedef EigenLDLTSolver CholmodSolver;
} // namespace COMISO

#endif // COMISO_SUITESPARSE_AVAILABLE
//=============================================================================
#endif // COMISO_CHOLMOD_SOLVER_HH defined
//=============================================================================
