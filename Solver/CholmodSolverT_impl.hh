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


//== COMPILE-TIME PACKAGE REQUIREMENTS ========================================
#include <CoMISo/Config/config.hh>
#if COMISO_SUITESPARSE_AVAILABLE
//=============================================================================


#define COMISO_CHOLMOD_SOLVER_TEMPLATES_C

#include <CoMISo/Solver/GMM_Tools.hh>
#include <CoMISo/Solver/Eigen_Tools.hh>
#include <CoMISo/Solver/CholmodSolver.hh>


namespace COMISO {

#if COMISO_GMM_AVAILABLE

template< class GMM_MatrixT>
bool CholmodSolver::calc_system_gmm( const GMM_MatrixT& _mat)
{
//   std::vector<int>    colptr;
//   std::vector<int>    rowind;
//   std::vector<double> values;


    if(show_timings_) sw_.start();

    COMISO_GMM::get_ccs_symmetric_data( _mat,
					 'u',
					 values_,
					 rowind_,
					 colptr_ );

    if(show_timings_)
    {
      std::cerr << "Cholmod Timing GMM convert: " << sw_.stop()/1000.0 << "s\n";
      std::cerr << "#nnz: " << values_.size() << std::endl;
    }

    return calc_system( colptr_, rowind_, values_);
}


//-----------------------------------------------------------------------------


template< class GMM_MatrixT>
bool CholmodSolver::update_system_gmm( const GMM_MatrixT& _mat)
{
//   std::vector<int>    colptr;
//   std::vector<int>    rowind;
//   std::vector<double> values;

  COMISO_GMM::get_ccs_symmetric_data( _mat,
				      'u',
				       values_,
				       rowind_,
				       colptr_ );

  return update_system( colptr_, rowind_, values_);
}

#endif // COMISO_GMM_AVAILABLE

//-----------------------------------------------------------------------------

template< class Eigen_MatrixT>
bool CholmodSolver::calc_system_eigen( const Eigen_MatrixT& _mat)
{
    if(show_timings_) sw_.start();

#if COMISO_EIGEN3_AVAILABLE
    COMISO_EIGEN::get_ccs_symmetric_data( _mat,
					 'u',
					 values_,
					 rowind_,
					 colptr_ );
#endif

    if(show_timings_)
    {
      std::cerr << "Cholmod Timing EIGEN convert: " << sw_.stop()/1000.0 << "s\n";
      std::cerr << "#nnz: " << values_.size() << std::endl;
    }

    return calc_system( colptr_, rowind_, values_);
}


//-----------------------------------------------------------------------------


template< class Eigen_MatrixT>
bool CholmodSolver::calc_system_eigen_prepare_pattern( const Eigen_MatrixT& _mat, const Eigen_MatrixT& _mat_pattern)
{
    if(show_timings_) sw_.start();

#if COMISO_EIGEN3_AVAILABLE
    COMISO_EIGEN::get_ccs_symmetric_data( _mat,
                                         'u',
                                         values_,
                                         rowind_,
                                         colptr_ );
#endif

    std::vector<double> values2;
    std::vector<int>    colptr2;
    std::vector<int>    rowind2;

#if COMISO_EIGEN3_AVAILABLE
    COMISO_EIGEN::get_ccs_symmetric_data( _mat_pattern,
                                         'u',
                                         values2,
                                         rowind2,
                                         colptr2 );
#endif

    if(show_timings_)
    {
      std::cerr << "Cholmod Timing EIGEN convert: " << sw_.stop()/1000.0 << "s\n";
      std::cerr << "#nnz: " << values_.size() << std::endl;
    }

    return calc_system_prepare_pattern( colptr_, rowind_, values_, colptr2, rowind2, values2);
}


//-----------------------------------------------------------------------------

template< class Eigen_MatrixT>
bool CholmodSolver::update_system_eigen( const Eigen_MatrixT& _mat)
{
#if COMISO_EIGEN3_AVAILABLE
  COMISO_EIGEN::get_ccs_symmetric_data( _mat,
				      'u',
				       values_,
				       rowind_,
				       colptr_ );
#endif
  return update_system( colptr_, rowind_, values_);
}

//-----------------------------------------------------------------------------


template< class Eigen_MatrixT>
bool CholmodSolver::update_downdate_factor_eigen( const Eigen_MatrixT& _mat, const bool _upd)
{
    if(show_timings_) sw_.start();

#if COMISO_EIGEN3_AVAILABLE
    COMISO_EIGEN::get_ccs_symmetric_data( _mat,
                                         'c',
                                         values_,
                                         rowind_,
                                         colptr_ );
#endif

    if(show_timings_)
    {
      std::cerr << "Cholmod Timing EIGEN convert: " << sw_.stop()/1000.0 << "s\n";
      std::cerr << "#nnz: " << values_.size() << std::endl;
    }

    return update_downdate_factor( colptr_, rowind_, values_, _upd);
}

//-----------------------------------------------------------------------------


template <class Eigen_VectorT>
bool CholmodSolver::solve(Eigen_VectorT& _x, const Eigen_VectorT& _b)
{
  // hopefully, cholmod_solve does not change b
  return solve(_x.data(), const_cast<Eigen_VectorT&>(_b).data());
}

}

//=============================================================================
#endif // COMISO_SUITESPARSE_AVAILABLE
//=============================================================================
