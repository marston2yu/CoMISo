//=============================================================================
//
//  CLASS LinearConstraintHandlerElimination
//
//=============================================================================


#ifndef COMISO_LINEARCONSTRAINTHANDLERPENALTY_HH
#define COMISO_LINEARCONSTRAINTHANDLERPENALTY_HH


//== INCLUDES =================================================================

#include <CoMISo/Config/CoMISoDefines.hh>
#include <CoMISo/Utils/gmm.hh>
#include <iostream>

//== COMPILE-TIME PACKAGE REQUIREMENTS ========================================
#include <CoMISo/Config/config.hh>
#if COMISO_GMM_AVAILABLE

//== FORWARDDECLARATIONS ======================================================

//== NAMESPACES ===============================================================

namespace COMISO {

//== CLASS DEFINITION =========================================================



	      
/** \class LinearConstraintHandler LinearConstraintHandler.hh <COMISO/.../LinearConstraintHandler.hh>

    Brief Description.
  
    A more elaborate description follows.
*/


class COMISODLLEXPORT LinearConstraintHandlerPenalty
{
public:
   
  typedef gmm::col_matrix< gmm::wsvector< double> > CMatrix;
  typedef gmm::row_matrix< gmm::wsvector< double> > RMatrix;

  // use c-arrays as vectors for gmm
  typedef gmm::array1D_reference<double*> VectorPT;

  /// Constructor
  LinearConstraintHandlerPenalty();

  // initialize Constructor
  template<class MatrixT, class VectorT>
  LinearConstraintHandlerPenalty( const MatrixT& _C, const VectorT& _c);
 
  // penalty weight
  double& penalty();

  // number of variables
  int n();

  // number of linearly independent constraints (n-n_reduced)
  int n_constraints();

  // initialize new constraints
  template<class MatrixT, class VectorT>
  void initialize( const MatrixT& _C, const VectorT& _c);

  // initialize new constraints rhs only
  void initialize( const std::vector<double>& _c);
  void initialize( double* _c);

  // transform energy
  double add_penalty_f( double* _x, const double _f);

  // transform gradient
  void add_penalty_gradient( const std::vector<double>& _x, std::vector<double>& _g);
  void add_penalty_gradient( double* _x, double* _g);

  // transform hessian
  void add_penalty_hessian( RMatrix& _H);

private:
  // penalty weight
  double penalty_;
  
  // Linear Constraints C_*x_ = b_
  RMatrix             C_;
  std::vector<double> b_;

  // precomputed penalty terms
  RMatrix penalty_H_;
  std::vector<double> penalty_grad_b_;

  // temp vector
  std::vector<double> temp_;

  // number of variables
  int n_;
  // number of constraints (linear independent)
  int m_;
};

//=============================================================================
} // namespace COMISO
//=============================================================================
#if defined(INCLUDE_TEMPLATES) && !defined(COMISO_LINEARCONSTRAINTHANDLERPENALTY_C)
#define COMISO_LINEARCONSTRAINTHANDLERPENALTY_TEMPLATES
#include "LinearConstraintHandlerPenaltyT_impl.hh"
#endif

//=============================================================================
#endif // COMISO_GMM_AVAILABLE
//=============================================================================
#endif // COMISO_LINEARCONSTRAINTHANDLERPENALTY_HH defined
//=============================================================================

