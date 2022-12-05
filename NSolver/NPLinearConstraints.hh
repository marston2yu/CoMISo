//=============================================================================
//
//  CLASS NPLinearConstraints
//
//=============================================================================


#ifndef COMISO_NPLINEARCONSTRAINTS_HH
#define COMISO_NPLINEARCONSTRAINTS_HH


//== COMPILE-TIME PACKAGE REQUIREMENTS ========================================
#include <CoMISo/Config/config.hh>
#if COMISO_GMM_AVAILABLE

//== INCLUDES =================================================================

#include <CoMISo/Config/CoMISoDefines.hh>
#include <CoMISo/Utils/gmm.hh>
#include "NProblemGmmInterface.hh"
#include "LinearConstraintHandlerElimination.hh"
#include "LinearConstraintHandlerPenalty.hh"
//#include "LinearConstraintHandlerLagrange.hh"


//== FORWARDDECLARATIONS ======================================================

//== NAMESPACES ===============================================================

namespace COMISO {

//== CLASS DEFINITION =========================================================

	      

/** \class NPLinearConstraints{ Elimination, Penalty, Lagrange}
providing different techniques to handle linear constraints within
non-linear optimization. All three techniques transfer the constraint
problem into an unconstrained one.
*/


// ConstraintHandler working with Elimination approach
class COMISODLLEXPORT NPLinearConstraintsElimination : public NProblemGmmInterface, public LinearConstraintHandlerElimination
{
public:

  // use c-arrays as vectors for gmm
  typedef gmm::array1D_reference<double*> VectorPT;

  /// Default constructor having a pointer to the main problem
  NPLinearConstraintsElimination( NProblemGmmInterface* _np);
 
  // initialize constraints
  template<class MatrixT, class VectorT>
  void initialize_constraints( const MatrixT& _C, const VectorT& _c);


  // NSolverInterface
  virtual int    n_unknowns   ();

  virtual void   initial_x    ( double* _x );

  virtual double eval_f( const double* _x );

  virtual void   eval_gradient( const double* _x, double*    _g);

  virtual void   eval_hessian ( const double* _x, SMatrixNP& _H);

  virtual void   store_result ( const double* _x );

protected:

  void update_x(const double* _xr);

  // convert vector into pointer
  double* P(const std::vector<double>& _v);

private:
  NProblemGmmInterface* base_;

  bool cache_initialized_;

  std::vector<double> x_;
  std::vector<double> xr_;
  std::vector<double> vtemp_;
};


//=============================================================================


// ConstraintHandler working with Penalty
class COMISODLLEXPORT NPLinearConstraintsPenalty : public NProblemGmmInterface, public LinearConstraintHandlerPenalty
{
public:

  // use c-arrays as vectors for gmm
  typedef gmm::array1D_reference<double*> VectorPT;

  /// Default constructor having a pointer to the main problem
  NPLinearConstraintsPenalty( NProblemGmmInterface* _np);

  // initialize constraints
  template<class MatrixT, class VectorT>
  void initialize_constraints( const MatrixT& _C, const VectorT& _c);

  // NSolverInterface
  virtual int    n_unknowns   ();

  virtual void   initial_x    ( double* _x );

  virtual double eval_f( const double* _x );

  virtual void   eval_gradient( const double* _x, double*    _g);

  virtual void   eval_hessian ( const double* _x, SMatrixNP& _H);

  virtual void   store_result ( const double* _x );

private:
  NProblemGmmInterface* base_;
};


//=============================================================================
} // namespace COMISO
//=============================================================================
#if defined(INCLUDE_TEMPLATES) && !defined(COMISO_NSOLVERGMMINTERFACE_C)
#define COMISO_NSOLVERGMMINTERFACE_TEMPLATES
#include "NPLinearConstraintsT_impl.hh"
#endif

//=============================================================================
#endif // COMISO_GMM_AVAILABLE
//=============================================================================
#endif // COMISO_NSOLVERGMMINTERFACE_HH defined
//=============================================================================

