//=============================================================================
//
//  CLASS NProblemInterface
//
//=============================================================================


#ifndef COMISO_NPROBLEMINTERFACE_HH
#define COMISO_NPROBLEMINTERFACE_HH


//== COMPILE-TIME PACKAGE REQUIREMENTS ========================================
#include <CoMISo/Config/config.hh>
#if COMISO_EIGEN3_AVAILABLE

//== INCLUDES =================================================================

#include <Base/Code/Quality.hh>
#include <iostream>
#include <cfloat>

LOW_CODE_QUALITY_SECTION_BEGIN
#include <Eigen/Eigen>
#if !(EIGEN_VERSION_AT_LEAST(3,1,0))
  #define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#endif
#include <Eigen/Sparse>
LOW_CODE_QUALITY_SECTION_END

#include <CoMISo/Config/CoMISoDefines.hh>

//== FORWARDDECLARATIONS ======================================================

//== NAMESPACES ===============================================================

namespace COMISO {

//== CLASS DEFINITION =========================================================



/** \class NProblemInterface NProblemInterface.hh <COMISO/NSolver/NProblemInterface.hh>

    Brief Description.

    A more elaborate description follows.
*/
class COMISODLLEXPORT NProblemInterface
{
public:

  // Sparse Matrix Type
 #if EIGEN_VERSION_AT_LEAST(3,1,0)
  typedef Eigen::SparseMatrix<double,Eigen::ColMajor> SMatrixNP;
 #else
  typedef Eigen::DynamicSparseMatrix<double,Eigen::ColMajor> SMatrixNP;
 #endif

  /// Default constructor
  NProblemInterface();

  /// Destructor
  virtual ~NProblemInterface();

  // problem definition
  virtual int    n_unknowns   (                                ) = 0;
  virtual void   initial_x    (       double* _x               ) = 0;
  virtual double eval_f       ( const double* _x               ) = 0;
  virtual void   eval_gradient( const double* _x, double*    _g) = 0;
  virtual void   eval_hessian ( const double* _x, SMatrixNP& _H) = 0;
  virtual void   store_result ( const double* _x               ) = 0;

  // advanced properties
  virtual bool   constant_gradient ()                                    const { return false; }
  virtual bool   constant_hessian  ()                                    const { return false; }
  virtual double max_feasible_step ( const double* _x, const double* _v)       { return DBL_MAX; }
};


//=============================================================================
} // namespace COMISO
//=============================================================================
#endif // COMISO_EIGEN3_AVAILABLE
//=============================================================================
#endif // COMISO_NPROBLEMINTERFACE_HH defined
//=============================================================================

