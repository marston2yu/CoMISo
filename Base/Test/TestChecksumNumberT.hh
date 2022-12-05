// (C) Copyright 2022 by Autodesk, Inc.

#ifndef BASE_TESTCHECKSUMNUMBERT_HH_INCLUDED
#define BASE_TESTCHECKSUMNUMBERT_HH_INCLUDED

#include <Base/Test/TestChecksum.hh>

#ifdef TEST_ON

#include <cmath>
#include <sstream>
#include <type_traits>

namespace Test
{
namespace Checksum
{

/*! Comparison class using opeartor== of ValueT.
*/
template <typename ValueT>
struct EqualityCompareT
{
  bool same(const ValueT& _a, const ValueT& _b) const { return _a == _b; }
};

/*! Comparison class for checksums that MUST only be logged.
*/
template <typename ValueT>
struct NoCompareT
{
  bool same(const ValueT& _a, const ValueT& _b) const { return true; }
};


enum class ToleranceType
{
  ABSOLUTE,
  RELATIVE
};

enum
{
  DUMMY_EXPONENT = 1234
};

/*!
Utility class to compare values with a tolerance. Tolerance is specified via the
template argument as 10^EXPONENT. If RELATIVE is true, the difference is divided
by the right (new) value before comparison.
*/
template <typename ScalarT, int EXPONENT,
    ToleranceType TOLERANCE_TYPE = ToleranceType::ABSOLUTE>
struct TolerantCompareT
{

  static_assert(EXPONENT != DUMMY_EXPONENT,
      "For floating point checksums, please explicitly define a comparison "
      "class, e.g. this TolerantCompareT with an exponent suitable for your "
      "checksum.");
  static_assert(
      EXPONENT <= std::numeric_limits<ScalarT>::max_exponent10 ||
          EXPONENT == DUMMY_EXPONENT, // disable assertion for default exponent
      "Exponent too large");
  static_assert(EXPONENT >= std::numeric_limits<ScalarT>::min_exponent10,
      "Exponent too small");

  bool same(const ScalarT& _a, const ScalarT& _b) const
  {
    static const auto TOLERANCE = std::pow(10, EXPONENT);

    // Prevent division by zero
    if constexpr (TOLERANCE_TYPE == ToleranceType::RELATIVE)
    {
      if (_b == 0)
        return _a == _b;
    }

    double diff = std::fabs(_a - _b);

    if constexpr (TOLERANCE_TYPE == ToleranceType::RELATIVE)
      diff /= _b;

    return diff <= TOLERANCE;
  }
};

// The default Compare type for ValueT. For integral types exact equality is
// checked. For floating point type TolerantCompare is used which actually
// generates a compile error and asks the developer to choose a meaningful
// threshold.
template <typename ValueT>
using DefaultCompareT = std::conditional_t<   // if
    std::is_floating_point_v<ValueT>,         // ValueT is floating point
    TolerantCompareT<ValueT, DUMMY_EXPONENT>, // use TolerantCompareT
    EqualityCompareT<ValueT>>;                // else use EqualityCompare

/*!
Generic checksum class to record and compare a value of a certain type.
*/
template <typename ValueT, class CompareT = DefaultCompareT<ValueT>>
class NumberT : public Object
{
public:
  NumberT(const char* const _name,       //!<[in] Checksum name
      const Level _lvl = L_ALL,          //!<[in] Checksum level
      const CompareT& _comp = CompareT() //!<[in] Comparison function
      )
      : Object(_name, _lvl), comp_(_comp)
  {
  }

protected:
  //! Generic implementation of number data comparison
  virtual Difference compare_data(const String& _old, const String& _new) const
  {
    std::istringstream strm_old(_old), strm_new(_new);
    ValueT val_old, val_new;
    strm_old >> val_old;
    strm_new >> val_new;

    // TODO: multiple comparisons, can we use just one?
    if (val_new == val_old) // bitwise comparison
      return Difference::EQUAL;

    Base::OStringStream diff;
    diff << (val_new - val_old);

    if (comp_.same(val_old, val_new)) // tolerance comparison
      return Difference(Difference::NEGLIGIBLE, diff.str);

    return Difference(Difference::UNKNOWN, diff.str);
  }

private:
  CompareT comp_; // Compare class.
};

template <typename ValueT, int EXPONENT,
    ToleranceType TOLERANCE_TYPE = ToleranceType::ABSOLUTE>
using TolerantNumberT =
    NumberT<ValueT, TolerantCompareT<ValueT, EXPONENT, TOLERANCE_TYPE>>;


}//namespace Checksum
}//namespace Test

#endif//TEST_ON
#endif//BASE_TESTCHECKSUMNUMBERT_HH_INCLUDED
