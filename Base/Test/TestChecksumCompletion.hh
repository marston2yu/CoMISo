// (C) Copyright 2021 by Autodesk, Inc.

#ifndef BASE_CHECKSUMCOMPLETION_HH_INCLUDE
#define BASE_CHECKSUMCOMPLETION_HH_INCLUDE

#ifdef TEST_ON

#include <Base/Security/Mandatory.hh>
#include <Base/Test/TestChecksum.hh>

namespace Test
{
namespace Checksum
{

// Writes a checksum when the test completes.
class Completion : public Object
{
public:
  Completion();

  //! Record the test "end" completion checksum
  void record_end();

  //! Get if the line contains the end completion checksum record
  static bool end(const std::string& _line);
};

// Register the checksum to check test completion.
extern Completion completion;

} // namespace Checksum
} // namespace Test

#endif // TEST_ON
#endif // BASE_CHECKSUMCOMPLETION_HH_INCLUDE
