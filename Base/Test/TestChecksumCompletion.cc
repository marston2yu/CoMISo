// (C) Copyright 2021 by Autodesk, Inc.

#ifdef TEST_ON

#include "TestChecksumCompletion.hh"
#include <Base/Test/TestChecksum.hh>

#include <sstream>
#include <string>

namespace Test
{
namespace Checksum
{

namespace
{
const char* const END = "END";
} // namespace

Completion::Completion() : Object("Completion", L_STABLE) {}

void Completion::record_end()
{
  std::stringstream mess;
  mess << "Success: " << END;

  add(Result::OK, mess.str());
}

bool Completion::end(const std::string& _line)
{
  return _line.find(Checksum::completion.name()) != std::string::npos &&
         _line.find(END) != std::string::npos;
}

// Register the checksum to check test completion.
Completion completion;

} // namespace Checksum
} // namespace Test

#endif // TEST_ON
