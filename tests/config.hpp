// ****************************************************************************
/// @file config.hpp
/// @author Kyle Webster
/// @date 20 Nov 2025
/// @details
/// Configuration for doctest tests
// ****************************************************************************
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include <print>
#include "numli.hpp"
// ****************************************************************************

// same-length array comparison
inline auto eq = [](const float* a, const float* b, int n) 
  { return std::equal(a, a + n, b); };

// ****************************************************************************