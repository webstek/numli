// ****************************************************************************
/// @file simd.cpp
/// @author Kyle Webster
/// @date 22 Nov 2025
/// @details
/// Tests for the simd module
// ****************************************************************************

// ** Includes ****************************************************************
#include <algorithm>

#include "config.hpp"
#include "simd.hpp"
// ****************************************************************************

TEST_CASE("AVX2 Tests") {
  alignas(32) float a5[5]  = {4,4,4,4,4};
  alignas(32) float b5[5]  = {2,2,2,2,2};
  alignas(32) float out5[5];

  float add5[5] = {6,6,6,6,6};
  float sub5[5] = {2,2,2,2,2};
  float mul5[5] = {8,8,8,8,8};
  float div5[5] = {2,2,2,2,2};

  nl::simd::add_avx2_float(out5,a5,b5,5); CHECK(eq(out5,add5,5));
  nl::simd::sub_avx2_float(out5,a5,b5,5); CHECK(eq(out5,sub5,5));
  nl::simd::mul_avx2_float(out5,a5,b5,5); CHECK(eq(out5,mul5,5));
  nl::simd::div_avx2_float(out5,a5,b5,5); CHECK(eq(out5,div5,5));

  alignas(32) float a8[8]  = {2,2,2,2,2,2,2,2};
  alignas(32) float b8[8]  = {2,2,2,2,2,2,2,2};
  alignas(32) float out8[8];

  float add8[8] = {4,4,4,4,4,4,4,4};
  float sub8[8] = {0,0,0,0,0,0,0,0};
  float mul8[8] = {4,4,4,4,4,4,4,4};
  float div8[8] = {1,1,1,1,1,1,1,1};

  nl::simd::add_avx2_float(out8,a8,b8,8); CHECK(eq(out8,add8,8));
  nl::simd::sub_avx2_float(out8,a8,b8,8); CHECK(eq(out8,sub8,8));
  nl::simd::mul_avx2_float(out8,a8,b8,8); CHECK(eq(out8,mul8,8));
  nl::simd::div_avx2_float(out8,a8,b8,8); CHECK(eq(out8,div8,8));

  alignas(32) float a17[17] = {3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3};
  alignas(32) float b17[17] = {3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3};
  alignas(32) float out17[17];

  float add17[17] = {6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6};
  float sub17[17] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  float mul17[17] = {9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9};
  float div17[17] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};

  nl::simd::add_avx2_float(out17,a17,b17,17); CHECK(eq(out17,add17,17));
  nl::simd::sub_avx2_float(out17,a17,b17,17); CHECK(eq(out17,sub17,17));
  nl::simd::mul_avx2_float(out17,a17,b17,17); CHECK(eq(out17,mul17,17));
  nl::simd::div_avx2_float(out17,a17,b17,17); CHECK(eq(out17,div17,17));
}