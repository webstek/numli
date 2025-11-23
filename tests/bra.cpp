// ****************************************************************************
/// @file bra.cpp
/// @author Kyle Webster
/// @date 22 Nov 2025
/// @details
/// Tests for the algebra module (bra.hpp)
// ****************************************************************************

// ** Includes ****************************************************************
#include "config.hpp"
#include "bra.hpp"
// ****************************************************************************

TEST_CASE("ℝn - constructors and element-wise operators")
{
  // array constructor
  float a3_arr[3] = {1.f, 2.f, 3.f};
  float b3_arr[3] = {2.f, 1.f, 4.f};

  nl::bra::ℝn<3,float> A(a3_arr);
  nl::bra::ℝn<3,float> B(b3_arr);

  // copy and value constructors
  nl::bra::ℝn<3,float> C(A);
  for (int i = 0; i < 3; ++i) CHECK(C.elem[i] == A.elem[i]);

  nl::bra::ℝn<3,float> V(1.f);
  for (int i = 0; i < 3; ++i) CHECK(V.elem[i] == 1.f);

  // addition
  auto S = A + B;
  CHECK(S.elem[0] == 3.f);
  CHECK(S.elem[1] == 3.f);
  CHECK(S.elem[2] == 7.f);

  // subtraction
  auto D = A - B;
  CHECK(D.elem[0] == -1.f);
  CHECK(D.elem[1] == 1.f);
  CHECK(D.elem[2] == -1.f);

  // multiplication
  auto M = A * B;
  CHECK(M.elem[0] == 2.f);
  CHECK(M.elem[1] == 2.f);
  CHECK(M.elem[2] == 12.f);

  // division
  auto R = A / B;
  CHECK(R.elem[0] == 0.5f);
  CHECK(R.elem[1] == 2.f);
  CHECK(R.elem[2] == 0.75f);

  // larger vector to exercise non-trivial sizes (serial path or SIMD depending on build)
  nl::bra::ℝn<64,float> L(1.f);
  nl::bra::ℝn<64,float> M1(2.f);
  auto Lsum = L + M1;
  for (int i = 0; i < 64; ++i) CHECK(Lsum.elem[i] == 3.f);
}
