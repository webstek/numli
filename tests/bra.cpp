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

TEST_CASE("ℝn - member functions: l2 (norm)")
{
  // test with simple values where norm is easy to compute
  // v = [3, 4] -> l2 = sqrt(9 + 16) = sqrt(25) = 5
  nl::bra::ℝn<2,float> v{3.f, 4.f};
  CHECK(v.l2() == doctest::Approx(5.f).epsilon(1e-6f));

  // v = [1, 1, 1] -> l2 = sqrt(1 + 1 + 1) = sqrt(3) ≈ 1.732
  nl::bra::ℝn<3,float> v2{1.f, 1.f, 1.f};
  CHECK(v2.l2() == doctest::Approx(std::sqrt(3.f)).epsilon(1e-6f));

  // v = [2, 0, 0] -> l2 = 2
  nl::bra::ℝn<3,float> v3{2.f, 0.f, 0.f};
  CHECK(v3.l2() == doctest::Approx(2.f).epsilon(1e-6f));

  // zero vector -> l2 = 0
  nl::bra::ℝn<3,float> zero{0.f, 0.f, 0.f};
  CHECK(zero.l2() == doctest::Approx(0.f).epsilon(1e-6f));
}

TEST_CASE("ℝn - member functions: normalize and normalized")
{
  // test normalize: v = [3, 4] -> normalized = [0.6, 0.8]
  nl::bra::ℝn<2,float> v{3.f, 4.f};
  v.normalize();
  CHECK(v.elem[0] == doctest::Approx(0.6f).epsilon(1e-6f));
  CHECK(v.elem[1] == doctest::Approx(0.8f).epsilon(1e-6f));
  // verify the normalized vector has l2 norm of 1
  CHECK(v.l2() == doctest::Approx(1.f).epsilon(1e-6f));

  // test normalized (const version): v = [1, 1, 1]
  nl::bra::ℝn<3,float> v2{1.f, 1.f, 1.f};
  auto v2_normalized = v2.normalized();
  float expected = 1.f / std::sqrt(3.f);
  CHECK(v2_normalized.elem[0] == doctest::Approx(expected).epsilon(1e-6f));
  CHECK(v2_normalized.elem[1] == doctest::Approx(expected).epsilon(1e-6f));
  CHECK(v2_normalized.elem[2] == doctest::Approx(expected).epsilon(1e-6f));
  // verify original is unchanged
  CHECK(v2.elem[0] == 1.f);
  CHECK(v2.elem[1] == 1.f);
  CHECK(v2.elem[2] == 1.f);
  // verify result has l2 norm of 1
  CHECK(v2_normalized.l2() == doctest::Approx(1.f).epsilon(1e-6f));

  // test normalized with unit vector (should remain unchanged)
  nl::bra::ℝn<3,float> unit{1.f, 0.f, 0.f};
  auto unit_normalized = unit.normalized();
  CHECK(unit_normalized.elem[0] == doctest::Approx(1.f).epsilon(1e-6f));
  CHECK(unit_normalized.elem[1] == doctest::Approx(0.f).epsilon(1e-6f));
  CHECK(unit_normalized.elem[2] == doctest::Approx(0.f).epsilon(1e-6f));
}

TEST_CASE("ℝn - member functions: negate and negated")
{
  // test negate: v = [1, -2, 3] -> negated = [-1, 2, -3]
  nl::bra::ℝn<3,float> v{1.f, -2.f, 3.f};
  v.negate();
  CHECK(v.elem[0] == -1.f);
  CHECK(v.elem[1] == 2.f);
  CHECK(v.elem[2] == -3.f);

  // test negated (const version)
  nl::bra::ℝn<3,float> v2{5.f, -3.f, 0.f};
  auto v2_negated = v2.negated();
  CHECK(v2_negated.elem[0] == -5.f);
  CHECK(v2_negated.elem[1] == 3.f);
  CHECK(v2_negated.elem[2] == 0.f);
  // verify original is unchanged
  CHECK(v2.elem[0] == 5.f);
  CHECK(v2.elem[1] == -3.f);
  CHECK(v2.elem[2] == 0.f);
}

TEST_CASE("ℝnxm - member functions: identity")
{
  // test 2x2 identity
  nl::bra::ℝnxm<2,2,float> I2;
  I2.identity();
  CHECK(I2(0,0) == 1.f);
  CHECK(I2(0,1) == 0.f);
  CHECK(I2(1,0) == 0.f);
  CHECK(I2(1,1) == 1.f);

  // test 3x3 identity
  nl::bra::ℝnxm<3,3,float> I3;
  I3.identity();
  for (int i = 0; i < 3; ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      if (i == j) CHECK(I3(i,j) == 1.f);
      else CHECK(I3(i,j) == 0.f);
    }
  }

  // test 4x4 identity
  nl::bra::ℝnxm<4,4,float> I4;
  I4.identity();
  for (int i = 0; i < 4; ++i)
  {
    for (int j = 0; j < 4; ++j)
    {
      if (i == j) CHECK(I4(i,j) == 1.f);
      else CHECK(I4(i,j) == 0.f);
    }
  }
}

TEST_CASE("ℝn - cross product")
{
  // cross product of [1,0,0] x [0,1,0] = [0,0,1]
  nl::ℝ3 x_axis{1.f, 0.f, 0.f};
  nl::ℝ3 y_axis{0.f, 1.f, 0.f};
  auto result = x_axis ^ y_axis;
  CHECK(result.elem[0] == doctest::Approx(0.f).epsilon(1e-6f));
  CHECK(result.elem[1] == doctest::Approx(0.f).epsilon(1e-6f));
  CHECK(result.elem[2] == doctest::Approx(1.f).epsilon(1e-6f));

  // cross product is anti-commutative: a x b = -(b x a)
  auto result2 = y_axis ^ x_axis;
  CHECK(result2.elem[0] == doctest::Approx(0.f).epsilon(1e-6f));
  CHECK(result2.elem[1] == doctest::Approx(0.f).epsilon(1e-6f));
  CHECK(result2.elem[2] == doctest::Approx(-1.f).epsilon(1e-6f));

  // cross product with parallel vectors is zero
  nl::ℝ3 v1{2.f, 4.f, 6.f};
  nl::ℝ3 v2{1.f, 2.f, 3.f};
  auto zero_result = v1 ^ v2;
  CHECK(zero_result.elem[0] == doctest::Approx(0.f).epsilon(1e-6f));
  CHECK(zero_result.elem[1] == doctest::Approx(0.f).epsilon(1e-6f));
  CHECK(zero_result.elem[2] == doctest::Approx(0.f).epsilon(1e-6f));
}
