// ****************************************************************************
/// @file bra.cpp
/// @author Kyle Webster
/// @date 20 Nov 2025
/// @details
/// Tests for the algebra module
// ****************************************************************************

// ** Includes ****************************************************************
#include "config.hpp"
#include "bra.hpp"
// ****************************************************************************

TEST_CASE("ℝn tests")
{ // test instantiation using all constructors
  nl::bra::ℝn<3,float> u;
  nl::bra::ℝn<3,float> v(1.f);
  for (int i=0; i<3; i++) { CHECK(v.elem[i] == 1.f); }
  nl::bra::ℝn<3,float> w(v);
  for (int i=0; i<3; i++) { CHECK(w.elem[i] == v.elem[i]); }

  float x[3] = {1.f, 2.f, 3.f};
  nl::bra::ℝn<3,float> y(x);
  for (int i=0; i<3; i++) { CHECK(y.elem[i] == i+1); }

  nl::bra::ℝn<3,float> z = v+y;
  for (int i=0; i<3; i++) { CHECK(z.elem[i] == i+2); }

  nl::bra::ℝn<64,double> s(1.);
  nl::bra::ℝn<64,double> t(1.);
  auto k = s+t;
  for (int i=0;i<64;i++) { CHECK(k.elem[i] == 2.); }
}

TEST_CASE("ℝnxm tests")
{
  nl::bra::ℝnxm<3,4,double> U;
  nl::bra::ℝnxm<5,3,int> M(1);
  for (int i=0;i<5;i++) for (int j=0;j<3;j++) { CHECK(M.elem[i*3+j] == 1); }
}