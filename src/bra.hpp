// ****************************************************************************
/// @file bra.hpp
/// @author Kyle Webster
/// @version 0.3
/// @date 20 Nov 2025
/// @brief Numerics Library - Algebra @ref bra
/// @details
/// Collection of algebraic structures and algorithms
// ****************************************************************************
#ifndef BRA_HPP
#define BRA_HPP
// ** Includes ****************************************************************
#include <concepts>
#include <cstdint>
#include <stdfloat>

#include "numli.hpp"
#include "simd.hpp"
// ****************************************************************************

namespace nl 
{ // ** nl ****************************

/// @namespace bra
/// @brief Algebraic data structures and algorithms
/// @details
/// Contents:
///  * @ref operators
///  * @ref vector_spaces
namespace bra
{ // ** nl::bra ***********************

// ** concepts ************************
template <typename T> concept arithmetic = 
  std::integral<T> || std::floating_point<T>;

template <typename Op, typename T> concept binary_operator = 
  requires(Op op, T a, T b) { {op(a, b)} -> std::same_as<T>; };

template <typename Op, typename T> concept unary_operator = 
  requires(Op op, T a) { {op(a)} -> std::same_as<T>; };
// ** end of concepts ***************

// ****************************************************************************
/// @name operators
/// @brief structures representing operators
/// @warning For SIMD, only AVX2 +,-,*,/ for floats are defined

template<arithmetic T, simd::SIMD ISA>
struct Operator {};

// ** Binary Operators ****************
template<arithmetic T> struct op_add : Operator<T, simd::SIMD_MODE>
{ 
  static constexpr T serial(T a, T b) { return a+b; }
  static void SIMD(
    T       * restrict x, 
    T const * restrict a, 
    T const * restrict b, 
    size_t             n)
  {
    if constexpr (simd::SIMD_MODE == simd::AVX2)
    {
      if constexpr (std::same_as<T,float>) simd::add_avx2_float(x,a,b,n);
    }
  }
  constexpr T operator()(T a, T b) { serial(a,b); }
};

template<arithmetic T> struct op_sub : Operator<T, simd::SIMD_MODE>
{ 
  static constexpr T serial(T a, T b) { return a-b; }
  static void SIMD(
    T       * restrict x, 
    T const * restrict a, 
    T const * restrict b, 
    size_t             n)
  {
    if constexpr (simd::SIMD_MODE == simd::AVX2)
    {
      if constexpr (std::same_as<T,float>) simd::sub_avx2_float(x,a,b,n);
    }
  }
  constexpr T operator()(T a, T b) { serial(a,b); }
};

template<arithmetic T> struct op_mul : Operator<T, simd::SIMD_MODE>
{ 
  static constexpr T serial(T a, T b) { return a*b; }
  static void SIMD(
    T       * restrict x, 
    T const * restrict a, 
    T const * restrict b, 
    size_t             n)
  {
    if constexpr (simd::SIMD_MODE == simd::AVX2)
    {
      if constexpr (std::same_as<T,float>) simd::mul_avx2_float(x,a,b,n);
    }
  }
  constexpr T operator()(T a, T b) { serial(a,b); }
};

template<arithmetic T> struct op_div : Operator<T, simd::SIMD_MODE>
{ 
  static constexpr T serial(T a, T b) { return a/b; }
  static void SIMD(
    T       * restrict x, 
    T const * restrict a, 
    T const * restrict b, 
    size_t             n)
  {
    if constexpr (simd::SIMD_MODE == simd::AVX2)
    {
      if constexpr (std::same_as<T,float>) simd::div_avx2_float(x,a,b,n);
    }
  }
  constexpr T operator()(T a, T b) { serial(a,b); }
};

// ** Reductions **********************
template<arithmetic T, size_t n> 
struct op_inner_prod : Operator<T, simd::SIMD_MODE>
{
  static constexpr T serial(T const (&a)[n], T const (&b)[n])
    {T sum = T(0); for (int i=0; i<n; i++) sum += a[i]*b[i]; return sum;}
  /// @todo write SIMD inner product  
  static void SIMD(T &x, T const* restrict a, T const* restrict b);
  constexpr T operator()(T const (&a)[n], T (&b)[n]) { serial(a,b); }
};

// ****************************************************************************
/// @name vector spaces
/// @details
/// Contents:
///  * @ref ℝn - number in ℝⁿ
///  * @ref ℝnxm - number in ℝⁿˣᵐ, stored in row-major order

// ************************************
/// @name array operators

/// @brief A binary operation performed element-wise between arrays
/// @attention a, b, and x must not alias each other
template <uint32_t n, arithmetic T, binary_operator<T> Op>
constexpr void opBinaryArray(
  T       * restrict x,
  T const * restrict a,
  T const * restrict b)
{
  if constexpr (n > simd::simd<T>::lanes)
  {
    Op::SIMD(x, a, b, n);
  } else { for (size_t i=0;i<n;i++) {x[i] = Op::serial(a[i],b[i]);} }
}

/// @brief A unary operation performed element-wise on an array
/// @attention a and x must not alias each other
template <uint32_t n, arithmetic T, unary_operator<T> Op>
constexpr void opUnaryArray(T * restrict x, T const * restrict a)
{
  if constexpr (n > simd::simd<T>::lanes)
  {
    /// @todo SIMD unary operation
  } else { for (size_t i=0;i<n;i++) {x[i] = Op{}(a[i]);} }
}
// ** end of Array Operators **********

/// @struct ℝn
/// @brief a number in ℝⁿ
/// @tparam n number of basis elements that square to 1
/// @tparam T storage type, default float64_t
/// @attention T must be arithmetic
/// @details
/// Element of ℝⁿ with type T. Stored aligned to SIMD width if available.
template<uint32_t n, arithmetic T = std::float64_t> struct ℝn
{
  /// @name members
  alignas(simd::simd<T>::alignment) T elem[n]; ///< basis coefficients

  /// @name constructors
  constexpr explicit ℝn( T const restrict (&a)[n] ) { copy(elem, a); }
  constexpr explicit ℝn( ℝn<n,T> const &x ) { copy(elem, x.elem); }
  constexpr explicit ℝn( ℝn<n,T> const *x ) { copy(elem, x->elem); }
  constexpr explicit ℝn( T v ) { for (uint32_t i=0;i<n;i++) elem[i]=v; }
  template <binary_operator<T> Op> constexpr explicit ℝn(
    ℝn<n,T> const &x, ℝn<n,T> const &y, std::type_identity<Op>)
    { opBinaryArray<n,T,Op>(elem, x.elem, y.elem); }
  constexpr ℝn() {}

  /// @name member functions
  constexpr T l2() const 
    {T sum=T(0); for(uint32_t i=0;i<n;i++)sum+=elem[i]; return std::sqrt(sum);}
  constexpr void normalize()
    { T len=l2(); for (uint32_t i=0;i<n;i++) elem[i]/=len; }
  constexpr ℝn normalized() const { ℝn x(this); x.normalize(); return x; }
  constexpr void negate() { for (uint32_t i=0;i<n;i++) elem[i]*=-1; }
  constexpr ℝn negated() { ℝn x(this); x.negate(); return x; }

  /// @name operators
  constexpr T& operator[](size_t i)       noexcept { return elem[i]; }
  constexpr T  operator[](size_t i) const noexcept { return elem[i]; }
};
// ** end ℝn **************************

/// @struct ℝnxm
/// @brief a number in ℝⁿˣᵐ
/// @tparam n number of rows
/// @tparam m number of columns
/// @tparam T storage type, default float64_t
/// @attention T must be arithmetic
/// @details
/// Element of ℝⁿˣᵐ with type T. Stored aligned to SIMD width if available.
template<uint32_t n, uint32_t m, arithmetic T = std::float64_t> 
struct ℝnxm
{
  /// @name members
  static constexpr uint64_t N = n*m;           ///< number of elements
  alignas(simd::simd<T>::alignment) T elem[N]; ///< basis coefficients

  /// @name constructors
  constexpr explicit ℝnxm( T const * restrict (&a)[n] ) {nl::copy(elem, a);}
  explicit ℝnxm( ℝn<n,T> const &x ) { nl::copy(elem,x.elem); }
  explicit ℝnxm( T v ) { for (uint64_t k=0;k<N;k++) elem[k]=v; }
  ℝnxm() {}

  /// @name member functions
  constexpr void identity() 
  {
    for (uint64_t i=0;i<n;i++) for (uint64_t j=0;j<m;j++) 
      if (i==j) {elem[i*m+j]=T(1);} else {elem[i*m+j]=T(0);}
  }

  /// @name operators
  constexpr T& operator()(size_t i, size_t j) noexcept 
    { return elem[i*m+j]; }
  constexpr T  operator()(size_t i, size_t j) const noexcept
    { return elem[i*m+j]; }
};
// ** end ℝnxm ************************

// ** Operators on ℝn and ℝnxm ********

/// @name binary operators
template <uint32_t n, arithmetic T>
constexpr ℝn<n,T> operator+(ℝn<n,T> const &x, ℝn<n,T> const &y)
  { return ℝn<n,T>(x,y,std::type_identity<op_add<T>>{}); }
template <uint32_t n, arithmetic T>
constexpr ℝn<n,T> operator-(ℝn<n,T> const &x, ℝn<n,T> const &y)
  { return ℝn<n,T>(x,y,std::type_identity<op_sub<T>>{}); }
template <uint32_t n, arithmetic T>
constexpr ℝn<n,T> operator*(ℝn<n,T> const &x, ℝn<n,T> const &y)
  { return ℝn<n,T>(x,y,std::type_identity<op_mul<T>>{}); }
template <uint32_t n, arithmetic T>
constexpr ℝn<n,T> operator/(ℝn<n,T> const &x, ℝn<n,T> const &y)
  { return ℝn<n,T>(x,y,std::type_identity<op_div<T>>{}); }

/// @brief cross product in ℝ3
template <arithmetic T>
constexpr ℝn<3,T> operator^(ℝn<3,T> const &x, ℝn<3,T> const &y)
  { return {x[1]*y[2]-x[2]*y[1], x[2]*y[0]-x[0]*y[2], x[0]*y[1]-x[1]*y[0]}; }

// ** end operators on ℝn and ℝnxm ****

// ** end of Vector Spaces ****************************************************
} // ** end of namespace bra **********

// ************************************
/// @name aliases
using ℝ3 = bra::ℝn<3,float>;
using ℝ4 = bra::ℝn<4,float>;
using ℝ3x4 = bra::ℝnxm<3,4,float>;

} // ** end of namespace nl ***********

// ****************************************************************************
#endif // #ifndef BRA_HPP