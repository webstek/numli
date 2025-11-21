// ****************************************************************************
/// @file bra.hpp
/// @author Kyle Webster
/// @version 0.2
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
  // ** forward declrations ***********

  
  // ** concepts **********************
  template <typename T> concept arithmetic = 
    std::integral<T> || std::floating_point<T>;
  
  template <typename Op, typename T> concept binary_operator = 
    requires(Op op, T a, T b) { {op(a, b)} -> std::same_as<T>; };
  
  template <typename Op, typename T> concept unary_operator = 
    requires(Op op, T a) { {op(a)} -> std::same_as<T>; };
  // ** end of concepts ***************

  // **************************************************************************
  /// @name operators
  /// @brief structures representing operators

  template<arithmetic T, SIMD ISA>
  struct Operator {};

  template<arithmetic T>
  struct op_add : Operator<T, SIMD_MODE>
    { constexpr T operator()(T const a, T const b) const {return a+b;} };
  
  // **************************************************************************
  /// @name vector spaces
  /// @details
  /// Contents:
  ///  * @ref ℝn - number in ℝⁿ
  ///  * @ref ℝnxm - number in ℝⁿˣᵐ, stored in row-major order

  // **********************************
  /// @name array operators
  
  /// @brief A binary operation performed element-wise between arrays
  /// @attention a, b, and x must not alias each other
  template <uint32_t n, arithmetic T, binary_operator<T> op>
  constexpr void opBinaryArray(
    T       * restrict x,
    T const * restrict a,
    T const * restrict b)
  {
    if constexpr (n > nl::simd<T>::lanes)
    {
      /// @todo SIMD binary operation
    } else { for (std::size_t i=0;i<n;i++) {x[i] = op(a[i],b[i]);} }
  }

  /// @brief A unary operation performed element-wise on an array
  /// @attention a and x must not alias each other
  template <uint32_t n, arithmetic T, unary_operator<T> op>
  constexpr void opUnaryArray(T * restrict x, T const * restrict a)
  {
    if constexpr (n > nl::simd<T>::lanes)
    {
      /// @todo SIMD unary operation
    } else { for (std::size_t i=0;i<n;i++) {x[i] = op(a[i]);} }
  }
  // ** end of Array Operators ********

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
    alignas(nl::simd<T>::alignment) T elem[n]; ///< basis element coefficients

    /// @name constructors
    constexpr explicit ℝn( T const restrict (&a)[n] ) { nl::copy(elem, a); }
    constexpr explicit ℝn( ℝn<n,T> const &x ) { nl::copy(elem, x.elem); }
    constexpr explicit ℝn( T v ) { for (uint32_t i=0;i<n;i++) elem[i]=v; }
    template <binary_operator<T> op> constexpr explicit ℝn(
      ℝn<n,T> const &x, ℝn<n,T> const &y)
      { opBinaryArray<n,T,op>(elem, x.elem, y.elem); }
      // { for (uint32_t i=0;i<n;i++) elem[i]=op(x.elem[i],y.elem[i]); }
    constexpr ℝn() {}
  };
  // ** end ℝn ************************

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
    static constexpr uint64_t N = n*m;         ///< number of elements
    alignas(nl::simd<T>::alignment) T elem[N]; ///< basis element coefficients

    /// @name constructors
    constexpr explicit ℝnxm( T const * restrict (&a)[n] ) {nl::copy(elem, a);}
    explicit ℝnxm( ℝn<n,T> const &x ) { nl::copy(elem,x.elem); }
    explicit ℝnxm( T v ) 
      { for(uint32_t i=0;i<n;i++) for(uint32_t j=0;j<m;j++) elem[i*m+j]=v; }
    ℝnxm() {}
  };
  // ** end ℝnxm **********************

  // ** Operators on ℝn and ℝnxm ******

  /// @brief Element-wise addition of two ℝn
  /// @bug figure out correct way to work in the op_add<T> operator
  template <uint32_t n, std::floating_point T>
  constexpr ℝn<n,T> operator+(ℝn<n,T> const &x, ℝn<n,T> const &y)
    { return ℝn<n,T>(x,y); }
  // ** end operators on ℝn and ℝnxm **
  // ** end of Vector Spaces **************************************************
} // ** end of namespace bra **********
} // ** end of namespace nl ***********

// ****************************************************************************
#endif // #ifndef BRA_HPP