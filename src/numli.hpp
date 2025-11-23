// ****************************************************************************
/// @file numli.hpp
/// @author Kyle Webster
/// @version 0.3
/// @date 22 Nov 2025
/// @brief Numerics Library - @ref nl
/// @details
/// Collection of constants and core numerics
// ****************************************************************************
#ifndef NUMLI_HPP
#define NUMLI_HPP
// ** Includes ****************************************************************
#include <numbers>
#include <stdfloat>
#include <cstring>
// ****************************************************************************

// ************************************
/// @namespace nl
/// @brief Numerics Library main module
/// @details
/// Contents:
///  * @ref compiler - compatibility for MSC and GCC
///  * @ref utilities - language utilities
///  * @ref constants - mathematical constants
namespace nl
{
  // **************************************************************************
  /// @name compiler
  /// @brief defines macros depending on which compiler is being used

  // ** restrict **********************
  #if defined(__GNUC__)
  #define restrict __restrict__
  #elif defined(_MSC_VER)
  #define restrict __restrict
  #else
  #define restrict
  #endif
  // ** end of compiler *******************************************************

  // **************************************************************************
  /// @name utilities
  /// @brief common shorthand functions and misc items

  // ** copy *********************
  template <std::size_t n, typename T>
  constexpr void copy_impl(T (&dst)[n], T const (&src)[n], std::true_type)
    { for (std::size_t i=0; i<n; ++i) {dst[i] = src[i];} }
  template <std::size_t n, typename T>
  void copy_impl(T (&dst)[n], T const (&src)[n], std::false_type)
    { std::memcpy(dst, src, n * sizeof(T)); }
  /// @brief copies from src to dst optimized for constexpr or runtime call
  /// @tparam T type to copy
  /// @tparam n number of elements of T to copy
  /// @param dst destination array
  /// @param src source array
  template <std::size_t n, typename T>
  constexpr void copy(T (&dst)[n], T const (&src)[n])
  {
    copy_impl(dst, src, 
      std::integral_constant<bool, std::is_constant_evaluated()>{});
  }
  // ** end of utilities ******************************************************

  // **************************************************************************
  /// @name constants
  /// @brief Mathematical and computational constants
  template<typename T> constexpr T π     = std::numbers::pi_v<T>;
  template<typename T> constexpr T π²    = 
    std::numbers::pi_v<T>*std::numbers::pi_v<T>;
  template<typename T> constexpr T π_2   = std::numbers::pi_v<T>/2;
  template<typename T> constexpr T inv_π = 1/std::numbers::pi_v<T>; 
  template<typename T> constexpr T e     = std::numbers::e_v<T>;
  template<typename T> constexpr T ε     = 
    std::numeric_limits<T>::epsilon();
  // ** end of constants ******************************************************
} // ** end of namespace nl ***********

// ****************************************************************************
#endif // #ifndef NUMLI_HPP