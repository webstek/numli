// ****************************************************************************
/// @file numli.hpp
/// @author Kyle Webster
/// @version 0.1
/// @date 27 Sep 2025
/// @brief Numerics Library - @ref nl
/// @details
/// Collection of constants and core numerics
// ****************************************************************************
#ifndef NUMLI_HPP
#define NUMLI_HPP
// ** Includes ****************************************************************
#include <numbers>
#include <cstdint>
#include <stdfloat>
// ****************************************************************************

// ************************************
/// @namespace nl
/// @brief Numerics Library main module
/// @details
/// Contents:
///  * @ref simd - simd details
///  * @ref constants - mathematical constants
namespace nl
{
  // ** SIMD ******************************************************************
  /// @name SIMD
  /// @brief simd support information/configuration
  #if defined(__AVX2__)
    constexpr bool avx2 = true;
  #else
    constexpr bool avx2 = false;
  #endif
  // **********************************
  /// @struct simd
  /// @brief simd information
  /// @tparam T storage type for alignment if avx2 is unavailable
  template<typename T> struct simd
  {
    static constexpr size_t alignment = avx2 ? 32B : sizeof(T);
  };
  // ** end of SIMD ***********************************************************

  // **************************************************************************
  /// @name constants
  /// @brief Mathematical and computational constants
  template<typename T> inline constexpr T π     = std::numbers::pi_v(T);
  template<typename T> inline constexpr T π²    = 
    std::numbers::pi_v(T)*std::numbers::pi_v(T);
  template<typename T> inline constexpr T π_2   = std::numbers::pi_v(T)/2;
  template<typename T> inline constexpr T inv_π = 1/std::numbers::pi_v(T); 
  template<typename T> inline constexpr T e     = std::numbers::e_v(T);
  template<typename T> inline constexpr T ε     = 
    std::numeric_limits<T>::epsilon();
  // ** end of constants **************
} // namespace nl


// ****************************************************************************
#endif // #ifndef NUMI_HPP