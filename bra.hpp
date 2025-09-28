// ****************************************************************************
/// @file bra.hpp
/// @author Kyle Webster
/// @version 0.1
/// @date 27 Sep 2025
/// @brief Numerics Library - Algebra @ref bra
/// @details
/// Collection of algebraic structures and algorithms
// ****************************************************************************
#ifndef BRA_HPP
#define BRA_HPP
// ** Includes ****************************************************************
#include "numli.hpp"
#include <concepts>
#include <cstdint>
#include <stdfloat>
// ****************************************************************************

namespace nl {
// ************************************
/// @namespace bra
/// @brief Algebraic data structures and algorithms
/// @details
/// Contents:
///  * @ref ℝn - number in ℝⁿ, specifiable floating point storage type
///  * @ref ℝnxm - number in ℝⁿˣᵐ, column-major order floating point storage
namespace bra
{ // nl::bra
  // **********************************
  /// @struct ℝn
  /// @brief a number in ℝⁿ
  /// @tparam n number of basis elements that square to 1
  /// @tparam T storage type, default float64_t
  /// @attention T must be floating point
  /// @details
  /// Element of ℝⁿ stored as type T. If AVX2 is available it is stored aligned
  /// to 32 Bytes (AVX2 register size).
  template<uint32_t n, std::floating_point T = std::float64_t> 
  struct ℝn
  {
    /// @name members
    alignas(nl::simd<T>::alignment) T x[n]; ///< basis element coefficients

    /// @name constructors
    explicit ℝn( T const * restrict a ) { memcpy(x, a.x, n*sizeof(T)); }
    explicit ℝn( T v ) { for (uint8_t i=0; i<n; ++i) x[i]=v; }
    explicit ℝn( ℝn<n,T> const &a ) { memcpy(x, a.x, n*sizeof(T)); }

    /// @name operators
  };
  // ** ℝn ****************************

  // **********************************
  /// @struct ℝnxm
  /// @brief a number in ℝⁿˣᵐ
  /// @tparam n number of rows
  /// @tparam m number of columns
  /// @tparam T floating point storage type
  /// @attention T must be one of the std::floatXX_t types
  /// @details
  /// Element of ℝⁿˣᵐ stored as type T. If AVX2 is available it is stored 
  /// aligned to 32 Bytes.
  template<uint32_t n, uint32_t m, std::floating_point T = std::float64_t>
  struct ℝnxm
  {
    /// @name members
    constexpr uint64_t N = n*m;            ///< number of elements
    alignas(nl::simd<T>::aligment) T x[N]; ///< basis element coefficients
  };
  
  // ** Old Implementation ****************************************************
  template<uint32_t n, typename T = double>
  class NdReal {
    static_assert(std::is_arithmetic<T>::value, "T must arithmetic.");
    public:
    // Members:
      T x[n];

    // Functions:
      /// @brief Initializes all components to 0. 
      NdReal(void) {
        for (uint32_t i = 0; i < n; i++) {
          x[i] = 0.0;
        }
      }

      /// @brief Initializes to the values in y
      /// @param y 
      NdReal(T y[n]) {
        for (uint32_t i = 0; i < n; i++) {
          x[i] = y[i];
        }
      }

      // TODO: arithmetic operators
  };

  /// @brief A dynamically allocated vector containing elements of type T.
  /// @tparam T The type of the contents of the vector.
  template<typename T = double>
  class Vec {
  public:
  // Functions:
    /// @brief Default constructor setting all components to zero.
    Vec(uint32_t n) : n_(n) {
      pelems_ = new T[n_];
      for (uint32_t i = 0; i < n_; i++) {
        pelems_[i] = 0.0;
      }
    }

    /// @brief Sets components equal to the provided values.
    /// @param u n element array of objects of type T.
    Vec(uint32_t n, T *p_u) : n_(n) {
      pelems_ = new T[n_];
      for (uint32_t i = 0; i < n_; i++) {
        pelems_[i] = p_u[i];
      }
    };

    /// @brief Deallocate memory for the elements
    ~Vec(void) {
      delete[] pelems_;
    }

    // TODO: arithmetic operators

    /// @brief Gets the ith component.
    /// @param i 0-indexed component to get.
    /// @return The value of the ith component.
    T& operator [](const uint32_t &i) const {
      return pelems_[i];
    }

    /// @brief Gives the size of the vector.
    /// @return The number of elements of the vector.
    uint32_t size(void) {return n_;}

  private:
  // Members:
    uint32_t n_;            /**< Number of elements */
    T *pelems_ = nullptr;   /**< Pointer to array of elements */

  // Functions:
    // TODO: resizing, pushing elements, etc.
  };


  /// @brief Dynamical matrix of elements of type T in column-major format.
  /// @tparam T Type of the elements;
  template<typename T = double> class Mat {
    public:
    // Functions:
      /// @brief Empty size 1 by 1 matrix
      Mat() : cols_(0), rows_(0), col_cap_(1), row_cap_(1) {
        pflat_a_ = new T[1];
      }

      /// @brief Empty m by n matrix
      /// @param m The number of rows to put data in.
      /// @param n The number of columns to put data in.
      /// @throws @c runtime_error if @p m or @p n is zero.
      Mat(uint32_t m, uint32_t n) : col_cap_(n), row_cap_(m) {
        if (m == 0 || n == 0) throw std::runtime_error("m nor n can be 0!");
        pflat_a_ = new T[col_cap_ * row_cap_];
      }

      ~Mat(void) {
        delete[] pflat_a_;
      }

      /// @brief Appends a column to the right side of the matrix;
      /// @param column Values to extend the matrix with.
      void pushCol(T *p_column) {
        if (cols_ == col_cap_) upsizeCols();
        for (uint32_t j = 0; j < rows_; j++) {
          pflat_a_[cols_+j*col_cap_] = p_column[j];
        }
        cols_++;
      }

      // TODO: pushRow, popCol, popRow
      
      /// @brief Gets the number of rows.
      /// @return The number of rows of the matrix.
      uint32_t rows(void) const {return rows_;}

      /// @brief Gets the number of columns.
      /// @return The number of columns of the matrix.
      uint32_t cols(void) const {return cols_;}

    private:
    // Members:
      uint32_t cols_;       /**< Number of columns */
      uint32_t rows_;       /**< Number of rows */
      uint32_t col_cap_;    /**< Allocated column size */
      uint32_t row_cap_;    /**< Allocated row size */
      T *pflat_a_;          /**< Pointer to  */

    // Functions
      /// @brief Up-size the matrix over the columns.
      void upsizeCols(void) {
        col_cap_ *= 2;

        // copy elements to new array
        T *pflat_a_dest = new T[col_cap_ * row_cap_];
        for (uint32_t i = 0; i < cols_; i++) {
          for (uint32_t j = 0; j < rows_; j++) {
            pflat_a_dest[i+col_cap_*j] = pflat_a_[i+col_cap_*j];
          }
        }
        delete[] pflat_a_;
        pflat_a_ = pflat_a_dest;
      }

      // TODO upsizeRow, downsizeCol, downsizeRow
  };


// Functions:
  /// @brief The square of the l2 norm of a vector.
  /// @tparam T Arithmetic type of elements.
  /// @param u Vector to take the norm of.
  /// @return The l2 norm of u squared.
  template<typename T>
  T l2Squared(Vec<T> u) {
    T val = 0;
    for (size_t i = 0; i < u.size(); i++) {
      val += u[i]*u[i];
    }
    return val;
  }

  /// @brief The l2 norm of a vector.
  /// @tparam T Arithmetic type of elements.
  /// @param u Vector to take the norm of.
  /// @return The l2 norm of u.
  template<typename T>
  T l2(Vec<T> u) {return sqrt(l2Squared(u));}


} // namespace bra




} // namespace nl

#endif // ifndef BRA_HPP