# Numli Coding Style

```c++
// ** Preamble ****************************************************************
/// @file file.hpp
/// @author Kyle Webster
/// @version 1.0
/// @date dd Month year
/// @brief Brief description
/// @details
/// A longer description
// ****************************************************************************
#ifndef PATH_TO_FILE_HPP
#define PATH_TO_FILE_HPP
// ** Includes ****************************************************************
#include <>
#include ""
// ****************************************************************************

// ************************************
/// @namespace style
/// @brief description
/// @see also, perhaps dependencies
/// @details
/// Contents:
///  * @ref DataCollection
namespace style {

// ************************************
/// @struct DataCollection
/// @brief brief description
/// @tparam T buffer type
/// @see related data defintions
template<std::arithmetic T> struct DataCollection  // CamelCase
{
  std::vector<color> data_vector;  ///< data stored on the heap
  const char*        name;         ///< for identification

  /// @brief resize the buffer
  /// @param size the new number of elements
  /// Details or references of implementation here...
  void resizeBuffer(int size);
}
// ** DataCollection ******************
.
.
.
} // namespace style

// ****************************************************************************
#endif // #ifndef PATH_TO_FILE_HPP
```

Generally follows the LLVM C++ coding standard in terms of practice, but with
minor formatting differences. The use of unicode characters is encouraged to
improve code clarity - restrict to use of latin, greek, double-struck, and APL 
glyphs only.

```
Struct CamelCase   -> collection of data for an interface
Struct lowercase   -> data with behaviour
Class  CamelCase   -> interface
var_x  snake_case  -> variable
void   camelCase() -> function
const  ALL_CAPS    -> constant
```

## Structure

- Include all implementation in a single `.hpp` file per module
- Use concepts to enforce template requirements heavily
- Use anonymous namespaces to encapsulate unexported functions in an interface
- Prefer stuctures of arrays over arrays of structures
- Flatten trees to avoid pointer chasing