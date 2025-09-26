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
#include ""
#include <>
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
/// @see related data defintions
struct DataCollection  // CamelCase
{
  std::vector<color> data_vector;  ///< data stored on the heap
  const char*        name;         ///< for identification

  /// @brief resize the buffer
  /// @param size the new number of elements
  /// Details or references of implementation here...
  void resizeBuffer(int size);
}
.
.
.
} // namespace style

//===---------------------------------------------------------------------===//
#endif // #ifndef PATH_TO_FILE_HPP
```

Generally follows the LLVM C++ coding standard in terms of practice, but with
minor formatting differences.

```
Struct CamelCase   -> collection of data for an interface
Struct lowercase   -> data with behaviour
Class  CamelCase   -> interface
var_x  snake_case  -> variable
void   camelCase() -> function
```

## Structure

- Separate files into `.hpp` header files and `.cpp` implementation files
- Use anonymous namespaces to encapsulate functions not exported by the interface
- Prefer stuctures of arrays over arrays of structures
- Flatten trees to avoid pointer chasing