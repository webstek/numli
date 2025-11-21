# numli
Numerics Library containing various useful C++ structures and functions. Numli is a header-only library, but has separated implementation into modules. See [Table 1](#table-1-namespaces-and-dependencies) for dependencies.

## Structure

`numli` is divided into namespaces, with the top level being `nl`. Each module
is contained in its associated `.hpp` file and is accessible by its namespace.

#### *Table 1: Namespaces and Dependencies*
| Module            | namespace | required modules |
|-------------------|-----------|------------------|
| Numli             | `nl`      | None
| Algebra           | `nl::bra` | Numli
| Computer Graphics | `nl::cg`  | Numli, Algebra

## Tests

Unit tests using [doctest](https://github.com/doctest/doctest?tab=readme-ov-file) are defined in associated `.cpp` files in the test directory. Use `make` to build the test executables and `make clean` to remove test materials.

## Status

#### *Current:*
* Working on implementing SIMD compatible $\mathbb{R}^n$ and
$\mathbb{R}^{n\times m}$ data types.

#### *Complete:*
* `restrict` support for MSVC and GCC compilers
* Mathematical constants - π, e, ε (floating point epsilon)

#### *Future:*
* Minimal working state of the Algebra module
* Significant work on Computer Graphics module