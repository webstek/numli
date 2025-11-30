# numli
Numerics Library containing various useful C++ structures and functions. Numli is a header-only library, but has separated implementation into modules. See [Table 1](#table-1-namespaces-and-dependencies) for dependencies.

## Structure

`numli` is divided into namespaces, with the top level being `nl`. Each module
is contained in its associated `.hpp` file and is accessible by its namespace.

#### *Table 1: Namespaces and Dependencies*
| Module            | namespace | required modules |
|-------------------|-----------|------------------|
| Numli             | `nl`      | None
| SIMD              | `nl::simd`| Numli
| Algebra           | `nl::bra` | Numli, SIMD
| Computer Graphics | `nl::cg`  | Algebra

## Tests

Unit tests using [doctest](https://github.com/doctest/doctest?tab=readme-ov-file) are defined in associated `.cpp` files in the test directory. Use `make` to build the test executables and `make clean` to remove test materials.

## Status

#### *Recently Complete:*
* Minimal working CG module
* Minimal `.nls` file loading
* Minimal working ℝnxm

#### *Future:*
* nl::stoch library for stochastic process related work
* Significant work on Computer Graphics module