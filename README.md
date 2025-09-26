# numli
Numerics Library containing various useful C++ structures and functions. Numli
is a header-only library, but has separated implementation into modules. See [Table 1](#table-1-namespaces-and-dependencies) for dependencies.

## Structure

`numli` is divided into namespaces, with the top level being `nl`. Each module
is contained in its associated `.hpp` file and is accessible by its namespace.

#### *Table 1: Namespaces and Dependencies*
| Module            | namespace | required modules |
|-------------------|-----------|------------------|
| Algebra           | `nl::bra` | None
| Computer Graphics | `nl::cg`  | Algebra

