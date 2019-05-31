## Pull Request Guidelines

- Squash the commit if there are too many small ones.

- Follow the [code style](#code-style).

- Make sure `HYPHYGTEST` and HBL tests in `tests/hbltests/libv3` passes. 

- If adding new feature:
    - Add accompanying test case.

- If fixing a bug:
    - Provide detailed description of the bug in the pull request.
    - Add appropriate test coverage if applicable.

## Code Style

### Variables

All lower case with underscores separating words, as in `loop_indexer`

#### Global variables 

Prefix global variables with `hy_` in C++ (e.g. `hy_global_list`), and with namespace id in HBL (e.g. `io.global_variable`)

#### Private/aux variables

Prefix variables that are meant to be private or auxiliary with an underscore, as in `_hy_ignore_me`

### Functions, classes  

Capital case, as in `ComputeExpression()` or `VariableContainer` (class). Prefix the names of *auxiliary* classes and types with an `_`.

### Simple types

Types such as enums, structs, and typedefs should prefixed with `hy` and followed by capital case identifiers, e.g. `hyFloat`

#### Cheap functions

If a function is *cheap* to call (i.e. don't worry about caching the result), name it as a variable, e.g. `member_getter()`. 

### Constants, enums (C++ only)

Capital case prefixed with a lower case `k`, as in `kHarvestFrequenciesResolveAmbiguities`

### Bit-flags 

Capital case prefixed with `f`, as in `fSearchCaseSensitive`

### Environment variables (HBL only)

All capitals with underscore-separated words, as in `LIKELIHOOD_FUNCTION_EXPORT_FORMAT`
