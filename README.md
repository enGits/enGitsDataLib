# enGitsDataLib

enGitsDataLib (EDL) is a C++ library providing data structures and utilities for geometry, interpolation, and numerical software. It is developed by [enGits](https://www.engits.com/en/).

EDL supports Windows, Linux, and macOS. It requires a C++17 compiler and CMake 3.22 or newer.

## Build and install

```sh
cmake -S . -B build -DBUILD_TESTING=OFF
cmake --build build --config Release
cmake --install build --config Release --prefix <install-prefix>
```

Consume the installed CMake package with:

```cmake
find_package(enGitsDataLib 1.0 CONFIG REQUIRED)
target_link_libraries(my_target PRIVATE enGitsDataLib::engitsdatalib)
```

Pass the installation location when configuring the consumer if necessary:

```sh
cmake -S . -B build -DCMAKE_PREFIX_PATH=<install-prefix>
```

## Legacy installation prefixes

EDL no longer assumes a fixed installation prefix. Existing software that expects `/local` or `C:\local` can continue using that layout by specifying it explicitly with `cmake --install --prefix`. New consumers should use the exported CMake package targets instead of hard-coded include and library paths.

```sh
cmake --install build --config Release --prefix "C:\local"
# or on Linux/macOS:
cmake --install build --config Release --prefix "/local"
```

## vcpkg

A vcpkg port is prepared for v1.0.0. It will be finalized and submitted to the upstream registry after the release tag is published. Until then, use the normal CMake installation.

## Licence

EDL is released under the [MIT License](LICENSE).
