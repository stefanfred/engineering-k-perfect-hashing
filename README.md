# Engineering Minimal k-Perfect Hash Functions

A k-perfect hash function is a hash function that maps a set of keys to a set of integers.
Each output value might collide with up to k other keys.

This repository contains a selection of k-perfect hash function implementations.
It also includes benchmarks comparing them to other approaches from the literature.

### Library Usage

Clone this repository (with submodules) and add the following to your `CMakeLists.txt`.

```
add_subdirectory(path/to/kPerfectHashing)
target_link_libraries(YourTarget PRIVATE Kphf::<hash function to use>)
```

You can link with one of the following targets.
Each target is a different k-perfect hash function with its own set of dependencies.

- `Kphf::PaCHash`
- `Kphf::RecSplit`
- `Kphf::HashDisplace`
- `Kphf::ThresholdBasedBumping`
