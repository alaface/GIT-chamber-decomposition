# README

## Description

This repository contains a collection of functions for working with cones, grading matrices, and GIT fans in the context of algebraic geometry. These functions are implemented in the Magma programming language and are designed to handle computations related to effective cones, moving cones, orbit cones, and GIT chambers.

## Functions Overview

1. **`Ffaces`**: Computes the F-faces of an ideal or the subsets of integers.
2. **`Eff`**: Computes the effective cone of a grading matrix.
3. **`Mov`**: Computes the moving cone of a grading matrix.
4. **`OrbitCones`**: Constructs orbit cones based on F-faces and a grading matrix.
5. **`GitChamber`**: Finds the GIT chamber containing a given class.
6. **`BunchCones`**: Returns the bunch of cones containing a class.
7. **`SameSbl`**: Checks if two classes have the same stable base locus.
8. **`GitFan`**: Constructs the GIT fan from orbit cones.

## Requirements

- Magma computational algebra system.

## Usage

To use these functions, load them in your Magma session. Below is an example demonstrating their usage.

## Example

### Input

Consider the following grading matrix `Q` and an integer `I`:

```magma
Q := Matrix(Rationals(), 2, 4, [1, 1, 0, 3, 0, 0, 1, 1]);
I := 4;
```

### Compute Effective Cone

```magma
eff_cone := Eff(Q);
print "Effective Cone:", eff_cone;
```

### Compute Moving Cone

```magma
mov_cone := Mov(Q);
print "Moving Cone:", mov_cone;
```

### Compute F-Faces

```magma
f_faces := Ffaces(I);
print "F-Faces:", f_faces;
```

### Orbit Cones

```magma
orbit_cones := OrbitCones(f_faces, Q);
print "Orbit Cones:", orbit_cones;
```

### GIT Fan

```magma
git_fan := GitFan(orbit_cones);
print "GIT Fan:", git_fan;
```

## Notes

- The `Ffaces` function can accept both integers and ideals as input.
- Ensure all inputs, such as matrices or ideals, are correctly defined before calling the functions.

## License

This project is licensed under the MIT License. See the LICENSE file for details.

