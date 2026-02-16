# 3D Exact Cover Solver & Generator (Rust)

A high-performance **3D exact cover puzzle solver and generator** written in **Rust**, featuring both a state-of-the-art research-based solver (**Dancing Links / Algorithm X**) and a fully original brute-force solver developed independently from first principles.

This project explores the algorithmic, mathematical, and systems-level challenges behind solving and generating **NP-complete** 3D packing problems, with a strong emphasis on performance, memory efficiency, and clean low-level implementation.

---

##  Problem Overview

The core problem solved here is a **3D exact cover problem**:

> Given a 3D grid with obstacles and a set of polycube blocks, determine whether the blocks can be placed (with arbitrary rotations) such that every empty cell is covered *exactly once*.

This class of problems is:
- **NP-complete**
- Highly sensitive to branching factor
- Dominated by combinatorial explosion
- Extremely reliant on pruning, symmetry reduction, and cache-friendly data structures

---

##  Key Features

- **Two independent solvers**
  - A research-backed **DLX / Algorithm X** implementation
  - A fully original, hand-optimized brute-force solver
- **3D block rotation system** with symmetry elimination
- **Exact cover matrix construction** for 3D placement constraints
- **Procedural puzzle generator** using flood-fill and randomized growth
- Written entirely in **safe Rust**, emphasizing performance and correctness

---

##  Solvers

### 1. Dancing Links Solver (`dlx.rs`)

This solver implements **Donald Knuth’s Algorithm X** using the **Dancing Links (DLX)** data structure—the gold standard for exact cover problems.

#### Highlights:
- Sparse exact cover matrix encoded as a **toroidal doubly-linked list**
- O(1) column and row removal / restoration
- **Minimum column heuristic** to minimize branching factor
- Depth-first recursive search with backtracking
- Exact modeling of:
  - One block used exactly once
  - One block occupies each empty cell exactly once

#### Concepts & Techniques:
- Dancing Links (DLX)
- Algorithm X
- Exact cover constraint encoding
- Backtracking with structural undo
- Constraint satisfaction via sparse matrix manipulation

This implementation focuses heavily on **memory locality** and **low-level efficiency**, storing the DLX structure in flat vectors rather than heap-allocated nodes.

---

### 2. Original Brute-Force Solver (`solver.rs`)

This solver was developed **entirely from scratch**, without reference to existing literature, and refined through multiple full rewrites.

While theoretically exponential, it performs surprisingly well due to aggressive pruning and lazy evaluation strategies.

#### Highlights:
- Explicit **depth-first search** over block placements
- Lazy generation of:
  - Rotational states
  - Placement coordinates
- Incremental collision and bounds checking
- Early termination on invalid partial placements
- Extensive symmetry and redundancy elimination

#### Concepts & Techniques:
- DFS backtracking
- Lazy evaluation
- Combinatorial pruning
- Rotation group enumeration
- Constraint propagation through partial state validation

This is the latest version of the very first program that started the entire project. It went through an iterative development process that began with a simple Python prototype and evolved through multiple Rust rewrites. Each iteration focused on improving clarity, performance, and architectural structure, with successive optimizations informed by profiling, refactoring, and growing problem understanding.

---

##  Puzzle Generator (`generator.rs`)

The generator creates **new, valid, non-trivial 3D exact cover instances** from scratch.

Rather than relying on predefined shapes, it procedurally constructs blocks using spatial growth algorithms.

#### Highlights:
- Randomized **flood-fill–style block growth**
- Breadth-first expansion using frontier queues
- Connectivity guarantees for all generated blocks
- Automatic volume partitioning
- Field-aware block generation that respects obstacles

#### Concepts & Techniques:
- Flood fill / BFS
- Randomized growth processes
- Spatial connectivity constraints
- Discrete volume partitioning
- Procedural generation

The generator ensures that produced puzzles are both solvable and structurally interesting, making it suitable for benchmarking and experimentation.

---

## 🛠 Architecture Overview

- **Field abstraction** for 3D grid operations
- **Rotator system** for enumerating unique block orientations
- Shared geometric primitives for coordinates and blocks
- Separation of:
  - Geometry
  - Search logic
  - Constraint modeling

The project is intentionally modular, allowing different solvers or generators to be swapped or compared.

---

##  Why This Project Matters

This repository demonstrates:
- Mastery of **advanced algorithms and data structures**
- Comfort with **NP-complete problems**
- Ability to implement **theoretical research algorithms** and *develop* **original solutions** from first-principles
- Strong understanding of **performance-oriented Rust**
- Persistence in iterative optimization and refactoring

It reflects real-world problem solving: starting from brute force, discovering bottlenecks, inventing optimizations, and ultimately implementing a state-of-the-art solution.

---

##  Future Work

- Cleaning up the code (I was fairly new to idiomatic Rust at the time of writing, and optimized for development speed instead of idiomatic code. I could significantly improve the code with what I've learned since).
- Parallelized DLX search (Would be really cool to learn more about threads in Rust!)
- Heuristic-driven block ordering
- Difficulty metrics for generated puzzles
- Visualization tooling (I've already played around with this using Tauri and ThreeJS)
- SAT / SMT solver comparison

---

##  Notes

This project was built as a learning exercise in algorithmic problem solving, low-level optimization, 3D computational geometry, algorithm-heavy systems work, constraint solving and performance-critical Rust. No external solver libraries were used.
