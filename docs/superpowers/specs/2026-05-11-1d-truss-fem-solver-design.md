# 1D Truss FEM Solver — Design Spec

## Overview

A 1D planar truss finite element solver that reads a JSON model file, assembles and solves the global stiffness equations using Gaussian elimination, and prints nodal displacements, element stresses, and support reactions to the console.

## Technology Stack

- **Language:** C++20
- **Build:** CMake
- **Test:** Google Test
- **Input:** JSON file
- **Output:** Console (stdout)

## Pipeline

```
JSON Input → Parse → Assemble K,F → Apply BCs → Solve K·u=F → Post-process → Console Output
```

1. **Parse** — Read JSON, build `TrussModel` with nodes, elements, BCs, loads
2. **Assemble** — Compute element stiffness matrices, sum into global `K` and `F`
3. **Apply BCs** — Row/column elimination on constrained DOFs
4. **Solve** — Gaussian elimination with back-substitution → displacement vector `u`
5. **Post-process** — Compute stresses, reactions, print formatted results

## Input Format (JSON)

```json
{
  "nodes": [
    { "id": 0, "x": 0.0,  "y": 0.0 },
    { "id": 1, "x": 1.0,  "y": 0.0 },
    { "id": 2, "x": 0.5,  "y": 1.0 }
  ],
  "elements": [
    { "id": 0, "n1": 0, "n2": 2, "E": 210e9, "A": 0.001 },
    { "id": 1, "n1": 1, "n2": 2, "E": 210e9, "A": 0.001 }
  ],
  "boundary_conditions": [
    { "node": 0, "dof": "ux", "value": 0.0 },
    { "node": 0, "dof": "uy", "value": 0.0 },
    { "node": 1, "dof": "uy", "value": 0.0 }
  ],
  "loads": [
    { "node": 2, "dof": "uy", "value": -10000.0 }
  ]
}
```

- `dof`: `"ux"` or `"uy"`
- Units: SI (meters, Newtons, Pascals) — convention, not enforced

## Output Format

```
=== Nodal Displacements ===
Node 0: ux =  0.000000e+00, uy =  0.000000e+00
Node 1: ux =  1.234567e-05, uy =  0.000000e+00
...

=== Element Stresses ===
Element 0: σ =  1.500000e+07 Pa
...

=== Support Reactions ===
Node 0: Fx =  5.000000e+03, Fy =  8.660000e+03
...
```

## Architecture

### Class Hierarchy

**Value types (POD-like structs):**
- `Node` — id, x, y, computed displacement (ux, uy)
- `Element` — id, node refs (n1, n2), E, A, L; methods: `computeLocalStiffness()`, `computeGlobalStiffness()`, `computeStress()`, `computeStrain()`
- `BoundaryCondition` — nodeId, dof (enum `DoF { UX, UY }`), value
- `NodalLoad` — nodeId, dof, value

**Aggregate:**
- `TrussModel` — owns vectors of all domain objects; `totalDofs()`, `dofIndex(nodeId, dof)`

**Algorithm classes (stateless or single-purpose):**
- `Assembler` — `assembleStiffness(model) → Matrix`, `assembleForceVector(model) → Vector`
- `BoundaryConditionHandler` — `apply(K, F, model)` using row/column elimination
- `Solver` (abstract) — `solve(K, F) → Vector`
  - `GaussianSolver` (concrete) — Gaussian elimination + back-substitution
- `JsonParser` — `parse(filepath) → TrussModel`
- `PostProcessor` — `printResults(model, displacements, reactions)`

### DOF Indexing

Node `i` → DOF `2i` (ux), DOF `2i+1` (uy). Total DOFs = `2 × |nodes|`.

Global stiffness matrix `K` is `n×n` dense, symmetric. Force and displacement vectors are length `n`.

### Main Orchestration

```cpp
TrussModel model = JsonParser::parse(inputFile);
Matrix K = Assembler::assembleStiffness(model);
Vector F = Assembler::assembleForceVector(model);
BoundaryConditionHandler::apply(K, F, model);
Solver* solver = new GaussianSolver();
Vector u = solver->solve(K, F);
PostProcessor::printResults(model, u);
```

## Row/Column Elimination (BC Application)

For each constrained DOF `d` with prescribed value `v`:
1. For all `j`: `K[d][j] = 0`, `K[j][d] = 0`
2. `K[d][d] = 1`
3. `F[d] = v`
4. For all `j ≠ d`: `F[j] -= K[j][d] * v`

## Error Handling

**Tier 1 — Hard errors (throw → catch in main, stderr, non-zero exit):**
- JSON parse failure
- Missing required fields
- Node ID out of range in element/BC/load
- Element with n1 == n2 (zero-length)
- E ≤ 0 or A ≤ 0
- Duplicate BC on same DOF

**Tier 2 — Warnings (cerr, continue if possible):**
- Near-zero pivot during Gauss elimination (mechanism detected)
- Unconstrained DOFs (rank deficiency)

## Testing Strategy

| Component | Tests |
|---|---|
| `Element` | Local stiffness at 0° and 45°, stress/strain from displacement |
| `Assembler` | 2-bar model: K dimensions, symmetry, correct assembly |
| `BoundaryConditionHandler` | Verify row/col zeroed, diagonal=1, F adjusted |
| `GaussianSolver` | Known 2×2 and 4×4 systems, zero-diagonal detection |
| `JsonParser` | Valid file → model, missing field → throw, bad ref → throw |
| **Integration** | 2-bar truss vs. analytic solution, 3-bar indeterminate truss |

Test fixtures: JSON files in `tests/fixtures/`.

## Project Structure

```
truss-fem-solver/
├── CMakeLists.txt
├── src/
│   ├── Node.h, Node.cpp
│   ├── Element.h, Element.cpp
│   ├── BoundaryCondition.h
│   ├── NodalLoad.h
│   ├── TrussModel.h, TrussModel.cpp
│   ├── Matrix.h, Matrix.cpp
│   ├── Assembler.h, Assembler.cpp
│   ├── BoundaryConditionHandler.h, BoundaryConditionHandler.cpp
│   ├── Solver.h
│   ├── GaussianSolver.h, GaussianSolver.cpp
│   ├── JsonParser.h, JsonParser.cpp
│   └── PostProcessor.h, PostProcessor.cpp
├── src/main.cpp
├── tests/
│   ├── CMakeLists.txt
│   ├── test_element.cpp
│   ├── test_assembler.cpp
│   ├── test_bc_handler.cpp
│   ├── test_gaussian_solver.cpp
│   ├── test_json_parser.cpp
│   ├── test_integration.cpp
│   └── fixtures/
│       ├── two_bar_truss.json
│       └── three_bar_truss.json
└── docs/superpowers/specs/
```

## Out of Scope (v1)

- Nonlinear geometry or material
- Dynamic analysis (no mass matrix, no eigenvalues)
- Distributed loads (point loads only)
- Multiple load cases
- Sparse matrix storage
- File output (console only)
- CLI flags beyond required input path
