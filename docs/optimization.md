# HyPhy Optimization Architecture

This document describes the optimization logic and architecture in HyPhy, centered around `_LikelihoodFunction::Optimize`.

## Entry Point: `_LikelihoodFunction::Optimize`

The `Optimize` method is the primary driver for likelihood-based parameter estimation in HyPhy. It supports several optimization modes, configurable via environment variables or an associative list of options.

### Supported Methods
- **Coordinate-Wise (0)**: Iterates through parameters one by one, performing a 1D line search for each.
- **Nedler-Mead Simplex (3)**: A derivative-free search algorithm.
- **Gradient Descent / Conjugate Gradient (7)**: Uses numerical gradients and the Polak-Ribière Conjugate Gradient method.
- **Hybrid (4, Default)**: A multi-stage approach that typically combines Gradient Descent (to jump-start or refine) with Coordinate-Wise refinement.

---

## 1. Stage-Based Hybrid Optimization

For most high-dimensional problems, HyPhy uses a **Hybrid** approach:

1.  **Initial Grid Search (Optional)**: If `LF_INITIAL_GRID` is provided, it evaluates the likelihood at specific points to find a good starting maximum.
2.  **Gradient Pass**: Performs a Conjugate Gradient pass to find the general vicinity of the maximum.
3.  **Coordinate-Wise Refinement**: Switches to coordinate-wise iteration to polish the estimates, especially for parameters at boundaries or those with complex constraints.
4.  **Final Gradient Pass**: A final CG pass at high precision to ensure convergence of all parameters simultaneously.

---

## 2. Coordinate-Wise Optimization Logic

The coordinate-wise optimizer is highly sophisticated and employs several heuristics to reduce the number of expensive likelihood evaluations:

### Adaptive Step and Precision
- **Precision Scheduling**: The algorithm starts with a coarse precision and refines it as it approaches the maximum.
- **Adaptive Bracketing**: Initial steps for the 1D search (`Bracket`) are estimated from the history of parameter changes in previous passes.
- **Shrink Factor**: The search interval is shrunk or expanded based on whether convergence is accelerating or slowing down.

### Variable Scheduling (Heuristics)
- **Large Change Only**: If a subset of variables is driving most of the likelihood improvement, the optimizer may skip unchanging variables for several passes to focus on the "active" parameters.
- **Proactive CG Trigger (Refined 2026)**: If convergence slows down (`convergenceMode >= 2`), the optimizer triggers a Conjugate Gradient pass. If the "Large Change Only" mode is active and the number of active variables is small (<= 32), CG is performed *only* on those variables. This targeted directional update is significantly more efficient than full gradient descent for models with many local variables.
- **Core Variables**: Identifies "core" variables that frequently drive changes and prioritizes them.
- **Variable Shuffling**: Optionally shuffles the order of variable updates to avoid getting stuck in narrow "valleys" on the likelihood surface.

### 1.D Line Search: `LocateTheBump`
- Performs a bracketing step followed by a Brent search.
- **Optimization (2026)**: Now utilizes all three bracket points (`left`, `middle`, `right`) to initialize the Brent search, allowing for an immediate parabolic fit and fewer evaluations.

---

## 3. Gradient-Based Optimization: `ConjugateGradientDescent`

The CG implementation handles multi-parameter updates simultaneously:

- **Numerical Gradients**: Gradients are computed using central or forward differences (`ComputeGradient`).
- **Parameter Partitioning (Gradient Blocks)**: Parameters can be grouped into independent blocks (e.g., branch lengths vs. substitution model parameters) to perform targeted CG passes on smaller subsets of variables.
- **Polak-Ribière (PR+)**: Uses the standard PR+ algorithm with periodic resets to the steepest descent direction for stability.
- **Line Search**: Uses `GradientLocateTheBump` along the conjugate direction.

---

## 4. Convergence Criteria

Optimization terminates when:
1.  The improvement in log-likelihood falls below the target `OPTIMIZATION_PRECISION`.
2.  The maximum iterations per variable (`MAXIMUM_ITERATIONS_PER_VARIABLE`) is reached.
3.  A hard time limit (`OPTIMIZATION_HARD_LIMIT`) is exceeded.
4.  A custom convergence callback (HBL function) returns true.

## 5. Implementation Details

- **Numerical Stability**: Includes checks for `NaN` and `-Infinity` likelihoods, with retry/reset logic.
- **MPI Parallelism**: Supports distributing likelihood evaluations across multiple nodes when compiled with MPI.
- **Parameter Mapping**: Maps bounded parameters (e.g., $[0, 1]$ or $[0, \infty)$) to an internal unbounded space (using expit or log transforms) for smoother gradient behavior.
