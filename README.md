# Genetic Perturbation Simulator

A Python tool for simulating the downstream effects of a single gene perturbation (overexpression or knockout) across a weighted, directed biological network. Built on a **recursive, depth-limited path summation model** using the `igraph` library.

- **Additive influence** from multiple paths converging on a target
- **Feedback loop amplification** via cyclic structures in the graph
- **Tunable depth** to distinguish short-term response from long-term steady-state

---

## The Model

### Biological Motivation

Gene regulatory networks can be modelled as directed graphs where nodes are genes and weighted edges represent regulatory relationships (activation or repression). When a gene's expression level is artificially perturbed — through CRISPR knockout, RNAi knockdown, or overexpression — the effect propagates downstream through its regulatory targets, their targets, and so on.

The question this tool answers is: **given a perturbation at a single source gene, what is the predicted change in expression for every other gene in the network?**

### Recursive Path Summation

The core idea is that the total influence a source gene $S$ exerts on a target gene $T$ is the **sum of contributions from all directed walks** from $S$ to $T$, up to a maximum depth $d$.

Each walk's contribution is the **product of edge weights** along that walk:

$$\text{Influence}(S \to T) = \sum_{\text{walks } p \,:\, S \rightsquigarrow T,\; |p| \leq d} \;\prod_{e \,\in\, p} w_e$$

where $w_e$ is the signed weight of edge $e$, and $|p|$ is the number of edges in walk $p$.

The final expression change at gene $T$ is then:

$$\Delta E_T = \Delta E_S \cdot \text{Influence}(S \to T)$$

where $\Delta E_S$ is the initial perturbation (e.g., $+0.5$ for an overexpression from baseline $1.0$ to $1.5$).

### Why Walks, Not Paths?

This model enumerates **walks** (vertices may repeat), not simple paths (vertices visited at most once). This is deliberate — in a biological network, a signal *does* cycle through feedback loops multiple times. A walk that traverses the loop $4 \to 3 \to 4 \to 3 \to \ldots$ represents the biological reality of repeated mutual activation between those genes.

### Feedback Loops and Convergence

Consider a loop between genes $A$ and $B$ with edge weights $w_{A \to B}$ and $w_{B \to A}$. The **loop gain** is:

$$g = w_{A \to B} \cdot w_{B \to A}$$

Each time the signal traverses the full loop, its magnitude is multiplied by $g$. Over many iterations the total amplification forms a geometric series:

$$1 + g + g^2 + g^3 + \cdots = \frac{1}{1 - g} \quad \text{(for } |g| < 1\text{)}$$

In the example network, the loop between Gene 4 and Gene 3 has gain $g = 0.6 \times 0.4 = 0.24$. The theoretical amplification factor is $\frac{1}{1 - 0.24} \approx 1.316$. Setting `max_depth=100` allows the simulation to approach this steady-state; setting `max_depth=2` captures only the immediate, single-pass response.

When $|g| \geq 1$, the loop is unstable and expression levels diverge — the `max_depth` parameter acts as a natural cutoff preventing infinite recursion in such cases.

### The Algorithm

The implementation uses depth-first recursion from the source gene. For each target gene $T$:

1. Start at the source vertex $S$ with `current_path_product = 1.0`
2. For each successor $N$ of the current vertex:
   - Multiply the path product by $w_{\text{current} \to N}$
   - If $N = T$, add the accumulated product to the running total
   - Regardless, recurse into $N$ (incrementing depth) to find longer walks through $T$
3. Stop when `current_depth > max_depth`

This yields the total influence factor, which is multiplied by $\Delta E_S$ to produce $\Delta E_T$.

---

## The Example Network

The included example uses the following 8-gene network:

```
          +-----+         +-----+
          |  2  |         |  5  |
          +-----+         +-----+
             |               ^ \
             | 1.0      0.6 /   \ 1.0
             v             /     v
+-----+   +-----+ <---+ +-----+   +-----+
|  3  |<--|  4  |      | |  6  |-->|  7  |
+-----+   +-----+      | +-----+   +-----+
  |  ^      ^  |        |   ^         
  |  |      |  | 0.3    |   | 0.7
  |  | 0.4  |  v        |   |
  |  +------+ +-----+   |   |
  |   0.6     |  8  |   |   |
  |            +-----+   |   |
  |                      |   |
  +-----------[1]--------+---+
              source
```

**Edge list:**

| Source | Target | Weight |
|:---:|:---:|:---:|
| 1 | 4 | 1.0 |
| 1 | 5 | 0.6 |
| 1 | 6 | 0.7 |
| 2 | 4 | 1.0 |
| 3 | 4 | 0.4 |
| 4 | 3 | 0.6 |
| 4 | 8 | 0.3 |
| 5 | 6 | 1.0 |
| 6 | 7 | 0.8 |

The perturbation is an overexpression of Gene 1 from $1.0$ to $1.5$ ($\Delta E_1 = +0.5$).

---

## How to Run

#### 1. Install Dependencies

```bash
pip install igraph
```

If you are using the included Dev Container (e.g. via GitHub Codespaces), all dependencies are pre-installed.

#### 2. Run the Simulation

```bash
python gene_cascade.py
```

#### 3. Example Output

```
--- Simulation (Short-Term Response, Depth=2) ---
Gene     Initial  New Level  Change (ΔE)
----------------------------------------
1           1.00       1.50        +0.50
2           1.00       1.00        +0.00
3           1.00       1.30        +0.30
4           1.00       1.50        +0.50
5           1.00       1.30        +0.30
6           1.00       1.65        +0.65
7           1.00       1.28        +0.28
8           1.00       1.15        +0.15

--- Simulation (Long-Term/Steady-State, Depth=100) ---
(Note the larger change in Genes 3, 4, and 8 due to the loop)
Gene     Initial  New Level  Change (ΔE)
----------------------------------------
1           1.00       1.50        +0.50
2           1.00       1.00        +0.00
3           1.00       1.39        +0.39
4           1.00       1.66        +0.66
5           1.00       1.30        +0.30
6           1.00       1.65        +0.65
7           1.00       1.52        +0.52
8           1.00       1.20        +0.20
```

#### Understanding the Results

In the **short-term** simulation (`max_depth=2`), only direct and two-hop walks are considered. Gene 4 receives influence solely from the direct edge $1 \to 4$ with weight $1.0$, giving $\Delta E_4 = 0.5 \times 1.0 = 0.50$.

In the **long-term** simulation (`max_depth=100`), the signal cycles through the $4 \leftrightarrow 3$ feedback loop many times. The geometric amplification factor of $\frac{1}{1 - 0.24} \approx 1.316$ increases Gene 4's change to $\approx 0.66$ and Gene 3's to $\approx 0.39$. Gene 8, downstream of Gene 4, also increases proportionally.

---

## Usage as a Module

The script can be imported directly into your own projects or Jupyter notebooks.

- **`compute_gene_effect(...)`** — Main function. Builds the graph, runs the recursive simulation, and returns final expression levels.
- **`get_total_influence(...)`** — Recursive helper that performs depth-limited walk enumeration from source to target.

#### Example

```python
from gene_cascade import compute_gene_effect

# Define the network as (source, target, weight) tuples
my_edges = [
    ('1', '4', 1.0),
    ('1', '5', 0.6),
    ('4', '3', 0.6),
    ('5', '6', 1.0),
]

# Baseline expression levels
my_initial_levels = {
    '1': 1.0,
    '2': 1.0,
    '3': 1.0,
    '4': 1.0,
    '5': 1.0,
    '6': 1.0,
}

# Simulate overexpression of Gene 1 (long-term steady-state)
new_levels = compute_gene_effect(
    edge_list=my_edges,
    initial_levels=my_initial_levels,
    source_gene='1',
    modified_level=1.5,
    max_depth=100
)

print(new_levels)
# {'1': 1.5, '2': 1.0, '3': 1.39..., '4': 1.66..., ...}
```

---

## Limitations

- **Computational complexity** grows exponentially with `max_depth` in densely connected graphs, since all walks (not just simple paths) are enumerated. For large networks, consider capping depth or switching to a matrix-based approach.
- **Edge weights are static** — the model does not account for nonlinear saturation, stochastic noise, or time-dependent regulation.
- **No inhibition semantics** — negative weights represent repression numerically, but the model does not enforce biological constraints (e.g., expression cannot go below zero).
