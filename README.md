# Genetic Perturbation Simulator

## Overview
This repository contains a Python tool for simulating the effects of a single gene perturbation (e.g., overexpression or knockout) across a weighted, directed biological network.

The simulation is based on a **recursive, depth-limited path summation model**. It uses the `igraph` library to calculate how an initial "shock" (e.g., a `+0.5` change in expression) propagates through the network.

The model correctly accounts for:
- **Additive Influence:** The cumulative effect from multiple paths (e.g., `1 -> 5 -> 6` and `1 -> 6`).
- **Feedback Loops:** The amplifying (or dampening) effects of loops (e.g., `4 <-> 3`).

By adjusting the `max_depth` parameter, the model can simulate both the **short-term response** (immediate effects from simple paths) and the **long-term steady-state** (the full effect after influence has cycled through loops).

---

## How to Run

The script is self-contained and includes an example simulation.

#### 1. Dependencies

The model requires the `igraph` library. You can install it via pip:
```bash
pip install igraph
```
(If you are using the included Dev Container (e.g. via Codespaces), all dependencies are already installed).

#### 2. Run the Example

Simply execute the Python script:
```bash
python gene_cascade.py
```

#### Example Output

Running the script will produce the following output, demonstrating the difference between a short-term and long-term simulation.
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

Notice the difference between the two simulations.

- In the **Short-Term (Depth=2)** run, the influence from `1 -> 4` (1.0) is the only major path affecting Gene 4.
- In the **Long-Term (Depth=100)** run, the simulation allows the influence to cycle through the `4 -> 3 -> 4` feedback loop many times. The loop's gain is `0.6 * 0.4 = 0.24`, which acts as an amplifier.
- This amplification causes the final change in **Genes 4, 3 and 8** to be significantly higher as they reach their "steady-state" expression levels.

This amplification causes the final change in Genes 4, 3, and 8 to be significantly higher as they reach their "steady-state" expression levels.

---

## Usage as a Module

The script is structured so its functions can be imported into your own projects or a Jupyter Notebook.

1. `compute_gene_effect(...)`: The main function that sets up the igraph object and calculates the final expression levels.
2. `get_total_influence(...)`: A recursive helper function that traverses the graph.

To use this in your own project, import the `compute_gene_effect` function and pass it your own `edge_list` and `initial_levels`.

#### Example:

```Python
# Import the function from your script
from gene_cascade import compute_gene_effect

# 1. Define your network as a list of (source, target, weight)
my_edges = [
    ('1', '4', 1.0),
    ('1', '5', 0.6),
    ('4', '3', 0.6), 
    ('5', '6', 1.0),
    # ... etc
]

# 2. Define the baseline levels for all genes
my_initial_levels = {
    '1': 1.0,
    '2': 1.0,
    '3': 1.0,
    '4': 1.0,
    # ... etc
}

# 3. Run the simulation
# Let's check the long-term effect of modifying Gene '1' to 1.5
new_levels = compute_gene_effect(
    edge_list=my_edges,
    initial_levels=my_initial_levels,
    source_gene='1',
    modified_level=1.5,
    max_depth=100
)

# The result is a dictionary:
# {'1': 1.5, '2': 1.0, '3': 1.39..., '4': 1.66..., ...}
print(new_levels)
```
