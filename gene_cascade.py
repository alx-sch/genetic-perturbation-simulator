import igraph as ig
from typing import Dict, List, Tuple

# -----------------
# --- Functions ---
# -----------------

def get_total_influence(
	G: ig.Graph,					# The network graph
	current_gene: int,				# Current gene vertex index
	target_gene: int,				# Target gene vertex index
	current_path_product: float,	# Product of weights along the current path
	max_depth: int,					# Maximum depth to traverse
	current_depth: int				# Current depth in recursion
) -> float:
	"""
	Recursively calculates the cumulative influence factor (sum of weighted "walks")
	from the current gene to the target gene.

	It relies purely on max_depth to stop the recursion.
	"""

	if current_depth > max_depth:
		return 0.0
	
	total_influence = 0.0	# Initialize total influence for this path
	
	# Iterate over all genes directly downstream of the current gene
	for neighbor_gene in G.successors(current_gene):
		# Get the weight of this edge
		edge_id = G.get_eid(current_gene, neighbor_gene)
		weight = G.es[edge_id]['weight']

		# Calculate the new path product (W1 * W2 * W3...)
		new_path_product = current_path_product * weight

		if neighbor_gene == target_gene:
			# Base Case: Found the target. Add this path's total influence.
			total_influence += new_path_product

		total_influence += get_total_influence(
			G,
			neighbor_gene,		# The "new" current gene
			target_gene,
			new_path_product,	# The updated path product
			max_depth,
			current_depth + 1	# Increment depth
		)

	return total_influence


def compute_gene_effect(
	edge_list: List[Tuple[str, str, float]],
	initial_levels: Dict[str, float],
	source_gene: str,
	modified_level: float,
	max_depth: int = 5
) -> Dict[str, float]:
	"""
	Computes the effect of a gene perturbation using a recursive path summation model.

	Parameters:
		edge_list: A list of tuples defining the network, e.g., [('1', '4', 1.0), ...]
		initial_levels: A dict of baseline expression levels, e.g., {'1': 1.0, ...}
		source_gene: The name of the gene being modified (e.g., '1').
		modified_level: The new expression level of the source gene (e.g., 1.5).
		max_depth: The max number of steps (edges) to follow a path.

	Returns:
		A dictionary of all genes and their new calculated expression levels.
	"""

	# --- Setup Graph and Calculate Initial Shock ---

	# Get all unique gene names
	all_gene_names = set(initial_levels.keys())
	for source, target, weight in edge_list:
		all_gene_names.add(source)
		all_gene_names.add(target)

	# Create the igraph object using DictList
	# see here: https://python.igraph.org/en/main/api/igraph.Graph.html#DictList
	vertex_dicts = [{"name": name} for name in all_gene_names]
	edge_dicts = [
		{"source": source, "target": target, "weight": weight} 
		for source, target, weight in edge_list
	]

	G = ig.Graph.DictList(
		vertices=vertex_dicts,
		edges=edge_dicts, 
		directed=True
	)

	# Calculate the initial "shock" (Delta E_source)
	initial_shock = modified_level - initial_levels.get(source_gene, 0.0)

	# Initialize the results dictionary
	new_levels = initial_levels.copy()
	new_levels[source_gene] = modified_level

	source_vid = G.vs.find(name=source_gene).index # use vertex index for igraph

	# --- Iterate and Calculate for all Target Genes ---

	for target_gene in all_gene_names:
		if target_gene == source_gene:
			continue

		try:
			target_vid = G.vs.find(name=target_gene).index
		except ig.InternalError:
			# This gene is not in the graph (like Gene '2'), so it can't be reached.
			continue 

		# Call the recursive helper function
		total_influence_factor = get_total_influence(
			G,
			source_vid,
			target_vid,
			current_path_product=1.0,
			max_depth=max_depth,
			current_depth=1,  # Start depth count at 1 (first edge)
	)

		# --- Apply the Final Change ---
		initial_level_target = initial_levels.get(target_gene, 0.0)

		# Change = (Initial Shock) * (Sum of all Path Influences)
		change = initial_shock * total_influence_factor

		new_levels[target_gene] = initial_level_target + change

	return new_levels

# ---------------------------------
# --- Example Usage and Testing ---
# ---------------------------------

if __name__ == "__main__":

	# Define the gene network edges (source, target, weight)
	EDGES = [
		('1', '4', 1.0),
		('1', '5', 0.6),
		('1', '6', 0.7),
		('2', '4', 1.0),
		('3', '4', 0.4),
		('4', '3', 0.6),
		('4', '8', 0.3), 
		('5', '6', 1.0),
		('6', '7', 0.8),
	]

	# Define the initial levels for genes
	INITIAL_LEVELS = {
		'1': 1.0,
		'2': 1.0,
		'3': 1.0,
		'4': 1.0,
		'5': 1.0,
		'6': 1.0,
		'7': 1.0,
		'8': 1.0
	}

	# --- Run Simulation (Short-Term) ---
	DEPTH = 2 # Using a smaller depth to limit loop effects
	print(f"--- Simulation (Short-Term Response, Depth={DEPTH}) ---")

	# Call the new recursive function
	new_levels_short = compute_gene_effect(
		EDGES,
		INITIAL_LEVELS,
		source_gene='1',
		modified_level=1.5,
		max_depth=DEPTH
	)

	print(f"{'Gene':<7} {'Initial':>8} {'New Level':>10} {'Change (ΔE)':>12}")
	print("-" * 40)
	for gene in sorted(INITIAL_LEVELS.keys(), key=int):
		initial = INITIAL_LEVELS[gene]
		new = new_levels_short.get(gene, initial)
		change = new - initial
		print(f"{gene:<7} {initial:>8.2f} {new:>10.2f} {change:>+12.2f}")

	# --- Run Simulation (Long-Term) ---
	DEPTH = 100  # Using a high depth to capture the loop's steady-state effect
	print(f"\n--- Simulation (Long-Term/Steady-State, Depth={DEPTH}) ---")
	print("(Note the larger change in Genes 3, 4, and 8 due to the loop)")

	# Call the new recursive function
	new_levels_long = compute_gene_effect(
		EDGES,
		INITIAL_LEVELS,
		source_gene='1',
		modified_level=1.5,
		max_depth=DEPTH
	)

	print(f"{'Gene':<7} {'Initial':>8} {'New Level':>10} {'Change (ΔE)':>12}")
	print("-" * 40)
	for gene in sorted(INITIAL_LEVELS.keys(), key=int):
		initial = INITIAL_LEVELS[gene]
		new = new_levels_long.get(gene, initial)
		change = new - initial
		print(f"{gene:<7} {initial:>8.2f} {new:>10.2f} {change:>+12.2f}")
