use std::collections::HashMap;

// Let's say that our graphs are don't frequently change once constructed
// They are weighted, but not directed
// For Eades84 and  Fruchterman/Reingold91 we'll need to, in each iteration,
// iterate through all pairs of vertices; for each pair, check the weight of edge
// so perhaps an adjancency matrix is a good way to store it, considered as spare i.e. a hashmap
// of (x, y) => weight pairs, where x < y
pub struct Graph {
    // one more than the largest value of a node in the graph
    pub nodes: usize,
    // map (lower node, higher node ) => weight of that edge
    pub weights: HashMap<(usize, usize), usize>,
}

impl Graph {
    pub fn new() -> Graph {
        Graph {
            nodes: 0,
            weights: HashMap::new(),
        }
    }

    pub fn edge_weight(&self, node1: usize, node2: usize) -> Option<usize> {
        self.weights
            .get(&(node1.min(node2), node1.max(node2)))
            .copied()
    }

    pub fn edges(&'_ self) -> impl Iterator<Item = ((usize, usize), usize)> + '_ {
        self.weights.iter().map(|((a, b), c)| ((*a, *b), *c))
    }

    pub fn all_pairs_shortest_paths(&self) -> Vec<Vec<Option<usize>>> {
        // We follow Floyd-Warshall (simplified for undirected graphs):
        // Incrementally solve and store the subproblems APSP(from_vertex, to_vertex)
        // which may be Some(value), or None if the vertices aren't connected yet.
        // Note that for us (as opposed to Wikipedia), we number vertices from 0.
        let mut subproblem_cache: Vec<Vec<Option<usize>>> =
            vec![vec![None; self.nodes]; self.nodes];

        // as a base case, for vertices that have an actual edge beteween them, there are no intermediate
        // vertices needed
        for ((i, j), w) in self.edges() {
            subproblem_cache[i][j] = Some(w);
            // This is the the part where undirected is used
            subproblem_cache[j][i] = Some(w);
        }
        for i in 0..self.nodes {
            subproblem_cache[i][i] = Some(0);
        }

        // We'll frame this iteratively (as opposed to recursively).
        for intermediate in 0..self.nodes {
            let mut next_subproblem_cache: Vec<Vec<Option<usize>>> =
                vec![vec![None; self.nodes]; self.nodes];
            for from_vertex in 0..self.nodes {
                for to_vertex in 0..self.nodes {
                    let previous = subproblem_cache[from_vertex][to_vertex];
                    let using_k: Option<usize> = subproblem_cache[from_vertex][intermediate]
                        .zip(subproblem_cache[intermediate][to_vertex])
                        .map(|(x, y)| x + y);
                    match (previous, using_k) {
                        (None, None) => {}
                        (Some(x), None) | (None, Some(x)) => {
                            next_subproblem_cache[from_vertex][to_vertex] = Some(x);
                        }
                        (Some(x), Some(y)) => {
                            next_subproblem_cache[from_vertex][to_vertex] = Some(x.min(y));
                        }
                    }
                }
            }
            subproblem_cache = next_subproblem_cache;
        }
        subproblem_cache
    }
}

impl Default for Graph {
    fn default() -> Self {
        Self::new()
    }
}

/// A square grid graph, with each node connected to the left, right, above, and below nodes if those
/// nodes aren't past the boundary.
///
/// For example, for nodes_across=3, the graph looks like:
/// ```
///  . - . - .
///  |   |   |
///  . - . - .
///  |   |   |
///  . - . - .
/// ```
/// The nodes are numbered consistent with left-to-right, top-to-bottom
pub fn grid_graph(nodes_across: usize) -> Graph {
    let mut g = Graph::new();
    g.nodes = nodes_across * nodes_across;
    // we number nodes left to right top to bottom
    for x in 0..nodes_across {
        for y in 0..nodes_across {
            let node_num = y * nodes_across + x;
            if x + 1 < nodes_across {
                g.weights.insert((node_num, node_num + 1), 1);
            }
            if y + 1 < nodes_across {
                g.weights.insert((node_num, node_num + nodes_across), 1);
            }
        }
    }
    g
}

/// A graph with each node connect to its two neighbors, forming a cycle.
pub fn cycle_graph(nodes: usize) -> Graph {
    let mut g = Graph::new();
    g.nodes = nodes;
    for i in 0..nodes {
        g.weights.insert((i, (i + 1) % nodes), 1);
    }
    g
}

/// A graph representation of a torus: a grid with the given dimensions, but with the
/// nodes at the edge connected "around" to the opposite edge
pub fn torus_graph(width: usize, height: usize) -> Graph {
    let mut g = Graph::new();
    g.nodes = width * height;
    let node_num = |x, y| y * width + x;

    for x in 0..width {
        for y in 0..height {
            let node = node_num(x, y);
            let right_node = node_num((x + 1) % width, y);
            let down_node = node_num(x, (y + 1) % height);
            g.weights.insert((node, right_node), 1);
            g.weights.insert((node, down_node), 1);
        }
    }
    g
}

// other graphs to test
// cubes
// serpinski triangle
// petersen graph
