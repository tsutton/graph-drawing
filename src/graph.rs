use std::collections::BinaryHeap;

// Let's say that our graphs are don't frequently change once constructed
// They are weighted, directed
pub struct Graph {
    // one more than the largest value of a node in the graph
    pub nodes: usize,
    // edges[i] = list of pairs (j, weight) where there's an edge of that weight from i to j
    edge_lists: Vec<Vec<(usize, usize)>>,
}

impl Graph {
    pub fn new() -> Graph {
        Graph {
            nodes: 0,
            edge_lists: Vec::new(),
        }
    }

    pub fn from_edges_unchecked(nodes: usize, edges: Vec<Vec<(usize, usize)>>) -> Self {
        Self {
            nodes,
            edge_lists: edges,
        }
    }

    pub fn edge_weight(&self, node1: usize, node2: usize) -> Option<usize> {
        self.edge_lists.get(node1).and_then(|edges| {
            edges
                .iter()
                .find(|(j, _)| *j == node2)
                .map(|(_, weight)| *weight)
        })
    }

    pub fn edges(&'_ self) -> impl Iterator<Item = ((usize, usize), usize)> + '_ {
        self.edge_lists
            .iter()
            .enumerate()
            .flat_map(|(source, edges)| {
                edges
                    .iter()
                    .map(move |(target, weight)| ((source, *target), *weight))
            })
    }

    pub fn edges_from(&'_ self, from: usize) -> impl Iterator<Item = (usize, usize)> + '_ {
        self.edge_lists[from].iter().copied()
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
            for from_vertex in 0..self.nodes {
                for to_vertex in 0..self.nodes {
                    let previous = subproblem_cache[from_vertex][to_vertex];
                    let using_k: Option<usize> = subproblem_cache[from_vertex][intermediate]
                        .zip(subproblem_cache[intermediate][to_vertex])
                        .map(|(x, y)| x + y);
                    match (previous, using_k) {
                        (None, None) => {}
                        (Some(x), None) | (None, Some(x)) => {
                            subproblem_cache[from_vertex][to_vertex] = Some(x);
                        }
                        (Some(x), Some(y)) => {
                            subproblem_cache[from_vertex][to_vertex] = Some(x.min(y));
                        }
                    }
                }
            }
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
    let nodes = nodes_across * nodes_across;
    let mut edges = Vec::with_capacity(nodes);
    // we number nodes left to right top to bottom
    for y in 0..nodes_across {
        for x in 0..nodes_across {
            let node_num = y * nodes_across + x;
            let mut edges_from_node = Vec::new();
            if x + 1 < nodes_across {
                edges_from_node.push((node_num + 1, 1));
            }
            if y + 1 < nodes_across {
                edges_from_node.push((node_num + nodes_across, 1));
            }
            edges.push(edges_from_node);
        }
    }
    Graph::from_edges_unchecked(nodes, edges)
}

/// A graph with each node connect to its two neighbors, forming a cycle.
pub fn cycle_graph(nodes: usize) -> Graph {
    let mut edges = Vec::with_capacity(nodes);
    for i in 0..nodes {
        let edges_from_node = if i == 0 {
            vec![(1, 1), (nodes - 1, 1)]
        } else {
            vec![(nodes + 1, 1)]
        };
        edges.push(edges_from_node)
    }
    Graph::from_edges_unchecked(nodes, edges)
}

/// A graph representation of a torus: a grid with the given dimensions, but with the
/// nodes at the edge connected "around" to the opposite edge
pub fn torus_graph(width: usize, height: usize) -> Graph {
    let nodes = width * height;
    let node_num = |x, y| y * width + x;
    let mut edges = vec![Vec::new(); nodes];

    for y in 0..height {
        for x in 0..width {
            let node = node_num(x, y);
            let right_node = node_num((x + 1) % width, y);
            let down_node = node_num(x, (y + 1) % height);
            // TODO use min max
            edges[node.min(right_node)].push((node.max(right_node), 1));
            edges[node.min(down_node)].push((node.max(down_node), 1));
        }
    }
    Graph::from_edges_unchecked(nodes, edges)
}

// other graphs to test
// cubes
// serpinski triangle
// petersen graph

// DijskstraState and associated impls adapted from std::collections::binary_heap example.
#[derive(Copy, Clone, Eq, PartialEq)]
struct DijkstraState {
    cost: usize,
    position: usize,
}

impl std::cmp::Ord for DijkstraState {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        // Notice that the we flip the ordering on costs.
        // In case of a tie we compare positions - this step is necessary
        // to make implementations of `PartialEq` and `Ord` consistent.
        other
            .cost
            .cmp(&self.cost)
            .then_with(|| self.position.cmp(&other.position))
    }
}

impl PartialOrd for DijkstraState {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
