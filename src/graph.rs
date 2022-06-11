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
        self.dijkstra()
    }

    pub fn dijkstra(&self) -> Vec<Vec<Option<usize>>> {
        (0..self.nodes)
            .map(|i| self.dijkstra_one_source(i))
            .collect()
    }

    // Returns the length of the shortest path from starting node to each node,
    // or None if the graph is not connected.
    // Very lightly adapted from the implementation example in std::collections::binary_heap.
    pub fn dijkstra_one_source(&self, starting_node: usize) -> Vec<Option<usize>> {
        // This is a medium-optimized version of Dijskstra's algorithm.
        // (T)his comment assumes basic familiarity with Dijsktra's alg.)
        // The best asymptotic complexity versions use a priority queue mapping keys to priorities,
        // that also includes a fast decrease_priority_of_key(key, new_priority) function,
        // or generally allow looking up the priority of a given key, rather than just the standard
        // insert_key_with_priority() and pop_min() functions.
        // In place of that, we'll use a strategy that might have a key appear twice in the queue, but
        // maintain the mapping of keys to priorities separatenly in a Vec (since keys are just usizes).
        // Then, during our loop, if we pop_min and find that that (key, priority) pair is actually
        // worse than some other path we've found since that pair was pushed to the queue, we skip the loop body.
        // let mut queue = BinaryHeap::new();
        let mut distances = vec![None; self.nodes];
        distances[starting_node] = Some(0);
        let mut pqueue = BinaryHeap::new();
        pqueue.push(DijkstraState {
            cost: 0,
            position: starting_node,
        });
        while let Some(DijkstraState { cost, position }) = pqueue.pop() {
            if distances[position]
                .map(|known_cost| known_cost < cost)
                .unwrap_or(false)
            {
                continue;
            }
            for (target, weight) in self.edges_from(position) {
                let next = DijkstraState {
                    cost: cost + weight,
                    position: target,
                };

                // If so, add it to the frontier and continue
                if distances[next.position]
                    .map(|known_cost| known_cost > next.cost)
                    .unwrap_or(true)
                {
                    pqueue.push(next);
                    // Relaxation, we have now found a better way
                    distances[next.position] = Some(next.cost);
                }
            }
        }

        distances
    }

    pub fn floyd_warshall(&self) -> Vec<Vec<Option<usize>>> {
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
/// ```text
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
            if x > 0 {
                edges_from_node.push((node_num - 1, 1));
            }
            if y + 1 < nodes_across {
                edges_from_node.push((node_num + nodes_across, 1));
            }
            if y > 0 {
                edges_from_node.push((node_num - nodes_across, 1));
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
        let edges_from_node = vec![((i + 1) % nodes, 1), ((nodes + i - 1) % nodes, 1)];
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
            edges[node].push((right_node, 1));
            edges[node].push((down_node, 1));
            edges[right_node].push((node, 1));
            edges[down_node].push((node, 1));
        }
    }
    Graph::from_edges_unchecked(nodes, edges)
}

// other graphs to test
// cubes
// serpinski triangle
// petersen graph

// DijskstraState and associated impls adapted from std::collections::binary_heap example.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
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

#[cfg(test)]
mod test {
    use super::grid_graph;

    #[test]
    fn test_grid() {
        let graph = grid_graph(3);
        assert_eq!(graph.edges().count(), 24);
        assert_eq!(graph.edge_weight(0, 1).unwrap(), 1);
        assert_eq!(graph.edge_weight(0, 3).unwrap(), 1);
        assert_eq!(graph.edge_weight(1, 0).unwrap(), 1);
        assert_eq!(graph.edge_weight(3, 0).unwrap(), 1);
        assert_eq!(
            graph.edges_from(0).collect::<Vec<_>>(),
            vec![(1, 1), (3, 1),]
        );
        assert!(graph.edge_weight(0, 2).is_none());
    }

    #[test]
    fn test_grid_dijkstra() {
        let graph = grid_graph(3);
        let paths_from_one = graph.dijkstra_one_source(0);
        assert_eq!(
            paths_from_one,
            vec![
                Some(0),
                Some(1),
                Some(2),
                Some(1),
                Some(2),
                Some(3),
                Some(2),
                Some(3),
                Some(4),
            ]
        )
    }
}
