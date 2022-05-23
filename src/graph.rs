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
}

impl Default for Graph {
    fn default() -> Self {
        Self::new()
    }
}

pub fn grid_graph(nodes_across: usize) -> Graph {
    let mut g = Graph::new();
    g.nodes = nodes_across * nodes_across;
    // for n=2, we have the graph
    //  . - .
    //  |   |
    //  . - .
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

pub fn cycle_graph(nodes: usize) -> Graph {
    let mut g = Graph::new();
    g.nodes = nodes;
    for i in 0..nodes {
        g.weights.insert((i, (i + 1) % nodes), 1);
    }
    g
}

// other graphs to test
// cubes
// serpinski triangle
// petersen graph
