use std::collections::HashMap;
use svg::node::element::Circle;
use svg::{node, Document};

fn main() {
    // let graph = Graph {
    //     nodes: 3,
    //     weights: HashMap::from([((0, 1), 1), ((0, 2), 1), ((1, 2), 3)]),
    // };
    let graph = grid_graph(3);
    let positions = eades(&graph);
    let svg_document = to_svg(&graph, &positions);
    svg::save("out.svg", &svg_document).unwrap();
}

// Let's say that our graphs are don't frequently change once constructed
// They are weighted, but not directed
// For Eades84 and  Fruchterman/Reingold91 we'll need to, in each iteration,
// iterate through all pairs of vertices; for each pair, check the weight of edge
// so perhaps an adjancency matrix is a good way to store it, considered as spare i.e. a hashmap
// of (x, y) => weight pairs, where x < y
struct Graph {
    // one more than the largest value of a node in the graph
    nodes: usize,
    // map (lower node, higher node ) => weight of that edge
    weights: HashMap<(usize, usize), usize>,
}

impl Graph {
    fn new() -> Graph {
        Graph {
            nodes: 0,
            weights: HashMap::new(),
        }
    }

    fn edge_weight(&self, node1: usize, node2: usize) -> Option<usize> {
        self.weights
            .get(&(node1.min(node2), node1.max(node2)))
            .copied()
    }
}

fn eades(graph: &Graph) -> Vec<Vector> {
    let hooke_constant: f64 = 2.0; // c1, 2.0 rec
    let ideal_length: f64 = 5.0; // c2, 1.0 rec
    let repel_constant: f64 = 1.0; // c3, 1.0 rec
    let force_multiplier: f64 = 0.15; // c4, 0.1 rec
    let iterations = 500;
    let viewport_size = 100.0;

    let mut positions: Vec<Vector> = Vec::with_capacity(graph.nodes);

    // initialize
    for _ in 0..graph.nodes {
        let x_rand: f64 = rand::random();
        let y_rand: f64 = rand::random();
        positions.push((x_rand * viewport_size, y_rand * viewport_size).into())
    }

    for _ in 0..iterations {
        let mut forces: Vec<Vector> = vec![(0.0, 0.0).into(); graph.nodes];
        for first_node in 0..graph.nodes {
            for second_node in (first_node + 1)..graph.nodes {
                let first_pt = positions[first_node];
                let second_pt = positions[second_node];
                let current_distance = first_pt.distance_to(&second_pt);
                let unit_vec = first_pt.unit_vector_to(&second_pt);
                match graph.edge_weight(first_node, second_node) {
                    Some(weight) => {
                        let factor = hooke_constant
                            * (current_distance / (weight as f64) / ideal_length).ln();
                        forces[first_node] = forces[first_node] + unit_vec.scaled(factor);
                        forces[second_node] = forces[second_node] + unit_vec.scaled(-factor);
                    }
                    None => {
                        let factor = repel_constant / (current_distance * current_distance);
                        forces[first_node] = forces[first_node] + unit_vec.scaled(-factor);
                        forces[second_node] = forces[second_node] + unit_vec.scaled(factor);
                    }
                }
            }
        }

        for i in 0..graph.nodes {
            positions[i] = positions[i] + forces[i].scaled(force_multiplier);
        }
    }
    positions
}

#[derive(Clone, Copy, PartialEq, Debug)]
struct Vector {
    x: f64,
    y: f64,
}

impl Vector {
    fn distance_to(&self, v: &Vector) -> f64 {
        ((self.x - v.x) * (self.x - v.x) + (self.y - v.y) * (self.y - v.y)).sqrt()
    }

    fn unit_vector_to(&self, v: &Vector) -> Vector {
        let distance = self.distance_to(v);
        (*v - *self).scaled(1.0 / distance)
    }

    fn scaled(&self, factor: f64) -> Vector {
        Vector {
            x: self.x * factor,
            y: self.y * factor,
        }
    }
}

impl std::ops::Add for Vector {
    type Output = Vector;

    fn add(self, rhs: Self) -> Self::Output {
        Vector {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
        }
    }
}

impl std::ops::Sub for Vector {
    type Output = Vector;

    fn sub(self, rhs: Self) -> Self::Output {
        Vector {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
        }
    }
}

impl From<(f64, f64)> for Vector {
    fn from(pt: (f64, f64)) -> Self {
        Vector { x: pt.0, y: pt.1 }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn vector_test() {
        struct TestCase {
            first: (f64, f64),
            second: (f64, f64),
            expected_distance: f64,
            expected_unit: (f64, f64),
        }
        fn case(first: (f64, f64), second: (f64, f64), dist: f64, unit: (f64, f64)) -> TestCase {
            TestCase {
                first,
                second,
                expected_distance: dist,
                expected_unit: unit,
            }
        }
        let cases = vec![
            case((0.0, 0.0), (1.0, 0.0), 1.0, (1.0, 0.0)),
            case((0.0, 0.0), (-1.0, 0.0), 1.0, (-1.0, 0.0)),
            case((0.0, 0.0), (2.0, 0.0), 2.0, (1.0, 0.0)),
            // case((0.0, 0.0), (3.0, 4.0), 5.0, (0.6, 0.8)),
        ];
        for (i, case) in cases.iter().enumerate() {
            let left: Vector = case.first.into();
            let right: Vector = case.second.into();
            assert_eq!(left.distance_to(&left), 0.0, "case {i}");
            assert_eq!(right.distance_to(&right), 0.0, "case {i}");
            assert_eq!(left.distance_to(&right), case.expected_distance, "case {i}");
            assert_eq!(right.distance_to(&left), case.expected_distance, "case {i}");
            assert_eq!(
                left.unit_vector_to(&right),
                case.expected_unit.into(),
                "case {i}"
            );
        }
    }
}

// TODO draw edges
fn to_svg(graph: &Graph, positions: &[Vector]) -> Document {
    debug_assert!(graph.nodes == positions.len());
    let mut document = Document::new();
    for (i, position) in positions.iter().enumerate() {
        let circle = Circle::new()
            .set("cx", position.x)
            .set("cy", position.y)
            .set("r", 2)
            .set("fill", "white")
            .set("stroke", "black")
            .set("stroke-width", 1)
            .add(node::element::Title::new().add(node::Text::new(format!("{}", i))));
        document = document.add(circle);
    }
    document
}

fn grid_graph(nodes_across: usize) -> Graph {
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
    println!("{:?}", g.weights);
    g
}
