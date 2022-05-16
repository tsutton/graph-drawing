pub mod graph;
pub mod layout;

use crate::graph::Graph;

use layout::Vector;

pub struct EadesDrawer {
    hooke_constant: f64,
    ideal_length: f64,
    repel_constant: f64,
    force_multiplier: f64,
    iterations: usize,
    viewport_size: f64,
}

impl EadesDrawer {
    pub fn paper_recommended() -> Self {
        Self {
            hooke_constant: 2.0,
            ideal_length: 1.0,
            repel_constant: 1.0,
            force_multiplier: 0.1,
            iterations: 100,
            // paper does not have any particular viewport size, which is confusing, since surely it matters
            viewport_size: 100.0,
        }
    }

    pub fn my_best_guess() -> Self {
        Self {
            hooke_constant: 2.0,
            ideal_length: 5.0,
            repel_constant: 1.0,
            force_multiplier: 0.15,
            iterations: 500,
            // paper does not have any particular viewport size, which is confusing, since surely it matters
            viewport_size: 100.0,
        }
    }

    pub fn draw(&self, graph: &Graph) -> Vec<Vector> {
        let mut positions: Vec<Vector> = Vec::with_capacity(graph.nodes);

        // initialize
        for _ in 0..graph.nodes {
            let x_rand: f64 = rand::random();
            let y_rand: f64 = rand::random();
            positions.push((x_rand * self.viewport_size, y_rand * self.viewport_size).into())
        }

        for _ in 0..self.iterations {
            let mut forces: Vec<Vector> = vec![(0.0, 0.0).into(); graph.nodes];
            for first_node in 0..graph.nodes {
                for second_node in (first_node + 1)..graph.nodes {
                    let first_pt = positions[first_node];
                    let second_pt = positions[second_node];
                    let current_distance = first_pt.distance_to(&second_pt);
                    let unit_vec = first_pt.unit_vector_to(&second_pt);
                    match graph.edge_weight(first_node, second_node) {
                        Some(weight) => {
                            let factor = self.hooke_constant
                                * (current_distance / (weight as f64) / self.ideal_length).ln();
                            forces[first_node] = forces[first_node] + unit_vec.scaled(factor);
                            forces[second_node] = forces[second_node] + unit_vec.scaled(-factor);
                        }
                        None => {
                            let factor =
                                self.repel_constant / (current_distance * current_distance);
                            forces[first_node] = forces[first_node] + unit_vec.scaled(-factor);
                            forces[second_node] = forces[second_node] + unit_vec.scaled(factor);
                        }
                    }
                }
            }

            for i in 0..graph.nodes {
                positions[i] = positions[i] + forces[i].scaled(self.force_multiplier);
            }
        }
        positions
    }
}

pub struct FruchtermanReingoldDrawer {
    iterations: usize,
    width: f64,
    height: f64,
    initial_temperature: f64,
}

impl FruchtermanReingoldDrawer {
    pub fn my_best_guess() -> Self {
        Self {
            iterations: 100,
            width: 400.0,
            height: 400.0,
            initial_temperature: 40.0,
        }
    }

    pub fn draw(&self, graph: &Graph) -> Vec<Vector> {
        let mut positions: Vec<Vector> = Vec::with_capacity(graph.nodes);
        let ideal_distance = (self.width * self.height / graph.nodes as f64).sqrt();

        // initialize
        for _ in 0..graph.nodes {
            let x_rand: f64 = rand::random();
            let y_rand: f64 = rand::random();
            positions.push((x_rand * self.width, y_rand * self.height).into())
        }

        let mut temperature = self.initial_temperature;
        for _ in 0..self.iterations {
            let mut displacements: Vec<Vector> = vec![(0.0, 0.0).into(); graph.nodes];
            for first_node in 0..graph.nodes {
                for second_node in (first_node + 1)..graph.nodes {
                    let first_pt = positions[first_node];
                    let second_pt = positions[second_node];
                    let current_distance = first_pt.distance_to(&second_pt);
                    let unit_vec = first_pt.unit_vector_to(&second_pt);

                    // all pairs repel each other
                    let repel_factor = ideal_distance.powi(2) / current_distance;
                    displacements[first_node] =
                        displacements[first_node] + unit_vec.scaled(-repel_factor);
                    displacements[second_node] =
                        displacements[second_node] + unit_vec.scaled(repel_factor);

                    if let Some(weight) = graph.edge_weight(first_node, second_node) {
                        let attractive_factor = current_distance.powi(2) / ideal_distance;
                        displacements[first_node] =
                            displacements[first_node] + unit_vec.scaled(attractive_factor);
                        displacements[second_node] =
                            displacements[second_node] + unit_vec.scaled(-attractive_factor);
                    }
                }
            }

            for i in 0..graph.nodes {
                let displacement = if displacements[i].len() > temperature {
                    displacements[i].scaled(temperature / displacements[i].len())
                } else {
                    displacements[i]
                };
                debug_assert!(temperature * 1.01 >= displacement.len());
                positions[i] = positions[i] + displacement;
                positions[i].x = positions[i].x.max(0.0).min(self.width);
                positions[i].y = positions[i].y.max(0.0).min(self.height);
            }
            temperature -= self.initial_temperature / self.iterations as f64;
        }
        positions
    }
}
