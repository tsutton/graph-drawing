pub mod graph;
pub mod layout;

use std::f64::consts::PI;

use crate::graph::Graph;

use layout::Vector;
use rand::prelude::SliceRandom;

/// Implements the algorithm from \[Ead84\], following description in \[Kou13\].
///
/// For each pair of nodes, if the nodes have an edge, we give them an attract force proportional
/// to log(distance between them/desired length).
/// Then all non-adjacent nodes repel proportional to 1/distance^2.
/// These forces are then simulated for a fixed number of iteration.
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

/// Implements the algorithm from \[FR91\], following description in \[K13\].
///
/// For each pair of nodes, if the nodes have an edge, we give them an attract force proportional
/// to (distance between them)^2. Then all nodes (even pairs with an edge) repel proportional to 1/distance.
/// These two forces balance each other out exactly when connected nodes are at a chosen ideal distance.
/// Then we iterate: calculate all the forces on all the nodes, apply those forces with some multiplier, and repeat.
/// The multiplier starts out high and then reduces over time so that the nodes settle into place.
pub struct FruchtermanReingoldDrawer {
    iterations: usize,
    width: f64,
    height: f64,
    initial_temperature: f64,
    init_strategy: InitializationStrategy,
}

impl FruchtermanReingoldDrawer {
    pub fn my_best_guess() -> Self {
        Self {
            iterations: 500,
            width: 400.0,
            height: 400.0,
            initial_temperature: 40.0,
            init_strategy: InitializationStrategy::RegularPolygon,
        }
    }

    pub fn draw(&self, graph: &Graph) -> Vec<Vector> {
        let mut positions: Vec<Vector> = Vec::with_capacity(graph.nodes);
        let ideal_distance = (self.width * self.height / graph.nodes as f64).sqrt();

        self.init_strategy.initialize(
            graph.nodes,
            &mut positions,
            0.0,
            self.width,
            0.0,
            self.height,
        );

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

                    if graph.edge_weight(first_node, second_node).is_some() {
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
                // TODO maybe instead of forcing inside the box, do a post-processing to translate/scale back to within the box?
                positions[i].x = positions[i].x.max(0.0).min(self.width);
                positions[i].y = positions[i].y.max(0.0).min(self.height);
            }
            temperature -= self.initial_temperature / self.iterations as f64;
        }
        positions
    }
}

mod kamada_kawai;

pub use kamada_kawai::{KamadaKawaiDrawer, KamadaKawaiFastDrawer};

pub enum InitializationStrategy {
    Random,
    RegularPolygon,
    RegularPolygonRandomizedOrder,
}

impl InitializationStrategy {
    fn initialize(
        &self,
        node_count: usize,
        positions: &mut Vec<Vector>,
        min_x: f64,
        max_x: f64,
        min_y: f64,
        max_y: f64,
    ) {
        match self {
            InitializationStrategy::Random => {
                for _ in 0..node_count {
                    let x_rand: f64 = rand::random();
                    let y_rand: f64 = rand::random();
                    positions.push(
                        (
                            x_rand * (max_x - min_x) + min_x,
                            y_rand * (max_y - min_y) + min_y,
                        )
                            .into(),
                    )
                }
            }
            InitializationStrategy::RegularPolygon => {
                let center_x = (max_x + min_x) / 2.0;
                let center_y = (max_y + min_y) / 2.0;
                let scale_factor_x = (max_x - min_x) / 2.0;
                let scale_factor_y = (max_y - min_y) / 2.0;
                for i in 0..node_count {
                    let angle = 2.0 * PI / node_count as f64 * i as f64;
                    let unit_circle_x = angle.cos();
                    let unit_circle_y = angle.sin();
                    positions.push(
                        (
                            unit_circle_x * scale_factor_x + center_x,
                            unit_circle_y * scale_factor_y + center_y,
                        )
                            .into(),
                    );
                }
            }
            InitializationStrategy::RegularPolygonRandomizedOrder => {
                let start_idx = positions.len();
                InitializationStrategy::RegularPolygon
                    .initialize(node_count, positions, min_x, max_x, min_y, max_y);
                positions[start_idx..].shuffle(&mut rand::thread_rng());
            }
        }
    }
}
