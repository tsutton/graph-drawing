pub mod graph;
pub mod layout;

use std::{cmp::Ordering, f64::consts::PI};

use crate::graph::Graph;

use layout::Vector;
use rand::prelude::SliceRandom;

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

pub struct KamadaKawaiDrawer {
    maximum_allowed_energy_derivative: f64,
    width: f64,
    height: f64,
    init_strategy: InitializationStrategy,
}

impl KamadaKawaiDrawer {
    // K scales energy, and all derivs, etc. So double it should double max energy, etc
    const K: f64 = 1.0;

    pub fn my_best_guess() -> Self {
        Self {
            width: 400.0,
            height: 400.0,
            maximum_allowed_energy_derivative: 1.0,
            init_strategy: InitializationStrategy::RegularPolygonRandomizedOrder,
        }
    }

    pub fn for_benchmarking() -> Self {
        Self {
            width: 400.0,
            height: 400.0,
            maximum_allowed_energy_derivative: 50.0,
            init_strategy: InitializationStrategy::RegularPolygon,
        }
    }

    pub fn draw(&self, graph: &Graph) -> Vec<Vector> {
        if graph.nodes == 0 {
            return Vec::new();
        }
        let mut positions: Vec<Vector> = Vec::with_capacity(graph.nodes);

        self.init_strategy.initialize(
            graph.nodes,
            &mut positions,
            0.0,
            self.width,
            0.0,
            self.height,
        );

        let distances = graph.all_pairs_shortest_paths();
        let max_distance = distances
            .iter()
            .flatten()
            .map(|x| x.expect("graph to be connected"))
            .max()
            .expect("graph to be nonempty");
        let draw_distance = self.width.min(self.height) / max_distance as f64;

        let mut deltas: Vec<f64> = (0..graph.nodes)
            .map(|i| {
                let (dx, dy) =
                    Self::jacobian(i, graph.nodes, &positions, &distances, draw_distance);
                (dx.powi(2) + dy.powi(2)).sqrt()
            })
            .collect();
        loop {
            // Possible optimizations:
            // - share subexpressions between jacobians and hessian
            // - don't recalc jacobian when solving hessian system
            let (m, max_delta) = deltas
                .iter()
                .enumerate()
                .max_by(|(_, f), (_, g)| match f.partial_cmp(g) {
                    Some(i) => i,
                    // if we have NaN or inf, ???
                    None => Ordering::Less,
                })
                .expect("graph to be nonempty");
            if *max_delta < self.maximum_allowed_energy_derivative {
                break;
            }
            let hessian = Self::hessian_m(m, graph.nodes, &positions, &distances, draw_distance);
            let jacobian = Self::jacobian(m, graph.nodes, &positions, &distances, draw_distance);
            let hessian_det = det(hessian);
            let delta_x =
                det([[-jacobian.0, hessian[0][1]], [-jacobian.1, hessian[1][1]]]) / hessian_det;
            let delta_y =
                det([[hessian[0][0], -jacobian.0], [hessian[1][0], -jacobian.1]]) / hessian_det;
            debug_assert!(
                (delta_x * hessian[0][0] + delta_y * hessian[0][1] + jacobian.0) / jacobian.0
                    < 0.01,
            );
            debug_assert!(
                (delta_x * hessian[1][0] + delta_y * hessian[1][1] + jacobian.1) / jacobian.1
                    < 0.01,
            );
            // TODO this is supposed to be +, not -, I must have flipped a sign some where, but I can't find where...
            positions[m].x -= delta_x;
            positions[m].y -= delta_y;
            let (new_dx, new_dy) =
                Self::jacobian(m, graph.nodes, &positions, &distances, draw_distance);
            deltas[m] = (new_dx.powi(2) + new_dy.powi(2)).sqrt();
        }
        positions
    }

    fn jacobian(
        m: usize,
        nodes: usize,
        positions: &[Vector],
        distances: &[Vec<Option<usize>>],
        draw_distance: f64,
    ) -> (f64, f64) {
        (0..nodes)
            .map(|i| {
                if i == m {
                    (0.0, 0.0)
                } else {
                    let distance = distances[m][i].unwrap() as f64;
                    let k_mi = Self::K / distance.powi(2);
                    let l_mi = draw_distance * distance;
                    let denominator = positions[m].distance_to(&positions[i]);
                    (
                        k_mi * (positions[m].x - positions[i].x) * (1.0 - l_mi / denominator),
                        k_mi * (positions[m].y - positions[i].y) * (1.0 - l_mi / denominator),
                    )
                }
            })
            .reduce(|(x1, y1), (x2, y2)| (x1 + x2, y1 + y2))
            .unwrap()
    }

    // Matrix:
    // [ [dxdx dxdy]
    //   [dydx dydy] ]
    fn hessian_m(
        m: usize,
        nodes: usize,
        positions: &[Vector],
        distances: &[Vec<Option<usize>>],
        draw_distance: f64,
    ) -> [[f64; 2]; 2] {
        let (dx_dx, dx_dy, dy_dx, dy_dy) = (0..nodes)
            .map(|i| {
                if i == m {
                    (0.0, 0.0, 0.0, 0.0)
                } else {
                    let distance = distances[m][i].unwrap() as f64;
                    let k_mi = Self::K / distance.powi(2);
                    let l_mi = draw_distance * distance;
                    let denominator = positions[m].distance_to(&positions[i]);
                    (
                        k_mi * (1.0
                            - l_mi * (positions[m].y - positions[i].y).powi(2) / denominator),
                        k_mi * (l_mi
                            * (positions[m].y - positions[i].y)
                            * (positions[m].x - positions[i].x)
                            / denominator),
                        k_mi * (l_mi
                            * (positions[m].y - positions[i].y)
                            * (positions[m].x - positions[i].x)
                            / denominator),
                        k_mi * (1.0
                            - l_mi * (positions[m].x - positions[i].x).powi(2) / denominator),
                    )
                }
            })
            .reduce(|(a1, b1, c1, d1), (a2, b2, c2, d2)| (a1 + a2, b1 + b2, c1 + c2, d1 + d2))
            .unwrap();
        [[dx_dx, dx_dy], [dy_dx, dy_dy]]
    }
}

fn det(matrix: [[f64; 2]; 2]) -> f64 {
    matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]
}

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
