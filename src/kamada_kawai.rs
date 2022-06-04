use super::InitializationStrategy;
use crate::graph::Graph;
use crate::layout::Vector;

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

        let nodes = graph.nodes;
        let distances_raw = graph.all_pairs_shortest_paths();
        let mut distances_flat = vec![0.0; graph.nodes * graph.nodes];
        for i in 0..nodes {
            for j in 0..nodes {
                distances_flat[i * nodes + j] =
                    distances_raw[i][j].expect("graph to be connected") as f64;
            }
        }

        let max_distance = distances_raw
            .iter()
            .flatten()
            .max()
            .expect("graph to be nonempty")
            .expect("graph to be nonempty") as f64;
        let draw_distance = self.width.min(self.height) / max_distance;

        let mut deltas: Vec<f64> = (0..graph.nodes)
            .map(|i| {
                let (dx, dy) =
                    Self::jacobian(i, graph.nodes, &positions, &distances_flat, draw_distance);
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
                .max_by(|(_, f), (_, g)| {
                    f.partial_cmp(g)
                        .expect("no non-comparable floats in deltas")
                })
                .expect("graph to be nonempty");
            if *max_delta < self.maximum_allowed_energy_derivative {
                break;
            }
            let hessian =
                Self::hessian_m(m, graph.nodes, &positions, &distances_flat, draw_distance);
            let jacobian =
                Self::jacobian(m, graph.nodes, &positions, &distances_flat, draw_distance);
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
                Self::jacobian(m, graph.nodes, &positions, &distances_flat, draw_distance);
            deltas[m] = (new_dx.powi(2) + new_dy.powi(2)).sqrt();
        }
        positions
    }

    fn jacobian(
        m: usize,
        nodes: usize,
        positions: &[Vector],
        distances: &[f64],
        draw_distance: f64,
    ) -> (f64, f64) {
        (0..nodes)
            .map(|i| {
                if i == m {
                    (0.0, 0.0)
                } else {
                    let distance = distances[m * nodes + i];
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
        distances: &[f64],
        draw_distance: f64,
    ) -> [[f64; 2]; 2] {
        // For reasons I do not understand, local benchmarking indicates that keeping all four floats,
        // instead of eliminating the redudant dy_dx (which is always the same as dx_dy), seems to be
        // better for performance
        let (dx_dx, dx_dy, dy_dx, dy_dy) = (0..nodes)
            .map(|i| {
                if i == m {
                    (0.0, 0.0, 0.0, 0.0)
                } else {
                    let distance = distances[m * nodes + i];
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
