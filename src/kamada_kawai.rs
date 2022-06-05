use super::InitializationStrategy;
use crate::graph::Graph;
use crate::layout::Vector;

/// Implements draw based on the original Kamada Kawai 1989 paper \[KK89\].
///
/// Assigns an energy function to the layout based on having all edges act as springs.
/// The length of the spring is equal to the path-length along the graph between the nodes.
/// The goal is to minimize this energy - equivalently, find a zero of the derivatives of it.
/// Those zeroes are found by Newton's method: we loop, picking the node whose contribution to the
/// derivative is highest, then apply Netwon's method (similar to gradient descent) to move that one node.
/// This is repeated until the derivatives are sufficiently small.
pub struct KamadaKawaiDrawer {
    maximum_allowed_energy_derivative: f64,
    width: f64,
    height: f64,
    init_strategy: InitializationStrategy,
}

// K scales energy, and all derivs, etc. So double it should double max energy, etc
const K: f64 = 1.0;

impl KamadaKawaiDrawer {
    pub fn my_best_guess() -> Self {
        Self {
            width: 400.0,
            height: 400.0,
            maximum_allowed_energy_derivative: 1e-7,
            init_strategy: InitializationStrategy::RegularPolygonRandomizedOrder,
        }
    }

    pub fn for_benchmarking() -> Self {
        Self {
            width: 400.0,
            height: 400.0,
            maximum_allowed_energy_derivative: 1e-7,
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
                let (dx, dy) = jacobian(i, graph.nodes, &positions, &distances_flat, draw_distance);
                (dx.powi(2) + dy.powi(2)).sqrt()
            })
            .collect();
        loop {
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
            let hessian = hessian_m(m, graph.nodes, &positions, &distances_flat, draw_distance);
            let jacobian_value =
                jacobian(m, graph.nodes, &positions, &distances_flat, draw_distance);
            let hessian_det = det(hessian);
            let delta_x = det([
                [-jacobian_value.0, hessian[0][1]],
                [-jacobian_value.1, hessian[1][1]],
            ]) / hessian_det;
            let delta_y = det([
                [hessian[0][0], -jacobian_value.0],
                [hessian[1][0], -jacobian_value.1],
            ]) / hessian_det;
            debug_assert!(
                (delta_x * hessian[0][0] + delta_y * hessian[0][1] + jacobian_value.0)
                    / jacobian_value.0
                    < 0.01,
            );
            debug_assert!(
                (delta_x * hessian[1][0] + delta_y * hessian[1][1] + jacobian_value.1)
                    / jacobian_value.1
                    < 0.01,
            );
            positions[m].x += delta_x;
            positions[m].y += delta_y;
            let (new_dx, new_dy) =
                jacobian(m, graph.nodes, &positions, &distances_flat, draw_distance);
            deltas[m] = (new_dx.powi(2) + new_dy.powi(2)).sqrt();
        }
        positions
    }
}

fn det(matrix: [[f64; 2]; 2]) -> f64 {
    matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]
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
                let k_mi = K / distance.powi(2);
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

fn jacobian_local(
    m: usize,
    positions: &[Vector],
    distances: &[(usize, f64)],
    draw_distance: f64,
) -> (f64, f64) {
    distances
        .iter()
        .map(|&(i, distance)| {
            if i == m {
                (0.0, 0.0)
            } else {
                let k_mi = K / distance.powi(2);
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
                let k_mi = K / distance.powi(2);
                let l_mi = draw_distance * distance;
                let denominator = positions[m].distance_to(&positions[i]).powi(3);
                (
                    k_mi * (1.0 - l_mi * (positions[m].y - positions[i].y).powi(2) / denominator),
                    k_mi * (l_mi
                        * (positions[m].y - positions[i].y)
                        * (positions[m].x - positions[i].x)
                        / denominator),
                    k_mi * (l_mi
                        * (positions[m].y - positions[i].y)
                        * (positions[m].x - positions[i].x)
                        / denominator),
                    k_mi * (1.0 - l_mi * (positions[m].x - positions[i].x).powi(2) / denominator),
                )
            }
        })
        .reduce(|(a1, b1, c1, d1), (a2, b2, c2, d2)| (a1 + a2, b1 + b2, c1 + c2, d1 + d2))
        .unwrap();
    [[dx_dx, dx_dy], [dy_dx, dy_dy]]
}

fn hessian_m_local(
    m: usize,
    positions: &[Vector],
    distances: &[(usize, f64)],
    draw_distance: f64,
) -> [[f64; 2]; 2] {
    let (dx_dx, dx_dy, dy_dx, dy_dy) = distances
        .iter()
        .map(|&(i, distance)| {
            if i == m {
                (0.0, 0.0, 0.0, 0.0)
            } else {
                let k_mi = K / distance.powi(2);
                let l_mi = draw_distance * distance;
                let denominator = positions[m].distance_to(&positions[i]).powi(3);
                (
                    k_mi * (1.0 - l_mi * (positions[m].y - positions[i].y).powi(2) / denominator),
                    k_mi * (l_mi
                        * (positions[m].y - positions[i].y)
                        * (positions[m].x - positions[i].x)
                        / denominator),
                    k_mi * (l_mi
                        * (positions[m].y - positions[i].y)
                        * (positions[m].x - positions[i].x)
                        / denominator),
                    k_mi * (1.0 - l_mi * (positions[m].x - positions[i].x).powi(2) / denominator),
                )
            }
        })
        .reduce(|(a1, b1, c1, d1), (a2, b2, c2, d2)| (a1 + a2, b1 + b2, c1 + c2, d1 + d2))
        .unwrap();
    [[dx_dx, dx_dy], [dy_dx, dy_dy]]
}

/// Implements a more optimized variant of KK based on \[HK01\].
///
/// The overall plan (minimize energy function by Newton's method) stays the same, but
/// We use a more efficient implementation. We also choose to limit how many iterations to run for,
/// although it means the final output might not be so great for large graphs - that's why \[HK01\] only uses
/// this as a "local beautification" as part of a multi-scale strategy.
pub struct KamadaKawaiFastDrawer {
    width: f64,
    height: f64,
    iteration_multiple: usize,
    radius: usize,
}

impl KamadaKawaiFastDrawer {
    pub fn fast() -> Self {
        Self {
            width: 400.0,
            height: 400.0,
            iteration_multiple: 100,
            radius: 7,
        }
    }

    pub fn good() -> Self {
        let mut r = Self::fast();
        r.iteration_multiple *= 10;
        r
    }

    pub fn draw(&self, graph: &Graph) -> Vec<Vector> {
        if graph.nodes == 0 {
            return Vec::new();
        }

        let nodes = graph.nodes;
        let distances_raw = graph.all_pairs_shortest_paths();
        let mut distances_flat = vec![0.0; graph.nodes * graph.nodes];
        for i in 0..nodes {
            for j in 0..nodes {
                distances_flat[i * nodes + j] =
                    distances_raw[i][j].expect("graph to be connected") as f64;
            }
        }

        // k_neighborhoods[i] = Vec of (node, distance) pairs of its k-neighborhood, in arbitrary order
        let mut k_neighborhoods: Vec<Vec<(usize, f64)>> = Vec::with_capacity(nodes);
        for distances in distances_raw.iter() {
            let mut nbd = Vec::new();
            for (j, distance_option) in distances.iter().enumerate() {
                let distance = distance_option.expect("graph to be connected");
                if distance < self.radius {
                    nbd.push((j, distance as f64))
                }
            }
            k_neighborhoods.push(nbd);
        }

        let max_distance = distances_raw
            .iter()
            .flatten()
            .max()
            .expect("graph to be nonempty")
            .expect("graph to be nonempty") as f64;
        let draw_distance = self.width.min(self.height) / max_distance;

        let mut positions: Vec<Vector> = Vec::with_capacity(graph.nodes);
        InitializationStrategy::Random.initialize(
            graph.nodes,
            &mut positions,
            0.0,
            self.width,
            0.0,
            self.height,
        );

        // we want to store the jacobians/deltas in a way such that:
        // - we can easily find which vertex has the max delta
        // - we can quickly access, for updating, the deltas and jacobians for a given vertex
        // One option is to use a binary min-heap-map, sorted by deltas and mapping to jacobians, plus
        // an auxiliary array which maps vertex-index to heap-index
        // If the k-neighborhoods are small in comparison to the number of vertices, this is nice:
        // O(log(n)*size of k-nbds) per loop
        // Another option is just to use an array of jacobians, and have linear-time min-finding but
        // super easy updates
        // let's start with the array since it's simpler
        let mut jacobians: Vec<(f64, f64)> = (0..graph.nodes)
            .map(|i| jacobian_local(i, &positions, &k_neighborhoods[i], draw_distance))
            .collect();

        for _ in 0..self.iteration_multiple * nodes {
            let m = jacobians
                .iter()
                .map(|(dx, dy)| dx.powi(2) + dy.powi(2))
                .enumerate()
                .max_by(|(_, f), (_, g)| {
                    f.partial_cmp(g)
                        .expect("no non-comparable floats in deltas")
                })
                .expect("graph to be nonempty")
                .0;
            // println!("on interation {i}, had max delta index {m}");
            // println!("jacobians: {jacobians:?}");
            let hessian = hessian_m_local(m, &positions, &k_neighborhoods[m], draw_distance);
            let jacobian_value = jacobians[m];
            let hessian_det = det(hessian);
            let delta_x = det([
                [-jacobian_value.0, hessian[0][1]],
                [-jacobian_value.1, hessian[1][1]],
            ]) / hessian_det;
            let delta_y = det([
                [hessian[0][0], -jacobian_value.0],
                [hessian[1][0], -jacobian_value.1],
            ]) / hessian_det;
            debug_assert!(
                (delta_x * hessian[0][0] + delta_y * hessian[0][1] + jacobian_value.0)
                    / jacobian_value.0
                    < 0.01,
            );
            debug_assert!(
                (delta_x * hessian[1][0] + delta_y * hessian[1][1] + jacobian_value.1)
                    / jacobian_value.1
                    < 0.01,
            );
            positions[m].x += delta_x;
            positions[m].y += delta_y;
            for &(k, _) in k_neighborhoods[m].iter() {
                jacobians[k] = jacobian_local(k, &positions, &k_neighborhoods[k], draw_distance);
            }
        }
        positions
    }
}
