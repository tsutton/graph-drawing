use nalgebra::{Const, DMatrix, Dynamic, Matrix, MatrixXx2, OMatrix};

use crate::{graph::Graph, layout::Vector};

/// Draws graphs using the Stress Majorization method from \[GKN05\]
///
/// This method uses finds a minimum of an energy function that tries to make each pair of nodes
/// have physical distance equal to the shortest-path distance (similar to Kamada Kawai), putting higher
/// value on having nodes with smaller path distance have the correct physical distance.
/// It uses a clever mathematical trick to turn this optimization into repeatedly solving matrix equations,
/// for which there are highly optimized algorithms and implementations.
pub struct StressMajorization {
    // Paper recommends 10e-4
    pub tolerance: f64,
}

// TODO benchmark conjugate gradient instead of cholesky
// TODO implement weighted edge lengths (paper section 3)
// TODO implement sparse energy functions (paper section 4)
// TODO (maybe?) implement restricted subspace and smart initialization (paper section 5)
impl StressMajorization {
    pub fn draw(&self, graph: &Graph, width: f64, height: f64) -> Vec<Vector> {
        // strategy:
        // - compute APSP and L^w
        // - initialize X(0)
        // - Iterate until X(t+1) and X(t) are within tolerance:
        //   - Compute L^X(t)
        //   - Solve for X(t+1): L^w X(t+1)^(a) = L^X(t) X(t)^a, a = 1,2 (for 2-d embedding)
        // X(t) is n-by-2 matrix, whose rows are the positions.
        let distances_raw = graph.all_pairs_shortest_paths();
        let distances = distances_raw
            .into_iter()
            .map(|distances_single_row| {
                distances_single_row
                    .into_iter()
                    .map(|w| w.expect("graph to be connected") as f64)
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();

        let nodes = graph.nodes;
        // Note: w_ij = 1/d_ij^2
        // Note: it's a bit more convenient to code the *last* vertex at (0,0), rather than the first like the paper
        let weighted_laplacan = DMatrix::from_fn(nodes - 1, nodes - 1, |row, col| {
            if row == col {
                distances[row]
                    .iter()
                    .enumerate()
                    .filter(|(i, _)| *i != row)
                    .map(|(_, w)| 1.0 / w.powi(2))
                    .sum::<f64>()
            } else {
                -1.0 / distances[row][col].powi(2)
            }
        });

        // Used for a debug assertion below.
        let original_weighted_laplacian = weighted_laplacan.clone();

        let decomposition = weighted_laplacan
            .cholesky()
            .expect("laplacian to be positive definite");

        // TODO as written, the alg sets width to d_ij, so we should start with them around max(d_ij)
        let mut layout: MatrixXx2<f64> = MatrixXx2::from_fn(nodes - 1, |_, col| {
            let scale = if col == 0 {
                width as f64
            } else {
                height as f64
            };
            rand::random::<f64>() * scale
        });
        let deltas = DMatrix::from_fn(nodes, nodes, |row, col| {
            if row == col {
                0.0
            } else {
                1.0 / distances[row][col]
            }
        });
        let mut current_stress = stress(&distances, &layout);
        loop {
            let l_x_prev = lz(&layout, &deltas);
            let sliced_lx = l_x_prev.slice((0, 0), (nodes - 1, nodes - 1));
            let x_coord_rhs = sliced_lx * layout.column(0);
            let new_x_coords = decomposition.solve(&x_coord_rhs);
            let y_coord_rhs = sliced_lx * layout.column(1);
            let new_y_coords: OMatrix<f64, Dynamic, Const<1>> = decomposition.solve(&y_coord_rhs);
            debug_assert!(
                ((original_weighted_laplacian.clone())
                    * MatrixXx2::from_columns(&[new_x_coords.clone(), new_y_coords.clone()])
                    - sliced_lx * layout)
                    .norm()
                    < 10e-4
            );
            layout = MatrixXx2::from_columns(&[new_x_coords, new_y_coords]);
            let new_stress = stress(&distances, &layout);
            let relative_diff = (current_stress - new_stress) / current_stress;
            if relative_diff.abs() < self.tolerance {
                break;
            }
            current_stress = new_stress;
        }
        let min_x = layout.column(0).min().min(0.0);
        let min_y = layout.column(1).min().min(0.0);
        let max_x = layout.column(0).max().max(0.0);
        let max_y = layout.column(1).max().max(0.0);
        let x_scale = width / (max_x - min_x);
        let y_scale = height / (max_y - min_y);
        let mut positions: Vec<Vector> = layout
            .row_iter()
            .map(|row| ((row[0] - min_x) * x_scale, (row[1] - min_y) * y_scale).into())
            .collect();
        positions.push((-min_x * x_scale, -min_y * y_scale).into());
        positions
    }
}

// assumptions: distances.len() == distances[i].len() == layout.shape().0 + 1, layout has a magic "0" at the end
fn stress(distances: &[Vec<f64>], layout: &MatrixXx2<f64>) -> f64 {
    let mut sum = 0.0;
    for i in 0..distances.len() {
        for j in 0..i {
            let distance = distances[i][j];
            let w = 1.0 / (distance * distance);
            let stress_ij = if i == distances.len() - 1 {
                w * (layout.row(j).norm() - distance * distance)
            } else {
                w * ((layout.row(j) - layout.row(i)).norm() - distance * distance)
            };
            sum += stress_ij;
        }
    }
    sum
}

// L^Z in the paper
// This function pretends that z has an extra (0,0) at the end, and returns a matching L^Z,
// so the caller should slice off the last row/column
// Forcing the caller to slice it off simplifies ownership: this function needs to return an owned matrix,
// which means it would need to allocate the non-truncated matrix, then allocate another matrix for the truncated matrix
// to return. But if we force the caller to slice off the last bit, the truncated matrix doesn't require an allocation since
// it can just refer to the allocation of the non-truncated one via borrowing.
fn lz(
    z: &MatrixXx2<f64>,    // (n-1) x 2
    deltas: &DMatrix<f64>, // n x n matrix of deltas_ij = delta_ij as in paper
) -> Matrix<f64, Dynamic, Dynamic, nalgebra::VecStorage<f64, nalgebra::Dynamic, nalgebra::Dynamic>>
{
    let nodes = z.shape().0 + 1;
    let mut ret = DMatrix::from_fn(nodes, nodes, |row, col| {
        // manifest a magical row of (0.0, 0.0) at z.row(nodes-1)
        let norm = if row == col {
            0.0
        } else if row == nodes - 1 {
            z.row(col).norm()
        } else if col == nodes - 1 {
            z.row(row).norm()
        } else {
            (z.row(row) - z.row(col)).norm()
        };
        if norm == 0.0 {
            0.0
        } else {
            -deltas[(row, col)] / norm
        }
    });
    for i in 0..z.shape().0 {
        ret[(i, i)] = -ret.column(i).sum()
    }
    ret
}
