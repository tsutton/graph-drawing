use nalgebra::{Const, DMatrix, Dynamic, Matrix, MatrixXx2, OMatrix};

use crate::{graph::Graph, layout::Vector};

pub struct StressMajorization {
    // Paper recommends 10e-4
    pub tolerance: f64,
}

impl StressMajorization {
    pub fn draw(&self, graph: &Graph, width: f64, height: f64) -> Vec<Vector> {
        // strategy:
        // - compute APSP and L^w
        // - initialize X(0)
        // - Iterate until X(t+1) and X(t) are within tolerance:
        //   - Compute L^X(t)
        //   - Solve for X(t+1): L^w X(t+1)^(a) = L^X(t) X(t)^a, a = 1,2 (for 2-d embedding)
        // X(t) is n-by-2 matrix,
        let distances_raw = graph.all_pairs_shortest_paths();

        let nodes = graph.nodes;
        // Note: w_ij = 1/d_ij^2
        // Note: it's a bit more convenient to code the *last* vertex at (0,0), rather than the first like the paper
        let weighted_laplacan = DMatrix::from_fn(nodes - 1, nodes - 1, |row, col| {
            if row == col {
                distances_raw[row]
                    .iter()
                    .enumerate()
                    .filter(|(i, _)| *i != row)
                    .map(|(_, w)| 1.0 / (w.expect("graph to be connected") as f64).powi(2))
                    .sum::<f64>()
            } else {
                -1.0 / (distances_raw[row][col].expect("graph to be connected") as f64).powi(2)
            }
        });
        // used for a debug assertion below
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
                1.0 / distances_raw[row][col].unwrap() as f64
            }
        });
        for _ in 0..1000 {
            let l_x_prev = lz(&layout, &deltas);
            let sliced_lx = l_x_prev.slice((0, 0), (nodes - 1, nodes - 1));
            // .clone() explanation: seems like Mul takes ownership
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

// L^Z in the paper
// This function pretends that z has an extra (0,0) at the end, and returns a matching L^Z,
// so the caller should slice off the last row/column
fn lz(
    z: &MatrixXx2<f64>,
    deltas: &DMatrix<f64>, // n x n matrix of deltas_ij = delta_ij as in paper
) -> Matrix<f64, Dynamic, Dynamic, nalgebra::VecStorage<f64, nalgebra::Dynamic, nalgebra::Dynamic>>
{
    let nodes = z.shape().0 + 1;
    let inverses = DMatrix::from_fn(nodes, nodes, |row, col| {
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
            -1.0 / norm
        }
    });
    let mut ret = inverses.component_mul(deltas);
    for i in 0..z.shape().0 {
        ret[(i, i)] = -ret.column(i).sum()
    }
    ret
}
