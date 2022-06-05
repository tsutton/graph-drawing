use layout::Vector;
use rand::prelude::SliceRandom;
use std::f64::consts::PI;

pub mod graph;
pub mod layout;

mod eades;
pub use eades::EadesDrawer;

mod fruchterman_reingold;
pub use fruchterman_reingold::FruchtermanReingoldDrawer;

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
