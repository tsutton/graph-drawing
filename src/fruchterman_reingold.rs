use super::InitializationStrategy;
use crate::graph::Graph;
use crate::layout::Vector;

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
