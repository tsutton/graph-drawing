use svg::node::element::{Circle, Line};
use svg::{node, Document};

use graph_drawing::graph::{cycle_graph, grid_graph, Graph};
use graph_drawing::layout::Vector;
use graph_drawing::KamadaKawaiDrawer;

fn main() {
    let drawer = KamadaKawaiDrawer::my_best_guess();

    // It seems like regular polygon is almost pathological for grid, because
    // it ends up *twisted*, as a result it's really slow to get to the "real" answer
    let grid_size = 5;
    let graph = grid_graph(grid_size);
    // for ((i, j), w) in graph.all_pairs_shortest_paths() {
    //     println!(
    //         "({}, {})->({}, {}): {}",
    //         i % grid_size,
    //         (i - i % grid_size) / grid_size,
    //         j % grid_size,
    //         (j - j % grid_size) / grid_size,
    //         w
    //     )
    // }
    let positions = drawer.draw(&graph);
    let svg_document = to_svg(&graph, &positions);
    svg::save("grid.svg", &svg_document).unwrap();

    let graph = cycle_graph(10);
    // for ((i, j), w) in graph.all_pairs_shortest_paths() {
    //     println!("{}->{}: {}", i, j, w)
    // }
    let positions = drawer.draw(&graph);
    let svg_document = to_svg(&graph, &positions);
    svg::save("cycle.svg", &svg_document).unwrap();
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
    for ((source, target), _) in graph.edges() {
        let path = Line::new()
            .set("x1", positions[source].x)
            .set("x2", positions[target].x)
            .set("y1", positions[source].y)
            .set("y2", positions[target].y)
            .set("stroke", "black")
            .set("stroke-width", 1);
        document = document.add(path);
    }
    document
}
