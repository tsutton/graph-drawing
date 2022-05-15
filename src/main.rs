use svg::node::element::Circle;
use svg::{node, Document};

use graph_drawing::graph::{cycle_graph, grid_graph, Graph};
use graph_drawing::layout::Vector;
use graph_drawing::EadesDrawer;

fn main() {
    let drawer = EadesDrawer::my_best_guess();

    let graph = grid_graph(3);
    let positions = drawer.draw(&graph);
    let svg_document = to_svg(&graph, &positions);
    svg::save("grid.svg", &svg_document).unwrap();

    let graph = cycle_graph(10);
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
    document
}
