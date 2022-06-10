use svg::node::element::{Circle, Line};
use svg::{node, Document};

use graph_drawing::graph::{grid_graph, torus_graph, Graph};
use graph_drawing::layout::Vector;
use graph_drawing::StressMajorization;

fn main() {
    // let drawer = KamadaKawaiFastDrawer::good();
    let drawer = StressMajorization { tolerance: 10e-4 };

    let grid_size = 3;
    let graph = grid_graph(grid_size);
    println!("{:?}", graph.all_pairs_shortest_paths());
    let positions = drawer.draw(&graph, 400.0, 400.0);
    let svg_document = to_svg(&graph, &positions);
    svg::save("out/grid.svg", &svg_document).unwrap();

    // let graph = torus_graph(10, 40);
    // let positions = drawer.draw(&graph, 400.0, 400.0);
    // let svg_document = to_svg(&graph, &positions);
    // svg::save("out/torus.svg", &svg_document).unwrap();

    // let graph = torus_graph(8, 8);
    // let positions = drawer.draw(&graph, 400.0, 400.0);
    // let svg_document = to_svg(&graph, &positions);
    // svg::save("out/sq_torus.svg", &svg_document).unwrap();
}

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
