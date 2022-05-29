use criterion::{criterion_group, criterion_main, Criterion};
use graph_drawing::graph::grid_graph;
use graph_drawing::KamadaKawaiDrawer;

pub fn criterion_benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("grid graph");
    group.bench_function("grid graph 5", |b| {
        b.iter(|| {
            let drawer = KamadaKawaiDrawer::for_benchmarking();
            let grid_size = 5;
            let graph = grid_graph(grid_size);
            let _positions = drawer.draw(&graph);
        })
    });
    group.finish()
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
