extern crate criterion;

use criterion::*;
use frolling::deal;

fn bench(c: &mut Criterion) {
    let mut group = c.benchmark_group("create proof");
    group.sample_size(10);
    group.bench_function("with_setup", move |b| {
        b.iter(|| {
            deal(16, 4);
        })
    });
}

criterion_group!(benches, bench);
criterion_main!(benches);
