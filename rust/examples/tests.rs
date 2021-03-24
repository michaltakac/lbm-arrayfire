use arrayfire::*;
use std::time::Instant;

fn first(tm: &Array<f32>, a: &Array<f32>) {
    let mut out = matmul(&tm, &a, MatProp::NONE, MatProp::NONE);
    out = out + 1.0f32;
    eval!(&out);
}

fn second(tm: &Array<f32>, a: &Array<f32>) {
    let mut out = matmul(&a, &tm, MatProp::NONE, MatProp::NONE);
    out = out + 1.0f32;
    eval!(&out);
}

fn main() {
    set_device(0);
    set_backend(Backend::CUDA);

    let nx: u64 = 200;
    let ny: u64 = 200;
    let nz: u64 = 200;
    let total_nodes = nx * ny * nz;

    let tm = randu(dim4!(27,27));
    let tm_t = transpose(&tm, false);
    let a = randu(dim4!(27, total_nodes as u64));
    let a_t = transpose(&a, false);

    let timer = Instant::now();

    for i in 0..10000 {
        first(&tm, &a);
        // second(&tm_t, &a_t); // TODO: IMPLEMENT THIS!!!
    }

    let total_time = timer.elapsed().as_secs_f32();
    println!("{}s elapsed.", total_time);
}