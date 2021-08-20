use arrayfire::*;
use std::time::Instant;
use std::error::Error;
use csv::Writer;

type FloatNum = f32;

fn normalize(a: &Array<FloatNum>) -> Array<FloatNum> {
    (a / (max_all(&abs(a)).0 * 1.1 as FloatNum)) + 0.1 as FloatNum
}

fn stream(f: &Array<FloatNum>) -> Array<FloatNum> {
    let mut pdf = f.clone();
    eval!(pdf[1:1:0, 1:1:0, 1:1:1] = shift(&view!(f[1:1:0, 1:1:0, 1:1:1]), &[1, 0, 0, 0]));
    eval!(pdf[1:1:0, 1:1:0, 2:2:1] = shift(&view!(f[1:1:0, 1:1:0, 2:2:1]), &[0, 1, 0, 0]));
    eval!(pdf[1:1:0, 1:1:0, 3:3:1] = shift(&view!(f[1:1:0, 1:1:0, 3:3:1]), &[-1, 0, 0, 0]));
    eval!(pdf[1:1:0, 1:1:0, 4:4:1] = shift(&view!(f[1:1:0, 1:1:0, 4:4:1]), &[0, -1, 0, 0]));
    eval!(pdf[1:1:0, 1:1:0, 5:5:1] = shift(&view!(f[1:1:0, 1:1:0, 5:5:1]), &[1, 1, 0, 0]));
    eval!(pdf[1:1:0, 1:1:0, 6:6:1] = shift(&view!(f[1:1:0, 1:1:0, 6:6:1]), &[-1, 1, 0, 0]));
    eval!(pdf[1:1:0, 1:1:0, 7:7:1] = shift(&view!(f[1:1:0, 1:1:0, 7:7:1]), &[-1, -1, 0, 0]));
    eval!(pdf[1:1:0, 1:1:0, 8:8:1] = shift(&view!(f[1:1:0, 1:1:0, 8:8:1]), &[1, -1, 0, 0]));
    pdf
}

fn output_csv(mlups: Vec<FloatNum>, size: u64) -> Result<(), Box<dyn Error>> {
  let mut wtr = Writer::from_path(format!("benchmarks/GPU_NAME_d2q9_bgk_lid_mlups_{}.csv",size))?;

  wtr.write_record(&["Iterations", "MLUPS"])?;
  for (i, item) in mlups.iter().enumerate() {
    wtr.write_record(&[i.to_string(), item.to_string()])?;
  }

  wtr.flush()?;
  Ok(())
}

fn lbm(write_csv: bool) {
    // Grid length, number and spacing
    let nx: u64 = 256;
    let ny: u64 = 256;

    let total_nodes = nx * ny;

    // Physical parameters.
    let ux_lid: FloatNum = 0.05; // horizontal lid velocity
    let uy_lid: FloatNum = 0.0; // vertical lid velocity
    let rho0: FloatNum = 1.0;

    // Reynolds number
    let re: FloatNum = 100.0;
    // Kinematic viscosity
    let nu: FloatNum = ux_lid * 2.0 * nx as FloatNum / re;
    // Relaxation time
    let tau: FloatNum = (3.0 as FloatNum) * nu + (0.5 as FloatNum);
    // Relaxation parameter
    let omega: FloatNum = (1.0 as FloatNum) / tau;

    println!("Horizontal lid velocity ux_lid: {}", ux_lid);
    println!("Vertical lid velocity uy_lid: {}", uy_lid);
    println!("Reynolds number: {}", re);
    println!("Lattice viscosity: {}", nu);
    println!("Relaxation time: {}", tau);
    println!("Relaxation parameter: {}", omega);

    let t1: FloatNum = 4. / 9.;
    let t2: FloatNum = 1. / 9.;
    let t3: FloatNum = 1. / 36.;

    let x: Array<FloatNum> = tile(&range(dim4!(nx), 0), dim4!(1, ny));
    let y: Array<FloatNum> = tile(&range(dim4!(1, ny), 1), dim4!(nx, 1));

    let dims = dim4!(nx, ny);

    let ux_lid_af = constant::<FloatNum>(ux_lid, dims);
    let uy_lid_af = constant::<FloatNum>(uy_lid, dims);

    let lid = seq!(1, nx as i32 - 2, 1);
    let end_y = seq!(nx as i32 - 1, ny as i32 - 1, 1);

    //  c6  c2   c5
    //    \  |  /
    //  c3 -c0 - c1
    //    /  |  \
    //  c7  c4   c8
    // Discrete velocities
    let ex = Array::<FloatNum>::new(&[0., 1., 0., -1., 0., 1., -1., -1., 1.], dim4!(9));
    let ey = Array::<FloatNum>::new(&[0., 0., 1., 0., -1., 1., 1., -1., -1.], dim4!(9));

    // weights
    let w = Array::new(&[t1, t2, t2, t2, t2, t3, t3, t3, t3], dim4!(9));

    let ci: Array<u64> = (range::<u64>(dim4!(1, 8), 1) + 1) * total_nodes;
    let nbidx = Array::new(&[2, 3, 0, 1, 6, 7, 4, 5], dim4!(8));
    let span = seq!();
    let nbi: Array<u64> = view!(ci[span, nbidx]);

    let main_index = moddims(&range(dim4!(total_nodes * 9), 0), dim4!(nx, ny, 9));
    let nb_index = flat(&stream(&main_index));

    // Open lid
    let mut bound = constant::<FloatNum>(1.0, dims);
    let zeros = constant::<FloatNum>(0.0, dims);
    let all_except_top_lid = seq!(1, ny as i32 - 1, 1);
    assign_seq(
        &mut bound,
        &[lid, all_except_top_lid],
        &index(&zeros, &[lid, all_except_top_lid]),
    );

    // matrix offset of each Occupied Node
    let on = locate(&bound);

    // Bounceback indexes
    let to_reflect = flat(&tile(&on, dim4!(ci.elements() as u64)))
        + flat(&tile(&ci, dim4!(on.elements() as u64)));
    let reflected = flat(&tile(&on, dim4!(nbi.elements() as u64)))
        + flat(&tile(&nbi, dim4!(on.elements() as u64)));

    let mut density = constant::<FloatNum>(rho0, dims);
    let mut ux = constant::<FloatNum>(0.0, dims);
    let mut uy = constant::<FloatNum>(0.0, dims);

    let zeroed_on = constant::<FloatNum>(0.0, on.dims());

    // Start in equilibrium state
    let mut u_sq: Array<FloatNum> =
        flat(&(pow(&ux, &(2.0 as FloatNum), false) + pow(&uy, &(2.0 as FloatNum), false)));
    let mut eu: Array<FloatNum> = flat(
        &(&mul(&transpose(&ex, false), &flat(&ux), true)
            + &mul(&transpose(&ey, false), &flat(&uy), true)),
    );
    let mut f: Array<FloatNum> = flat(&mul(&transpose(&w, false), &flat(&density), true))
        * ((1.0 as FloatNum)
            + (3.0 as FloatNum) * &eu
            + (4.5 as FloatNum) * (&pow(&eu, &(2.0 as FloatNum), false))
            - (1.5 as FloatNum) * (&tile(&flat(&u_sq), dim4!(9))));

    // Create a window to show the waves.
    // let mut win = Window::new(1536, 768, "LBM solver using ArrayFire".to_string());
    // win.grid(1, 2);

    let mut iter: u64 = 0;
    let maxiter: u64 = 5000;
    let mut mlups: Vec<FloatNum> = Vec::with_capacity(5000);

    sync(0);
    let timer = Instant::now();

    mem_info!("Before benchmark");

    while iter < maxiter {
        // Streaming by reading from neighbors (with pre-built index) - pull scheme
        let f_streamed = view!(f[nb_index]);

        let bouncedback = view!(f_streamed[to_reflect]); // Densities bouncing back at next timestep

        let f_2d = moddims(&f_streamed, dim4!(total_nodes, 9));

        // Compute macroscopic variables
        let rho = sum(&f_2d, 1);
        density = moddims(&rho, dims);

        let fex = mul(&transpose(&ex, false), &f_2d, true);
        let fey = mul(&transpose(&ey, false), &f_2d, true);

        ux = moddims(&(sum(&fex, 1) / &rho), dims);
        uy = moddims(&(sum(&fey, 1) / &rho), dims);

        // Macroscopic (Dirichlet) boundary conditions
        eval!(ux[lid, end_y] = view!(ux_lid_af[lid, end_y]));
        eval!(uy[lid, end_y] = view!(uy_lid_af[lid, end_y]));

        eval!(ux[on] = zeroed_on);
        eval!(uy[on] = zeroed_on);
        eval!(density[on] = zeroed_on);

        // Collision
        u_sq = flat(&(pow(&ux, &(2.0 as FloatNum), false) + pow(&uy, &(2.0 as FloatNum), false)));
        eu = flat(
            &(&mul(&transpose(&ex, false), &flat(&ux), true)
                + &mul(&transpose(&ey, false), &flat(&uy), true)),
        );
        let feq = flat(&mul(&transpose(&w, false), &flat(&density), true))
            * ((1.0 as FloatNum)
                + (3.0 as FloatNum) * &eu
                + (4.5 as FloatNum) * (&pow(&eu, &(2.0 as FloatNum), false))
                - (1.5 as FloatNum) * (&tile(&flat(&u_sq), dim4!(9))));

        f = omega * feq + (1.0 - omega) * f_streamed;

        eval!(f[reflected] = bouncedback);

        // Visualization
        if iter % 10 == 0 {
            let mut uu = moddims(&sqrt(&u_sq), dims);
            eval!(uu[on] = constant::<FloatNum>(FloatNum::NAN, on.dims()));

            let filter = seq!(0, nx as i32 - 1, nx as i32 / 30);

            win.set_view(0, 0);
            win.set_colormap(ColorMap::SPECTRUM);
            win.draw_image(
                &flip(&transpose(&normalize(&uu), false), 0),
                Some(format!("XY domain in iteration {}", &iter).to_string()),
            );

            win.set_view(0, 1);
            win.set_axes_limits_2d(0.0, nx as f32, 0.0, ny as f32, true);
            win.draw_vector_field2(
                &flat(&view!(x[filter,filter])),
                &flat(&view!(y[filter,filter])),
                &flat(&view!(ux[filter,filter])),
                &flat(&view!(uy[filter,filter])),
                Some(format!("Velocity field in iteration {}", &iter).to_string()),
            );

            win.show();
        }

        let time = timer.elapsed().as_secs_f32();
        let updates = (total_nodes as FloatNum * iter as FloatNum * 10e-6) / time;

        if !updates.is_nan() && !updates.is_infinite() {
          mlups.push(updates);
        }

        if iter % 100 == 0 {
            println!(
                "{} iterations completed, {}s elapsed ({} MLUPS).",
                iter, time, updates
            );
        }

        iter += 1;
        sync(0);
    }

    mem_info!("After benchmark");
    sync(0);

    // output CSV of MLUPS data
    if write_csv {
      output_csv(mlups, nx);
    }
}

fn main() {
    set_backend(Backend::OPENCL);
    set_device(0);
    info();
    println!("LBM D2Q9 simulation\n");
    let write_csv = false;
    lbm(write_csv);
}
