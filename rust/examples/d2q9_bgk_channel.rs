use arrayfire::*;
use std::time::Instant;

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

fn lbm() {
    // Grid length, number and spacing
    let nx: u64 = 700;
    let ny: u64 = 300;

    let total_nodes = nx * ny;

    // Physical parameters.
    let rho0: FloatNum = 1.0;

    let obstacle_x: u64 = nx / 5 + 1; // x location of the cylinder
    let obstacle_y: u64 = ny / 2 + ny / 30; // y location of the cylinder
    let obstacle_r: u64 = ny / 10 + 1; // radius of the cylinder

    // Reynolds number
    let re: FloatNum = 220.0;
    // Lattice speed
    let u_max: FloatNum = 0.1;
    let u_max_af = constant::<FloatNum>(u_max, dim4!(ny));
    // Kinematic viscosity
    let nu: FloatNum = u_max * 2.0 * obstacle_r as FloatNum / re;
    // Relaxation time
    let tau: FloatNum = (3.0 as FloatNum) * nu + (0.5 as FloatNum);
    // Relaxation parameter
    let omega: FloatNum = (1.0 as FloatNum) / tau;

    println!("Reynolds number: {}", re);
    println!("Lattice speed: {}", u_max);
    println!("Lattice viscosity: {}", nu);
    println!("Relaxation time: {}", tau);
    println!("Relaxation parameter: {}", omega);

    let t1: FloatNum = 4. / 9.;
    let t2: FloatNum = 1. / 9.;
    let t3: FloatNum = 1. / 36.;

    let x: Array<FloatNum> = tile(&range(dim4!(nx), 0), dim4!(1, ny));
    let y: Array<FloatNum> = tile(&range(dim4!(1, ny), 1), dim4!(nx, 1));

    let dims = dim4!(nx, ny);

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

    // Flow around obstacle
    // circle
    let mut bound = constant::<FloatNum>(1.0, dims);
    let r = constant::<FloatNum>(obstacle_r as FloatNum, dims);
    let r_sq = &r * &r;
    let circle = moddims(
        &le(
            &(pow(
                &(flat(&x) - obstacle_x as FloatNum),
                &(2.0 as FloatNum),
                false,
            ) + pow(
                &(flat(&y) - obstacle_y as FloatNum),
                &(2.0 as FloatNum),
                false,
            )),
            &flat(&r_sq),
            false,
        ),
        dims,
    );
    bound = selectr(&bound, &circle, 0.0 as f64);
    set_col(&mut bound, &constant::<FloatNum>(1.0, dim4!(nx)), 0); //top
    set_col(
        &mut bound,
        &constant::<FloatNum>(1.0, dim4!(nx)),
        ny as i64 - 1,
    ); //bottom

    // matrix offset of each Occupied Node
    let on = locate(&bound);

    // Bounceback indexes
    let to_reflect = flat(&tile(&on, dim4!(ci.elements() as u64)))
        + flat(&tile(&ci, dim4!(on.elements() as u64)));
    let reflected = flat(&tile(&on, dim4!(nbi.elements() as u64)))
        + flat(&tile(&nbi, dim4!(on.elements() as u64)));

    let mut density = constant::<FloatNum>(rho0, dims);
    let mut ux = constant::<FloatNum>(u_max, dims);
    let mut uy = constant::<FloatNum>(0.0, dims);

    let zeroed_on = constant::<FloatNum>(0.0, on.dims());

    eval!(ux[on] = zeroed_on);
    eval!(uy[on] = zeroed_on);
    eval!(density[on] = zeroed_on);

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
    let mut win = Window::new(1536, 768, "LBM solver using ArrayFire".to_string());
    win.grid(2, 1);

    let mut iter: u64 = 0;
    let maxiter: u64 = 10000;
    let mut mlups: Vec<FloatNum> = Vec::with_capacity(10000);

    sync(0);
    let timer = Instant::now();

    mem_info!("Before benchmark");

    while !win.is_closed() && iter < maxiter {
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

        // inlet speed
        eval!(ux[0:0:1, 1:1:0] = u_max_af);

        eval!(ux[on] = zeroed_on);
        eval!(uy[on] = zeroed_on);
        eval!(density[on] = zeroed_on);
        set_row(&mut density, &constant::<FloatNum>(1.0, dim4!(ny)), 0);
        set_row(
            &mut density,
            &constant::<FloatNum>(1.0, dim4!(ny)),
            nx as i64 - 1,
        );

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

            let filter_x = seq!(0, nx as i32 - 1, nx as i32 / 25);
            let filter_y = seq!(0, ny as i32 - 1, ny as i32 / 20);

            win.set_view(0, 0);
            win.set_colormap(ColorMap::COLORS);
            win.draw_image(
                &transpose(&normalize(&uu), false),
                Some(format!("XY domain in iteration {}", &iter).to_string()),
            );

            win.set_view(1, 0);
            win.set_axes_limits_2d(0.0, nx as f32, 0.0, ny as f32, true);
            win.draw_vector_field2(
                &flat(&view!(x[filter_x,filter_y])),
                &flat(&view!(y[filter_x,filter_y])),
                &flat(&view!(ux[filter_x,filter_y])),
                &flat(&view!(uy[filter_x,filter_y])),
                Some(format!("Velocity field in iteration {}", &iter).to_string()),
            );

            win.show();
        }

        sync(0);
        let time = timer.elapsed().as_secs() as FloatNum;
        mlups.push((total_nodes as FloatNum * iter as FloatNum * 10e-6) / time);

        if iter % 100 == 0 {
            println!(
                "{} iterations completed, {}s elapsed ({} MLUPS).",
                iter, time, mlups[iter as usize]
            );
        }

        iter += 1;
    }

    sync(0);
    mem_info!("After benchmark");
}

fn main() {
    set_device(0);
    set_backend(Backend::OPENCL);
    info();
    println!("LBM D2Q9 simulation\n");
    lbm();
}
