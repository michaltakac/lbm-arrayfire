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
    // nearest-neighbours
    eval!(pdf[1:1:0, 1:1:0, 1:1:1] = shift(&view!(f[1:1:0, 1:1:0, 1:1:1]), &[ 1, 0, 0, 0]));
    eval!(pdf[1:1:0, 1:1:0, 2:2:1] = shift(&view!(f[1:1:0, 1:1:0, 2:2:1]), &[ -1, 0, 0, 0]));
    eval!(pdf[1:1:0, 1:1:0, 3:3:1] = shift(&view!(f[1:1:0, 1:1:0, 3:3:1]), &[ 0, 1, 0, 0]));
    eval!(pdf[1:1:0, 1:1:0, 4:4:1] = shift(&view!(f[1:1:0, 1:1:0, 4:4:1]), &[ 0, -1, 0, 0]));
    eval!(pdf[1:1:0, 1:1:0, 5:5:1] = shift(&view!(f[1:1:0, 1:1:0, 5:5:1]), &[ 0, 0, 1, 0]));
    eval!(pdf[1:1:0, 1:1:0, 6:6:1] = shift(&view!(f[1:1:0, 1:1:0, 6:6:1]), &[ 0, 0, -1, 0]));
    // next-nearest neighbours
    // xy plane
    eval!(pdf[1:1:0, 1:1:0, 7:7:1] = shift(&view!(f[1:1:0, 1:1:0, 7:7:1]), &[ 1, 1, 0, 0]));
    eval!(pdf[1:1:0, 1:1:0, 8:8:1] = shift(&view!(f[1:1:0, 1:1:0, 8:8:1]), &[ -1, 1, 0, 0]));
    eval!(pdf[1:1:0, 1:1:0, 9:9:1] = shift(&view!(f[1:1:0, 1:1:0, 9:9:1]), &[ 1, -1, 0, 0]));
    eval!(pdf[1:1:0, 1:1:0, 10:10:1] = shift(&view!(f[1:1:0, 1:1:0, 10:10:1]), &[ -1, -1, 0, 0]));
    // xz plane
    eval!(pdf[1:1:0, 1:1:0, 11:11:1] = shift(&view!(f[1:1:0, 1:1:0, 11:11:1]), &[ 1, 0, 1, 0]));
    eval!(pdf[1:1:0, 1:1:0, 12:12:1] = shift(&view!(f[1:1:0, 1:1:0, 12:12:1]), &[ -1, 0, 1, 0]));
    eval!(pdf[1:1:0, 1:1:0, 13:13:1] = shift(&view!(f[1:1:0, 1:1:0, 13:13:1]), &[ 1, 0, -1, 0]));
    eval!(pdf[1:1:0, 1:1:0, 14:14:1] = shift(&view!(f[1:1:0, 1:1:0, 14:14:1]), &[ -1, 0, -1, 0]));
    // yz plane
    eval!(pdf[1:1:0, 1:1:0, 15:15:1] = shift(&view!(f[1:1:0, 1:1:0, 15:15:1]), &[ 0, 1, 1, 0]));
    eval!(pdf[1:1:0, 1:1:0, 16:16:1] = shift(&view!(f[1:1:0, 1:1:0, 16:16:1]), &[ 0, -1, 1, 0]));
    eval!(pdf[1:1:0, 1:1:0, 17:17:1] = shift(&view!(f[1:1:0, 1:1:0, 17:17:1]), &[ 0, 1, -1, 0]));
    eval!(pdf[1:1:0, 1:1:0, 18:18:1] = shift(&view!(f[1:1:0, 1:1:0, 18:18:1]), &[ 0, -1, -1, 0]));
    // next next-nearest neighbours
    eval!(pdf[1:1:0, 1:1:0, 19:19:1] = shift(&view!(f[1:1:0, 1:1:0, 19:19:1]), &[ 1, 1, 1, 0]));
    eval!(pdf[1:1:0, 1:1:0, 20:20:1] = shift(&view!(f[1:1:0, 1:1:0, 20:20:1]), &[ -1, 1, 1, 0]));
    eval!(pdf[1:1:0, 1:1:0, 21:21:1] = shift(&view!(f[1:1:0, 1:1:0, 21:21:1]), &[ 1, -1, 1, 0]));
    eval!(pdf[1:1:0, 1:1:0, 22:22:1] = shift(&view!(f[1:1:0, 1:1:0, 22:22:1]), &[ -1, -1, 1, 0]));
    eval!(pdf[1:1:0, 1:1:0, 23:23:1] = shift(&view!(f[1:1:0, 1:1:0, 23:23:1]), &[ 1, 1, -1, 0]));
    eval!(pdf[1:1:0, 1:1:0, 24:24:1] = shift(&view!(f[1:1:0, 1:1:0, 24:24:1]), &[ -1, 1, -1, 0]));
    eval!(pdf[1:1:0, 1:1:0, 25:25:1] = shift(&view!(f[1:1:0, 1:1:0, 25:25:1]), &[ 1, -1, -1, 0]));
    eval!(pdf[1:1:0, 1:1:0, 26:26:1] = shift(&view!(f[1:1:0, 1:1:0, 26:26:1]), &[ -1, -1, -1, 0]));

    pdf
}

fn output_csv(mlups: Vec<f32>) -> Result<(), Box<dyn Error>> {
  let mut wtr = Writer::from_path("d3q27_bgk_lid_mlups.csv")?;

  wtr.write_record(&["Iterations", "MLUPS"])?;
  for (i, item) in mlups.iter().enumerate() {
    wtr.write_record(&[i.to_string(), item.to_string()])?;
  }

  wtr.flush()?;
  Ok(())
}

fn lbm(write_csv: bool) {
    // Grid length, number and spacing
    let nx: u64 = 100;
    let ny: u64 = 100;
    let nz: u64 = 100;

    let total_nodes = nx * ny * nz;

    // Physical parameters.
    let ux_lid: FloatNum = 0.05; // horizontal lid velocity
    let uy_lid: FloatNum = 0.0; // vertical lid velocity
    let uz_lid: FloatNum = 0.0; // lid velocity in z-direction
    let rho0: FloatNum = 1.0;

    // Reynolds number
    let re: FloatNum = 100.0;
    // Kinematic viscosity
    let nu: FloatNum = ux_lid * nx as FloatNum / re;
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

    let t1: FloatNum = 8. / 27.;
    let t2: FloatNum = 2. / 27.;
    let t3: FloatNum = 1. / 54.;
    let t4: FloatNum = 1. / 216.;

    let x: Array<FloatNum> = tile(&range(dim4!(nx), 1), dim4!(ny * nz));
    let y: Array<FloatNum> = tile(&range(dim4!(1, ny), 1), dim4!(nx, nz));
    let z: Array<FloatNum> = tile(&range(dim4!(1, nz), 1), dim4!(nx, ny));

    let dims = dim4!(nx, ny, nz);
    let span = seq!();

    let ux_lid_af = constant::<FloatNum>(ux_lid, dims);
    let uy_lid_af = constant::<FloatNum>(uy_lid, dims);
    let uz_lid_af = constant::<FloatNum>(uz_lid, dims);

    let lidx = seq!(1, nx as i32 - 2, 1);
    let lidy = seq!(1, ny as i32 - 2, 1);
    let end_y = seq!(nx as i32 - 1, ny as i32 - 1, 1);

    // Discrete velocities
    let ex = Array::<FloatNum>::new(&[0., 1., 0., -1., 0., 1., -1., -1., 1.], dim4!(27));
    let ey = Array::<FloatNum>::new(&[0., 0., 1., 0., -1., 1., 1., -1., -1.], dim4!(27));
    let ez = Array::<FloatNum>::new(&[0., 0., 1., 0., -1., 1., 1., -1., -1.], dim4!(27));

    // weights
    let w = Array::new(&[t1,t2,t2,t2,t2,t2,t2,t3,t3,t3,t3,t3,t3,t3,t3,t3,t3,t3,t3,t4,t4,t4,t4,t4,t4,t4,t4], dim4!(27));

    let ci: Array<u64> = (range::<u64>(dim4!(1, 26), 1) + 1) * total_nodes;
    let nbidx = Array::new(&[1,0,3,2,5,4,9,8,7,6,13,12,11,10,17,16,15,14,25,24,23,22,21,20,19,18], dim4!(26));
    let nbi: Array<u64> = view!(ci[span, nbidx]);

    let main_index = moddims(&range(dim4!(total_nodes * 27), 0), dim4!(nx, ny, nz, 27));
    let nb_index = flat(&stream(&main_index));

    // Open lid 3D
    let mut bound = constant::<FloatNum>(1.0, dims);
    let zeros = constant::<FloatNum>(0.0, dims);
    let all_except_top_lid = seq!(1, nz as i32 - 1, 1);
    assign_seq(
        &mut bound,
        &[lidx, lidy, all_except_top_lid],
        &index(&zeros, &[lidx, lidy, all_except_top_lid]),
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
    let mut uz = constant::<FloatNum>(0.0, dims);

    let zeroed_on = constant::<FloatNum>(0.0, on.dims());

    // Start in equilibrium state
    let mut u_sq: Array<FloatNum> =
        flat(&(pow(&ux, &(2.0 as FloatNum), false) + pow(&uy, &(2.0 as FloatNum), false) + pow(&uz, &(2.0 as FloatNum), false)));
    let mut eu: Array<FloatNum> = flat(
        &(&mul(&transpose(&ex, false), &flat(&ux), true)
            + &mul(&transpose(&ey, false), &flat(&uy), true)
            + &mul(&transpose(&ez, false), &flat(&uy), true)),
    );
    let mut f: Array<FloatNum> = flat(&mul(&transpose(&w, false), &flat(&density), true))
        * ((1.0 as FloatNum)
            + (3.0 as FloatNum) * &eu
            + (4.5 as FloatNum) * (&pow(&eu, &(2.0 as FloatNum), false))
            - (1.5 as FloatNum) * (&tile(&flat(&u_sq), dim4!(27))));

    // Create a window to show the waves.
    // let mut win = Window::new(1536, 768, "LBM solver using ArrayFire".to_string());
    // win.grid(2, 2);

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

        let f_2d = moddims(&f_streamed, dim4!(total_nodes, 27));

        // Compute macroscopic variables
        let rho = sum(&f_2d, 1);
        density = moddims(&rho, dims);

        let fex = mul(&transpose(&ex, false), &f_2d, true);
        let fey = mul(&transpose(&ey, false), &f_2d, true);
        let fez = mul(&transpose(&ez, false), &f_2d, true);

        ux = moddims(&(sum(&fex, 1) / &rho), dims);
        uy = moddims(&(sum(&fey, 1) / &rho), dims);
        uz = moddims(&(sum(&fez, 1) / &rho), dims);

        // Macroscopic (Dirichlet) boundary conditions
        eval!(ux[lidx, lidy, end_y] = view!(ux_lid_af[lidx, lidy, end_y]));
        eval!(uy[lidx, lidy, end_y] = view!(uy_lid_af[lidx, lidy, end_y]));
        eval!(uz[lidx, lidy, end_y] = view!(uz_lid_af[lidx, lidy, end_y]));

        eval!(ux[on] = zeroed_on);
        eval!(uy[on] = zeroed_on);
        eval!(uz[on] = zeroed_on);
        eval!(density[on] = zeroed_on);

        // Collision
        u_sq = flat(&(pow(&ux, &(2.0 as FloatNum), false) + pow(&uy, &(2.0 as FloatNum), false) + pow(&uz, &(2.0 as FloatNum), false)));
        eu = flat(
            &(&mul(&transpose(&ex, false), &flat(&ux), true)
                + &mul(&transpose(&ey, false), &flat(&uy), true)
                + &mul(&transpose(&ez, false), &flat(&uy), true)),
        );
        let feq = flat(&mul(&transpose(&w, false), &flat(&density), true))
            * ((1.0 as FloatNum)
                + (3.0 as FloatNum) * &eu
                + (4.5 as FloatNum) * (&pow(&eu, &(2.0 as FloatNum), false))
                - (1.5 as FloatNum) * (&tile(&flat(&u_sq), dim4!(27))));

        f = omega * feq + (1.0 - omega) * f_streamed;

        eval!(f[reflected] = bouncedback);

        // Visualization
        // if iter % 10 == 0 {
        //     let mut uu = moddims(&sqrt(&u_sq), dims);
        //     eval!(uu[on] = constant::<FloatNum>(FloatNum::NAN, on.dims()));

        //     let filter = seq!(0, nx as i32 - 1, nx as i32 / 30);
        //     let z_section = seq!(nz as i32 / 2, nz as i32 / 2, 1);

        //     let xy_view = index(&reorder_v2(&normalize(&uu), 1, 0, Some(vec![2])), &[span, span, z_section]);

        //     win.set_view(0, 0);
        //     win.set_colormap(ColorMap::SPECTRUM);
        //     win.draw_image(
        //       &flip(&transpose(&xy_view, false), 0),
        //       Some(format!("XY domain in iteration {}", &iter).to_string()),
        //     );

        //     // win.set_view(0, 1);
        //     // win.set_axes_limits_2d(0.0, nx as f32, 0.0, ny as f32, true);
        //     // win.draw_vector_field2(
        //     //     &flat(&view!(x[filter,filter])),
        //     //     &flat(&view!(y[filter,filter])),
        //     //     &flat(&view!(ux[filter,filter,z_section])),
        //     //     &flat(&view!(uy[filter,filter,z_section])),
        //     //     Some(format!("Velocity field in iteration {}", &iter).to_string()),
        //     // );

        //     win.set_view(1, 0);
        //     win.set_colormap(ColorMap::SPECTRUM);
        //     win.draw_image(
        //       &normalize(&index(&uu, &[span, span, z_section])),
        //       Some(format!("XY domain in iteration {}", &iter).to_string()),
        //     );

        //     win.show();
        // }

        let time = timer.elapsed().as_secs() as FloatNum;
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
      output_csv(mlups);
    }
}

fn main() {
    set_device(0);
    set_backend(Backend::OPENCL);
    info();
    println!("LBM D2Q9 simulation\n");
    let write_csv = true;
    lbm(write_csv);
}
