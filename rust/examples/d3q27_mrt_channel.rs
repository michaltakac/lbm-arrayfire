use arrayfire::*;
use std::time::Instant;
use std::error::Error;
use csv::Writer;

// TODO: WORK IN PROGRESS!!!

type FloatNum = f32;

fn normalize(a: &Array<FloatNum>) -> Array<FloatNum> {
    (a / (max_all(&abs(a)).0 * 1.1 as FloatNum)) + 0.1 as FloatNum
}

fn stream(f: &Array<FloatNum>) -> Array<FloatNum> {
    let mut pdf = f.clone();
    // nearest-neighbours
    eval!(pdf[1:1:0, 1:1:0, 1:1:0, 1:1:1] = shift(&view!(f[1:1:0, 1:1:0, 1:1:0, 1:1:1]), &[ 1, 0, 0, 0]));
    eval!(pdf[1:1:0, 1:1:0, 1:1:0, 2:2:1] = shift(&view!(f[1:1:0, 1:1:0, 1:1:0, 2:2:1]), &[ -1, 0, 0, 0]));
    eval!(pdf[1:1:0, 1:1:0, 1:1:0, 3:3:1] = shift(&view!(f[1:1:0, 1:1:0, 1:1:0, 3:3:1]), &[ 0, 1, 0, 0]));
    eval!(pdf[1:1:0, 1:1:0, 1:1:0, 4:4:1] = shift(&view!(f[1:1:0, 1:1:0, 1:1:0, 4:4:1]), &[ 0, -1, 0, 0]));
    eval!(pdf[1:1:0, 1:1:0, 1:1:0, 5:5:1] = shift(&view!(f[1:1:0, 1:1:0, 1:1:0, 5:5:1]), &[ 0, 0, 1, 0]));
    eval!(pdf[1:1:0, 1:1:0, 1:1:0, 6:6:1] = shift(&view!(f[1:1:0, 1:1:0, 1:1:0, 6:6:1]), &[ 0, 0, -1, 0]));
    // next-nearest neighbours
    // xy plane
    eval!(pdf[1:1:0, 1:1:0, 1:1:0, 7:7:1] = shift(&view!(f[1:1:0, 1:1:0, 1:1:0, 7:7:1]), &[ 1, 1, 0, 0]));
    eval!(pdf[1:1:0, 1:1:0, 1:1:0, 8:8:1] = shift(&view!(f[1:1:0, 1:1:0, 1:1:0, 8:8:1]), &[ -1, 1, 0, 0]));
    eval!(pdf[1:1:0, 1:1:0, 1:1:0, 9:9:1] = shift(&view!(f[1:1:0, 1:1:0, 1:1:0, 9:9:1]), &[ 1, -1, 0, 0]));
    eval!(pdf[1:1:0, 1:1:0, 1:1:0, 10:10:1] = shift(&view!(f[1:1:0, 1:1:0, 1:1:0, 10:10:1]), &[ -1, -1, 0, 0]));
    // xz plane
    eval!(pdf[1:1:0, 1:1:0, 1:1:0, 11:11:1] = shift(&view!(f[1:1:0, 1:1:0, 1:1:0, 11:11:1]), &[ 1, 0, 1, 0]));
    eval!(pdf[1:1:0, 1:1:0, 1:1:0, 12:12:1] = shift(&view!(f[1:1:0, 1:1:0, 1:1:0, 12:12:1]), &[ -1, 0, 1, 0]));
    eval!(pdf[1:1:0, 1:1:0, 1:1:0, 13:13:1] = shift(&view!(f[1:1:0, 1:1:0, 1:1:0, 13:13:1]), &[ 1, 0, -1, 0]));
    eval!(pdf[1:1:0, 1:1:0, 1:1:0, 14:14:1] = shift(&view!(f[1:1:0, 1:1:0, 1:1:0, 14:14:1]), &[ -1, 0, -1, 0]));
    // yz plane
    eval!(pdf[1:1:0, 1:1:0, 1:1:0, 15:15:1] = shift(&view!(f[1:1:0, 1:1:0, 1:1:0, 15:15:1]), &[ 0, 1, 1, 0]));
    eval!(pdf[1:1:0, 1:1:0, 1:1:0, 16:16:1] = shift(&view!(f[1:1:0, 1:1:0, 1:1:0, 16:16:1]), &[ 0, -1, 1, 0]));
    eval!(pdf[1:1:0, 1:1:0, 1:1:0, 17:17:1] = shift(&view!(f[1:1:0, 1:1:0, 1:1:0, 17:17:1]), &[ 0, 1, -1, 0]));
    eval!(pdf[1:1:0, 1:1:0, 1:1:0, 18:18:1] = shift(&view!(f[1:1:0, 1:1:0, 1:1:0, 18:18:1]), &[ 0, -1, -1, 0]));
    // next next-nearest neighbours
    eval!(pdf[1:1:0, 1:1:0, 1:1:0, 19:19:1] = shift(&view!(f[1:1:0, 1:1:0, 1:1:0, 19:19:1]), &[ 1, 1, 1, 0]));
    eval!(pdf[1:1:0, 1:1:0, 1:1:0, 20:20:1] = shift(&view!(f[1:1:0, 1:1:0, 1:1:0, 20:20:1]), &[ -1, 1, 1, 0]));
    eval!(pdf[1:1:0, 1:1:0, 1:1:0, 21:21:1] = shift(&view!(f[1:1:0, 1:1:0, 1:1:0, 21:21:1]), &[ 1, -1, 1, 0]));
    eval!(pdf[1:1:0, 1:1:0, 1:1:0, 22:22:1] = shift(&view!(f[1:1:0, 1:1:0, 1:1:0, 22:22:1]), &[ -1, -1, 1, 0]));
    eval!(pdf[1:1:0, 1:1:0, 1:1:0, 23:23:1] = shift(&view!(f[1:1:0, 1:1:0, 1:1:0, 23:23:1]), &[ 1, 1, -1, 0]));
    eval!(pdf[1:1:0, 1:1:0, 1:1:0, 24:24:1] = shift(&view!(f[1:1:0, 1:1:0, 1:1:0, 24:24:1]), &[ -1, 1, -1, 0]));
    eval!(pdf[1:1:0, 1:1:0, 1:1:0, 25:25:1] = shift(&view!(f[1:1:0, 1:1:0, 1:1:0, 25:25:1]), &[ 1, -1, -1, 0]));
    eval!(pdf[1:1:0, 1:1:0, 1:1:0, 26:26:1] = shift(&view!(f[1:1:0, 1:1:0, 1:1:0, 26:26:1]), &[ -1, -1, -1, 0]));

    pdf
}

fn output_csv(mlups: Vec<f32>, nx: u64, ny: u64, nz: u64) -> Result<(), Box<dyn Error>> {
  let mut wtr = Writer::from_path(format!("benchmarks/GPU_NAME_d3q27_mrt_channel_mlups_{}_{}_{}.csv", nx, ny, nz))?;

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
    let ny: u64 = 30;
    let nz: u64 = 30;

    let total_nodes = nx * ny * nz;

    // Physical parameters
    let u_max: FloatNum = 0.1; // Inlet speed
    let u_max_af = constant::<FloatNum>(u_max, dim4!(ny));
    let rho0: FloatNum = 1.0;

    let dt: FloatNum = 1.0;
    let c_squ: FloatNum = 1. / 3.;

    let obstacle_x: u64 = nx / 5 + 1; // x location of the cylinder
    let obstacle_y: u64 = ny / 2 + ny / 15; // y location of the cylinder
    let obstacle_z: u64 = nz / 2; // z location of the cylinder
    let obstacle_r: u64 = ny / 10 + 1; // radius of the cylinder

    let x: Array<FloatNum> = tile(&range(dim4!(nx), 1), dim4!(1, ny * nz));
    let y: Array<FloatNum> = tile(&range(dim4!(1, ny), 1), dim4!(nx, nz));
    let z: Array<FloatNum> = tile(&range(dim4!(1, nz), 1), dim4!(nx, ny));

    let dims = dim4!(nx, ny, nz);
    let span = seq!();

    // Discrete velocities
    let ex = Array::<FloatNum>::new(&[0., 1.,-1., 0., 0., 0., 0., 1.,-1., 1.,-1., 1.,-1., 1.,-1., 0., 0., 0., 0., 1.,-1., 1.,-1., 1.,-1., 1.,-1.], dim4!(27));
    let ey = Array::<FloatNum>::new(&[0., 0., 0., 1.,-1., 0., 0., 1., 1.,-1.,-1., 0., 0., 0., 0., 1.,-1., 1.,-1., 1., 1.,-1.,-1., 1., 1.,-1.,-1.], dim4!(27));
    let ez = Array::<FloatNum>::new(&[0., 0., 0., 0., 0., 1.,-1., 0., 0., 0., 0., 1., 1.,-1.,-1., 1., 1.,-1.,-1., 1., 1., 1., 1.,-1.,-1.,-1.,-1.], dim4!(27));

    let ci: Array<u64> = (range::<u64>(dim4!(1, 26), 1) + 1) * total_nodes;
                                          // 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26
    let nbidx = Array::new(&[1,0,3,2,5,4,9,8,7, 6,13,12,11,10,17,16,15,14,25,24,23,22,21,20,19,18], dim4!(26));
    let nbi: Array<u64> = view!(ci[span, nbidx]);

    let main_index = moddims(&range(dim4!(total_nodes * 27), 0), dim4!(nx, ny, nz, 27));
    let nb_index = flat(&stream(&main_index));

    // Multiple relaxation parameters for MRT scheme
    let s0: FloatNum = 1.3;
    let s1: FloatNum = 1.3;
    let sv: FloatNum = 1.0;
    let sb: FloatNum = 0.55;
    let s3: FloatNum = 1.2;
    let s3b: FloatNum = 0.55;
    let s4: FloatNum = 1.2;
    let s4b: FloatNum = 1.2;
    let s5: FloatNum = 1.2;
    let s6: FloatNum = 1.2;
    let relax_params_arr: [FloatNum; 27] = [
        s0,
        s1,
        s1,
        s1,
        sv,
        sv,
        sv,
        sb,
        sv,
        sv,
        s3,
        s3,
        s3,
        s3,
        s3,
        s3,
        s3b,
        s4,
        s4,
        s4,
        s4b,
        s4b,
        s4b,
        s5,
        s5,
        s5,
        s6,
    ];
    let relax_params = Array::new(&relax_params_arr, dim4!(27));
    let s = diag_create(&relax_params, 0);
    // let s = diag_create(&constant(omega, dim4!(27)), 0); // = BGK, or
    // let s = diag_create(&constant(1.2f32, dim4!(27)), 0); // = SRT set to constant
    let s_t = transpose(&s, false);

    // Kinetic viscosity
    let nu = (1.0 as FloatNum / sv - 0.5 as FloatNum) * c_squ * dt;
    // Bulk viscosity
    let bv = (2. as FloatNum/3. as FloatNum) * (1.0 as FloatNum / sb - 0.5 as FloatNum) * c_squ * dt;

    println!("Inlet speed: {}", &u_max);
    println!("Lattice viscosity: {}", &nu);
    println!("Bulk viscosity: {}", &bv);

    // Transformation matrix
    let mut tm = constant(1.0 as FloatNum, dim4!(27, 27));
    set_row(&mut tm, &ex, 1);
    set_row(&mut tm, &ey, 2);
    set_row(&mut tm, &ez, 3);
    set_row(&mut tm, &(&ex * &ey), 4);
    set_row(&mut tm, &(&ex * &ez), 5);
    set_row(&mut tm, &(&ey * &ez), 6);
    set_row(&mut tm, &(&ex * &ex + &ey * &ey + &ez * &ez), 7);
    set_row(&mut tm, &(&ex * &ex - &ey * &ey), 8);
    set_row(&mut tm, &(&ex * &ex - &ez * &ez), 9);
    set_row(&mut tm, &(&ex * &ey * &ey), 10);
    set_row(&mut tm, &(&ex * &ez * &ez), 11);
    set_row(&mut tm, &(&ey * &ex * &ex), 12);
    set_row(&mut tm, &(&ez * &ex * &ex), 13);
    set_row(&mut tm, &(&ey * &ez * &ez), 14);
    set_row(&mut tm, &(&ez * &ey * &ey), 15);
    set_row(&mut tm, &(&ex * &ey * &ez), 16);
    set_row(&mut tm, &(&ex * &ex * &ey * &ey), 17);
    set_row(&mut tm, &(&ex * &ex * &ez * &ez), 18);
    set_row(&mut tm, &(&ey * &ey * &ez * &ez), 19);
    set_row(&mut tm, &(&ex * &ex * &ey * &ez), 20);
    set_row(&mut tm, &(&ex * &ey * &ey * &ez), 21);
    set_row(&mut tm, &(&ex * &ey * &ez * &ez), 22);
    set_row(&mut tm, &(&ex * &ey * &ey * &ez * &ez), 23);
    set_row(&mut tm, &(&ex * &ex * &ey * &ez * &ez), 24);
    set_row(&mut tm, &(&ex * &ex * &ey * &ey * &ez), 25);
    set_row(&mut tm, &(&ex * &ex * &ey * &ey * &ez * &ez), 26);

    // IM = M^-1
    let tm_inv = inverse(&tm, MatProp::NONE);

    let tm_t = transpose(&tm, false);
    let tm_inv_t = transpose(&tm_inv, false);

    // Flow around obstacle
    // circle
    let mut bound = constant::<FloatNum>(1.0, dims);
    let r = constant::<FloatNum>(obstacle_r as FloatNum, dims);
    let r_sq = &r * &r;
    let pipe = moddims(
        &ge(
            &(pow(&(flat(&y) - ny as FloatNum / 2.0 as FloatNum), &(2.0 as FloatNum), false) + pow(&(flat(&z) - nz as FloatNum / 2.0 as FloatNum), &(2.0 as FloatNum), false)),
            &flat(&constant::<FloatNum>(ny as FloatNum / 2.0 as FloatNum, dims)),
            false
        ), dims);
    let sphere = moddims(
        &le(
            &(pow(&(flat(&x) - obstacle_x as FloatNum), &(2.0 as FloatNum), false) + pow(&(flat(&y) - obstacle_y as FloatNum), &(2.0 as FloatNum), false) + pow(&(flat(&z) - obstacle_z as FloatNum), &(2.0 as FloatNum), false)),
            &flat(&r_sq),
            false
        ), dims);
    bound = selectr(&bound, &(pipe | sphere), 0.0 as f64);

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
    let mut uz = constant::<FloatNum>(0.0, dims);

    let zeroed_on = constant::<FloatNum>(0.0, on.dims());

    eval!(ux[0:0:1, 1:1:0] = u_max_af);
    eval!(ux[on] = zeroed_on);
    eval!(uy[on] = zeroed_on);
    eval!(uz[on] = zeroed_on);
    eval!(density[on] = zeroed_on);

    // Particle distribution function in initial equilibrium state
    let mut f = constant(rho0/27. as FloatNum, dim4!(total_nodes, 27));
    /*
    * MOMENTS
    */
    let mut m = matmul(&f, &tm_t, MatProp::NONE, MatProp::NONE);

    let cs_2 = c_squ;
    let cs_4 = c_squ * c_squ;
    let cs_6 = c_squ * c_squ * c_squ;
    let mut ux_1d = flat(&ux);
    let mut uy_1d = flat(&uy);
    let mut uz_1d = flat(&uz);
    let mut rho_1d = flat(&density);
    // Update equilibrium moments
    let mut rho_vx = &rho_1d * &ux_1d;
    let mut rho_vy = &rho_1d * &uy_1d;
    let mut rho_vz = &rho_1d * &uz_1d;
    let mut rho_cs_2 = cs_2 * &rho_1d;
    let mut rho_cs_4 = cs_4 * &rho_1d;
    let mut rho_cs_6 = cs_6 * &rho_1d;
    let mut v_x_sq = &ux_1d * &ux_1d;
    let mut v_y_sq = &uy_1d * &uy_1d;
    let mut v_z_sq = &uz_1d * &uz_1d;
    let mut v_sq = &v_x_sq + &v_y_sq + &v_z_sq;

    /*
    * EQUILIBRIUM MOMENTS
    */
    let mut meq = constant(0.0 as FloatNum, dim4!(total_nodes, 27));
    set_col(&mut meq, &rho_1d, 0);
    set_col(&mut meq, &rho_vx, 1);
    set_col(&mut meq, &rho_vy, 2);
    set_col(&mut meq, &rho_vz, 3);
    set_col(&mut meq, &(&rho_vx * &uy_1d), 4);
    set_col(&mut meq, &(&rho_vx * &uz_1d), 5);
    set_col(&mut meq, &(&rho_vy * &uz_1d), 6);
    set_col(&mut meq, &(&rho_1d + &rho_vx * &ux_1d + &rho_vy * &uy_1d + &rho_vz * &uz_1d), 7);
    set_col(&mut meq, &(&rho_1d * (&v_x_sq - &v_y_sq)), 8);
    set_col(&mut meq, &(&rho_1d * (&v_x_sq - &v_z_sq)), 9);
    set_col(&mut meq, &(&rho_cs_2 * &ux_1d), 10);
    set_col(&mut meq, &(&rho_cs_2 * &ux_1d), 11);
    set_col(&mut meq, &(&rho_cs_2 * &uy_1d), 12);
    set_col(&mut meq, &(&rho_cs_2 * &uz_1d), 13);
    set_col(&mut meq, &(&rho_cs_2 * &uy_1d), 14);
    set_col(&mut meq, &(&rho_cs_2 * &uz_1d), 15);
    // col 16 already initialized to zero
    set_col(&mut meq, &(&rho_cs_2 * (cs_2 + &v_x_sq + &v_y_sq)), 17);
    set_col(&mut meq, &(&rho_cs_2 * (cs_2 + &v_x_sq + &v_z_sq)), 18);
    set_col(&mut meq, &(&rho_cs_2 * (cs_2 + &v_y_sq + &v_z_sq)), 19);
    set_col(&mut meq, &(&rho_cs_2 * &uy_1d * &uz_1d), 20);
    set_col(&mut meq, &(&rho_cs_2 * &ux_1d * &uz_1d), 21);
    set_col(&mut meq, &(&rho_cs_2 * &ux_1d * &uy_1d), 22);
    set_col(&mut meq, &(&rho_cs_4 * &ux_1d), 23);
    set_col(&mut meq, &(&rho_cs_4 * &uy_1d), 24);
    set_col(&mut meq, &(&rho_cs_4 * &uz_1d), 25);
    set_col(&mut meq, &(&v_sq + &rho_cs_6), 26);

    let relaxation_part = matmul(&(&m - &meq), &s_t, MatProp::NONE, MatProp::NONE);
    let collided_moments = &m - &relaxation_part;
    f = matmul(&collided_moments, &tm_inv_t, MatProp::NONE, MatProp::NONE);

    // Create a window to show the waves.
    let mut win = Window::new(1500, 1500, "LBM solver using ArrayFire".to_string());
    win.grid(2, 2);

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

        // Boundary conditions
        eval!(ux[0:0:1, 1:1:0] = u_max_af); // inlet speed
        eval!(ux[on] = zeroed_on);
        eval!(uy[on] = zeroed_on);
        eval!(density[on] = zeroed_on);
        set_row(&mut density, &constant::<FloatNum>(1.0, dim4!(ny)), 0);
        set_row(&mut density, &constant::<FloatNum>(1.0, dim4!(ny)), nx as i64 - 1);

        // Collision
        /*
         * MOMENTS
         */
         m = matmul(&f_2d, &tm_t, MatProp::NONE, MatProp::NONE);

         ux_1d = flat(&ux);
         uy_1d = flat(&uy);
         uz_1d = flat(&uz);
         rho_1d = flat(&density);
         // Update equilibrium moments
         rho_vx = &rho_1d * &ux_1d;
         rho_vy = &rho_1d * &uy_1d;
         rho_vz = &rho_1d * &uz_1d;
         rho_cs_2 = cs_2 * &rho_1d;
         rho_cs_4 = cs_4 * &rho_1d;
         rho_cs_6 = cs_6 * &rho_1d;
         v_x_sq = &ux_1d * &ux_1d;
         v_y_sq = &uy_1d * &uy_1d;
         v_z_sq = &uz_1d * &uz_1d;
         v_sq = &v_x_sq + &v_y_sq + &v_z_sq;

         /*
         * EQUILIBRIUM MOMENTS
         */
         set_col(&mut meq, &rho_1d, 0);
         set_col(&mut meq, &rho_vx, 1);
         set_col(&mut meq, &rho_vy, 2);
         set_col(&mut meq, &rho_vz, 3);
         set_col(&mut meq, &(&rho_vx * &uy_1d), 4);
         set_col(&mut meq, &(&rho_vx * &uz_1d), 5);
         set_col(&mut meq, &(&rho_vy * &uz_1d), 6);
         set_col(&mut meq, &(&rho_1d + &rho_vx * &ux_1d + &rho_vy * &uy_1d + &rho_vz * &uz_1d), 7);
         set_col(&mut meq, &(&rho_1d * (&v_x_sq - &v_y_sq)), 8);
         set_col(&mut meq, &(&rho_1d * (&v_x_sq - &v_z_sq)), 9);
         set_col(&mut meq, &(&rho_cs_2 * &ux_1d), 10);
         set_col(&mut meq, &(&rho_cs_2 * &ux_1d), 11);
         set_col(&mut meq, &(&rho_cs_2 * &uy_1d), 12);
         set_col(&mut meq, &(&rho_cs_2 * &uz_1d), 13);
         set_col(&mut meq, &(&rho_cs_2 * &uy_1d), 14);
         set_col(&mut meq, &(&rho_cs_2 * &uz_1d), 15);
         // col 16 already initialized to zero
         set_col(&mut meq, &(&rho_cs_2 * (cs_2 + &v_x_sq + &v_y_sq)), 17);
         set_col(&mut meq, &(&rho_cs_2 * (cs_2 + &v_x_sq + &v_z_sq)), 18);
         set_col(&mut meq, &(&rho_cs_2 * (cs_2 + &v_y_sq + &v_z_sq)), 19);
         set_col(&mut meq, &(&rho_cs_2 * &uy_1d * &uz_1d), 20);
         set_col(&mut meq, &(&rho_cs_2 * &ux_1d * &uz_1d), 21);
         set_col(&mut meq, &(&rho_cs_2 * &ux_1d * &uy_1d), 22);
         set_col(&mut meq, &(&rho_cs_4 * &ux_1d), 23);
         set_col(&mut meq, &(&rho_cs_4 * &uy_1d), 24);
         set_col(&mut meq, &(&rho_cs_4 * &uz_1d), 25);
         set_col(&mut meq, &(&v_sq + &rho_cs_6), 26);

         let relaxation_part = matmul(&(&m - &meq), &s_t, MatProp::NONE, MatProp::NONE);
         let collided_moments = &m - &relaxation_part;
         f = matmul(&collided_moments, &tm_inv_t, MatProp::NONE, MatProp::NONE);

         eval!(f[reflected] = bouncedback);

        // Visualization
        if iter % 10 == 0 {
          let mut uu = moddims(&sqrt(&v_sq), dims);
          eval!(uu[on] = constant::<FloatNum>(FloatNum::NAN, on.dims()));

          let filter_x = seq!(0, nx as i32 - 1, nx as i32 / 15);
          let filter_y = seq!(0, ny as i32 - 1, ny as i32 / 30);
          let z_section = seq!(nz as i32 / 2, nz as i32 / 2, 1);

          win.set_view(0, 0);
          win.set_colormap(ColorMap::SPECTRUM);
          win.draw_image(&index(&transpose(&normalize(&uu), false), &[span, span, z_section]), Some("XY domain".to_string()));

          win.set_view(1, 0);
          win.set_axes_limits_2d(0.0, nx as f32, 0.0, ny as f32, true);
          win.draw_vector_field2(&flat(&view!(x[filter_x,filter_y])), &flat(&view!(y[filter_x,filter_y])), &flat(&view!(ux[filter_x,filter_y,z_section])), &flat(&view!(uy[filter_x,filter_y,z_section])), Some("Velocity field".to_string()));

          win.show();
        }

        let time = timer.elapsed().as_secs_f32();
        mlups.push((total_nodes as FloatNum * iter as FloatNum * 10e-6) / time);

        if iter % 100 == 0 {
            println!("{} iterations completed, {}s elapsed ({} MLUPS).", iter, time, mlups[iter as usize]);
        }

        sync(0);
        iter += 1;
    }

    mem_info!("After benchmark");

    // output CSV of MLUPS data
    if write_csv {
      output_csv(mlups, nx, ny, nz);
    }
}

fn main() {
    set_device(0);
    set_backend(Backend::OPENCL);
    info();
    println!("LBM D3Q27 simulation\n");
    let write_csv = false;
    lbm(write_csv);
}
