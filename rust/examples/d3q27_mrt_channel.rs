use arrayfire::*;
use std::time::Instant;

// WORK IN PROGRESS

type FloatNum = f32;

fn normalize(a: &Array<FloatNum>) -> Array<FloatNum> {
  (a / (max_all(&abs(a)).0 * 1.1 as FloatNum)) + 0.1 as FloatNum
}

fn lbm() {
    // Grid length, number and spacing
    let nx: u64 = 160;
    let ny: u64 = 80;
    let nz: u64 = 80;

    let total_nodes = nx * ny * nz;

    // Physical parameters.
    let rho0: FloatNum = 1.0;

    let obstacle_x: u64 = nx / 5 + 1; // x location of the cylinder
    let obstacle_y: u64 = ny / 2 + ny / 15; // y location of the cylinder
    let obstacle_z: u64 = nz / 2; // z location of the cylinder
    let obstacle_r: u64 = ny / 10 + 1; // radius of the cylinder

    // Reynolds number
    let re: FloatNum = 150.0;
    // Lattice speed
    let u_max: FloatNum = 0.05;
    let u_max_af = constant::<FloatNum>(u_max, dim4!(ny));
    // Kinematic viscosity
    let nu: FloatNum = u_max * 2.0 as FloatNum * obstacle_r as FloatNum / re;
    // Relaxation time
    let tau: FloatNum = (3.0 as FloatNum) * nu + (0.5 as FloatNum);
    // Relaxation parameter
    let omega: FloatNum = (1.0 as FloatNum) / tau;

    println!("Reynolds number: {}", re);
    println!("Lattice speed: {}", u_max);
    println!("Lattice viscosity: {}", nu);
    println!("Relaxation time: {}", tau);
    println!("Relaxation parameter: {}", omega);

    let t1: FloatNum = 8. / 27.;
    let t2: FloatNum = 2. / 27.;
    let t3: FloatNum = 1. / 54.;
    let t4: FloatNum = 1. / 216.;
    let c_squ: FloatNum = 1. / 3.;

    let x = tile(&range(dim4!(nx), 0), dim4!(1, ny * nz));
    let y = tile(&range(dim4!(1, ny), 1), dim4!(nx, nz));
    let z = tile(&range(dim4!(1, nz), 1), dim4!(nx * nz));

    let dims = dim4!(nx, ny, nz);

    // Discrete velocities
    let ex = Array::<FloatNum>::new(&[0., 1.,-1., 0., 0., 0., 0., 1.,-1., 1.,-1., 1.,-1., 1.,-1., 0., 0., 0., 0., 1.,-1., 1.,-1., 1.,-1., 1.,-1.], dim4!(27));
    let ey = Array::<FloatNum>::new(&[0., 0., 0., 1.,-1., 0., 0., 1., 1.,-1.,-1., 0., 0., 0., 0., 1.,-1., 1.,-1., 1., 1.,-1.,-1., 1., 1.,-1.,-1.], dim4!(27));
    let ez = Array::<FloatNum>::new(&[0., 0., 0., 0., 0., 1.,-1., 0., 0., 0., 0., 1., 1.,-1.,-1., 1., 1.,-1.,-1., 1., 1., 1., 1.,-1.,-1.,-1.,-1.], dim4!(27));

    // weights
    let w = Array::new(&[t1,t2,t2,t2,t2,t2,t2,t3,t3,t3,t3,t3,t3,t3,t3,t3,t3,t3,t3,t4,t4,t4,t4,t4,t4,t4,t4], dim4!(27));

    let ci: Array<u64> = (range::<u64>(dim4!(1, 26),1) + 1) * total_nodes;
                          // 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26
    let nbidx = Array::new(&[1,0,3,2,5,4,9,8,7, 6,13,12,11,10,17,16,15,14,25,24,23,22,21,20,19,18], dim4!(8));
    let span = seq!();
    let nbi: Array<u64> = view!(ci[span, nbidx]);

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
    let relax_params = Array::<FloatNum>::new(&[
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
    ], dim4!(27));
    let S = diag_create(&relax_params, 0);
    // let S = diag_create(constant(&omega,27), 0); // = BGK, or
    // let S = diag_create(constant(1.8 as FloatNum, dim4!(27)), 0); // = SRT set to constant

    // Transformation matrix
    let mut tm = constant::<FloatNum>(0.0, dim4!(27,27));
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
    let to_reflect = flat(&tile(&on,dim4!(ci.elements() as u64))) + flat(&tile(&ci,dim4!(on.elements() as u64)));
    let reflected = flat(&tile(&on,dim4!(nbi.elements() as u64))) + flat(&tile(&nbi,dim4!(on.elements() as u64)));

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

    // Start in equilibrium state
    let mut f = constant::<FloatNum>(0.0, dim4!(total_nodes, 27));

    /*
     * MOMENTS
     */
    let moments = matmul(&tm, &f, MatProp::NONE, MatProp::NONE);

    let ux_1d = flat(&ux);
    let uy_1d = flat(&uy);
    let uz_1d = flat(&uz);
    let rho_1d = flat(&density);
    // Update equilibrium moments
    let rho_vx = &rho_1d * &ux_1d;
    let rho_vy = &rho_1d * &uy_1d;
    let rho_vz = &rho_1d * &uz_1d;

    let cs_2 = constant::<FloatNum>(c_squ, dim4!(total_nodes));
    let cs_4 = &cs_2 * &cs_2;
    let cs_6 = &cs_2 * &cs_2 * &cs_2;

    let rho_cs_2 = &rho_1d * &cs_2;
    let rho_cs_4 = &rho_1d * &cs_4;
    let rho_cs_6 = &rho_1d * &cs_6;

    let v_x_sq = &ux_1d * &ux_1d;
    let v_y_sq = &uy_1d * &uy_1d;
    let v_z_sq = &uz_1d * &uz_1d;

    let v_sq = &v_x_sq + &v_y_sq + &v_z_sq;
    /*
     * EQUILIBRIUM MOMENTS
     */
    let mut meq = constant::<FloatNum>(
        0.0 as FloatNum,
        dim4!(27, total_nodes),
    );

    set_row(&mut meq, &rho_1d, 0);
    set_row(&mut meq, &rho_vx, 1);
    set_row(&mut meq, &rho_vy, 2);
    set_row(&mut meq, &rho_vz, 3);
    set_row(&mut meq, &(&rho_vx * &uy_1d), 4);
    set_row(&mut meq, &(&rho_vx * &uz_1d), 5);
    set_row(&mut meq, &(&rho_vy * &uz_1d), 6);
    set_row(
        &mut meq,
        &(&rho_1d + &(&rho_vx * &ux_1d) + &(&rho_vy * &uy_1d) + &(&rho_vz * &uz_1d)),
        7,
    );
    set_row(&mut meq, &(&rho_1d * &(&v_x_sq - &v_y_sq)), 8);
    set_row(&mut meq, &(&rho_1d * &(&v_x_sq - &v_z_sq)), 9);
    set_row(&mut meq, &(&rho_cs_2 * &ux_1d), 10);
    set_row(&mut meq, &(&rho_cs_2 * &ux_1d), 11);
    set_row(&mut meq, &(&rho_cs_2 * &uy_1d), 12);
    set_row(&mut meq, &(&rho_cs_2 * &uz_1d), 13);
    set_row(&mut meq, &(&rho_cs_2 * &uy_1d), 14);
    set_row(&mut meq, &(&rho_cs_2 * &uz_1d), 15);
    // meq(16,span) Already initialized, zeroed
    set_row(&mut meq, &(&rho_cs_2 * &(&cs_2 + &v_x_sq + &v_y_sq)), 17);
    set_row(&mut meq, &(&rho_cs_2 * &(&cs_2 + &v_x_sq + &v_z_sq)), 18);
    set_row(&mut meq, &(&rho_cs_2 * &(&cs_2 + &v_y_sq + &v_z_sq)), 19);
    set_row(&mut meq, &(&rho_cs_2 * &uy_1d * &uz_1d), 20);
    set_row(&mut meq, &(&rho_cs_2 * &ux_1d * &uz_1d), 21);
    set_row(&mut meq, &(&rho_cs_2 * &ux_1d * &uy_1d), 22);
    set_row(&mut meq, &(&rho_cs_4 * &ux_1d), 23);
    set_row(&mut meq, &(&rho_cs_4 * &uy_1d), 24);
    set_row(&mut meq, &(&rho_cs_4 * &uz_1d), 25);
    set_row(&mut meq, &(&v_sq + &rho_cs_6), 26);

    let mut relaxation_part = matmul(&S, &(&moments - &meq), MatProp::NONE, MatProp::NONE);
    let mut collided_moments = moments - relaxation_part;
    let mut F_postcol = matmul(&tm_inv, &collided_moments, MatProp::NONE, MatProp::NONE);
    f = transpose(&F_postcol, false);

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
        f = moddims(&f, dim4!(nx, ny, nz, 27));
        // Streaming
        let pre_stream = &view!(f[span, span, span]);
        // nearest-neighbours
        assign_seq(&mut f, &[span, span, seq!(1:1:1)], &shift(&index(pre_stream, &[span, span, seq!(1:1:1)]), &[1, 0, 0, 0]));
        assign_seq(&mut f, &[span, span, seq!(2:2:1)], &shift(&index(pre_stream, &[span, span, seq!(2:2:1)]), &[-1, 0, 0, 0]));
        assign_seq(&mut f, &[span, span, seq!(3:3:1)], &shift(&index(pre_stream, &[span, span, seq!(3:3:1)]), &[0, 1, 0, 0]));
        assign_seq(&mut f, &[span, span, seq!(4:4:1)], &shift(&index(pre_stream, &[span, span, seq!(4:4:1)]), &[0, -1, 0, 0]));
        assign_seq(&mut f, &[span, span, seq!(5:5:1)], &shift(&index(pre_stream, &[span, span, seq!(5:5:1)]), &[0, 0, 1, 0]));
        assign_seq(&mut f, &[span, span, seq!(6:6:1)], &shift(&index(pre_stream, &[span, span, seq!(6:6:1)]), &[0, 0, -1, 0]));
        // next-nearest neighbours
        // xy plane
        assign_seq(&mut f, &[span, span, seq!(7:7:1)], &shift(&index(pre_stream, &[span, span, seq!(7:7:1)]), &[1, 1, 0, 0]));
        assign_seq(&mut f, &[span, span, seq!(8:8:1)], &shift(&index(pre_stream, &[span, span, seq!(8:8:1)]), &[-1, 1, 0, 0]));
        assign_seq(&mut f, &[span, span, seq!(9:9:1)], &shift(&index(pre_stream, &[span, span, seq!(9:9:1)]), &[1, -1, 0, 0]));
        assign_seq(&mut f, &[span, span, seq!(10:10:1)], &shift(&index(pre_stream, &[span, span, seq!(10:10:1)]), &[-1, -1, 0, 0]));
        // xz plane
        assign_seq(&mut f, &[span, span, seq!(11:11:1)], &shift(&index(pre_stream, &[span, span, seq!(11:11:1)]), &[1, 0, 1, 0]));
        assign_seq(&mut f, &[span, span, seq!(12:12:1)], &shift(&index(pre_stream, &[span, span, seq!(12:12:1)]), &[-1, 0, 1, 0]));
        assign_seq(&mut f, &[span, span, seq!(13:13:1)], &shift(&index(pre_stream, &[span, span, seq!(13:13:1)]), &[1, 0, -1, 0]));
        assign_seq(&mut f, &[span, span, seq!(14:14:1)], &shift(&index(pre_stream, &[span, span, seq!(14:14:1)]), &[-1, 0, -1, 0]));
        // yz plane
        assign_seq(&mut f, &[span, span, seq!(15:15:1)], &shift(&index(pre_stream, &[span, span, seq!(15:15:1)]), &[0, 1, 1, 0]));
        assign_seq(&mut f, &[span, span, seq!(16:16:1)], &shift(&index(pre_stream, &[span, span, seq!(16:16:1)]), &[0, -1, 1, 0]));
        assign_seq(&mut f, &[span, span, seq!(17:17:1)], &shift(&index(pre_stream, &[span, span, seq!(17:17:1)]), &[0, 1, -1, 0]));
        assign_seq(&mut f, &[span, span, seq!(18:18:1)], &shift(&index(pre_stream, &[span, span, seq!(18:18:1)]), &[0, -1, -1, 0]));
        // next next-nearest neighbours
        assign_seq(&mut f, &[span, span, seq!(19:19:1)], &shift(&index(pre_stream, &[span, span, seq!(19:19:1)]), &[1, 1, 1, 0]));
        assign_seq(&mut f, &[span, span, seq!(20:20:1)], &shift(&index(pre_stream, &[span, span, seq!(20:20:1)]), &[-1, 1, 1, 0]));
        assign_seq(&mut f, &[span, span, seq!(21:21:1)], &shift(&index(pre_stream, &[span, span, seq!(21:21:1)]), &[1, -1, 1, 0]));
        assign_seq(&mut f, &[span, span, seq!(22:22:1)], &shift(&index(pre_stream, &[span, span, seq!(22:22:1)]), &[-1, -1, 1, 0]));
        assign_seq(&mut f, &[span, span, seq!(23:23:1)], &shift(&index(pre_stream, &[span, span, seq!(23:23:1)]), &[1, 1, -1, 0]));
        assign_seq(&mut f, &[span, span, seq!(24:24:1)], &shift(&index(pre_stream, &[span, span, seq!(24:24:1)]), &[-1, 1, -1, 0]));
        assign_seq(&mut f, &[span, span, seq!(25:25:1)], &shift(&index(pre_stream, &[span, span, seq!(25:25:1)]), &[1, -1, -1, 0]));
        assign_seq(&mut f, &[span, span, seq!(26:26:1)], &shift(&index(pre_stream, &[span, span, seq!(26:26:1)]), &[-1, -1, -1, 0]));

        let bouncedback = view!(f[to_reflect]); // Densities bouncing back at next timestep

        // Compute macroscopic variables
        density = sum(&f, 2);
        let rho = flat(&density);

        let f_2d = moddims(&f, dim4!(total_nodes, 27));

        let fex = tile(&transpose(&ex, false), dim4!(total_nodes)) * &f_2d;
        let fey = tile(&transpose(&ey, false), dim4!(total_nodes)) * &f_2d;
        let fez = tile(&transpose(&ez, false), dim4!(total_nodes)) * &f_2d;

        ux = moddims(&(sum(&fex, 1) / &rho), dim4!(nx,ny,nz));
        uy = moddims(&(sum(&fey, 1) / &rho), dim4!(nx,ny,nz));
        uz = moddims(&(sum(&fez, 1) / &rho), dim4!(nx,ny,nz));

        // inlet speed
        eval!(ux[0:0:1, 1:1:0] = u_max_af);

        eval!(ux[on] = zeroed_on);
        eval!(uy[on] = zeroed_on);
        eval!(density[on] = zeroed_on);
        set_row(&mut density, &constant::<FloatNum>(1.0, dim4!(ny)), 0);
        set_row(&mut density, &constant::<FloatNum>(1.0, dim4!(ny)), nx as i64 - 1);

        // Collision
        u_sq = pow(&ux, &(2.0 as FloatNum), false) + pow(&uy, &(2.0 as FloatNum), false);
        eu = (flat(&tile(&transpose(&ex, false), dim4!(total_nodes))) * tile(&flat(&ux),dim4!(9))) + (flat(&tile(&transpose(&ey, false), dim4!(total_nodes))) * &tile(&flat(&uy), dim4!(9)));
        let mut feq = flat(&tile(&transpose(&w, false), dim4!(total_nodes))) * &tile(&flat(&density), dim4!(9)) * ((1.0 as FloatNum) + (3.0 as FloatNum) * &eu + (4.5 as FloatNum) * (&pow(&eu, &(2.0 as FloatNum), false)) - (1.5 as FloatNum) * (&tile(&flat(&u_sq), dim4!(9))));
        feq = moddims(&feq,dim4!(nx,ny,9));

        f = omega * feq + (1.0 - omega) * f;

        eval!(f[reflected] = bouncedback);

        // Visualization
        if iter % 10 == 0 {
            let mut uu = sqrt(&u_sq);
            eval!(uu[on] = constant::<FloatNum>(FloatNum::NAN, on.dims()));

            let filter_x = seq!(0, nx as i32 - 1, nx as i32 / 15);
            let filter_y = seq!(0, ny as i32 - 1, ny as i32 / 30);

            win.set_view(0, 0);
            win.set_colormap(ColorMap::SPECTRUM);
            win.draw_image(&transpose(&normalize(&uu), false), Some("XY domain".to_string()));

            win.set_view(1, 0);
            win.set_axes_limits_2d(0.0, nx as f32, 0.0, ny as f32, true);
            win.draw_vector_field2(&flat(&view!(x[filter_x,filter_y])), &flat(&view!(y[filter_x,filter_y])), &flat(&view!(ux[filter_x,filter_y])), &flat(&view!(uy[filter_x,filter_y])), Some("Velocity field".to_string()));

            win.show();
        }

        let time = timer.elapsed().as_secs() as FloatNum;
        mlups.push((total_nodes as FloatNum * iter as FloatNum * 10e-6) / time);

        if iter % 100 == 0 {
            println!("{} iterations completed, {}s elapsed ({} MLUPS).", iter, time, mlups[iter as usize]);
        }

        sync(0);
        iter += 1;
    }

    mem_info!("After benchmark");
}

fn main() {
    set_device(0);
    set_backend(Backend::OPENCL);
    info();
    println!("LBM D2Q9 simulation\n");
    lbm();
}
