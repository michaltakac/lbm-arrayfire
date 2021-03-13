#include <arrayfire.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

using namespace af;

Window *win;

array normalize(array a)
{
  return (a / (max<float>(abs(a)) * 1.1)) + 0.1;
}

array stream(array f) {
    // nearest-neighbours
    f(span, span, span, 1) = shift(f, 1, 0, 0)(span, span, span, 1);
    f(span, span, span, 2) = shift(f, -1, 0, 0)(span, span, span, 2);
    f(span, span, span, 3) = shift(f, 0, 1, 0)(span, span, span, 3);
    f(span, span, span, 4) = shift(f, 0, -1, 0)(span, span, span, 4);
    f(span, span, span, 5) = shift(f, 0, 0, 1)(span, span, span, 5);
    f(span, span, span, 6) = shift(f, 0, 0, -1)(span, span, span, 6);
    // next-nearest neighbours
    // xy plane
    f(span, span, span, 7) = shift(f, 1, 1, 0)(span, span, span, 7);
    f(span, span, span, 8) = shift(f, -1, 1, 0)(span, span, span, 8);
    f(span, span, span, 9) = shift(f, 1, -1, 0)(span, span, span, 9);
    f(span, span, span, 10) = shift(f, -1, -1, 0)(span, span, span, 10);
    // xz plane
    f(span, span, span, 11) = shift(f, 1, 0, 1)(span, span, span, 11);
    f(span, span, span, 12) = shift(f, -1, 0, 1)(span, span, span, 12);
    f(span, span, span, 13) = shift(f, 1, 0, -1)(span, span, span, 13);
    f(span, span, span, 14) = shift(f, -1, 0, -1)(span, span, span, 14);
    // yz plane
    f(span, span, span, 15) = shift(f, 0, 1, 1)(span, span, span, 15);
    f(span, span, span, 16) = shift(f, 0, -1, 1)(span, span, span, 16);
    f(span, span, span, 17) = shift(f, 0, 1, -1)(span, span, span, 17);
    f(span, span, span, 18) = shift(f, 0, -1, -1)(span, span, span, 18);
    // next next-nearest neighbours
    f(span, span, span, 19) = shift(f, 1, 1, 1)(span, span, span, 19);
    f(span, span, span, 20) = shift(f, -1, 1, 1)(span, span, span, 20);
    f(span, span, span, 21) = shift(f, 1, -1, 1)(span, span, span, 21);
    f(span, span, span, 22) = shift(f, -1, -1, 1)(span, span, span, 22);
    f(span, span, span, 23) = shift(f, 1, 1, -1)(span, span, span, 23);
    f(span, span, span, 24) = shift(f, -1, 1, -1)(span, span, span, 24);
    f(span, span, span, 25) = shift(f, 1, -1, -1)(span, span, span, 25);
    f(span, span, span, 26) = shift(f, -1, -1, -1)(span, span, span, 26);
    return f;
}

static void lbm(bool console)
{
  // Grid length, number and spacing
  const unsigned nx = 128;
  const unsigned ny = 128;
  const unsigned nz = 128;

  const unsigned total_nodes = nx * ny * nz;

  // Physical parameters.
  const float rho0 = 1.0;

  const int obstacle_x = nx / 4; // x location of the cylinder
  const int obstacle_y = ny / 2; // y location of the cylinder
  const int obstacle_z = nz / 2; // z location of the cylinder
  const int obstacle_r = ny / 9; // radius of the cylinder

  // Reynolds number
  float Re = 150.0;
  // Lattice speed
  float u_max = 0.01;
  // Kinematic viscosity
  float nu = u_max * 2 * obstacle_r / Re; // dt / dh_sq / Re;
  // Relaxation time
  float tau = 3 * nu + 0.5;
  // Relaxation parameter
  float omega = 1.0 / tau;

  printf("Reynolds number: %f\n", Re);
  printf("Lattice speed: %f\n", u_max);
  printf("Lattice viscosity: %f\n", nu);
  printf("Relaxation time: %f\n", tau);
  printf("Relaxation parameter: %f\n", omega);

  const float t1 = 8. / 27.;
  const float t2 = 2. / 27.;
  const float t3 = 1. / 54.;
  const float t4 = 1. / 216.;
  const float c_squ = 1. / 3.;

  array x = tile(range(nx), 1, ny * nz);
  array y = tile(range(dim4(1, ny), 1), nx, nz);
  array z = tile(range(dim4(1, nz), 1), nx * ny);

  array coords = join(1, flat(x), flat(y), flat(z));
  // print("coords",coords);
  array codes = flat(x) + flat(y) * nx + flat(z) * nx * ny;
  // print("codes",codes);

  // Discrete velocities
  float cx[27] = {0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1};
  float cy[27] = {0, 0, 0, 1,-1, 0, 0, 1, 1,-1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1, 1,-1,-1, 1, 1,-1,-1};
  float cz[27] = {0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1, 1,-1,-1, 1, 1,-1,-1, 1, 1, 1, 1,-1,-1,-1,-1};

  array ex(27, cx);
  array ey(27, cy);
  array ez(27, cz);

  // weights
  float weights[27] = {t1,t2,t2,t2,t2,t2,t2,t3,t3,t3,t3,t3,t3,t3,t3,t3,t3,t3,t3,t4,t4,t4,t4,t4,t4,t4,t4};
  array w(27, weights);

  array CI = (range(dim4(1, 26), 1) + 1) * total_nodes;
  int nbindex[] = {1, 2, 3, 4, 5, 6, 8, 9, 6, 7, 12, 13, 10, 11, 16, 17, 14, 15, 22, 23, 24, 25, 18, 19, 20, 21};
  array nbidx(26, nbindex);
  array NBI = CI(span, nbidx);

  // Multiple relaxation parameters for MRT scheme
  float s0 = 1.0;
  float s1 = 1.0;
  float sv = 1.1;
  float sb = 0.8;
  float s3 = 1.1;
  float s3b = 0.8;
  float s4 = 1.1;
  float s4b = 1.1;
  float s5 = 2.1;
  float s6 = 2.1;
  float relax_params[27] = {
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
  };
  // float relax_params[27] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
  array relax_params_af(27, relax_params);
  array S = diag(relax_params_af, 0, false);
  // array S = diag(constant(omega,27), 0, false); // = SRT, or

  // Particle distribution function
  // array F = constant(rho0/27, nx,ny,nz,27);

  // Transformation matrix
  array M = constant(0, 27, 27);
  M(0, span) = 1;
  M(1, span) = ex;
  M(2, span) = ey;
  M(3, span) = ez;
  M(5, span) = ex * ey;
  M(5, span) = ex * ez;
  M(6, span) = ey * ez;
  M(7, span) = ex * ex + ey * ey + ez * ez;
  M(8, span) = ex * ex - ey * ey;
  M(9, span) = ex * ex - ez * ez;
  M(10, span) = ex * ey * ey;
  M(11, span) = ex * ez * ez;
  M(12, span) = ey * ex * ex;
  M(13, span) = ez * ex * ex;
  M(14, span) = ey * ez * ez;
  M(15, span) = ez * ey * ey;
  M(16, span) = ex * ey * ez;
  M(17, span) = ex * ex * ey * ey;
  M(18, span) = ex * ex * ez * ez;
  M(19, span) = ey * ey * ez * ez;
  M(20, span) = ex * ex * ey * ez;
  M(21, span) = ex * ey * ey * ez;
  M(22, span) = ex * ey * ez * ez;
  M(23, span) = ex * ey * ey * ez * ez;
  M(24, span) = ex * ex * ey * ez * ez;
  M(25, span) = ex * ex * ey * ey * ez;
  M(26, span) = ex * ex * ey * ey * ez * ez;

  // IM = M^-1
  array IM = pinverse(M);

  // Flow around obstacle
  // array BOUND = constant(0, nx, ny, nz);
  // // circle
  // BOUND(span,span,span) = moddims((af::pow(flat(x) - obstacle_x, 2) + af::pow(flat(y) - obstacle_y, 2) + af::pow(flat(z) - obstacle_z, 2)) <= pow(obstacle_r,2), nx, ny, nz);

  // BOUND(span, span, end) = 1; // top
  // BOUND(span, span, 0) = 1;   // bottom
  // BOUND(span, end, span) = 1; // front
  // BOUND(span, 0, span) = 1;   // back

  // Porous media
  array BOUND = randu(nx, ny, nz) < 0.8;

  // matrix offset of each Occupied Node
  array ON = where(BOUND);

  // Bounceback indexes
  array TO_REFLECT = flat(tile(ON, CI.elements())) + flat(tile(CI, ON.elements()));
  array REFLECTED = flat(tile(ON, NBI.elements())) + flat(tile(NBI, ON.elements()));

  array DENSITY = constant(rho0, nx, ny, nz);
  array UX = constant(u_max, nx, ny, nz);
  array UY = constant(0, nx, ny, nz);
  array UZ = constant(0, nx, ny, nz);

  UX(0,span,span) = u_max;
  UX(ON) = 0;
  UY(ON) = 0;
  UZ(ON) = 0;
  DENSITY(ON) = 0;

  // Particle distribution function in initial equilibrium state
  array u_sq = af::pow(UX, 2) + af::pow(UY, 2) + af::pow(UZ, 2);
  array eu = (flat(tile(transpose(ex), total_nodes)) * tile(flat(UX),27)) + (flat(tile(transpose(ey), total_nodes)) * tile(flat(UY),27)) + (flat(tile(transpose(ez), total_nodes)) * tile(flat(UZ),27));
  array F = flat(tile(transpose(w), total_nodes)) * tile(flat(DENSITY),27) * (1.0f + 3.0f*eu + 4.5f*(af::pow(eu,2)) - 1.5f*(tile(flat(u_sq),27)));

  array uu;

  if (!console)
  {
    win = new Window(1536, 768, "LBM solver using ArrayFire");
    win->grid(2, 2);
  }

  unsigned iter = 0;
  unsigned maxiter = 1000;

  sync(0);
  timer::start();

  // while (iter < maxiter)
  while (!win->close())
  {
    F = moddims(F, nx, ny, nz, 27);
    F = stream(F);

    array BOUNCEDBACK = F(TO_REFLECT); // Densities bouncing back at next timestep

    array F_2D = moddims(F, total_nodes, 27);
    array F_t = transpose(F_2D);

    // Compute macroscopic variables
    array rho = sum(F_2D, 1);
    DENSITY = moddims(rho,nx,ny,nz);

    array fex = tile(transpose(ex), total_nodes) * F_2D;
    array fey = tile(transpose(ey), total_nodes) * F_2D;
    array fez = tile(transpose(ez), total_nodes) * F_2D;

    UX = moddims((sum(fex, 1) / rho),nx,ny,nz);
    UY = moddims((sum(fey, 1) / rho),nx,ny,nz);
    UZ = moddims((sum(fez, 1) / rho),nx,ny,nz);

    // printf("UX dims = [%lld %lld %lld]\n", UX.dims(0), UX.dims(1), UX.dims(2));

    // int uxleft[] = {1, 7, 9, 11, 13, 19, 21, 23, 25};
    // int uxright[] = {2, 8, 10, 12, 14, 20, 22, 24, 26};
    // int uyleft[] = {3, 7, 8, 15, 17, 19, 20, 23, 24};
    // int uyright[] = {4, 9, 10, 16, 18, 21, 22, 25, 26};
    // int uzleft[] = {5, 11, 12, 15, 16, 19, 20, 21, 22};
    // int uzright[] = {6, 13, 14, 17, 18, 23, 24, 25, 26};
    // array idxleft(9, uxleft);
    // array idxright(9, uxright);
    // array idyleft(9, uyleft);
    // array idyright(9, uyright);
    // array idzleft(9, uyleft);
    // array idzright(9, uyright);
    // UX = (sum(F(span, span, span, idxleft), 3) - sum(F(span, span, span, idxright), 3)) / DENSITY;
    // // printf("UX dims = [%lld %lld %lld]\n", UX.dims(0), UX.dims(1), UX.dims(2));
    // UY = (sum(F(span, span, span, idyleft), 3) - sum(F(span, span, span, idyright), 3)) / DENSITY;
    // // printf("UY dims = [%lld %lld %lld]\n", UY.dims(0), UY.dims(1), UY.dims(2));
    // UZ = (sum(F(span, span, span, idzleft), 3) - sum(F(span, span, span, idzright), 3)) / DENSITY;
    // // printf("UZ dims = [%lld %lld %lld]\n", UZ.dims(0), UZ.dims(1), UZ.dims(2));
    // // UX(0, seq(1,end-1)) = UX(0, seq(1,end-1)) + u_lb;
    UX(0, span, span) = u_max;
    UX(ON) = 0;
    UY(ON) = 0;
    UZ(ON) = 0;
    DENSITY(ON) = 0;
    // DENSITY(0,span,span) = 1;
    // DENSITY(end,span,span) = 1;

    /*
     * MOMENTS
     */
    array moments = matmul(M, F_t);
    // printf("moments dims = [%lld %lld %lld]\n", moments.dims(0), moments.dims(1), moments.dims(2));

    array UX_1d = flat(UX);
    array UY_1d = flat(UY);
    array UZ_1d = flat(UZ);
    array DENSITY_1d = flat(DENSITY);
    // Update equilibrium moments
    array rho_vx = DENSITY_1d * UX_1d;
    array rho_vy = DENSITY_1d * UY_1d;
    array rho_vz = DENSITY_1d * UZ_1d;
    float cs_2 = c_squ;
    float cs_4 = c_squ * c_squ;
    float cs_6 = c_squ * c_squ * c_squ;
    array rho_cs_2 = cs_2 * DENSITY_1d;
    array rho_cs_4 = cs_4 * DENSITY_1d;
    array rho_cs_6 = cs_6 * DENSITY_1d;
    array v_x_sq = UX_1d * UX_1d;
    array v_y_sq = UY_1d * UY_1d;
    array v_z_sq = UZ_1d * UZ_1d;
    array v_sq = v_x_sq + v_y_sq + v_z_sq;


    /*
      * EQUILIBRIUM MOMENTS
      */
    array meq = constant(0, 27, total_nodes);

    meq(0, span) = DENSITY_1d;
    meq(1, span) = rho_vx;
    meq(2, span) = rho_vy;
    meq(3, span) = rho_vz;
    meq(4, span) = rho_vx * UY_1d;
    meq(5, span) = rho_vx * UZ_1d;
    meq(6, span) = rho_vy * UZ_1d;
    meq(7, span) = DENSITY_1d + rho_vx * UX_1d + rho_vy * UY_1d + rho_vz * UZ_1d;
    meq(8, span) = DENSITY_1d * (v_x_sq - v_y_sq);
    meq(9, span) = DENSITY_1d * (v_x_sq - v_z_sq);
    meq(10, span) = rho_cs_2 * UX_1d;
    meq(11, span) = rho_cs_2 * UX_1d;
    meq(12, span) = rho_cs_2 * UY_1d;
    meq(13, span) = rho_cs_2 * UZ_1d;
    meq(14, span) = rho_cs_2 * UY_1d;
    meq(15, span) = rho_cs_2 * UZ_1d;
    // meq(16,span) = 0; // Already initialized, zeroed
    meq(17, span) = rho_cs_2 * (cs_2 + v_x_sq + v_y_sq);
    meq(18, span) = rho_cs_2 * (cs_2 + v_x_sq + v_z_sq);
    meq(19, span) = rho_cs_2 * (cs_2 + v_y_sq + v_z_sq);
    meq(20, span) = rho_cs_2 * UY_1d * UZ_1d;
    meq(21, span) = rho_cs_2 * UX_1d * UZ_1d;
    meq(22, span) = rho_cs_2 * UX_1d * UY_1d;
    meq(23, span) = rho_cs_4 * UX_1d;
    meq(24, span) = rho_cs_4 * UY_1d;
    meq(25, span) = rho_cs_4 * UZ_1d;
    meq(26, span) = v_sq + rho_cs_6;
    // printf("meq dims = [%lld %lld %lld]\n", meq.dims(0), meq.dims(1), meq.dims(2));


    // array u_sq = af::pow(UX, 2) + af::pow(UY, 2) + af::pow(UZ, 2);
    // // printf("u_sq dims = [%lld %lld %lld]\n", u_sq.dims(0), u_sq.dims(1), u_sq.dims(2));
    // array eu = (flat(tile(transpose(ex), total_nodes)) * tile(flat(UX),27)) + (flat(tile(transpose(ey), total_nodes)) * tile(flat(UY),27)) + (flat(tile(transpose(ez), total_nodes)) * tile(flat(UZ),27));
    // // printf("eu dims = [%lld %lld %lld]\n", eu.dims(0), eu.dims(1), eu.dims(2));
    // array FEQ = flat(tile(transpose(w), total_nodes)) * tile(flat(DENSITY),27) * (1.0f + 3.0f*eu + 4.5f*(af::pow(eu,2)) - 1.5f*(tile(flat(u_sq),27)));
    // // printf("F dims = [%lld %lld %lld]\n", F.dims(0), F.dims(1), F.dims(2));

    // FEQ = transpose(moddims(FEQ,total_nodes,27));
    // array meq = matmul(M, FEQ);
    // array moments = matmul(M, F_t);

     /*
      * COLLISION STEP
      */
    // Partials
    array relaxation_part = matmul(S, (moments - meq));
    // printf("relaxation_part dims = [%lld %lld %lld]\n", relaxation_part.dims(0), relaxation_part.dims(1), relaxation_part.dims(2));
    array collided_moments = moments - relaxation_part;
    // printf("collided_moments dims = [%lld %lld %lld]\n", collided_moments.dims(0), collided_moments.dims(1), collided_moments.dims(2));

    array F_postcol = matmul(IM, collided_moments);
    // println!("Post-collision F from moments (dims): {:?}", &updated_f.dims());
    // printf("Post-collision F from moments (dims) = [%lld %lld %lld]\n", F_postcol.dims(0), F_postcol.dims(1), F_postcol.dims(2));
    // Shuffle into 1D array representation of the domain
    array F_postcol_2D = transpose(F_postcol);

    F = moddims(F, total_nodes, 27);
    F = F_postcol_2D;

    // FEQ = moddims(FEQ,nx,ny,nz,27);

    // F = omega * FEQ + (1 - omega) * F;

    F(REFLECTED) = BOUNCEDBACK;

    if (iter % 10 == 0) {
      uu = sqrt(moddims(u_sq,nx,ny,nz));
      uu(ON) = af::NaN;

      std::stringstream titleUXY;
      std::stringstream titleUXZ;
      titleUXY << "Velocity field XY, iteration " << iter;
      titleUXZ << "Velocity field XZ, iteration " << iter;
      (*win)(0, 0).setColorMap(AF_COLORMAP_SPECTRUM);
      (*win)(0, 0).image(transpose(normalize(uu))(span, span, nz / 2));
      (*win)(0, 1).vectorField(flat(x), flat(y), flat(UX), flat(UY), std::move(titleUXY).str().c_str());
      (*win)(1, 0).image(transpose(reorder(normalize(uu)(span, ny / 2, span), 0, 2, 1)));
      (*win)(1, 1).vectorField(flat(x), flat(z), flat(UX), flat(UZ), std::move(titleUXZ).str().c_str());
      win->show();
    }

    if (iter % 10 == 0) {
      float time = timer::stop();
      float mlups = (total_nodes * iter * 10e-6) / time;
      printf("%u iterations completed, %fs elapsed (%f MLUPS).\n", iter, time, mlups);
    }

    iter++;
  }

  sync(0);

  float end = timer::stop();
  float mlups = (total_nodes * iter * 10e-6) / end;
  printf("Iterations: %d\n", iter);
  printf("Time: %fs\n", end);
  printf("MLUPS: %f\n", mlups);
  af::info();

  // while (!win->close())
  // {
  //   std::stringstream titleUXY;
  //   std::stringstream titleUXZ;
  //   titleUXY << "Velocity field XY, iteration " << iter;
  //   titleUXZ << "Velocity field XZ, iteration " << iter;
  //   (*win)(0, 0).setColorMap(AF_COLORMAP_SPECTRUM);
  //   (*win)(0, 0).image(transpose(normalize(uu)(span, span, nz / 2)));
  //   (*win)(0, 1).vectorField(flat(x), flat(y), flat(UX), flat(UY), std::move(titleUXY).str().c_str());
  //   (*win)(1, 0).image(transpose(reorder(normalize(uu)(span, ny / 2, span), 0, 2, 1)));
  //   (*win)(1, 1).vectorField(flat(x), flat(z), flat(UX), flat(UZ), std::move(titleUXZ).str().c_str());
  //   win->show();
  // }
}

int main(int argc, char *argv[])
{
  int device = argc > 1 ? atoi(argv[1]) : 0;
  bool console = argc > 2 ? argv[2][0] == '-' : false;
  try
  {
    af::setDevice(device);
    af::info();
    printf("LBM D3Q27 simulation\n");
    lbm(console);
  }
  catch (af::exception &e)
  {
    fprintf(stderr, "%s\n", e.what());
    throw;
  }

  return 0;
}
