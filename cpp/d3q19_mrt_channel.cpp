#include <arrayfire.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// WORK IN PROGRESS !!!

// Suga et. al (2015)
// https://www.sciencedirect.com/science/article/pii/S0898122115000346

using namespace af;

float cx[27] = {0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1};
float cy[27] = {0, 0, 0, 1,-1, 0, 0, 1, 1,-1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1, 1,-1,-1, 1, 1,-1,-1};
float cz[27] = {0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1, 1,-1,-1, 1, 1,-1,-1, 1, 1, 1, 1,-1,-1,-1,-1};

Window *win;

array normalize(array a)
{
  return (a / (max<float>(abs(a)) * 1.1)) + 0.1;
}

array stream(array f) {
    for (int i = 1; i < 27; i++) {
      f(span, span, span, i) = shift(f, cx[i], cy[i], cz[i])(span, span, span, i);
    }
    return f;
}

static void lbm()
{
  // Grid length, number and spacing
  const unsigned nx = 320;
  const unsigned ny = 80;
  const unsigned nz = 80;

  const unsigned total_nodes = nx * ny * nz;

  const float dt = 1.0;

  // Physical parameters.
  const float rho0 = 1.0;

  const int obstacle_x = nx / 4; // x location of the cylinder
  const int obstacle_y = ny / 2 + ny / 30; // y location of the cylinder
  const int obstacle_z = nz / 2 + nz / 30; // z location of the cylinder
  const int obstacle_r = ny / 9; // radius of the cylinder

  const float t1 = 8. / 27.;
  const float t2 = 2. / 27.;
  const float t3 = 1. / 54.;
  const float t4 = 1. / 216.;
  const float c_squ = 1. / 3.;

  array x = tile(range(nx), 1, ny * nz);
  array y = tile(range(dim4(1, ny), 1), nx, nz);
  array z = tile(range(dim4(1, nz), 1), nx * ny);

  // Discrete velocities
  array ex = {0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1};
  array ey = {0, 0, 0, 1,-1, 0, 0, 1, 1,-1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1, 1,-1,-1, 1, 1,-1,-1};
  array ez = {0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1, 1,-1,-1, 1, 1,-1,-1, 1, 1, 1, 1,-1,-1,-1,-1};

  array CI = (range(dim4(1, 26), 1) + 1) * total_nodes;
              // 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26
  array nbidx = {1,0,3,2,5,4,9,8,7, 6,13,12,11,10,17,16,15,14,25,24,23,22,21,20,19,20};
  array NBI = CI(span, nbidx);

  // Multiple relaxation parameters for MRT scheme
  float s0 = 1.0;
  float s1 = 1.0;
  float sv = 1.2;
  float sb = 0.6;
  float s3 = 1.2;
  float s3b = 0.6;
  float s4 = 1.2;
  float s4b = 1.2;
  float s5 = 1.2;
  float s6 = 1.2;
  array relax_params = {
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
  array S = diag(relax_params, 0, false);
  // array S = diag(constant(omega,27), 0, false); // = BGK, or
  // array S = diag(constant(0.8,27), 0, false); // = SRT set to constant

  // Lattice speed
  float u_max = 0.1;
  // Kinetic viscosity
  float nu = (1.0 / sv - 0.5) * c_squ * dt;
  // Bulk viscosity
  float bv = (2./3.) * (1.0 / sb - 0.5) * c_squ * dt;

  // printf("Reynolds number: %f\n", Re);
  printf("Lattice speed: %f\n", u_max);
  printf("Lattice viscosity: %f\n", nu);
  printf("Bulk viscosity: %f\n", bv);
  // printf("Relaxation time: %f\n", tau);
  // printf("Relaxation parameter: %f\n", omega);

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
  array BOUND = constant(0, nx, ny, nz);
  // circle
  BOUND(span,span,span) = moddims((af::pow(flat(x) - obstacle_x, 2) + af::pow(flat(y) - obstacle_y, 2) + af::pow(flat(z) - obstacle_z, 2)) <= pow(obstacle_r,2), nx, ny, nz);
  BOUND(span, span, end) = 1; // top
  BOUND(span, span, 0) = 1;   // bottom
  BOUND(span, end, span) = 1; // front
  BOUND(span, 0, span) = 1;   // back

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

  DENSITY(0,span,span) = 1;
  DENSITY(end,span,span) = 1;

  // Particle distribution function in initial equilibrium state
  array F = constant(rho0/27., 27, total_nodes);
  /*
   * MOMENTS
   */
  array moments = matmul(M, F);

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

  array u_sq = af::pow(flat(UX), 2) + af::pow(flat(UY), 2) + af::pow(flat(UZ), 2);
  array eu = (flat(tile(transpose(ex), total_nodes)) * tile(flat(UX),27)) + (flat(tile(transpose(ey), total_nodes)) * tile(flat(UY),27)) + (flat(tile(transpose(ez), total_nodes)) * tile(flat(UZ),27));
  array meq = flat(tile(transpose(w), total_nodes)) * tile(flat(DENSITY),27) * (1.0f + 3.0f*eu + 4.5f*(af::pow(eu,2)) - 1.5f*(tile(u_sq,27)));
  meq = moddims(F,nx,ny,nz,27);

  array relaxation_part = matmul(S, (moments - meq));
  array collided_moments = moments - relaxation_part;
  array F_postcol = matmul(IM, collided_moments);
  F = transpose(F_postcol);

  array uu;

  // Setup Window
  win = new Window(1500, 1500, "LBM solver using ArrayFire");
  win->grid(2, 2);

  unsigned iter = 0;
  unsigned maxiter = 15000;

  sync(0);
  timer::start();

  while (!win->close() && iter < maxiter)
  {
    F = moddims(F, nx, ny, nz, 27);
    F = stream(F);

    array BOUNCEDBACK = F(TO_REFLECT); // Densities bouncing back at next timestep

    array F_2D = moddims(F, total_nodes, 27);

    // Compute macroscopic variables
    array rho = sum(F_2D, 1);
    DENSITY = moddims(rho,nx,ny,nz);

    array fex = tile(transpose(ex), total_nodes) * F_2D;
    array fey = tile(transpose(ey), total_nodes) * F_2D;
    array fez = tile(transpose(ez), total_nodes) * F_2D;

    UX = moddims((sum(fex, 1) / rho),nx,ny,nz);
    UY = moddims((sum(fey, 1) / rho),nx,ny,nz);
    UZ = moddims((sum(fez, 1) / rho),nx,ny,nz);

    UX(0, span, span) = u_max;
    UX(ON) = 0;
    UY(ON) = 0;
    UZ(ON) = 0;
    DENSITY(ON) = 0;
    DENSITY(0,span,span) = 1;
    DENSITY(end,span,span) = DENSITY(end-1,span,span);

    F(end, span, span, span) = F(end-1, span, span, span);
    array F_t = transpose(moddims(F, total_nodes, 27));

    /*
     * MOMENTS
     */
    moments = matmul(M, F_t);

    UX_1d = flat(UX);
    UY_1d = flat(UY);
    UZ_1d = flat(UZ);
    DENSITY_1d = flat(DENSITY);
    rho_vx = DENSITY_1d * UX_1d;
    rho_vy = DENSITY_1d * UY_1d;
    rho_vz = DENSITY_1d * UZ_1d;
    rho_cs_2 = cs_2 * DENSITY_1d;
    rho_cs_4 = cs_4 * DENSITY_1d;
    rho_cs_6 = cs_6 * DENSITY_1d;
    v_x_sq = UX_1d * UX_1d;
    v_y_sq = UY_1d * UY_1d;
    v_z_sq = UZ_1d * UZ_1d;
    v_sq = v_x_sq + v_y_sq + v_z_sq;

    /*
     * EQUILIBRIUM MOMENTS
     */
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

     /*
      * COLLISION STEP
      */
    // Partials
    relaxation_part = matmul(S, (moments - meq));
    // printf("relaxation_part dims = [%lld %lld %lld]\n", relaxation_part.dims(0), relaxation_part.dims(1), relaxation_part.dims(2));
    collided_moments = moments - relaxation_part;
    // printf("collided_moments dims = [%lld %lld %lld]\n", collided_moments.dims(0), collided_moments.dims(1), collided_moments.dims(2));

    F_postcol = matmul(IM, collided_moments);
    // println!("Post-collision F from moments (dims): {:?}", &updated_f.dims());
    // printf("Post-collision F from moments (dims) = [%lld %lld %lld]\n", F_postcol.dims(0), F_postcol.dims(1), F_postcol.dims(2));
    // Shuffle into 1D array representation of the domain
    F = transpose(F_postcol);

    F(REFLECTED) = BOUNCEDBACK;

    if (iter % 10 == 0) {
      uu = sqrt(moddims(v_sq,nx,ny,nz));
      uu(ON) = af::NaN;

      seq filterX = seq(0,nx-1,(int)ceil(nx/15));
      seq filterY = seq(0,ny-1,(int)ceil(ny/30));
      seq filterZ = seq(0,nz-1,(int)ceil(nz/30));

      std::stringstream titleUXY;
      std::stringstream titleUXZ;
      titleUXY << "Velocity field XY, iteration " << iter;
      titleUXZ << "Velocity field XZ, iteration " << iter;
      (*win)(0, 0).setColorMap(AF_COLORMAP_SPECTRUM);
      (*win)(0, 0).image(transpose(normalize(uu))(span, span, (int)ceil(nz / 2)));
      (*win)(0, 1).vectorField(flat(x(filterX,filterY)), flat(y(filterX,filterY)), flat(UX(filterX,filterY,(int)ceil(nz / 2))), flat(UY(filterX,filterY,(int)ceil(nz / 2))), std::move(titleUXY).str().c_str());
      // (*win)(1, 0).image(transpose(reorder(normalize(uu)(span, ny / 2, span), 0, 2, 1)));
      // (*win)(1, 1).vectorField(flat(x(filterX,filterZ)), flat(z), flat(UX), flat(UZ), std::move(titleUXZ).str().c_str());
      win->show();
    }

    if (iter % 100 == 0) {
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
  try
  {
    af::setDevice(device);
    af::info();
    printf("LBM D3Q27 MRT simulation\n");
    lbm();
  }
  catch (af::exception &e)
  {
    fprintf(stderr, "%s\n", e.what());
    throw;
  }

  return 0;
}
