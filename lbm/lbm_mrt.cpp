#include <arrayfire.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

using namespace af;

Window *win;

array normalize(array a, float max)
{
  float mx = max * 1.1;
  float mn = -max * 1.1;
  return (a - mn) / (mx - mn);
}

static void lbmD3Q27(bool console)
{
  // Grid length, number and spacing
  const unsigned nx = 50;
  const unsigned ny = 50;
  const unsigned nz = 50;

  const unsigned total_nodes = nx * ny * nz;
  printf("total nodes: %i\n", total_nodes);

  // Physical parameters.
  const float L_p = 2.5;     //1.1; // Cavity dimension.
  const float U_p = 5.0;     //1.1; // Cavity lid velocity.
  const float nu_p = 1.2e-3; // 1.586e-5; // Physical kinematic viscosity.
  const float rho0 = 1.0;
  // Discrete/numerical parameters.
  const float dt = 0.2; //0.002;

  const int obstacle_x = nx / 4; // x location of the cylinder
  const int obstacle_y = ny / 2; // y location of the cylinder
  const int obstacle_z = nz / 2; // z location of the cylinder
  const int obstacle_r = ny / 9; // radius of the cylinder

  // Derived nondimensional parameters.
  float Re = 220.0; //L_p * U_p / nu_p;
  // Derived physical parameters.
  float t_p = L_p / U_p;
  // Derived discrete parameters.
  float dh = 1.0 / ((float)nx - 1.0);
  float dh_sq = dh * dh;
  // Lattice speed
  float u_lb = 1e-7; //dt / dh;
  // Lattice viscosity
  float nu_lb = u_lb * obstacle_r / Re; // dt / dh_sq / Re;
  // Relaxation time
  float tau = 3 * nu_lb + 0.5;
  float omega = 1.0 / tau; // 1.0;

  printf("Reynolds number: %f\n", Re);
  printf("Physical time scale: %fs\n", t_p);
  printf("dh: %f\n", dh);
  printf("Lattice speed: %f\n", u_lb);
  printf("Lattice viscosity: %f\n", nu_lb);
  printf("Relaxation time: %f\n", tau);
  printf("Relaxation parameter: %f\n", omega);

  const float t1 = 4. / 9.;
  const float t2 = 1. / 9.;
  const float t3 = 1. / 36.;
  const float c_squ = 1. / 3.;

  array x = tile(range(nx), 1, ny * nz);
  array y = tile(range(dim4(1, ny), 1), nx, nz);
  array z = tile(range(dim4(1, nz), 1), nx * ny);

  array coords = join(1, flat(x), flat(y), flat(z));
  // print("coords",coords);
  array codes = flat(x) + flat(y) * nx + flat(z) * nx * ny;
  // print("codes",codes);

  array CI = (range(dim4(1, 26), 1) + 1) * total_nodes;
  int nbindex[] = {1, 2, 3, 4, 5, 6, 8, 9, 6, 7, 12, 13, 10, 11, 16, 17, 14, 15, 22, 23, 24, 25, 18, 19, 20, 21};
  array nbidx(26, nbindex);
  array NBI = CI(span, nbidx);

  // Discrete velocities
  float cx[27] = {0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1};
  float cy[27] = {0, 0, 0, 1,-1, 0, 0, 1, 1,-1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1, 1,-1,-1, 1, 1,-1,-1};
  float cz[27] = {0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1, 1,-1,-1, 1, 1,-1,-1, 1, 1, 1, 1,-1,-1,-1,-1};

  array ex(27, cx);
  array ey(27, cy);
  array ez(27, cz);

  // Multiple relaxation parameters for MRT scheme
  float s0 = 1.0;
  float s1 = 1.0;
  float sv = 1.1;
  float sb = 0.8;
  float s3 = 1.1;
  float s3b = 0.8;
  float s4 = 1.1;
  float s4b = 1.1;
  float s5 = 1.1;
  float s6 = 1.1;
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
  array F = constant(rho0 / 27, total_nodes, 27);

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

  // print("CI: ", CI);

  // Flow around obstacle
  // array BOUND = constant(0, nx, ny, nz);
  // // circle
  // BOUND(span,span,span) = moddims((af::pow(flat(x) - obstacle_x, 2) + af::pow(flat(y) - obstacle_y, 2) + af::pow(flat(z) - obstacle_z, 2)) <= pow(obstacle_r,2), nx, ny, nz);
  // BOUND(span, span, end) = 1; // top
  // BOUND(span, span, 0) = 1;   // bottom
  // BOUND(span, end, span) = 1; // front
  // BOUND(span, 0, span) = 1;   // back

  array BOUND = randu(nx, ny, nz) < 0.8;

  // print("bound: ", BOUND);

  // matrix offset of each Occupied Node
  array ON = where(BOUND);

  // Bounceback indexes
  array TO_REFLECT = flat(tile(ON, CI.elements())) + flat(tile(CI, ON.elements()));
  array REFLECTED = flat(tile(ON, NBI.elements())) + flat(tile(NBI, ON.elements()));


  array DENSITY = constant(rho0, nx, ny, nz);
  array UX = constant(rho0, nx, ny, nz);
  array UY = constant(rho0, nx, ny, nz);
  array UZ = constant(rho0, nx, ny, nz);

  UX(ON) = 0;
  UY(ON) = 0;
  UZ(ON) = 0;
  DENSITY(ON) = 0;

  float avu = 1;
  float prevavu = 1;
  float numactivenodes = sum<float>(count(BOUND));
  array uu;

  if (!console)
  {
    win = new Window(1536, 768, "LBM solver using ArrayFire");
    win->grid(2, 2);
  }

  unsigned iter = 0;
  unsigned maxiter = 100;
  float m_f = max<float>(F);

  sync();
  timer::start();

  // while (iter < maxiter) 
  while (!win->close())
  {
    /*
    * STREAMING
    */
    // Propagate
    F = moddims(F, nx, ny, nz, 27);
    // x+ => [nx 1:nx-1] | y+ => [ny 1:ny-1] | z+ => [nz 1:nz-1]
    // x- => [2:nx 1]    | y- => [2:ny 1]    | z- => [2:nz 1]
    // nearest-neighbours
    F(span, span, span, 1) = shift(F, 1, 0, 0)(span, span, span, 1);
    F(span, span, span, 2) = shift(F, -1, 0, 0)(span, span, span, 2);
    F(span, span, span, 3) = shift(F, 0, 1, 0)(span, span, span, 3);
    F(span, span, span, 4) = shift(F, 0, -1, 0)(span, span, span, 4);
    F(span, span, span, 5) = shift(F, 0, 0, 1)(span, span, span, 5);
    F(span, span, span, 6) = shift(F, 0, 0, -1)(span, span, span, 6);
    // next-nearest neighbours
    // xy plane
    F(span, span, span, 7) = shift(F, 1, 1, 0)(span, span, span, 7);
    F(span, span, span, 8) = shift(F, -1, 1, 0)(span, span, span, 8);
    F(span, span, span, 9) = shift(F, 1, -1, 0)(span, span, span, 9);
    F(span, span, span, 10) = shift(F, -1, -1, 0)(span, span, span, 10);
    // xz plane
    F(span, span, span, 11) = shift(F, 1, 0, 1)(span, span, span, 11);
    F(span, span, span, 12) = shift(F, -1, 0, 1)(span, span, span, 12);
    F(span, span, span, 13) = shift(F, 1, 0, -1)(span, span, span, 13);
    F(span, span, span, 14) = shift(F, -1, 0, -1)(span, span, span, 14);
    // yz plane
    F(span, span, span, 15) = shift(F, 0, 1, 1)(span, span, span, 15);
    F(span, span, span, 16) = shift(F, 0, -1, 1)(span, span, span, 16);
    F(span, span, span, 17) = shift(F, 0, 1, -1)(span, span, span, 17);
    F(span, span, span, 18) = shift(F, 0, -1, -1)(span, span, span, 18);
    // next next-nearest neighbours
    F(span, span, span, 19) = shift(F, 1, 1, 1)(span, span, span, 19);
    F(span, span, span, 20) = shift(F, -1, 1, 1)(span, span, span, 20);
    F(span, span, span, 21) = shift(F, 1, -1, 1)(span, span, span, 21);
    F(span, span, span, 22) = shift(F, -1, -1, 1)(span, span, span, 22);
    F(span, span, span, 23) = shift(F, 1, 1, -1)(span, span, span, 23);
    F(span, span, span, 24) = shift(F, -1, 1, -1)(span, span, span, 24);
    F(span, span, span, 25) = shift(F, 1, -1, -1)(span, span, span, 25);
    F(span, span, span, 26) = shift(F, -1, -1, -1)(span, span, span, 26);

    array BOUNCEDBACK = F(TO_REFLECT); // Densities bouncing back at next timestep

    array F_2D = moddims(F, total_nodes, 27);
    array F_t = transpose(F_2D);

    // Compute macroscopic variables
    DENSITY = sum(F, 3);
    // printf("DENSITY dims = [%lld %lld %lld]\n", DENSITY.dims(0), DENSITY.dims(1), DENSITY.dims(2));
    array fex = moddims(tile(transpose(ex), total_nodes) * F_2D,nx,ny,nz,27);
    array fey = moddims(tile(transpose(ey), total_nodes) * F_2D,nx,ny,nz,27);
    array fez = moddims(tile(transpose(ez), total_nodes) * F_2D,nx,ny,nz,27);
    // printf("fex dims = [%lld %lld %lld]\n", fex.dims(0), fex.dims(1), fex.dims(2));

    UX = (sum(fex, 3) / DENSITY);
    UY = (sum(fey, 3) / DENSITY);
    UZ = (sum(fez, 3) / DENSITY);

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
    UX(0, span, span) = u_lb;
    UX(ON) = 0;
    UY(ON) = 0;
    UZ(ON) = 0;
    DENSITY(ON) = 0;

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

    F(REFLECTED) = BOUNCEDBACK;

    prevavu = avu;
    avu = sum<float>(sum(sum(UX))) / numactivenodes;
    uu = sqrt(moddims(v_sq, nx,ny,nz));// / avu;

    if (!console)
    {
      std::stringstream titleUXY;
      std::stringstream titleUXZ;
      titleUXY << "Velocity field XY, iteration " << iter;
      titleUXZ << "Velocity field XZ, iteration " << iter;
      (*win)(0, 0).setColorMap(AF_COLORMAP_SPECTRUM);
      (*win)(0, 0).image(transpose(uu(span, span, nz / 2)));
      (*win)(0, 1).vectorField(flat(x), flat(y), flat(UX), flat(UY), std::move(titleUXY).str().c_str());
      (*win)(1, 0).image(transpose(reorder(uu(span, ny / 2, span), 0, 2, 1)));
      (*win)(1, 1).vectorField(flat(x), flat(z), flat(UX), flat(UZ), std::move(titleUXZ).str().c_str());
      win->show();
    }
    else
    {
      // eval(uu);
      // eval(F);
      // eval(DENSITY);
    }
    iter++;

    if (iter % 100 == 0) {
      float time = timer::stop();
      float mlups = (total_nodes * iter * 10e-6) / time;
      printf("%u iterations completed, %fs elapsed (%f MLUPS).\n", iter, time, mlups);
    }
  }

  sync();

  float end = timer::stop();
  float mlups = (total_nodes * iter * 10e-6) / end;
  printf("Iterations: %d\n", iter);
  printf("Time: %fs\n", end);
  printf("MLUPS: %f\n", mlups);
  af::info();

  while (!win->close())
  {
    std::stringstream titleUXY;
    std::stringstream titleUXZ;
    titleUXY << "Velocity field XY, iteration " << iter;
    titleUXZ << "Velocity field XZ, iteration " << iter;
    (*win)(0, 0).setColorMap(AF_COLORMAP_SPECTRUM);
    (*win)(0, 0).image(transpose(uu(span, span, nz / 2)));
    (*win)(0, 1).vectorField(flat(x), flat(y), flat(UX), flat(UY), std::move(titleUXY).str().c_str());
    (*win)(1, 0).image(transpose(reorder(uu(span, ny / 2, span), 0, 2, 1)));
    (*win)(1, 1).vectorField(flat(x), flat(z), flat(UX), flat(UZ), std::move(titleUXZ).str().c_str());
    win->show();
  }
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
    lbmD3Q27(console);
  }
  catch (af::exception &e)
  {
    fprintf(stderr, "%s\n", e.what());
    throw;
  }

  return 0;
}
