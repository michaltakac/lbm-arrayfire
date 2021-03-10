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
  const unsigned nx = 100;
  const unsigned ny = 60;
  const unsigned nz = 60;

  const unsigned total_nodes = nx * ny * nz;

  // Physical parameters.
  const float rho0 = 1.0;

  const int obstacle_x = nx / 4 + 1;       // x location of the cylinder
  const int obstacle_y = ny / 2 + ny / 30; // y location of the cylinder
  const int obstacle_z = nz / 2 + nz / 30; // z location of the cylinder
  const int obstacle_r = ny / 10 + 1;      // radius of the cylinder

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

  array F = constant(rho0, nx, ny, nz, 27);
  array FEQ = F.copy();

  array CI = (range(dim4(1, 26), 1) + 1) * total_nodes;
  int nbindex[] = {1, 2, 3, 4, 5, 6, 8, 9, 6, 7, 12, 13, 10, 11, 16, 17, 14, 15, 22, 23, 24, 25, 18, 19, 20, 21};
  array nbidx(26, nbindex);
  array NBI = CI(span, nbidx);

  // Flow around obstacle
  array BOUND = constant(0, nx, ny, nz);
  // circle
  BOUND(span,span,span) = moddims((af::pow(flat(x) - obstacle_x, 2) + af::pow(flat(y) - obstacle_y, 2) + af::pow(flat(z) - obstacle_z, 2)) <= pow(obstacle_r,2), nx, ny, nz);
  BOUND(span, span, end) = 1; // top
  BOUND(span, span, 0) = 1;   // bottom
  BOUND(span, end, span) = 1; // front
  BOUND(span, 0, span) = 1;   // back

  // Porous media
  // array BOUND = randu(nx, ny, nz) < 0.8;

  // matrix offset of each Occupied Node
  array ON = where(BOUND);

  // Bounceback indexes
  array TO_REFLECT = flat(tile(ON, CI.elements())) + flat(tile(CI, ON.elements()));
  array REFLECTED = flat(tile(ON, NBI.elements())) + flat(tile(NBI, ON.elements()));

  array DENSITY = constant(rho0, nx, ny, nz);
  array UX = constant(u_max, nx, ny, nz);
  array UY = constant(0, nx, ny, nz);
  array UZ = constant(0, nx, ny, nz);

  UX(ON) = 0;
  UY(ON) = 0;
  UZ(ON) = 0;
  DENSITY(ON) = 0;

  // Particle distribution function in initial equilibrium state
  array u_sq = af::pow(flat(UX), 2) + af::pow(flat(UY), 2) + af::pow(flat(UZ), 2);
  array eu = (flat(tile(transpose(ex), total_nodes)) * tile(flat(UX),27)) + (flat(tile(transpose(ey), total_nodes)) * tile(flat(UY),27)) + (flat(tile(transpose(ez), total_nodes)) * tile(flat(UZ),27));
  F = flat(tile(transpose(w), total_nodes)) * tile(flat(DENSITY),27) * (1.0f + 3.0f*eu + 4.5f*(af::pow(eu,2)) - 1.5f*(tile(u_sq,27)));
  F = moddims(F,nx,ny,nz,27);

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

    UX(0, span, span) = u_max;
    UX(ON) = 0;
    UY(ON) = 0;
    UZ(ON) = 0;
    DENSITY(ON) = 0;
    DENSITY(0,span,span) = 1;
    DENSITY(end,span,span) = 1;

    u_sq = af::pow(flat(UX), 2) + af::pow(flat(UY), 2) + af::pow(flat(UZ), 2);
    eu = (flat(tile(transpose(ex), total_nodes)) * tile(flat(UX),27)) + (flat(tile(transpose(ey), total_nodes)) * tile(flat(UY),27)) + (flat(tile(transpose(ez), total_nodes)) * tile(flat(UZ),27));
    FEQ = flat(tile(transpose(w), total_nodes)) * tile(flat(DENSITY),27) * (1.0f + 3.0f*eu + 4.5f*(af::pow(eu,2)) - 1.5f*(tile(u_sq,27)));
    FEQ = moddims(FEQ,nx,ny,nz,27);

    F = omega * FEQ + (1 - omega) * F;

    F(REFLECTED) = BOUNCEDBACK;

    if (iter % 10 == 0) {
      uu = sqrt(moddims(u_sq,nx,ny,nz));
      // printf("uu dims = [%lld %lld %lld]\n", uu.dims(0), uu.dims(1), uu.dims(2));
      uu(ON) = af::NaN;
      seq filterX = seq(0,nx-1,nx/15);
      seq filterY = seq(0,ny-1,ny/30);
      seq filterZ = seq(0,nz-1,nz/30);

      std::stringstream titleUXY;
      std::stringstream titleUXZ;
      titleUXY << "Velocity field XY, iteration " << iter;
      titleUXZ << "Velocity field XZ, iteration " << iter;
      (*win)(0, 0).setColorMap(AF_COLORMAP_SPECTRUM);
      (*win)(0, 0).image(reorder(normalize(uu), 2, 1, 0)(span, span, (int)ceil(nz / 2)));
      // (*win)(0, 1).vectorField(flat(x(filterX,filterY,(int)ceil(nz / 2))), flat(y(filterX,filterY,(int)ceil(nz / 2))), flat(UX(filterX,filterY,(int)ceil(nz / 2))), flat(UY(filterX,filterY,(int)ceil(nz / 2))), std::move(titleUXY).str().c_str());
      (*win)(1, 0).image(transpose(reorder(normalize(uu)(span, (int)ceil(ny / 2), span), 0, 2, 1)));
      // (*win)(1, 1).vectorField(flat(x(filterX,filterY,filterZ)), flat(z(filterX,filterY,filterZ)), flat(UX(filterX,filterY,filterZ)), flat(UZ(filterX,filterY,filterZ)), std::move(titleUXZ).str().c_str());
      win->show();
    }
    iter++;

    if (iter % 100 == 0) {
      float time = timer::stop();
      float mlups = (total_nodes * iter * 10e-6) / time;
      printf("%u iterations completed, %fs elapsed (%f MLUPS).\n", iter, time, mlups);
    }
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
