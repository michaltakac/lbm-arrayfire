#include <arrayfire.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

using namespace af;

Window *win;

array mul(const array &a, const array &b) { return a * b; }

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

static void lbm()
{
  // Grid length, number and spacing
  const unsigned nx = 200;
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
  float Re = 220.0;
  // Lattice speed
  float u_max = 0.1;
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

  const float t1 = 4. / 9.;
  const float t2 = 1. / 9.;
  const float t3 = 1. / 36.;
  const float t4 = 1. / 144.;
  const float c_squ = 1. / 3.;

  array x = tile(range(nx), 1, ny * nz);
  array y = tile(range(dim4(1, ny), 1), nx, nz);
  array z = tile(range(dim4(1, nz), 1), nx * ny);

  // Discrete velocities
  int cx[27] = {0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1};
  int cy[27] = {0, 0, 0, 1,-1, 0, 0, 1, 1,-1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1, 1,-1,-1, 1, 1,-1,-1};
  int cz[27] = {0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1, 1,-1,-1, 1, 1,-1,-1, 1, 1, 1, 1,-1,-1,-1,-1};
  array ex(27, cx);
  array ey(27, cy);
  array ez(27, cz);

  // weights
  float weights[27] = {t1,t2,t2,t2,t2,t2,t2,t3,t3,t3,t3,t3,t3,t3,t3,t3,t3,t3,t3,t4,t4,t4,t4,t4,t4,t4,t4};
  array w(27, weights);

  array CI = (range(dim4(1, 26), 1) + 1) * total_nodes;
                         // 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26
  unsigned int nb_index_arr[26] = {1,0,3,2,5,4,9,8,7, 6,13,12,11,10,17,16,15,14,25,24,23,22,21,20,19,18};
  array nbidx(26, nb_index_arr);
  array NBI = CI(span, nbidx);

  array main_index = moddims(range(dim4(total_nodes*27)),nx,ny,nz,27);
  array nb_index = flat(stream(main_index));

  // Flow around obstacle
  array BOUND = constant(0, nx, ny, nz);
   // pipe
  BOUND(span,span,span) = moddims((pow(flat(y), 2) + pow(flat(z), 2)) >= pow(ny,2)-1, nx, ny, nz);
  // spherical obstacle
  BOUND(span,span,span) = moddims((pow(flat(x) - obstacle_x, 2) + pow(flat(y) - obstacle_y, 2) + pow(flat(z) - obstacle_z, 2)) <= pow(obstacle_r,2), nx, ny, nz);

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
  array u_sq = flat(pow(UX, 2) + pow(UY, 2) + pow(UZ, 2));
  array eu = flat(batchFunc(transpose(ex), flat(UX), mul) + batchFunc(transpose(ey), flat(UY), mul) + batchFunc(transpose(ez), flat(UZ), mul));
  array F = flat(batchFunc(transpose(w), flat(DENSITY), mul)) * (1.0f + 3.0f*eu + 4.5f*(pow(eu,2)) - 1.5f*(tile(u_sq,27)));

  array uu;

  // Setup Window
  win = new Window(1536, 768, "LBM solver using ArrayFire");
  win->grid(2, 2);

  unsigned iter = 0;
  unsigned maxiter = 10000;

  sync();
  timer::start();

  while (!win->close() && iter < maxiter)
  {
    array F_streamed = F(nb_index);

    array BOUNCEDBACK = F_streamed(TO_REFLECT); // Densities bouncing back at next timestep

    array F_2D = moddims(F_streamed, total_nodes, 27);
    array F_flat = flat(F_2D);

    // Compute macroscopic variables
    array rho = sum(F_2D, 1);
    DENSITY = moddims(rho,nx,ny,nz);

    array fex = batchFunc(transpose(ex), F_2D, mul);
    array fey = batchFunc(transpose(ey), F_2D, mul);
    array fez = batchFunc(transpose(ez), F_2D, mul);

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

    // Collision
    u_sq = flat(pow(UX, 2) + pow(UY, 2) + pow(UZ, 2));
    eu = flat(batchFunc(transpose(ex), flat(UX), mul) + batchFunc(transpose(ey), flat(UY), mul) + batchFunc(transpose(ez), flat(UZ), mul));
    array FEQ = flat(batchFunc(transpose(w), flat(DENSITY), mul)) * (1.0f + 3.0f*eu + 4.5f*(pow(eu,2)) - 1.5f*(tile(u_sq,27)));

    F = omega * FEQ + (1 - omega) * F_flat;

    F(REFLECTED) = BOUNCEDBACK;

    if (iter % 10 == 0) {
      uu = moddims(sqrt(u_sq),nx,ny,nz);
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
      (*win)(1, 0).image(reorder(normalize(uu)((int)ceil(nx / 4), span, span), 1, 2, 0));
      // (*win)(1, 1).vectorField(flat(y(filterY,filterZ)), flat(z(filterY,filterZ)), flat(reorder(UY((int)ceil(nx / 4),filterY,filterZ), 1, 2, 0)), flat(reorder(UZ((int)ceil(nx / 4),filterY,filterZ), 1, 2, 0)), std::move(titleUXZ).str().c_str());
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
  try
  {
    af::setDevice(device);
    af::info();
    printf("LBM D3Q27 simulation\n");
    lbm();
  }
  catch (af::exception &e)
  {
    fprintf(stderr, "%s\n", e.what());
    throw;
  }

  return 0;
}
