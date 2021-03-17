#include <arrayfire.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// WORK IN PROGRESS !!!

using namespace af;

Window *win;

array mul(const array &a, const array &b) { return a * b; }

array normalize(array a)
{
  return (a / (max<float>(abs(a)) * 1.1)) + 0.1;
}

array stream(array f)
{
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
  const unsigned ny = 100;
  const unsigned nz = 100;

  const unsigned total_nodes = nx * ny * nz;

  // Physical parameters.
  const float ux_lid = 0.1; // horizontal lid velocity
  const float uy_lid = 0;
  const float uz_lid = 0;
  const float rho0 = 1.0;
  // Reynolds number
  float Re = 150.0;
  // Kinematic viscosity
  float nu = ux_lid * nx / Re;
  // Relaxation time
  float tau = 3 * nu + 0.5;
  // Relaxation parameter
  float omega = 1.0 / tau;

  printf("Horizontal lid velocity ux_lid: %f\n", ux_lid);
  printf("Reynolds number: %f\n", Re);
  printf("Lattice viscosity: %f\n", nu);
  printf("Relaxation time: %f\n", tau);
  printf("Relaxation parameter: %f\n", omega);

  const float t1 = 8. / 27.;
  const float t2 = 2. / 27.;
  const float t3 = 1. / 54.;
  const float t4 = 1. / 216.;

  array x = tile(range(nx), 1, ny);
  array y = tile(range(dim4(1, ny), 1), nx, 1);
  array z = tile(range(dim4(1, nz), 1), nx * ny);
  seq lidx = seq(1,nx-2);
  seq lidy = seq(1,ny-2);

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

  array main_index_4d = moddims(range(dim4(total_nodes*27)),nx,ny,nz,27);
  array nb_index = flat(stream(main_index_4d));

  // Open lid
  array BOUND = constant(1, nx, ny, nz);
  BOUND(lidx,lidy,seq(1,end)) = 0; // all empty except top lid

  // matrix offset of each Occupied Node
  array ON = where(BOUND);

  // Bounceback indexes
  array TO_REFLECT = flat(tile(ON,CI.elements())) + flat(tile(CI,ON.elements()));
  array REFLECTED = flat(tile(ON,NBI.elements())) + flat(tile(NBI,ON.elements()));

  array DENSITY = constant(rho0, nx, ny, nz);
  array UX = constant(0, nx, ny, nz);
  array UY = constant(0, nx, ny, nz);
  array UZ = constant(0, nx, ny, nz);

  array w_tiled = flat(tile(transpose(w), total_nodes));
  array ex_tiled = flat(tile(transpose(ex), total_nodes));
  array ey_tiled = flat(tile(transpose(ey), total_nodes));
  array ez_tiled = flat(tile(transpose(ez), total_nodes));

  // Particle distribution function in initial equilibrium state
  array u_sq = flat(pow(UX, 2) + af::pow(UY, 2) + af::pow(UZ, 2));
  array eu = flat(batchFunc(transpose(ex), flat(UX), mul) + batchFunc(transpose(ey), flat(UY), mul) + batchFunc(transpose(ez), flat(UZ), mul));
  array F = flat(batchFunc(transpose(w), flat(DENSITY), mul)) * (1.0f + 3.0f*eu + 4.5f*(pow(eu,2)) - 1.5f*(tile(u_sq,27)));

  array uu = constant(0,nx,ny,nz);


  // Setup Window
  win = new Window(1500, 1500, "LBM solver using ArrayFire");
  win->grid(2, 2);

  unsigned iter = 0;
  unsigned maxiter = 15000;

  sync();
  timer::start();

  while (!win->close() && iter < maxiter)
  {
    // Streaming by reading from neighbors (with pre-built index) - pull scheme
    array F_streamed = F(nb_index);

    array BOUNCEDBACK = F_streamed(TO_REFLECT); // Densities bouncing back at next timestep

    array F_2D = moddims(F_streamed, total_nodes, 27);

    // Compute macroscopic variables
    array rho = sum(F_2D, 1);
    DENSITY = moddims(rho,nx,ny,nz);

    array fex = batchFunc(transpose(ex), F_2D, mul);
    array fey = batchFunc(transpose(ey), F_2D, mul);
    array fez = batchFunc(transpose(ez), F_2D, mul);

    UX = moddims((sum(fex, 1) / rho),nx,ny,nz);
    UY = moddims((sum(fey, 1) / rho),nx,ny,nz);
    UZ = moddims((sum(fez, 1) / rho),nx,ny,nz);

    // Macroscopic (Dirichlet) boundary conditions
    UX(lidx,lidy,end) = ux_lid; // lid x - velocity
    UY(lidx,lidy,end) = uy_lid; // lid y - velocity

    UX(ON) = 0;
    UY(ON) = 0;
    UZ(ON) = 0;
    DENSITY(ON) = 0;

    // Collision
    u_sq = flat(pow(UX, 2) + af::pow(UY, 2) + af::pow(UZ, 2));
    eu = flat(batchFunc(transpose(ex), flat(UX), mul) + batchFunc(transpose(ey), flat(UY), mul) + batchFunc(transpose(ez), flat(UZ), mul));
    array FEQ = flat(batchFunc(transpose(w), flat(DENSITY), mul)) * (1.0f + 3.0f*eu + 4.5f*(pow(eu,2)) - 1.5f*(tile(u_sq,27)));

    F = omega * FEQ + (1 - omega) * F_streamed;

    F(REFLECTED) = BOUNCEDBACK;


    if (iter % 10 == 0) {
      uu = moddims(sqrt(u_sq),nx,ny,nz);
      uu(ON) = af::NaN;

      seq filter =  seq(0,nx-1,(int)ceil(nx/30));

      const char *str = "Velocity field for iteration ";
      std::stringstream title;
      title << str << iter;
      (*win)(0, 0).setColorMap(AF_COLORMAP_SPECTRUM);
      (*win)(0, 0).image(flip(transpose(reorder(normalize(uu), 1, 0, 2)(span, span, (int)ceil(nz / 2))), 0));
      (*win)(0, 1).setAxesLimits(0.0f,(float)nx,0.0f,(float)ny,true);
      (*win)(0, 1).vectorField(flat(x(filter,filter)), flat(y(filter,filter)), flat(UX(filter,filter,(int)ceil(nz / 2))), flat(UY(filter,filter,(int)ceil(nz / 2))), std::move(title).str().c_str());
      (*win)(1, 0).setColorMap(AF_COLORMAP_SPECTRUM);
      (*win)(1, 0).image(normalize(uu)(span, span, (int)ceil(nz / 2)));
      (*win)(1, 1).setAxesLimits(0.0f,(float)nx,0.0f,(float)nz,true);
      (*win)(1, 1).vectorField(flat(x(filter,filter)), flat(z(filter,filter)), flat(UX(filter,filter,(int)ceil(nz / 2))), flat(UZ(filter,filter,(int)ceil(nz / 2))), std::move(title).str().c_str());
      win->show();
    }

    if (iter % 100 == 0) {
      float time = timer::stop();
      float mlups = (total_nodes * iter * 10e-6) / time;
      printf("%u iterations completed, %fs elapsed (%f MLUPS).\n", iter, time, mlups);
    }

    iter++;
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
    uu = sqrt(moddims(u_sq,nx,ny,nz));
    uu(ON) = af::NaN;

    seq filter = seq(0,nx-1,(int)ceil(nx/30));

    const char *str = "Velocity field for iteration ";
    std::stringstream title;
    title << str << iter;
    (*win)(0, 0).setColorMap(AF_COLORMAP_SPECTRUM);
    (*win)(0, 0).image(flip(transpose(normalize(uu)),0));
    (*win)(0, 1).setAxesLimits(0.0f,(float)nx,0.0f,(float)ny,true);
    (*win)(0, 1).vectorField(flat(x(filter,filter)), flat(y(filter,filter)), flat(UX(filter,filter,(int)ceil(nz / 2))), flat(UY(filter,filter,(int)ceil(nz / 2))), std::move(title).str().c_str());
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
    printf("LBM D2Q9 simulation\n");
    lbm(console);
  }
  catch (exception &e)
  {
    fprintf(stderr, "%s\n", e.what());
    throw;
  }

  return 0;
}
