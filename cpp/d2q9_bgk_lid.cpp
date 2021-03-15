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
  f(span, span, 1) = shift(f, 1, 0)(span, span, 1);
  f(span, span, 2) = shift(f, 0, 1)(span, span, 2);
  f(span, span, 3) = shift(f,-1, 0)(span, span, 3);
  f(span, span, 4) = shift(f, 0,-1)(span, span, 4);
  f(span, span, 5) = shift(f, 1, 1)(span, span, 5);
  f(span, span, 6) = shift(f,-1, 1)(span, span, 6);
  f(span, span, 7) = shift(f,-1,-1)(span, span, 7);
  f(span, span, 8) = shift(f, 1,-1)(span, span, 8);
  return f;
}

static void lbm()
{
  // Grid length, number and spacing
  const unsigned nx = 512;
  const unsigned ny = 512;

  const unsigned total_nodes = nx * ny;

  // Physical parameters.
  const float ux_lid  = 0.05; // horizontal lid velocity
  const float uy_lid = 0;     // vertical lid velocity
  const float rho0 = 1.0;
  // Reynolds number
  float Re = 200.0;
  // Kinematic viscosity
  float nu = ux_lid * nx / Re;
  // Relaxation time
  float tau = 3 * nu + 0.5;
  // Relaxation parameter
  float omega = 1.0 / tau;

  printf("Horizontal lid velocity ux_lid: %f\n", ux_lid);
  printf("Vertical lid velocity uy_lid: %f\n", uy_lid);
  printf("Reynolds number: %f\n", Re);
  printf("Lattice viscosity: %f\n", nu);
  printf("Relaxation time: %f\n", tau);
  printf("Relaxation parameter: %f\n", omega);

  const float t1 = 4. / 9.;
  const float t2 = 1. / 9.;
  const float t3 = 1. / 36.;

  array x = tile(range(nx), 1, ny);
  array y = tile(range(dim4(1, ny), 1), nx, 1);
  seq lid = seq(1,nx-2);

  //  c6  c2   c5
  //    \  |  /
  //  c3 -c0 - c1
  //    /  |  \
  //  c7  c4   c8
  // Discrete velocities
  array ex = {0, 1, 0,-1, 0, 1,-1,-1, 1};
  array ey = {0, 0, 1, 0,-1, 1, 1,-1,-1};

  // weights
  array w = {t1,t2,t2,t2,t2,t3,t3,t3,t3};

  array CI = (range(dim4(1,8),1)+1) * total_nodes;
  array nbidx = {2,3,0,1,6,7,4,5};
  array NBI = CI(span,nbidx);

  array main_index = moddims(range(dim4(total_nodes*9)),nx,ny,9);
  array nb_index = flat(stream(main_index));

  // Open lid
  array BOUND = constant(1, nx, ny);
  BOUND(lid,seq(1,end)) = 0; // all empty except top lid

  // matrix offset of each Occupied Node
  array ON = where(BOUND);

  // Bounceback indexes
  array TO_REFLECT = flat(tile(ON,CI.elements())) + flat(tile(CI,ON.elements()));
  array REFLECTED = flat(tile(ON,NBI.elements())) + flat(tile(NBI,ON.elements()));

  array DENSITY = constant(rho0, nx, ny);
  array UX = constant(0, nx, ny);
  array UY = constant(0, nx, ny);

  // Start in equilibrium state
  array u_sq = flat(pow(UX, 2) + pow(UY, 2));
  array eu = flat(batchFunc(transpose(ex), flat(UX), mul) + batchFunc(transpose(ey), flat(UY), mul));
  array F = flat(batchFunc(transpose(w), flat(DENSITY), mul)) * (1.0f + 3.0f*eu + 4.5f*(af::pow(eu,2)) - 1.5f*(tile(u_sq,9)));

  array uu = constant(0,nx,ny);

  // Setup Window
  win = new Window(1536, 768, "LBM solver using ArrayFire");
  win->grid(1, 2);

  unsigned iter = 0;
  unsigned maxiter = 15000;

  sync();
  timer::start();

  while (!win->close() && iter < maxiter)
  {
    array F_streamed = F(nb_index);

    array BOUNCEDBACK = F_streamed(TO_REFLECT); // Densities bouncing back at next timestep

    array F_2D = moddims(F_streamed, total_nodes, 9);
    array F_flat = flat(F_2D);

    // Compute macroscopic variables
    array rho = sum(F_2D, 1);
    DENSITY = moddims(rho,nx,ny);

    array fex = batchFunc(transpose(ex), F_2D, mul);
    array fey = batchFunc(transpose(ey), F_2D, mul);

    UX = moddims((sum(fex, 1) / rho),nx,ny);
    UY = moddims((sum(fey, 1) / rho),nx,ny);

    // MACROSCOPIC (DIRICHLET) BOUNDARY CONDITIONS
    UX(lid,end) = ux_lid; // lid x - velocity
    UY(lid,end) = uy_lid; // lid y - velocity

    UX(ON) = 0;
    UY(ON) = 0;
    DENSITY(ON) = 0;

    // Collision
    u_sq = flat(pow(UX, 2) + pow(UY, 2));
    eu = flat(batchFunc(transpose(ex), flat(UX), mul) + batchFunc(transpose(ey), flat(UY), mul));
    array FEQ = flat(batchFunc(transpose(w), flat(DENSITY), mul)) * (1.0f + 3.0f*eu + 4.5f*(af::pow(eu,2)) - 1.5f*(tile(u_sq,9)));

    F = omega * FEQ + (1 - omega) * F_flat;

    F(REFLECTED) = BOUNCEDBACK;

    if (iter % 10 == 0) {
      uu = moddims(sqrt(u_sq),nx,ny);
      uu(ON) = af::NaN;

      seq filter = seq(0,nx-1,nx/30);

      const char *str = "Velocity field for iteration ";
      std::stringstream title;
      title << str << iter;
      (*win)(0, 0).setColorMap(AF_COLORMAP_SPECTRUM);
      (*win)(0, 0).image(flip(transpose(normalize(uu)),0));
      (*win)(0, 1).setAxesLimits(0.0f,(float)nx,0.0f,(float)ny,true);
      (*win)(0, 1).vectorField(flat(x(filter,filter)), flat(y(filter,filter)), flat(UX(filter,filter)), flat(UY(filter,filter)), std::move(title).str().c_str());
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

  // while (!win->close())
  // {
  //   uu = sqrt(u_sq) / ux_lid;
  //   uu(ON) = af::NaN;

  //   seq filter = seq(0,nx-1,nx/30);

  //   const char *str = "Velocity field for iteration ";
  //   std::stringstream title;
  //   title << str << iter;
  //   (*win)(0, 0).setColorMap(AF_COLORMAP_SPECTRUM);
  //   (*win)(0, 0).image(flip(transpose(normalize(uu)),0));
  //   (*win)(0, 1).setAxesLimits(0.0f,(float)nx,0.0f,(float)ny,true);
  //   (*win)(0, 1).vectorField(flat(x(filter,filter)), flat(y(filter,filter)), flat(UX(filter,filter)), flat(UY(filter,filter)), std::move(title).str().c_str());
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
    printf("LBM D2Q9 simulation\n");
    lbm();
  }
  catch (af::exception &e)
  {
    fprintf(stderr, "%s\n", e.what());
    throw;
  }

  return 0;
}
