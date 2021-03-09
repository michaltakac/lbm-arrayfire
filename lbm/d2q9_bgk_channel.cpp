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

static void lbm(bool console)
{
  // Grid length, number and spacing
  const unsigned nx = 1000;
  const unsigned ny = 300;

  const unsigned total_nodes = nx * ny;

  // Physical parameters.
  const float rho0 = 1.0;

  const int obstacle_x = nx / 5 + 1; // x location of the cylinder
  const int obstacle_y = ny / 2 + ny / 30; // y location of the cylinder
  const int obstacle_r = ny / 10 + 1; // radius of the cylinder

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
  const float c_squ = 1. / 3.;

  array x = tile(range(nx), 1, ny);
  array y = tile(range(dim4(1, ny), 1), nx, 1);

  //  c6  c2   c5
  //    \  |  /
  //  c3 -c0 - c1
  //    /  |  \
  //  c7  c4   c8
  // Discrete velocities
  float cx[9] = {0, 1, 0,-1, 0, 1,-1,-1, 1};
  float cy[9] = {0, 0, 1, 0,-1, 1, 1,-1,-1};
  array ex(9, cx);
  array ey(9, cy);

  // weights
  float weights[9] = {t1,t2,t2,t2,t2,t3,t3,t3,t3};
  array w(9, weights);

  array F = constant(rho0/9, nx, ny, 9);
  array FEQ = F.copy();

  array CI = (range(dim4(1,8),1)+1) * total_nodes;
  int nbindex[8] = {2,3,0,1,6,7,4,5};
  array nbidx(8, nbindex);
  array NBI = CI(span,nbidx);

  // Flow around obstacle
  // circle
  array BOUND = constant(0,nx,ny);
  BOUND(span,span) = moddims((af::pow(flat(x) - obstacle_x, 2) + af::pow(flat(y) - obstacle_y, 2)) <= pow(obstacle_r,2), nx, ny);
  BOUND(span,0) = 1; //top
  BOUND(span,end) = 1; //bottom

  // matrix offset of each Occupied Node
  array ON = where(BOUND);

  // Bounceback indexes
  array TO_REFLECT = flat(tile(ON,CI.elements())) + flat(tile(CI,ON.elements()));
  array REFLECTED = flat(tile(ON,NBI.elements())) + flat(tile(NBI,ON.elements()));

  array DENSITY = constant(rho0, nx, ny);
  array UX = constant(u_max, nx, ny);
  array UY = constant(0, nx, ny);

  UX(ON) = 0;
  UY(ON) = 0;
  DENSITY(ON) = 0;

  // Start in equilibrium state
  array u_sq = pow(UX, 2) + pow(UY, 2);
  array eu = (flat(tile(transpose(ex), total_nodes)) * tile(flat(UX),9)) + (flat(tile(transpose(ey), total_nodes)) * tile(flat(UY),9));
  F = flat(tile(transpose(w), total_nodes)) * tile(flat(DENSITY),9) * (1.0f + 3.0f*eu + 4.5f*(af::pow(eu,2)) - 1.5f*(tile(flat(u_sq),9)));
  F = moddims(F,nx,ny,9);
  
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
  unsigned maxiter = 5000;

  sync(0);
  timer::start();

  // while (iter < maxiter)
  while (!win->close())// (iter < maxiter) // (!win->close())
  {
    F = stream(F);

    array BOUNCEDBACK = F(TO_REFLECT); // Densities bouncing back at next timestep

    // Compute macroscopic variables
    DENSITY = sum(F, 2);

    array F_2D = moddims(F, total_nodes, 9);
    array F_t = transpose(F_2D);

    array fex = moddims(tile(transpose(ex), total_nodes) * F_2D,nx,ny,9);
    array fey = moddims(tile(transpose(ey), total_nodes) * F_2D,nx,ny,9);

    UX = (sum(fex, 2) / DENSITY);
    UY = (sum(fey, 2) / DENSITY);

    UX(0,span) = u_max;
    UX(ON) = 0;
    UY(ON) = 0;
    DENSITY(ON) = 0;
    DENSITY(0,span) = 1;
    DENSITY(end,span) = 1;

    // Collision
    u_sq = pow(UX, 2) + pow(UY, 2);
    eu = (flat(tile(transpose(ex), total_nodes)) * tile(flat(UX),9)) + (flat(tile(transpose(ey), total_nodes)) * tile(flat(UY),9));
    FEQ = flat(tile(transpose(w), total_nodes)) * tile(flat(DENSITY),9) * (1.0f + 3.0f*eu + 4.5f*(af::pow(eu,2)) - 1.5f*(tile(flat(u_sq),9)));
    FEQ = moddims(FEQ,nx,ny,9);

    F = omega * FEQ + (1 - omega) * F;

    F(REFLECTED) = BOUNCEDBACK;

    if (!console)
    {
      if (iter % 10 == 0) {
        prevavu = avu;
        avu = sum<float>(sum(UX)) / numactivenodes;
        uu = sqrt(u_sq);
        uu(ON) = af::NaN;

        seq filterX = seq(0,nx-1,ny/15);
        seq filterY = seq(0,ny-1,ny/30);

        const char *str = "Velocity field for iteration ";
        std::stringstream title;
        title << str << iter;
        (*win)(0, 0).setColorMap(AF_COLORMAP_SPECTRUM);
        (*win)(0, 0).image(transpose(uu));
        (*win)(0, 1).setAxesLimits(0.0f,(float)nx,0.0f,(float)ny,true);
        (*win)(0, 1).vectorField(flat(x(filterX,filterY)), flat(y(filterX,filterY)), flat(UX(filterX,filterY)), flat(UY(filterX,filterY)), std::move(title).str().c_str());
        (*win)(1, 0).setColorMap(AF_COLORMAP_HEAT);
        (*win)(1, 0).image(normalize(transpose(DENSITY), max<float>(DENSITY)));
        win->show();
      }
    }
    else
    {
      // eval(uu, F, FEQ);
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

  while (!win->close())
  {
    (*win)(0, 0).setColorMap(AF_COLORMAP_SPECTRUM);
    (*win)(0, 0).image(transpose(uu));
    (*win)(0, 1).vectorField(flat(x), flat(y), flat(UX), flat(UY), "Velocity field");
    (*win)(1, 0).image(transpose(DENSITY));
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
  catch (af::exception &e)
  {
    fprintf(stderr, "%s\n", e.what());
    throw;
  }

  return 0;
}
