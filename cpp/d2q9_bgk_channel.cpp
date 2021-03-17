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

array stream(array f)
{
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
  const unsigned nx = 700;
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
  float nu = u_max * 2 * obstacle_r / Re;
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
  array ey(9, cx);

  // weights
  float weights[9] = {t1,t2,t2,t2,t2,t3,t3,t3,t3};
  array w(9, weights);

  array CI = (range(dim4(1,8),1)+1) * total_nodes;
  float nb_index_arr[8] = {2,3,0,1,6,7,4,5};
  array nbidx(8, nb_index_arr);
  array NBI = CI(span,nbidx);

  array main_index = moddims(range(dim4(total_nodes*9)),nx,ny,9);
  array nb_index = flat(stream(main_index));

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
  array u_sq = flat(pow(UX, 2) + pow(UY, 2));
  array eu = flat(batchFunc(transpose(ex), flat(UX), mul) + batchFunc(transpose(ey), flat(UY), mul));
  array F = flat(batchFunc(transpose(w), flat(DENSITY), mul)) * (1.0f + 3.0f*eu + 4.5f*(af::pow(eu,2)) - 1.5f*(tile(u_sq,9)));

  array uu;

  // Setup Window
  win = new Window(1536, 1024, "LBM solver using ArrayFire");
  win->grid(2, 1);

  unsigned iter = 0;
  unsigned maxiter = 10000;
  float MLUPS[10000];

  sync(0);
  timer::start();

  while (!win->close() && iter < maxiter)
  {
    // Streaming by reading from neighbors (with pre-built index) - pull scheme
    array F_streamed = F(nb_index);

    array BOUNCEDBACK = F_streamed(TO_REFLECT); // Densities bouncing back at next timestep

    array F_2D = moddims(F_streamed, total_nodes, 9);

    // Compute macroscopic variables
    array rho = sum(F_2D, 1);
    DENSITY = moddims(rho,nx,ny);

    array fex = batchFunc(transpose(ex), F_2D, mul);
    array fey = batchFunc(transpose(ey), F_2D, mul);

    UX = moddims((sum(fex, 1) / rho),nx,ny);
    UY = moddims((sum(fey, 1) / rho),nx,ny);

    UX(0,span) = u_max;
    UX(ON) = 0;
    UY(ON) = 0;
    DENSITY(ON) = 0;
    DENSITY(0,span) = 1;
    DENSITY(end,span) = 1;

    // Collision
    u_sq = flat(pow(UX, 2) + pow(UY, 2));
    eu = flat(batchFunc(transpose(ex), flat(UX), mul) + batchFunc(transpose(ey), flat(UY), mul));
    array FEQ = flat(batchFunc(transpose(w), flat(DENSITY), mul)) * (1.0f + 3.0f*eu + 4.5f*(af::pow(eu,2)) - 1.5f*(tile(u_sq,9)));

    F = omega * FEQ + (1 - omega) * F_streamed;

    F(REFLECTED) = BOUNCEDBACK;

    if (iter % 10 == 0) {
      uu = moddims(sqrt(u_sq),nx,ny);
      uu(ON) = af::NaN;

      seq filterX = seq(0,nx-1,nx/15);
      seq filterY = seq(0,ny-1,ny/30);

      const char *str = "Velocity field for iteration ";
      std::stringstream title;
      title << str << iter;
      (*win)(0, 0).setColorMap(AF_COLORMAP_SPECTRUM);
      (*win)(0, 0).image(transpose(normalize(uu)));
      (*win)(1, 0).setAxesLimits(0.0f,(float)nx,0.0f,(float)ny,true);
      (*win)(1, 0).vectorField(flat(x(filterX,filterY)), flat(y(filterX,filterY)), flat(UX(filterX,filterY)), flat(UY(filterX,filterY)), std::move(title).str().c_str());
      win->show();
    }

    float time = timer::stop();
    MLUPS[iter] = (total_nodes * iter * 10e-6) / time;

    if (iter % 100 == 0) {
      printf("%u iterations completed, %fs elapsed (%f MLUPS).\n", iter, time, MLUPS[iter]);
    }

    sync(0);
    iter++;
  }

  float end = timer::stop();
  float mlups = (total_nodes * iter * 10e-6) / end;
  printf("Iterations: %d\n", iter);
  printf("Time: %fs\n", end);
  printf("MLUPS: %f\n", mlups);
  af::info();

  // array mlups_y(maxiter, MLUPS);
  // array mlups_x = range(maxiter);

  // while (!win->close())
  // {
  //   uu = sqrt(u_sq);
  //   uu(ON) = af::NaN;

  //   seq filterX = seq(0,nx-1,ny/10);
  //   seq filterY = seq(0,ny-1,ny/20);

  //   const char *str = "Velocity field for iteration ";
  //   std::stringstream title;
  //   title << str << iter;
  //   (*win)(0, 0).setColorMap(AF_COLORMAP_SPECTRUM);
  //   (*win)(0, 0).image(transpose(normalize(uu)));
  //   (*win)(1, 0).setAxesLimits(0.0f,(float)nx,0.0f,(float)ny,true);
  //   (*win)(1, 0).vectorField(flat(x(filterX,filterY)), flat(y(filterX,filterY)), flat(UX(filterX,filterY)), flat(UY(filterX,filterY)), std::move(title).str().c_str());
  //   (*win)(2, 0).setAxesLabelFormat("%.1f","%.1f");
  //   (*win)(2, 0).setAxesTitles("N of Iterations", "MLUPS");
  //   (*win)(2, 0).setAxesLimits(mlups_x, mlups_y);
  //   (*win)(2, 0).plot(mlups_x, mlups_y, "MLUPS");
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
