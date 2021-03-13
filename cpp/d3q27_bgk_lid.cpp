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
  const unsigned nx = 50;
  const unsigned ny = 50;
  const unsigned nz = 50;

  const unsigned total_nodes = nx * ny * nz;

  // Physical parameters.
  const float ux_lid  = 0.05; // horizontal lid velocity
  const float uy_lid = 0;
  const float uz_lid = 0;
  const float rho0 = 1.0;
  // Reynolds number
  float Re = 100.0;
  // Kinematic viscosity
  float nu = ux_lid * nx * ny / Re;
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

  // Indexes for directions, corner case
  int dirs1[3] = {0, 1, 3};
  int dirs2[3] = {2, 5, 6};
  array idxdirs1(3, dirs1);
  array idxdirs2(3, dirs2);

  // Start in equilibrium state
  array u_sq = pow(UX, 2) + pow(UY, 2) + pow(UZ, 2);
  array eu = (flat(tile(transpose(ex), total_nodes)) * tile(flat(UX),27)) + (flat(tile(transpose(ey), total_nodes)) * tile(flat(UY),27)) + (flat(tile(transpose(ez), total_nodes)) * tile(flat(UZ),27));
  F = flat(tile(transpose(w), total_nodes)) * tile(flat(DENSITY),27) * (1.0f + 3.0f*eu + 4.5f*(af::pow(eu,2)) - 1.5f*(tile(flat(u_sq),27)));
  F = moddims(F,nx,ny,nz,27);

  array uu = constant(0,nx,ny,nz);

  if (!console)
  {
    win = new Window(1536, 768, "LBM solver using ArrayFire");
    win->grid(1, 2);
  }

  unsigned iter = 0;
  unsigned maxiter = 15000;

  sync();
  timer::start();

  // while (iter < maxiter)
  while (!win->close())// (iter < maxiter) // (!win->close())
  {
    F = stream(F);

    array BOUNCEDBACK = F(TO_REFLECT); // Densities bouncing back at next timestep

    // Compute macroscopic variables
    DENSITY = sum(F, 3);

    array F_2D = moddims(F, total_nodes, 27);
    array F_t = transpose(F_2D);

    array fex = moddims(tile(transpose(ex), total_nodes) * F_2D,nx,ny,nz,27);
    array fey = moddims(tile(transpose(ey), total_nodes) * F_2D,nx,ny,nz,27);
    array fey = moddims(tile(transpose(ey), total_nodes) * F_2D,nx,ny,nz,27);

    UX = (sum(fex, 2) / DENSITY);
    UY = (sum(fey, 2) / DENSITY);

    // MACROSCOPIC (DIRICHLET) BOUNDARY CONDITIONS
    UX(lid,end) = ux_lid; // lid x - velocity
    UY(lid,end) = uy_lid; // lid y - velocity
    DENSITY(lid,end) = 1 / (1+UY(lid,end)) * (sum(F(lid,end,idxdirs1), 2) + 2*sum(F(lid,end,idxdirs2), 2));

    // MICROSCOPIC BOUNDARY CONDITIONS: LID (Zou/He BC)
    F(lid,end,4) = F(lid,end,2) - 2./3.*DENSITY(lid,end)*UY(lid,end);
    F(lid,end,8) = F(lid,end,6) + 1./2.*(F(lid,end,3)-F(lid,end,1))+1./2.*DENSITY(lid,end)*UX(lid,end) - 1./6.*DENSITY(lid,end)*UY(lid,end);
    F(lid,end,7) = F(lid,end,5) + 1./2.*(F(lid,end,1)-F(lid,end,3))-1./2.*DENSITY(lid,end)*UX(lid,end) - 1./6.*DENSITY(lid,end)*UY(lid,end);

    UX(ON) = 0;
    UY(ON) = 0;
    DENSITY(ON) = 0;

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
        uu = sqrt(u_sq) / ux_lid;
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

  sync();

  float end = timer::stop();
  float mlups = (total_nodes * iter * 10e-6) / end;
  printf("Iterations: %d\n", iter);
  printf("Time: %fs\n", end);
  printf("MLUPS: %f\n", mlups);
  af::info();

  while (!win->close())
  {
    uu = sqrt(u_sq) / ux_lid;
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
