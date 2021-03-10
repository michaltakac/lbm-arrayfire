#include <arrayfire.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
using namespace af;
int main(int argc, char *argv[]) {
  af::info();
  const unsigned nx = 1000, ny = 300;
  const unsigned total_nodes = nx * ny;
  const float rho0 = 1.0;
  const int obstacle_x = nx / 5 + 1; // x location of the cylinder
  const int obstacle_y = ny / 2 + ny / 30; // y location of the cylinder
  const int obstacle_r = ny / 10 + 1; // radius of the cylinder
  float Re = 220.0, u_max = 0.1; // Reynolds number, Lattice speed
  float nu = u_max * 2 *obstacle_r / Re; // Kinematic viscosity
  float tau = 3 * nu + 0.5; // Relaxation time
  float omega = 1.0 / tau; // Relaxation parameter
  const float t1 = 4. / 9., t2 = 1. / 9., t3 = 1. / 36.;
  array x = tile(range(nx), 1, ny);
  array y = tile(range(dim4(1, ny), 1), nx, 1);
  float cx[9] = {0, 1, 0,-1, 0, 1,-1,-1, 1};
  float cy[9] = {0, 0, 1, 0,-1, 1, 1,-1,-1};
  array ex(9, cx);
  array ey(9, cy);
  float weights[9] = {t1,t2,t2,t2,t2,t3,t3,t3,t3};
  array w(9, weights);
  array F = constant(rho0/9, nx, ny, 9);
  array FEQ = F.copy();
  array CI = (range(dim4(1,8),1)+1) * total_nodes;
  int nbindex[8] = {2,3,0,1,6,7,4,5};
  array nbidx(8, nbindex);
  array NBI = CI(span,nbidx);
  array BOUND = constant(0,nx,ny);
  BOUND(span,span) = moddims((af::pow(flat(x) - obstacle_x, 2) + af::pow(flat(y) - obstacle_y, 2)) <= pow(obstacle_r,2), nx, ny);
  BOUND(span,0) = 1; //top
  BOUND(span,end) = 1; //bottom
  array ON = where(BOUND); // matrix offset of each Occupied Node
  array TO_REFLECT = flat(tile(ON,CI.elements())) + flat(tile(CI,ON.elements()));
  array REFLECTED = flat(tile(ON,NBI.elements())) + flat(tile(NBI,ON.elements()));
  array DENSITY = constant(rho0, nx, ny);
  array UX = constant(u_max, nx, ny);
  array UY = constant(0, nx, ny);
  UX(ON) = 0;
  DENSITY(ON) = 0;
  array u_sq = pow(UX, 2) + pow(UY, 2);
  array eu = (flat(tile(transpose(ex), total_nodes)) * tile(flat(UX),9)) + (flat(tile(transpose(ey), total_nodes)) * tile(flat(UY),9));
  F = flat(tile(transpose(w), total_nodes)) * tile(flat(DENSITY),9) * (1.0f + 3.0f*eu + 4.5f*(af::pow(eu,2)) - 1.5f*(tile(flat(u_sq),9)));
  F = moddims(F,nx,ny,9);
  Window *win = new Window(1536, 768, "LBM solver using ArrayFire"); win->grid(2, 2);
  unsigned iter = 0;
  timer::start();
  while (!win->close()) {
    for (int i=0;i<9;i++) { // Streaming
      F(span, span, i) = shift(F, cx[i], cy[i])(span, span, i);
    }
    array BOUNCEDBACK = F(TO_REFLECT); // Densities bouncing back at next timestep
    array F_2D = moddims(F, total_nodes, 9);
    array F_t = transpose(F_2D);
    array fex = moddims(tile(transpose(ex), total_nodes) * F_2D,nx,ny,9);
    array fey = moddims(tile(transpose(ey), total_nodes) * F_2D,nx,ny,9);
    DENSITY = sum(F, 2);
    UX = (sum(fex, 2) / DENSITY);
    UY = (sum(fey, 2) / DENSITY);
    UX(0,span) = u_max;
    UX(ON) = 0;
    UY(ON) = 0;
    DENSITY(ON) = 0;
    DENSITY(0,span) = 1;
    DENSITY(end,span) = 1;
    u_sq = pow(UX, 2) + pow(UY, 2);
    eu = (flat(tile(transpose(ex), total_nodes)) * tile(flat(UX),9)) + (flat(tile(transpose(ey), total_nodes)) * tile(flat(UY),9));
    FEQ = flat(tile(transpose(w), total_nodes)) * tile(flat(DENSITY),9) * (1.0f + 3.0f*eu + 4.5f*(af::pow(eu,2)) - 1.5f*(tile(flat(u_sq),9)));
    FEQ = moddims(FEQ,nx,ny,9);
    F = omega * FEQ + (1 - omega) * F;
    F(REFLECTED) = BOUNCEDBACK;
    if (iter % 10 == 0) { // Visualization
      array uu = sqrt(u_sq);
      uu(ON) = af::NaN;
      (*win)(0, 0).setColorMap(AF_COLORMAP_SPECTRUM);
      (*win)(0, 0).image(transpose(uu));
      (*win)(0, 1).setAxesLimits(0.0f,(float)nx,0.0f,(float)ny,true);
      (*win)(0, 1).vectorField(flat(x(seq(0,nx-1,ny/15),seq(0,ny-1,ny/30))), flat(y(seq(0,nx-1,ny/15),seq(0,ny-1,ny/30))), flat(UX(seq(0,nx-1,ny/15),seq(0,ny-1,ny/30))), flat(UY(seq(0,nx-1,ny/15),seq(0,ny-1,ny/30))), "Velocity field");
      win->show();
    }
    iter++;
  }
  af::sync(0);
  float time = timer::stop();
  float mlups = (total_nodes * iter * 10e-6) / time;
  printf("%u iterations completed, %fs elapsed (%f MLUPS).\n", iter, time, mlups);
  return 0;
}
