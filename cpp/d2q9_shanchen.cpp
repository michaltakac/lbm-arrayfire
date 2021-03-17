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
  const unsigned nx = 201;
  const unsigned ny = 201;

  const unsigned total_nodes = nx * ny;

  // Physical parameters.
  const float G = -1.2f; // Amplitude of the molecular interaction force
  const float omega1 = 1.f;  // Relaxation parameter for fluid 1
  const float omega2 = 1.f;  // Relaxation parameter for fluid 2
  const float drho = 0.001f;

  array delta_rho = -drho * (1.f - 2.f * randu(nx, ny));

  float Gomega1 = G/omega1;
  float Gomega2 = G/omega2;

  const float t1 = 4.f / 9.f;
  const float t2 = 1.f / 9.f;
  const float t3 = 1.f / 36.f;

  array x = tile(range(nx), 1, ny);
  array y = tile(range(dim4(1, ny), 1), nx, 1);

  //  c6  c2   c5
  //    \  |  /
  //  c3 -c0 - c1
  //    /  |  \
  //  c7  c4   c8
  // Discrete velocities
  int cx[9] = {0, 1, 0,-1, 0, 1,-1,-1, 1};
  int cy[9] = {0, 0, 1, 0,-1, 1, 1,-1,-1};
  array ex(9, cx);
  array ey(9, cy);

  // weights
  float weights[9] = {t1,t2,t2,t2,t2,t3,t3,t3,t3};
  array w(9, weights);

  array main_index = moddims(range(dim4(total_nodes*9)),nx,ny,9);
  array nb_index = flat(stream(main_index));

  // Initial condition for both distribution functions: (T=0) ==> TIn(i) = t(i)
  array f = batchFunc(transpose(w), flat(1.f + delta_rho), mul);
  array g = batchFunc(transpose(w), flat(1.f - delta_rho), mul);
  printf("f dims = [%lld %lld %lld]\n", f.dims(0), f.dims(1), f.dims(2));
  printf("g dims = [%lld %lld %lld]\n", g.dims(0), g.dims(1), g.dims(2));
  array f_out = f.copy();
  array g_out = g.copy();

  array feq = f.copy();
  array geq = g.copy();

  array rhoContrib1x = constant(0.f, total_nodes);
  array rhoContrib2x = constant(0.f, total_nodes);

  array rhoContrib1y = constant(0.f, total_nodes);
  array rhoContrib2y = constant(0.f, total_nodes);

  // Setup Window
  win = new Window(800, 800, "LBM solver using ArrayFire");
  win->grid(1, 1);

  unsigned iter = 0;
  unsigned maxiter = 15000;

  sync();
  timer::start();

  while (!win->close() && iter < maxiter)
  {
    // Streaming both fluids by reading from neighbors (with pre-built index) - pull scheme
    array F_streamed = f(nb_index);
    array G_streamed = g(nb_index);

    array F_2D = moddims(F_streamed, total_nodes, 9);
    array G_2D = moddims(G_streamed, total_nodes, 9);

    // Compute macroscopic variables
    array rho1 = sum(F_2D, 1);
    array rho2 = sum(G_2D, 1);
    printf("rho1 dims = [%lld %lld %lld]\n", rho1.dims(0), rho1.dims(1), rho1.dims(2));

    array fex = batchFunc(transpose(ex), F_2D, mul);
    array fey = batchFunc(transpose(ey), F_2D, mul);
    array gex = batchFunc(transpose(ex), G_2D, mul);
    array gey = batchFunc(transpose(ey), G_2D, mul);
    printf("fex dims = [%lld %lld %lld]\n", fex.dims(0), fex.dims(1), fex.dims(2));

    array JX1 = sum(fex, 1) / rho1;
    array JY1 = sum(fey, 1) / rho1;
    array JX2 = sum(gex, 1) / rho2;
    array JY2 = sum(gey, 1) / rho2;
    printf("JX1 dims = [%lld %lld %lld]\n", JX1.dims(0), JX1.dims(1), JX1.dims(2));

    array rho_to_t_omega = rho1 * omega1 + rho2 * omega2;
    printf("rho_to_t_omega dims = [%lld %lld %lld]\n", rho_to_t_omega.dims(0), rho_to_t_omega.dims(1), rho_to_t_omega.dims(2));
    array uTotX = (JX1*omega1+JX2*omega2) / rho_to_t_omega;
    array uTotY = (JY1*omega1+JY2*omega2) / rho_to_t_omega;
    printf("uTotX dims = [%lld %lld %lld]\n", uTotX.dims(0), uTotX.dims(1), uTotX.dims(2));

    for (int i=1; i < 9; i++) {
        rhoContrib1x += shift(rho1*weights[i], cx[i], cy[i]) * cx[i]; // TODO: nesedi!
        rhoContrib1y += shift(rho1*weights[i], cx[i], cy[i]) * cy[i];

        rhoContrib2x += shift(rho2*weights[i], cx[i], cy[i]) * cx[i];
        rhoContrib2y += shift(rho2*weights[i], cx[i], cy[i]) * cy[i];
    }
printf("rhoContrib2x dims = [%lld %lld %lld]\n", rhoContrib2x.dims(0), rhoContrib2x.dims(1), rhoContrib2x.dims(2));
    array uTotX1 = uTotX - Gomega1 * rhoContrib2x;
    array uTotY1 = uTotY - Gomega1 * rhoContrib2y;

    array uTotX2 = uTotX - Gomega2 * rhoContrib1x;
    array uTotY2 = uTotY - Gomega2 * rhoContrib1y;

    printf("uTotX1 dims = [%lld %lld %lld]\n", uTotX1.dims(0), uTotX1.dims(1), uTotX1.dims(2));

    // Collision for both fluids
    for (int i=0; i < 9; i++) {
      array cu1 = 3*(cx[i]*uTotX1+cy[i]*uTotY1);
      array cu2 = 3*(cx[i]*uTotX2+cy[i]*uTotY2);
      printf("cu1 dims = [%lld %lld %lld]\n", cu1.dims(0), cu1.dims(1), cu1.dims(2));

      feq(span,i) = rho1 * weights[i] * (1.f + cu1 + 0.5f *(cu1*cu1) - 1.5f * (pow(uTotX1,2)+pow(uTotY1,2)));
      geq(span,i) = rho2 * weights[i] * (1.f + cu2 + 0.5f *(cu2*cu2) - 1.5f * (pow(uTotX2,2)+pow(uTotY2,2)));
      printf("f dims = [%lld %lld %lld]\n", f.dims(0), f.dims(1), f.dims(2));

      f_out(span,i)  = f(span,i) - omega1 * (f(span,i)-feq(span,i));
      g_out(span,i)  = g(span,i) - omega2 * (g(span,i)-geq(span,i));
    }

    // STREAMING STEP FLUID 1 AND 2
    // for (int i=0; i < 9; i++) {
    //   f(i,span,span) = shift(f_out(i,span,span), 0,cx[i],cy[i]);
    //   g(i,span,span) = shift(g_out(i,span,span), 0,cx[i],cy[i]);
    // }

    if (iter % 10 == 0) {
      const char *str = "Fluid 1 density, iteration ";
      std::stringstream title;
      title << str << iter;
      (*win)(0, 0).setColorMap(AF_COLORMAP_SPECTRUM);
      (*win)(0, 0).image(moddims(rho1, nx, ny), std::move(title).str().c_str());
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
