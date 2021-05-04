#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef struct _PARAMS_{
  double AP;
  double AE;
  double AW;
  double AN;
  double AS;
  double AEE;
  double AWW;
  double ANN;
  double ASS;
}PARAMS, *PPARAMS;

typedef struct _GRID_{
  double x;
  double y;
}GRID, *PGRID;

typedef struct _U_{
  double value;
  GRID grid;
  PARAMS params;
  PARAMS pressure_params;
}U, *PU;

typedef struct _V_{
  double value;
  GRID grid;
  PARAMS params;
  PARAMS pressure_params
}V, *PV;

typedef struct _P_{
  double value;
  GRID grid;
  PARAMS pressure_params;   // TODO:
}P, *PP;

// allocate the memory
void AllocateMemory(U** u, V** v, P** p, int nx, int ny)
{
  *u = (U*)malloc(sizeof(U)*nx*ny);
  *v = (V*)malloc(sizeof(V)*nx*ny);
  *p = (P*)malloc(sizeof(P)*nx*ny);
  // set the memory to zero
  memset(*u, 0, sizeof(U)*nx*ny);
  memset(*v, 0, sizeof(V)*nx*ny);
  memset(*p, 0, sizeof(P)*nx*ny);
}

// free the memory
void FreeMemory(U* u, V* v, P* p)
{
  free(u);
  free(v);
  free(p);
  u = NULL;
  v = NULL;
  p = NULL;
}

// grid generation, only one grid system
// walk through the J direction first.
//          =========================== ->ue
//      j-> |                         |
//          |                         |
//          |                         |
//          O==========================
//                   /|\
//                    |
//                    i
void Grid(U* u, V* v, P* p, int nx, int ny, double X, double Y)
{
  int i, j;
  double dx = X / (double)nx;
  double dy = Y / (double)ny;
  for (i = 0; i < nx; i++){
    for (j = 0; j < ny; j++){
      u[i*ny + j].grid.x = i*dx;
      u[i*ny + j].grid.y = j*dy;
      v[i*ny + j].grid.x = i*dx;
      v[i*ny + j].grid.y = j*dy;
      p[i*ny + j].grid.x = i*dx;
      p[i*ny + j].grid.y = j*dy;
      printf("(%lf, %lf)\n", p[i*ny + j].grid.x, p[i*ny + j].grid.y);
    }
  }
}

// initialize the field
void Initialize(U *u, V* v, P* p, 
  int nx, int ny, 
  double init_u, double init_v, double init_p)
{
  int i, j;
  for (j = 0; j < ny; j++){
    for (i = 0; i < nx; i++){
      u[i*ny + j].value = init_u;
      v[i*ny + j].value = init_v;
      p[i*ny + j].value = init_p;
    }
  }
}

void DebugPrintGrid(U* u, V* v, P* p, int nx, int ny)
{
  int i, j;
  FILE *fp = NULL;
  fp = fopen("grid.txt", "w");
  if (fp == NULL){
    printf("[*] DebugPrintGrid : can not open the file....\n");
    return;
  }
  fputs("U location\n", fp);
  for (j = ny - 1; j >= 0; j--){
    for (i = 0; i < nx; i++){
      fprintf(fp, "(%lf, %lf)    ", u[i*ny + j].grid.x, u[i*ny + j].grid.y);
    }
    fputs("\n", fp);
  }
  fputs("V location\n", fp);
  for (j = ny - 1; j >= 0; j--){
    for (i = 0; i < nx; i++){
      fprintf(fp, "(%lf, %lf)    ", v[i*ny + j].grid.x, v[i*ny + j].grid.y);
    }
    fputs("\n", fp);
  }
  fputs("P location\n", fp);
  for (j = ny - 1; j >= 0; j--){
    for (i = 0; i < nx; i++){
      fprintf(fp, "(%lf, %lf)    ", p[i*ny + j].grid.x, p[i*ny + j].grid.y);
    }
    fputs("\n", fp);
  }
  fclose(fp);
}

// Boundary dynamic_viscosityition ( driven cavity flow )
void Boundary(U* u, V* v, P* p, int nx, int ny, double boundary_u)
{
  int i, j;
  // upper wall
  for (i = 0; i < nx; i++){
    u[i*ny + ny - 1].value = boundary_u;
    u[i*ny + ny - 2].value = boundary_u;
    v[i*ny + ny - 1].value = 0.0;
    v[i*ny + ny - 2].value = 0.0;
  }
  // bottom wall
  for (i = 0; i < nx; i++){
    u[i*ny + 0].value = 0.0;
    u[i*ny + 1].value = 0.0;
    v[i*ny + 0].value = 0.0;
    v[i*ny + 1].value = 0.0;
  }
  // left wall
  for (j = 0; j < ny; j++){
    u[0 * ny + j].value = 0.0;
    u[1 * ny + j].value = 0.0;
    v[0 * ny + j].value = 0.0;
    v[1 * ny + j].value = 0.0;
  }
  // right wall
  for (j = 0; j < ny; j++){
    u[(nx - 1)*ny + j].value = 0.0;
    u[(nx - 2)*ny + j].value = 0.0;
    v[(nx - 1)*ny + j].value = 0.0;
    v[(nx - 2)*ny + j].value = 0.0;
  }
}

// Debug the boundary
void DebugPrintBoundary(U* u, V* v, P* p, int nx, int ny)
{
  int i, j;
  FILE *fp = NULL;
  fp = fopen("boundary.txt", "w");
  if (fp == NULL){
    printf("[*] DebugPrintBoundary : can not open the file....\n");
    return;
  }
  fputs("U\n", fp);
  for (j = ny - 1; j >= 0; j--){
    for (i = 0; i < nx; i++){
      fprintf(fp, "%6.2f    ", u[i*ny + j].value);
    }
    fputs("\n", fp);
  }
  fputs("V\n", fp);
  for (j = ny - 1; j >= 0; j--){
    for (i = 0; i < nx; i++){
      fprintf(fp, "%6.2f    ", v[i*ny + j].value);
    }
    fputs("\n", fp);
  }
  fclose(fp);
}

// get the coefficient of momentum equation on X direction
void MomentumX(U* u, V* v, P* p, 
  int nx, int ny,
  U* last_u, V* last_v, P* last_p,
  double rho, double dynamic_viscosity, 
  double dx, double dy)
{
  int i, j;
  double flow_u = 0.0;
  double flow_v = 0.0;
  int x = nx - 4;
  int y = ny - 4;
  for (i = 2; i < nx-2; i++){
    for (j = 2; j < ny-2; j++){
      flow_u = rho*last_u[i*ny + j].value;
      flow_v = rho*last_v[i*ny + j].value;
      if (u[i*ny + j].value>0){
        // velocity parameters
        u[i*ny + j].params.AP = (3 * flow_u / (2 * dx)) +
          (3 * flow_v / (2 * dy)) +
          (2 * dynamic_viscosity / (dx*dx)) +
          (2 * dynamic_viscosity / (dy*dy));
        u[i*ny + j].params.AW = (4 * flow_u / (2 * dx)) + dynamic_viscosity / (dx*dx);
        u[i*ny + j].params.AWW = -2 * flow_u / (2 * dx);
        u[i*ny + j].params.AE = dynamic_viscosity / (dx*dx);
        u[i*ny + j].params.AS = (4 * flow_v / (2 * dy)) + dynamic_viscosity / (dy*dy);
        u[i*ny + j].params.ASS = -flow_v / (2 * dy);
        u[i*ny + j].params.AN = dynamic_viscosity / (dy*dy);
        u[i*ny + j].params.ANN = 0.0;
        u[i*ny + j].params.AEE = 0.0;
        // pressure parameters
        u[i*ny + j].pressure_params.AP = -3 / (2 * dx);
        u[i*ny + j].pressure_params.AW = 4 / (2 * dx);
        u[i*ny + j].pressure_params.AWW = -1 / (2 * dx);
        u[i*ny + j].pressure_params.AE = 0;
        u[i*ny + j].pressure_params.AEE = 0;
        u[i*ny + j].pressure_params.AS = 0;
        u[i*ny + j].pressure_params.ASS = 0;
        u[i*ny + j].pressure_params.AN = 0;
        u[i*ny + j].pressure_params.ANN = 0;
      }
      else if (u[i*ny + j].value < 0){
        // TODO:

      }
    }
  }
}

// get the coefficient of momentum equation on Y direction
void MomentumY(U*u, V* v, P* p, 
  int nx, int ny, 
  U* last_u, V* last_v, P* last_p,
  double rho, double dynamic_viscosity)
{
  int i, j;
  int x = nx - 4;
  int y = ny - 4;
  for (i = 2; i < nx - 2; i++){
    for (j = 2; j < ny - 2; j++){
      if (v[i*ny + j].value>0){
        // todo

      }
    }
  }
}

int main(int argc, char *argv[])
{
  double X = 10.0;
  double Y = 5.0;
  int nx = 20;      // col
  int ny = 20;      // row
  double init_u = 1.0;
  double init_v = 1.0;
  double init_p = 0.0;

  double boundary_u = 100.0;

  U* u = NULL;
  V* v = NULL;
  P* p = NULL;
  
  AllocateMemory(&u, &v, &p, nx, ny);
  printf("Allocate memory has done...\n");

  Grid(u, v, p, nx, ny, X, Y);
  printf("Grid generation has done...\n");

  DebugPrintGrid(u, v, p, nx, ny);
  printf("DebugPrintGrid has done...\n");

  Initialize(u, v, p, nx, ny, init_u, init_v, init_p);
  printf("Initialize has done...\n");

  Boundary(u, v, p, nx, ny, boundary_u);
  printf("Boundary has done...\n");

  DebugPrintBoundary(u, v, p, nx, ny);

  FreeMemory(u, v, p);
  printf("All the allocated memory has been freed...\n");
  printf("Done...\n");
  return 0;
}