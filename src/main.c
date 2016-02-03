/*************** 2D LBM-DEM Code **************************/
#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <time.h>

#define _FLUIDE_  // Switch on or off FLUID

// Maximum number of soil grains
#define nbgrainsMax 5000

// Minimum number of soil grains
#define nbgrains_bas 1

// Dimension of the LBM Fluid domain
#define lx 2000
#define ly 1000

#define pi 3.14159265358979

#define rhoS 2650  // Density of solids
#define rhoW 1000  // Density of water

#define duration 1.5 // Duration of simulation
//*********************   Data LBM    ************************
#define Q 9
int nbgrains;
// width of LBM grid size, time step, lattice speed
double dx, dtLB, c, c_squ;
double w[Q] = {4. / 9,  1. / 36, 1. / 9,  1. / 36, 1. / 9,
               1. / 36, 1. / 9,  1. / 36, 1. / 9};
double f[lx][ly][Q];

// ************************************
// *                                  *
// *         e2   e9   e8             *
// *          \   |   /               *
// *            \ | /                 *
// *        e3--- e1 ---e7            *
// *            / |  \                *
// *          /   |    \              *
// *        e4    e5   e6             *
// *                                  *
// ************************************

int ex[Q] = {0, -1, -1, -1, 0, 1, 1, 1, 0};
int ey[Q] = {0, 1, 0, -1, -1, -1, 0, 1, 1};

// average fluid density
double rho_moy = 1000;  // air density =1 or water =1000 or 999.7 at 20
double rho_outlet, q_outlet;

// relaxation parameter
double tau = 0.506;
double s2 = 1.5, s3 = 1.4, s5 = 1.5, s7 = 1.5, s8 = 1.8868,
       s9 = 1.8868;  // s8=1.6666667,s9=1.6666667;

// obstacle array
int obst[lx][ly];
double grains[lx][ly];
// obstacle activity  array
int act[lx][ly];
double delta[lx][ly][Q];

// LB diameter for the smallest disk (in nodes number)
double rMin_LB = 10.;

// Fluid kinematic viscosity
double nu = 1e-6;  // 15.5e-6 for air and 1e-6 for water at 293K or 20C
double press[lx][ly];
double reductionR = 0.9;  // LBM reduced grain diameter

//***********   Data DEM    ********************
double G = 9.81;
double angleG = 0.;
double xG, yG;
double dt;  //=5.e-8;
double dt2;

// Spring stiffness
double km = 3e+6, kg = 1.6e+6;  // 2e8 1.6e8
double kt = 1.0e+6, ktm = 2e+6;  /// changed to higher for lesser
                                 /// interpenetration in sample generation 1.3e8
double nug = 6.4e+1;  // 1.1e1
double num = 8.7e+1, numb = 8.7e+1;  // 1.5e1
double nuf = 1.5e-1, nugt = 5e-1;  // frictionless packing nugt
double mu = .5317;
double mug = 0.0;  // Mu for assembling
double mum = .466, mumb = .466;  // 0.466 //0.53 0.51 0.43
double murf = 0.01;  // 0.01
double r = 1e-3;  // 5e-4;v
double distVerlet = 5e-4;  // changed to 1e-6 from 1e-3 for error
long UpdateVerlet = 100.;
double dtt = 0.0;  // Time after which wall to prepare sample is removed
double iterDEM = 100.;  // number of DEM iterations per LBM

// Tracked stats
// EPE-Effective potential energy
// Total Wall Friction - WF,
// SE- Strain Energy ESE & IFR- Total Internal Friction
double xfront, height, energie_cin, energie_x, energie_y, energie_teta,
    energy_p, energy_EPE, zmean, SE, ESE, WF,
    IFR;
double TSE = 0.0, TBW = 0.0, INCE = 0.0, TSLIP = 0.0,
       TRW =
           0.0;  // Total Body Work and Strain Energy TRW_ Total Rotational Work
double pf = 0., pft = 0., pff = 0.;  // previous force
double ic = 0;

// ********   Control parameters   *************
// Number of DEM steps in LB
int npDEM;
double rLB[nbgrainsMax];

int stepView = 100;
int stepPrint = 200;
int stepConsole = 100;

int stepStrob = 4000;  //visualisation steps
int stepFilm = 2000;

FILE* s_stats;
char filename_sample[25];

int nFile = 0;  // Nth File saves LB

int cptFlash;
int cumul[nbgrainsMax];
int neighbours[nbgrainsMax * 10];
// NeighbourWall Bottom, Right, Left & Top
int neighbourWallB[nbgrainsMax];
int neighbourWallR[nbgrainsMax];
int neighbourWallL[nbgrainsMax];
int neighbourWallT[nbgrainsMax];

int nNeighWallb, nNeighWallt, nNeighWallL, nNeighWallR;

int start = 0;
long nbsteps = 0;
int vib = 0;
double freq = 5;
double amp = 4.e-4;
double t = 0;

// Luding Friction Model
struct contact {
  int i, j;
  double nx, ny;  // normal vector from j to i
  double fn, ft;  // force component in local frame
} c1[nbgrainsMax][nbgrainsMax], c2[nbgrainsMax];

struct force {
  double f1, f2, f3;
};
struct force fhf[nbgrainsMax];
double fhf1[nbgrainsMax], fhf2[nbgrainsMax], fhf3[nbgrainsMax];

struct grain {
  double x1, x2, x3;
  double v1, v2, v3;
  double a1, a2, a3;
  double r, m, mw, It;
  double p;  // Pressure on grain
  double s;  // shear
  double f1, f2;  // force
  double ifm, fm;  // friction mobility
  double fr, ifr;  // frictional energy wall and Internal
  double M11, M12, M21, M22;  // Moments M11, M12, M21, M22
  double ice, slip,
      rw;  // Inelastic Collisional Energy, slip, & Rotational work
  int z;   // number of contacts
  int zz;  // number of contacts sans the Walls
};
struct grain g[nbgrainsMax];
static double rMax, rMin;

// Wall
static double Mby = 0.;
static double Mgx = 0.;
static double Mhy = 0.;
static double Mdx = 0.;


double aff_transx = 2e-1, aff_transy = 5e-0, aff_zoom = 150.;
double ecart_ini = 1.0;

// ***************************************************************************
// *   utilities
// ***************************************************************************

double Min(double x, double y) {
  if (x < y)
    return x;
  else
    return y;
}

double Max(double x, double y) {
  if (x < y)
    return y;
  else
    return x;
}

double Maxt(double x, double y) {
  if (x < y)
    return 0.;
  else
    return y;
}

//---------------------------------------------------
void stat() {
  int i;

  rMax = g[0].r;
  rMin = g[0].r;

  for (i = 0; i <= nbgrains - 1; i++) {
    rMax = Max(rMax, g[i].r);
    rMin = Min(rMin, g[i].r);
  }
}
//----------------------------------------------------
void swap(double* a, double* b) {
  double tmp;
  tmp = *a;
  *a = *b;
  *b = tmp;
}

// *******************************************************************
// *   Output files                                             *
// *******************************************************************
void write_grains() {
  int x, y, i;
  char filename[25];
  double pasxyz;
  FILE* outfile;
  sprintf(filename, "LBM_grains%.6i.vtk", nFile);

  pasxyz = 1. / lx;
  outfile = fopen(filename, "w");
  fprintf(outfile, "# vtk DataFile Version 2.0\n");
  fprintf(outfile, "Domain LB+LINK\n");
  fprintf(outfile, "ASCII\n");
  fprintf(outfile, "DATASET RECTILINEAR_GRID\n");
  fprintf(outfile, "DIMENSIONS %d %d 1\n", lx, ly);
  fprintf(outfile, "X_COORDINATES %d float\n", lx);
  for (i = 0; i <= lx - 1; i++) {
    fprintf(outfile, "%e ", (float)i * pasxyz);
  }
  fprintf(outfile, "\n");
  fprintf(outfile, "Y_COORDINATES %d float\n", ly);
  for (i = 0; i <= ly - 1; i++) {
    fprintf(outfile, "%e ", (float)i * pasxyz);
  }
  fprintf(outfile, "\n");
  fprintf(outfile, "Z_COORDINATES 1 float\n");
  fprintf(outfile, "0\n");

  // For LB
  fprintf(outfile, "POINT_DATA %d\n", lx * ly);
  fprintf(outfile, "SCALARS Grains float 1\n");
  fprintf(outfile, "LOOKUP_TABLE default\n");

  for (y = 0; y < ly; y++) {
    for (x = 0; x < lx; x++) {
      fprintf(outfile, "%.4lf\n", grains[x][y]);
    }
  }

  fclose(outfile);
}

void write_DEM() {
  int i;
  char filename[25];
  double N0, N1, N2, N3, N4, N5;  // Percentage of particles in contact
  double xgrainmax;
  FILE* outfile;
  // Output file
  // sprintf(filename,"DEM_Grains%.6i.dat",nFile);
  sprintf(filename, "DEM%.6i.dat", nFile);
  outfile = fopen(filename, "w");
  xfront = g[0].x1 + g[0].r;
  height = g[0].x2 + g[0].r;
  energie_cin = 0.;
  energie_x = 0.;
  energie_y = 0.;
  energie_teta = 0.;
  energy_p = 0.;
  energy_EPE = 0.;
  SE = 0.;
  ESE = 0.;
  WF = 0.;
  IFR = 0.;
  INCE = 0.;
  TSLIP = 0.;
  TRW = 0.;
  zmean = 0;
  xgrainmax = g[0].x1;
  N0 = 0;
  N1 = 0;
  N2 = 0;
  N3 = 0;
  N4 = 0;
  N5 = 0;
  for (i = 0; i < nbgrains; i++) {
    zmean += g[i].z;
    if (g[i].z == 0) N0 += 1;
    if (g[i].z == 1) N1 += 1;
    if (g[i].z == 2) N2 += 1;
    if (g[i].z == 3) N3 += 1;
    if (g[i].z == 4) N4 += 1;
    if (g[i].z == 5) N5 += 1;
    energie_x += 0.5 * g[i].m * g[i].v1 * g[i].v1;
    energie_y += 0.5 * g[i].m * g[i].v2 * g[i].v2;
    energie_teta += 0.5 * g[i].It * g[i].v3 * g[i].v3;
    energy_p += g[i].m * G * g[i].x2;
    /* if (nbsteps*dt>=dtt) */
    SE += 0.5 * (((g[i].p * g[i].p) / kg) + ((g[i].s * g[i].s) / kt));
    WF += g[i].fr;
    g[i].ifr =
        fabs(((g[i].m * G + g[i].f2) * (dt * g[i].v2 + dt2 * g[i].a2 / 2.)) +
             (g[i].f1 * (dt * g[i].v1 + dt2 * g[i].a1 / 2.)));
    //		g[i].ifr=(g[i].f2*(dt*g[i].v2))+(g[i].f1*(dt*g[i].v1));
    IFR += g[i].ifr;
    TSLIP += g[i].slip;
    TRW += g[i].rw;
    INCE += g[i].ice;
    TBW += g[i].ifr;
    ESE = 0.5 * (((g[i].p * g[i].p) / kg) + ((g[i].s * g[i].s) / kt));
    TSE += ESE;
    // if(g[i].x2>(-1.632*g[i].x1+0.0408)){energy_EPE+=g[i].m*G*g[i].x2;}
    if (g[i].x1 + g[i].r > xgrainmax) {
      xgrainmax = g[i].x1 + g[i].r;
    }
    if (g[i].x2 + g[i].r > height) {
      height = g[i].x2 + g[i].r;
    }
    if (g[i].zz > 0 && g[i].x1 + g[i].r >= xfront) {
      xfront = g[i].x1 + g[i].r;
    }
    if (g[i].z == 0) {
      g[i].fm = 0;
    } else
      g[i].fm = g[i].ifm / g[i].z;
    fprintf(outfile,
            "%i\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%"
            "le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%"
            "le\t%i\n",
            i, g[i].r, g[i].x1, g[i].x2, g[i].x3, g[i].v1, g[i].v2, g[i].v3,
            g[i].a1, g[i].a2, g[i].a3, fhf1[i], fhf2[i], fhf3[i], g[i].p,
            g[i].s, ESE, g[i].fr, g[i].ifr, g[i].ice, g[i].slip, g[i].rw,
            g[i].fm, g[i].M11, g[i].M12, g[i].M21, g[i].M22, g[i].z);
  }

  energie_cin = energie_x + energie_y + energie_teta;
  zmean = zmean / nbgrains;

  s_stats = fopen("stats.data", "a");

  fprintf(s_stats,
          "%le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le "
          "%le %le %le %le %le\n",
          nbsteps * dt - dtt, xfront, xgrainmax, height, zmean, energie_x,
          energie_y, energie_teta, energie_cin, N0 / nbgrains, N1 / nbgrains,
          N2 / nbgrains, N3 / nbgrains, N4 / nbgrains, N5 / nbgrains, energy_p,
          SE, WF, IFR, INCE, TSLIP, TRW);

  fclose(s_stats);
  fclose(outfile);
}

void write_forces() {
  int i, j;
  double dn;
  char nomfile[25];
  FILE* outfile1;
  // Ouverture du fichier
  // sprintf(filename,"DEM_Grains%.6i.dat",nFile);
  sprintf(nomfile, "DEM%.6i.ps", nFile);
  outfile1 = fopen(nomfile, "w");
  float margin = 10 * g[0].r, hrx1 = 3000, hry2 = 500;
  fprintf(outfile1, "%%!PS-Adobe-3.0 EPSF-3.0 \n");
  fprintf(outfile1, "%%BoundingBox: %f %f %f %f \n", -margin, -margin,
          hrx1 + margin, hry2 + margin);
  fprintf(outfile1, "%%Creator: Krishna Kumar \n");
  fprintf(outfile1, "%%Title: DEM Grains & Forces \n");
  fprintf(outfile1, "0.1 setlinewidth 0.0 setgray \n");
  for (i = 0; i <= nbgrains; i++)
    fprintf(outfile1,
            "newpath %le %le %le 0.0 setlinewidth %.2f setgray 0 360 arc gsave "
            "fill grestore\n",
            g[i].x1 * 10000, g[i].x2 * 10000, g[i].r * 10000,
            (0.8 - g[i].fm / 2));
  for (i = 0; i <= nbgrains; i++) {
    for (j = 0; j <= nbgrains; j++) {
      dn = (sqrt((g[i].x1 - g[j].x1) * (g[i].x1 - g[j].x1) +
                 (g[i].x2 - g[j].x2) * (g[i].x2 - g[j].x2))) -
           g[i].r - g[j].r;
      if (dn < -1e-10 && i != j) {
        //  printf("dn for i %i and j %i are: %le \n",i,j,dn);
        fprintf(outfile1, "%le setlinewidth \n 0.0 setgray \n", c1[i][j].fn);
        fprintf(outfile1, "1 setlinecap \n newpath \n");
        fprintf(outfile1, "%le %le moveto \n %le %le lineto\n", g[i].x1 * 10000,
                g[i].x2 * 10000, g[j].x1 * 10000, g[j].x2 * 10000);
        fprintf(outfile1, "stroke \n");
      }
    }
  }
  fclose(outfile1);
}

// --------------------------

void write_densities() {

  int x, y, i;
  double pasxyz;
  double P, u_x, u_y;

  char filename[25];
  char filename_press[25];
  FILE* outfile;
  FILE* s_press;

  sprintf(filename, "densities%.6i.vtk", nFile);
  sprintf(filename_press, "pressure_base%.6i.dat", nFile);
  pasxyz = 1. / lx;
  outfile = fopen(filename, "w");
  s_press = fopen(filename_press, "w");
  fprintf(outfile, "# vtk DataFile Version 2.0\n");
  fprintf(outfile, "Outfile domain LB t: %e\n", t);
  fprintf(outfile, "ASCII\n");
  fprintf(outfile, "DATASET RECTILINEAR_GRID\n");
  fprintf(outfile, "DIMENSIONS %d %d 1\n", lx, ly);
  fprintf(outfile, "X_COORDINATES %d float\n", lx);
  for (i = 0; i <= lx - 1; i++) {
    fprintf(outfile, "%e ", (float)i * pasxyz);
  }
  fprintf(outfile, "\n");
  fprintf(outfile, "Y_COORDINATES %d float\n", ly);
  for (i = 0; i <= ly - 1; i++) {
    fprintf(outfile, "%e ", (float)i * pasxyz);
  }
  fprintf(outfile, "\n");

  fprintf(outfile, "Z_COORDINATES 1 float\n");
  fprintf(outfile, "0\n");

  // Pour LB
  fprintf(outfile, "POINT_DATA %d\n", lx * ly);
  fprintf(outfile, "SCALARS Pressure float 1\n");
  fprintf(outfile, "LOOKUP_TABLE default\n");

  for (y = 0; y < ly; y++) {
    for (x = 0; x < lx; x++) {
      P = 0.;
      for (i = 0; i < Q; i++) {
        P += f[x][y][i];
      }
      P = (1. / 3.) * rho_moy * (P - 1.);
      if (obst[x][y] < 0) {
        fprintf(outfile, "%.4lf\n", P);
        if (y == 2) {
          fprintf(s_press, "%le %le\n", x * pasxyz, P);
        }
      } else {
        fprintf(outfile, "%.4lf\n", 0.);
        if (y == 2) {
          fprintf(s_press, "%le %le\n", x * pasxyz, 0.0);
        }
      }
    }
  }

  fprintf(outfile, "VECTORS VecVelocity float\n");

  for (y = 0; y < ly; y++) {
    for (x = 0; x < lx; x++) {
      // P=rho_moy;
      u_x = 0.;
      u_y = 0.;
      for (i = 0; i < Q; i++) {
        u_x += f[x][y][i] * ex[i];
        u_y += f[x][y][i] * ey[i];
      }
      // P = (P-rho_moy)*1./3.;
      // P = (1./3.)*rho_moy*(P-1.);
      if (obst[x][y] < 0) {
        fprintf(outfile, "%.4lf %.4lf 0.\n", u_x, u_y);
      } else {
        fprintf(outfile, "%.4lf %.4lf 0.\n", 0., 0.);
      }
    }
  }

  fclose(s_press);
  fclose(outfile);
}

// *******************************************************************
// *   sample initial                                           *
// *******************************************************************
void temp_sample() {
  long i, j, k;
  j = 0;
  k = 0;
  for (i = 0; i < nbgrains; ++i) {
    g[i].r = r;  //*(float)(i+1)/nbgrains;
    g[i].m = rhoS * pi * g[i].r * g[i].r;
#ifdef _FLUIDE_
    g[i].mw = rhoW * pi * g[i].r * g[i].r;
#else
    g[i].mw = 0;
#endif
    // g[i].m=(4./3.)*rhoS*pi*g[i].r*g[i].r*g[i].r;
    g[i].It = g[i].m * g[i].r * g[i].r / 2;
    // g[i].It=(2./5.)*g[i].m*g[i].r*g[i].r;
    g[i].x1 = r * 1.5 + 2. * r * j;  // Pente r*(1.5+(k/10.))
    g[i].x2 = r + 2. * r * k;
    g[i].x3 = 0.;
    g[i].v1 = 0.;
    g[i].v2 = 0.;
    g[i].v3 = 0.;
    g[i].a1 = 0.;
    g[i].a2 = 0.;
    g[i].a3 = 0.;
    // if(j<=4) {j++;} else {j=0;k++;};
    if (j < nbgrains_bas - 1 && k == 0) {
      j++;
    } else {
      if (j <= 13) {
        j++;
      } else {
        j = 0;
        k++;
      };
    };
  }
}

void sample() {
  int i;
  FILE* sample_file;
  char com[256];
  double L0, H0, MassGrain, xMax, xMin, yMax, yMin;

  sample_file = fopen(filename_sample, "r");

  fgets(com, 256, sample_file);
  printf("%s\n", com);
  fscanf(sample_file, "%d\n", &nbgrains);
  printf("Nb grains %d\n", nbgrains);
  for (i = 0; i < nbgrains; ++i) {
    fscanf(sample_file, "%le %le %le;\n", &g[i].r, &g[i].x1, &g[i].x2);
    // printf("%le %le %le\n",g[i].r,g[i].x1,g[i].x2);
    g[i].r = g[i].r * r;
    g[i].m = rhoS * pi * g[i].r * g[i].r;
    g[i].It = g[i].m * g[i].r * g[i].r / 2;
    g[i].x1 = g[i].x1 * r;
    g[i].x2 = g[i].x2 * r;
    g[i].x3 = 0.;
    g[i].v1 = 0.;
    g[i].v2 = 0.;
    g[i].v3 = 0.;
    g[i].a1 = 0.;
    g[i].a2 = 0.;
    g[i].a3 = 0.;
  }

  xMax = g[0].x1;
  xMin = g[0].x1;
  yMax = g[0].x2;
  yMin = g[0].x2;
  MassGrain = 0.;

  for (i = 0; i < nbgrains; ++i) {
    MassGrain += g[i].m;
    xMax = Max(xMax, g[i].x1 + g[i].r);
    xMin = Min(xMin, g[i].x1 - g[i].r);
    yMax = Max(yMax, g[i].x2 + g[i].r);
    yMin = Min(yMin, g[i].x2 - g[i].r);
  }
  L0 = xMax - xMin;
  H0 = yMax - yMin;
  printf("L0=%le H0=%le Mass of Grains=%le Phi=%le\n", L0, H0, MassGrain,
         MassGrain / (rhoS * (L0 * H0)));
}

// *******************************************************************************************
// *   Initialise obstacle array *
// *******************************************************************************************
void init_obst() {
  int x, y, i, xi, yi, xf, yf;
  // c.d.g. sphere
  //  double xc,yc;
  double dist2, r2, R2, xc, yc, rbl0;
  //	l2=lx*lx;
  //#pragma omp parallel for
  for (y = 1; y < ly - 1; y++) {
    for (x = 1; x < lx - 1; x++) {
      obst[x][y] = -1;
      grains[x][y] = -1;  // JYD
    }
  }

  //#pragma omp parallel for
  for (x = 0; x < lx; x++) {
    obst[x][0] = obst[x][ly - 1] = nbgrains;
    act[x][0] = act[x][ly - 1] = 0;
  }

  //#pragma omp parallel for
  for (y = 1; y < ly - 1; y++) {
    obst[0][y] = obst[lx - 1][y] = nbgrains;
    act[0][y] = act[lx - 1][y] = 0;
  }

#pragma omp parallel for
  for (i = 0; i < nbgrains; i++) {
    xc = (g[i].x1 - Mgx) / dx;
    yc = (g[i].x2 - Mby) / dx;
    r2 = rLB[i] * rLB[i];
    // Unreduced grain radius
    rbl0 = g[i].r / dx;
    R2 = rbl0 * rbl0;

    xi = (int)(xc - rbl0);
    xf = (int)(xc + rbl0);
    if (xi < 1) xi = 1;
    if (xf >= lx - 1) xf = lx - 2;
    yi = (int)(yc - rbl0);
    yf = (int)(yc + rbl0);
    if (yi < 1) yi = 1;
    if (yf >= ly - 1) yf = ly - 2;



    for (y = yi; y <= yf; y++) {
      for (x = xi; x <= xf; x++) {
        dist2 = (x - xc) * (x - xc) + (y - yc) * (y - yc);
        if (dist2 <= R2) {
          grains[x][y] = g[i].p;
          if (dist2 <= r2) obst[x][y] = i;
        }
        // else {
        //	grains[x][y]=-1;
        //}
      }
    }
  }
}
// *******************************************************************************************
// *   Initialise density distribution function with equilibrium to zero density
// *
// *******************************************************************************************
void init_density() {
  int x, y, iLB;
  // double u_squ,eu;
  for (x = 0; x < lx; x++) {
    for (y = 0; y < ly; y++) {
      for (iLB = 0; iLB < Q; iLB++) {
        f[x][y][iLB] = w[iLB];
      }  //*rho_moy; } sous hyp. rho=1
         /*		u_squ=u_inlet*u_inlet;
          for (i=0; i<Q;i++)
          {
          eu=ex[i]*u_inlet;
          f[x][y][i]=w[i]*rho_moy*(1.+3*eu+4.5*eu*eu-1.5*u_squ);
          }
          */
    }
  }
}

// *******************************************************************
// *   Calculate the forces on grains                                *
// *******************************************************************
struct force force_grains(long i, long j) {
  //   distance normale
  double dn, xOiOj, yOiOj, OiOj;
  double xn, yn;
  double vn, vxOiOj, vyOiOj, vt;
  double ftest;
  struct force f;

  // distance relative
  xOiOj = g[i].x1 - g[j].x1;
  yOiOj = g[i].x2 - g[j].x2;
  OiOj = sqrt(xOiOj * xOiOj + yOiOj * yOiOj);
  dn = OiOj - g[i].r - g[j].r;
  // calculate the forces
  if (dn >= 0) {
    f.f1 = 0;
    f.f2 = 0;
    f.f3 = 0;
  } else {
    // relative normal velocity
    vxOiOj = g[i].v1 - g[j].v1;
    vyOiOj = g[i].v2 - g[j].v2;
    xn = xOiOj / OiOj;
    yn = yOiOj / OiOj;

    // Compute velocities at contact
    vn = vxOiOj * xn + vyOiOj * yn;
    // Tangential Velocity
    vt = -vxOiOj * yn + vyOiOj * xn - g[i].v3 * g[i].r - g[j].v3 * g[j].r;

    // calculate normal force
    c1[i][j].fn = -kg * dn - nug * vn;

    if (c1[i][j].fn < 0) c1[i][j].fn = 0.0;
    c1[i][j].ft = c1[i][j].ft - kt * vt * dt;
    ftest = mu * c1[i][j].fn;
    if (fabs(c1[i][j].ft) > ftest) {
      if (c1[i][j].ft < 0.0)
        c1[i][j].ft = ftest;
      else
        c1[i][j].ft = -ftest;
    }
    f.f1 = c1[i][j].fn * xn - c1[i][j].ft * yn;
    f.f2 = c1[i][j].fn * yn + c1[i][j].ft * xn;
    f.f3 = -Maxt(c1[i][j].ft * g[i].r, c1[i][j].fn * murf * g[i].r * g[j].r);

    g[i].p += c1[i][j].fn;
    g[j].p += c1[i][j].fn;
    g[i].f1 += f.f1;
    g[i].f2 += f.f2;
    g[i].s += c1[i][j].ft;
    g[j].s += c1[i][j].ft;
    g[i].slip +=
        fabs(c1[i][j].ft) * (fabs(vt * dt) + (fabs(c1[i][j].ft - pft)) / kt);
    pft = c1[i][j].ft;
    g[i].rw += fabs(f.f3) * (fabs(g[i].v3 * dt) + (fabs(f.f3 - pff)) / kt);
    pff = f.f3;
    g[i].z += 1;
    g[i].zz += 1;
    g[i].ice += ic;
    if (c1[i][j].fn == 0)
      g[i].ifm = 0;
    else
      g[i].ifm += fabs(c1[i][j].ft / (mu * c1[i][j].fn));

    // Stress computations
    g[i].M11 += f.f1 * xOiOj;
    g[i].M12 += f.f1 * yOiOj;
    g[i].M21 += f.f2 * xOiOj;
    g[i].M22 += f.f2 * yOiOj;
  }

  return f;
}

// *******************************************************************
// *   Calculation of forces between the grains and Walls            *
// *******************************************************************

struct force force_WallB(long i, double dn) {
  double vn, vt, ftest;
  struct force f;
  vn = g[i].v2;
  vt = g[i].v1;
  c2[i].fn = -km * dn - num * vn;
  if (c2[i].fn < 0) c2[i].fn = 0.;
  c2[i].ft += ktm * vt;  //*dt; //Krishna
  ftest = mumb * c2[i].fn;
  if (fabs(c2[i].ft) > ftest) {
    if (c2[i].ft < 0.0)
      c2[i].ft = ftest;
    else
      c2[i].ft = -ftest;
  }

  f.f1 = c2[i].ft;
  f.f2 = c2[i].fn;
  f.f3 = -(c2[i].ft * g[i].r * murf);

  g[i].p += c2[i].fn;
  g[i].s += c2[i].ft;
  g[i].f1 += f.f1;
  g[i].z += 1;
  // Stress computations
  g[i].M11 += 0;
  g[i].M12 += f.f1 * dt;
  g[i].M21 += 0;
  g[i].M22 += f.f2 * dt;

  g[i].rw += fabs(f.f3) * (fabs(g[i].v3 * dt) + (fabs(f.f3 - pff)) / kt);
  g[i].fr += fabs(c2[i].ft) * (fabs(vt * dt) + (fabs(c2[i].ft - pft)) / kt);
  pff = f.f3;
  pft = c2[i].ft;
  return f;
}
struct force force_WallT(long i, double dn) {
  double vn, vt, fn, ft, ftmax;
  struct force f;
  vn = g[i].v2;
  fn = km * dn - num * vn;
  ic += num * vn * vn * dt;
  if (fn > 0.) fn = 0.;
  // relative tangential velocity
  vt = g[i].v1 + g[i].v3 * g[i].r -
       amp * freq * cos(freq * t);
  ft = fabs(ktm * vt);
  
  //	if (nbsteps*dt<dtt){mumb=mug;nugt=0;}
  if (vt >= 0) {
    ftmax = mumb * fn - nugt * vt;
  }  // ic+=nugt*vt*vt*dt;}
  else {
    ftmax = mumb * fn + nugt * vt;
  }  // ic+=nug*vt*vt*dt;}
  // ftmax=mum*fn-num*vt;
  if (ft > ftmax) ft = ftmax;
  if (vt > 0) ft = -ft;
  f.f1 = ft;
  f.f2 = fn;
  // f.f3=ft*g[i].r-fabs(murf*g[i].r*g[i].v3*fn);
  f.f3 = ft * g[i].r * murf;
  // f.f3=(ft-fabs(vt*nuf))*g[i].r;

  // Stress computations
  g[i].M11 += 0;
  g[i].M12 += f.f1 * fabs(dt);
  g[i].M21 += 0;
  g[i].M22 += f.f2 * fabs(dt);

  g[i].p += fn;
  g[i].s += ft;
  g[i].z += 1;
  // g[i].rw+=fabs(f.f3)*(fabs(g[i].v3*dt)+(fabs(f.f3-pff))/kt);
  // g[i].fr+=fabs(ft)*(fabs(vt*dt)+(fabs(ft-pft))/kt);
  //	pff=f.f3;pft=ft;
  return f;
}
struct force force_WallL(long i, double dn) {
  double vn, fn, vt, ft;
  struct force f;
  vn = g[i].v1;
  fn = -km * dn + num * vn;
  ic += num * vn * vn * dt;
  if (fn < 0.) fn = 0.;
  vt = g[i].v2;
  if (vt > 0)
    ft = mum * fn;
  else
    ft = mum * fn;
  if (vt > 0) ft = -ft;
  f.f1 = fn;
  f.f2 = ft;
  // f.f2=ft*g[i].r-fabs(murf*g[i].r*g[i].v3*fn);
  f.f3 = ft * g[i].r * murf;
  // f.f3=(ft-fabs(vt*nuf))*g[i].r;
  // Stress computations
  g[i].M11 += f.f1 * fabs(dt);
  g[i].M12 += 0;
  g[i].M21 += f.f2 * fabs(dt);
  g[i].M22 += 0;

  g[i].p += fn;
  g[i].s += ft;
  g[i].f1 += f.f1;
  g[i].z += 1;
  g[i].ice += ic;
  g[i].rw += fabs(f.f3) * fabs(g[i].v3 * dt);
  g[i].fr += fabs(ft) * (fabs(vt * dt) + (fabs(ft - pft)) / kt);
  pft = ft;
  return f;
}

struct force force_WallR(long i, double dn) {
  double vn, fn, vt, ft;
  struct force f;
  vn = g[i].v1;
  fn = km * dn - num * vn;
  //	ic+=num*vn*vn*dt;
  vt = g[i].v2;  // tangential velcoty
  ft = mum * fn;  // ic+=nugt*vt*vt*dt;
  // ftmax=mum*fn-num*vt;
  if (vt > 0) ft = -ft;
  if (fn > 0.) fn = 0.;
  f.f1 = fn;
  f.f2 = -ft;
  f.f3 = ft * g[i].r * murf;

  g[i].p += fn;
  g[i].f1 += f.f1;
  //	g[i].ice +=ic;
  //	g[i].fr+=fabs(ft)*(fabs(vt*dt)+(fabs(f.f3-pft))/kt);
  pft = ft;
  // Stress computations
  g[i].M11 += f.f1 * fabs(dt);
  g[i].M12 += 0;
  g[i].M21 += f.f2 * fabs(dt);
  g[i].M22 += 0;
  //	g[i].s += ft;
  g[i].z += 1;
  return f;
}
// *******************************************************************
// *                                                                 *
// *                                                                 *
// *                                                                 *
// *   Calculate the hydrodynamic forces                             *
// *                                                                 *
// *                                                                 *
// *                                                                 *
// *                                                                 *
// *******************************************************************
// *******************************************************************************************
// *   Reinitialise density distributions for nodes that change state from solid
// to fluid   *
// *******************************************************************************************
void reinit_obst_density() {
  int x, y, iLB;
  double u_squ, eu;

  for (x = 1; x < lx - 1; x++) {
    // #pragma omp parallel for
    for (y = 1; y < ly - 1; y++) {
      if (obst[x][y] != -1) { 
       // Usqu: the standard (^2) the speed of the node from the portion
       // solid to the fluid portion
        u_squ = ((g[obst[x][y]].v1 -
                  (y * dx + Mby - g[obst[x][y]].x2) * g[obst[x][y]].v3) *
                     (g[obst[x][y]].v1 -
                      (y * dx + Mby - g[obst[x][y]].x2) * g[obst[x][y]].v3) +
                 (g[obst[x][y]].v2 +
                  (x * dx + Mgx - g[obst[x][y]].x1) * g[obst[x][y]].v3) *
                     (g[obst[x][y]].v2 +
                      (x * dx + Mgx - g[obst[x][y]].x1) * g[obst[x][y]].v3)) /
                (c * c);
        for (iLB = 0; iLB < Q; iLB++) {
          // eu : e.u in formula feq

          eu = (ex[iLB] *
                    (g[obst[x][y]].v1 -
                     (y * dx + Mby - g[obst[x][y]].x2) * g[obst[x][y]].v3) +
                ey[iLB] *
                    (g[obst[x][y]].v2 +
                     (x * dx + Mgx - g[obst[x][y]].x1) * g[obst[x][y]].v3)) /
               c;
          f[x][y][iLB] =
              w[iLB] * (1. + 3 * eu + 4.5 * eu * eu - 1.5 * u_squ);  //*rho_moy;
        }
      }
    }
  }
}

// *******************************************************************************************
// *   Obstacle array construction && nodes activity *
// *******************************************************************************************
void obst_construction() {
  int x, y, xp, next_x, next_y, i, iLB, xi, yi, xf, yf;
  // c.d.g. sphere
  //  double xc,yc;
  double dist2, aa, bb, cc, r2, xc, yc, R2, rbl0;
//	l2=lx*lx;
#pragma omp parallel for
  for (y = 1; y < ly - 1; y++) {
    //	    #pragma omp parallel for
    for (x = 1; x < lx - 1; x++) {
      obst[x][y] = -1;
      grains[x][y] = -1;  // JYD
      act[x][y] = 1;
      for (iLB = 1; iLB < Q; iLB++) {
        delta[x][y][iLB] = 0;
      }
    }
  }
//	#pragma acc kernels
#pragma omp parallel for
  for (i = 0; i < nbgrains; i++) {
    xc = (g[i].x1 - Mgx) / dx;
    yc = (g[i].x2 - Mby) / dx;
    r2 = rLB[i] * rLB[i];
    rbl0 = g[i].r / dx;  // JYD2
    R2 = rbl0 * rbl0;
    // xi=xc-rLB[i]; xf=xc+rLB[i]; if(xi<1) xi=1; if(xf>=lx-1) xf=lx-2;
    // yi=yc-rLB[i]; yf=yc+rLB[i]; if(yi<1) yi=1; if(yf>=ly-1) yf=ly-2;
    xi = (int)(xc - rbl0);
    xf = (int)(xc + rbl0);
    if (xi < 1) xi = 1;
    if (xf >= lx - 1) xf = lx - 2;
    yi = (int)(yc - rbl0);
    yf = (int)(yc + rbl0);
    if (yi < 1) yi = 1;
    if (yf >= ly - 1) yf = ly - 2;
    //		#pragma loop for independent
    for (y = yi; y <= yf; y++) {
      for (x = xi; x <= xf; x++) {
        dist2 = (x - xc) * (x - xc) + (y - yc) * (y - yc);
        if (dist2 <= R2) {
          grains[x][y] = g[i].p;
          // printf("%le\n",grains[x][y]);
          if (dist2 <= r2) obst[x][y] = i;
        }
        // if(dist2<=r2) obst[x][y]=i;
      }
    }

    // *   Obstacle in inteaction with fluid (active obstacles)

    //	#pragma loop for independent
    for (y = yi; y <= yf; y++) {
      for (x = xi; x <= xf; x++) {
        if (obst[x][y] == i) {
          act[x][y] = 0;

          // Search fluid node neighbourss
          for (iLB = 1; iLB < Q; iLB++) {
            next_x = x + ex[iLB];  // if (next_x<0) next_x=0; if (next_x>=lx)
                                   // next_x=lx-1;
            next_y = y + ey[iLB];  // if (next_y<0) next_y=0; if (next_y>=ly)
                                   // next_y=ly-1;
            if (obst[next_x][next_y] == -1) {

              // Calculating the distance between the node fluid and the wall of the particle
              // (Klaus-Nils-Ulrich)

              act[x][y] = 1;
              xp = x;
              aa = fabs(ex[iLB]) + fabs(ey[iLB]);
              bb = (xp + ex[iLB] - xc) * ex[iLB] + (y + ey[iLB] - yc) * ey[iLB];
              cc = (xp + ex[iLB] - xc) * (xp + ex[iLB] - xc) +
                   (y + ey[iLB] - yc) * (y + ey[iLB] - yc) - r2;
              delta[x][y][iLB] = (bb - sqrt(fabs(bb * bb - aa * cc))) / aa;
            }
          }
        }
      }
    }
  }
}

// ************************************************************************************************
// *   Principal LB: Collision - (Streaming + Boundary Conditions) *
// ************************************************************************************************

void collision_streaming() {
  int x, y, iLB, next_x, next_y, next_xx, next_yy;
  double rho, e, eps, j_x, j_y, q_x, q_y, p_xx, p_xy, eO, epsO, q_xO, q_yO,
      p_xxO, p_xyO, j_x2, j_y2;
  // double s2=1.944,s3=1.944,s5=1.944,s7=1.944,s8=1.944,s9=1.944;
  const int half = (Q - 1) / 2;
  double a = 1. / 36;

// Post-collision part computation
// (Yu-Mei-Luo-Shyy)
#pragma acc kernels
  for (x = 1; x < lx - 1; x++) {
    for (y = 1; y < ly - 1; y++) {
      if (obst[x][y] == -1) {

        rho = f[x][y][0] + f[x][y][1] + f[x][y][2] + f[x][y][3] + f[x][y][4] +
              f[x][y][5] + f[x][y][6] + f[x][y][7] + f[x][y][8];

        e = -4 * f[x][y][0] + 2 * f[x][y][1] - f[x][y][2] + 2 * f[x][y][3] -
            f[x][y][4] + 2 * f[x][y][5] - f[x][y][6] + 2 * f[x][y][7] -
            f[x][y][8];
        eps = 4 * f[x][y][0] + f[x][y][1] - 2 * f[x][y][2] + f[x][y][3] -
              2 * f[x][y][4] + f[x][y][5] - 2 * f[x][y][6] + f[x][y][7] -
              2 * f[x][y][8];
        j_x = f[x][y][5] + f[x][y][6] + f[x][y][7] - f[x][y][1] - f[x][y][2] -
              f[x][y][3];
        q_x = -f[x][y][1] + 2 * f[x][y][2] - f[x][y][3] + f[x][y][5] -
              2 * f[x][y][6] + f[x][y][7];
        j_y = f[x][y][1] + f[x][y][8] + f[x][y][7] - f[x][y][3] - f[x][y][4] -
              f[x][y][5];
        q_y = f[x][y][1] - f[x][y][3] + 2 * f[x][y][4] - f[x][y][5] +
              f[x][y][7] - 2 * f[x][y][8];
        p_xx = f[x][y][2] - f[x][y][4] + f[x][y][6] - f[x][y][8];
        p_xy = -f[x][y][1] + f[x][y][3] - f[x][y][5] + f[x][y][7];

        j_x2 = j_x * j_x;
        j_y2 = j_y * j_y;

        eO = e - s2 * (e + 2 * rho - 3 * (j_x2 + j_y2) / rho);
        epsO = eps - s3 * (eps - rho + 3 * (j_x2 + j_y2) / rho);
        q_xO = q_x - s5 * (q_x + j_x);
        q_yO = q_y - s7 * (q_y + j_y);
        p_xxO = p_xx - s8 * (p_xx - (j_x2 - j_y2) / rho);
        p_xyO = p_xy - s9 * (p_xy - j_x * j_y / rho);

        //! à revoir. Consulter : ../général/Yu-Mei-Luo-Shyy

        f[x][y][0] = 4 * a * (rho - eO + epsO);
        f[x][y][1] = a * (4 * rho + 2 * eO + epsO - 6 * j_x - 3 * q_xO +
                          6 * j_y + 3 * q_yO - 9 * p_xyO);
        f[x][y][2] =
            a * (4 * rho - eO - 2 * epsO - 6 * j_x + 6 * q_xO + 9 * p_xxO);
        f[x][y][3] = a * (4 * rho + 2 * eO + epsO - 6 * j_x - 3 * q_xO -
                          6 * j_y - 3 * q_yO + 9 * p_xyO);
        f[x][y][4] =
            a * (4 * rho - eO - 2 * epsO - 6 * j_y + 6 * q_yO - 9 * p_xxO);
        f[x][y][5] = a * (4 * rho + 2 * eO + epsO + 6 * j_x + 3 * q_xO -
                          6 * j_y - 3 * q_yO - 9 * p_xyO);
        f[x][y][6] =
            a * (4 * rho - eO - 2 * epsO + 6 * j_x - 6 * q_xO + 9 * p_xxO);
        f[x][y][7] = a * (4 * rho + 2 * eO + epsO + 6 * j_x + 3 * q_xO +
                          6 * j_y + 3 * q_yO + 9 * p_xyO);
        f[x][y][8] =
            a * (4 * rho - eO - 2 * epsO + 6 * j_y - 6 * q_yO - 9 * p_xxO);
      }
    }
  }

  // To calculate the edges; see the book Lattice Boltzmann Modeling
  // Bounce back for y=0 & y=ly-1
  //	#pragma omp parallel for
  for (x = 1; x < lx - 1; x++) {
    f[x][0][8] = f[x][1][4];
    f[x][0][7] = f[x + 1][1][3];  //;+uw_b/6;
    f[x][0][1] = f[x - 1][1][5];  //-uw_b/6;
    // Top plate
    f[x][ly - 1][4] = f[x][ly - 2][8];
    f[x][ly - 1][3] = f[x - 1][ly - 2][7];  //-uw_h/6;
    f[x][ly - 1][5] = f[x + 1][ly - 2][1];  //+uw_h/6;
  }
  //#pragma omp parallel for
  for (y = 1; y < ly - 1; y++) {
    f[0][y][6] = f[1][y][2];
    f[0][y][7] = f[1][y + 1][3];  //;+uw_b/6;
    f[0][y][5] = f[1][y - 1][1];  //-uw_b/6;
    f[lx - 1][y][2] = f[lx - 2][y][6];
    f[lx - 1][y][3] = f[lx - 2][y - 1][7];  //-uw_h/6;
    f[lx - 1][y][1] = f[lx - 2][y + 1][5];  //+uw_h/6;
  }

  // corner nodes
  f[0][0][7] = f[1][1][3];
  f[lx - 1][0][1] = f[lx - 2][1][5];  //+uw_b/6
  f[0][ly - 1][5] = f[1][ly - 2][1];  //-uw_b/6
  f[lx - 1][ly - 1][3] = f[lx - 2][ly - 2][7];

   // bounce back in obstacles
  /////////////////////////////////////////////////////////
  //  To calculate force f[][][])                        //
  //  1: articlel of JYD-Mouloud                     //
  //  2: article of Klaus-Nils-Ulrich                //
  //  3: article of Yu-Mei-Luo-Shyy                  //
  /////////////////////////////////////////////////////////
  //	#pragma acc kernels
  //#pragma omp parallel for
  for (x = 1; x < lx - 1; x++) {
    for (y = 1; y < ly - 1; y++) {
      if (obst[x][y] != -1 && act[x][y] == 1) {
        for (iLB = 1; iLB <= half; iLB++) {
          next_x = x + ex[iLB];  // if (next_x<0) next_x=lx-1; if (next_x>=lx)
                                 // next_x=0;
          next_y = y + ey[iLB];
          if (obst[next_x][next_y] != -1)
            f[x][y][iLB] = w[iLB];
          else  //(obst[next_x][next_y]==-1)
          {
            // Calculation is based on JYD-Mouloud (2.3.3)
            if (delta[x][y][iLB] >= 0.5) {
              f[x][y][iLB] =
                  f[next_x][next_y][iLB + half] / (2 * delta[x][y][iLB]) +
                  (2 * delta[x][y][iLB] - 1) * f[next_x][next_y][iLB] /
                      (2 * delta[x][y][iLB]) +
                  3 * (w[iLB] / c) *
                      (ex[iLB] * (g[obst[x][y]].v1 -
                                  (y * dx + Mby - g[obst[x][y]].x2) *
                                      g[obst[x][y]].v3) +
                       ey[iLB] * (g[obst[x][y]].v2 +
                                  (x * dx + Mgx - g[obst[x][y]].x1) *
                                      g[obst[x][y]].v3)) /
                      delta[x][y][iLB];
            }
            if (delta[x][y][iLB] > 0. && delta[x][y][iLB] < 0.5) {
              next_xx = next_x + ex[iLB];  // if (next_xx<0) next_xx=lx-1; if
                                           // (next_xx>=lx) next_xx=0;
              next_yy = next_y + ey[iLB];
              f[x][y][iLB] =
                  2 * delta[x][y][iLB] * f[next_x][next_y][iLB + half] +
                  (1 - 2 * delta[x][y][iLB]) * f[next_xx][next_yy][iLB + half] +
                  6 * (w[iLB] / c) *
                      (ex[iLB] * (g[obst[x][y]].v1 -
                                  (y * dx + Mby - g[obst[x][y]].x2) *
                                      g[obst[x][y]].v3) +
                       ey[iLB] * (g[obst[x][y]].v2 +
                                  (x * dx + Mgx - g[obst[x][y]].x1) *
                                      g[obst[x][y]].v3));
            }
          }
        }

        for (iLB = 1 + half; iLB < Q; iLB++) {
          next_x = x + ex[iLB];  // if (next_x<0) next_x=lx-1; if (next_x>=lx)
                                 // next_x=0;
          next_y = y + ey[iLB];
          if (obst[next_x][next_y] != -1)
            f[x][y][iLB] = w[iLB];

          else  //(obst[next_x][next_y]==-1)
          {
            // Calculation is based on JYD-Mouloud (2.3.3) 
            if (delta[x][y][iLB] >= 0.5) {
              f[x][y][iLB] =
                  f[next_x][next_y][iLB - half] / (2 * delta[x][y][iLB]) +
                  (2 * delta[x][y][iLB] - 1) * f[next_x][next_y][iLB] /
                      (2 * delta[x][y][iLB]) +
                  3 * (w[iLB] / c) *
                      (ex[iLB] * (g[obst[x][y]].v1 -
                                  (y * dx + Mby - g[obst[x][y]].x2) *
                                      g[obst[x][y]].v3) +
                       ey[iLB] * (g[obst[x][y]].v2 +
                                  (x * dx + Mgx - g[obst[x][y]].x1) *
                                      g[obst[x][y]].v3)) /
                      delta[x][y][iLB];
            }
            if (delta[x][y][iLB] > 0. && delta[x][y][iLB] < 0.5) {
              next_xx = next_x + ex[iLB];  // if (next_xx<0) next_xx=lx-1; if
                                           // (next_xx>=lx) next_xx=0;
              next_yy = next_y + ey[iLB];
              f[x][y][iLB] =
                  2 * delta[x][y][iLB] * f[next_x][next_y][iLB - half] +
                  (1 - 2 * delta[x][y][iLB]) * f[next_xx][next_yy][iLB - half] +
                  6 * (w[iLB] / c) *
                      (ex[iLB] * (g[obst[x][y]].v1 -
                                  (y * dx + Mby - g[obst[x][y]].x2) *
                                      g[obst[x][y]].v3) +
                       ey[iLB] * (g[obst[x][y]].v2 +
                                  (x * dx + Mgx - g[obst[x][y]].x1) *
                                      g[obst[x][y]].v3));
            }
          }
        }
      }
    }
  }

  for (x = 0; x < lx; x++) {
    for (y = 0; y < ly; y++) {
      for (iLB = 1; iLB <= half; iLB++) {
        swap(&f[x][y][iLB], &f[x][y][iLB + half]);
      }
    }
  }

  //*****************************************************************************************

  // Stream by swapping; moving the densities according to their velocities

  for (x = 0; x < lx; x++) {
    for (y = 0; y < ly; y++) {
      for (iLB = 1; iLB <= half; iLB++) {
        next_x = x + ex[iLB];  // if(next_x<0) next_x=lx-1;
        next_y = y + ey[iLB];
        if (next_x >= 0 && next_y >= 0 && next_x < lx && next_y < ly) {
          swap(&f[x][y][iLB + half], &f[next_x][next_y][iLB]);
        }
      }
    }
  }
}

// ************************************************************************************************
// *   Compute total density to verify of the system  that no divergence occurs
// when iterating    *
// ************************************************************************************************
void check_density() {
  int x, y, iLB;
  double sum = 0;
  //#pragma omp parallel for
  for (x = 0; x < lx; x++) {
    for (y = 0; y < ly; y++) {
      for (iLB = 0; iLB < Q; iLB++) {
        sum = sum + f[x][y][iLB];
      }
    }
  }
  printf("Iteration Number %ld, Total density in the system %f\n", nbsteps,
         sum);
}
// ****************************************************************************
// *   Compute hydrodynamic forces                                          *
// ****************************************************************************
void forces_fluid() {
  int x, y, next_x, next_y, i, iLB;
  double fnx, fny;  //,sum_f1,sum_f2,sum_f3;
  //	double sum2_f1, sum2_f2,sum2_f3;
  int halfq;

  //	sum_f1=sum_f2=sum_f3=0.0;
  //	sum2_f1=sum2_f2=sum2_f3=0.0;

  const int half = (Q - 1) / 2;

#pragma omp parallel
  {
#pragma omp sections nowait
    {
#pragma omp section
      {

        reinit_obst_density();
        obst_construction();
        collision_streaming();
      }
    }
  }

  if (nbsteps % stepConsole == 0) check_density();
#pragma acc kernels copyout(                         \
    fhf1[0 : nbgrains],                              \
         fhf2[0 : nbgrains], fhf3[0 : nbgrains])     \
             copyin(obst[0 :][0 :], g[0 : nbgrains], \
                                      ey[0 :], f[0 :][0 :][0 :], ex[0 :])
  {
    for (i = 0; i < nbgrains; i++) {
      // sum_f1=sum_f2=sum_f3=0.0;
      // sum2_f1=sum2_f2=sum2_f3=0.0;
      fhf1[i] = fhf2[i] = fhf3[i] = 0.;
      //#pragma acc loop
      for (y = 0; y < ly; y++) {
        //	  #pragma acc loop
        for (x = 0; x < lx; x++) {
          if (obst[x][y] == i)  //&& act[x][y][z]==1
          {
// Search fluid nodes neighbors
#pragma acc for independent
            for (iLB = 1; iLB < Q; iLB++) {
              next_x = x + ex[iLB];
              next_y = y + ey[iLB];

              if (iLB <= half)
                halfq = half;
              else
                halfq = -half;
              if (obst[next_x][next_y] != i) {
                fnx = (f[x][y][iLB + halfq] + f[next_x][next_y][iLB]) *
                      ex[iLB + halfq];
                fny = (f[x][y][iLB + halfq] + f[next_x][next_y][iLB]) *
                      ey[iLB + halfq];
                fhf1[i] = fhf1[i] + fnx;
                fhf2[i] = fhf2[i] + fny;
                fhf3[i] = fhf3[i] - fnx * (y - (g[i].x2 - Mby) / dx) +
                          fny * (x - (g[i].x1 - Mgx) / dx);
              }
            }
          }
        }
      }

      fhf1[i] =
          fhf1[i] * rho_moy * 9 * nu * nu / (dx * (tau - 0.5) * (tau - 0.5));
      fhf2[i] =
          fhf2[i] * rho_moy * 9 * nu * nu / (dx * (tau - 0.5) * (tau - 0.5));
      fhf3[i] = fhf3[i] * dx * rho_moy * 9 * nu * nu /
                (dx * (tau - 0.5) * (tau - 0.5));
    }
  }
}

//**********************************************************
void acceleration_grains() {
  long i, j;
  int jdep;
  double dn, ftest;
  struct force fji; 
  if (nbsteps % stepFilm == 0 && start == 1) {  // Outfile MGPost
    //   distance normale
    double xOiOj, yOiOj, OiOj;
    double fn, xn, yn;
    double vn, vxOiOj, vyOiOj;
    double ft, ftmax, vt;

    for (i = 0; i <= nbgrains - 1; i++) {
      g[i].a1 = fhf1[i];
      g[i].a2 = fhf2[i];
      g[i].a3 = fhf3[i];
    }
    // Summation of forces on the grains
    for (i = 0; i <= nbgrains - 1; i++) {
       if (i == 0)
        jdep = 0;
      else
        jdep = cumul[i - 1];
      for (j = jdep; j < cumul[i]; j++) {
        //*
        // fji=force_grains(i,neighbours[j]);
        //*
        // forces_grains
        xOiOj = g[i].x1 - g[neighbours[j]].x1;
        yOiOj = g[i].x2 - g[neighbours[j]].x2;
        OiOj = sqrt(xOiOj * xOiOj + yOiOj * yOiOj);
        dn = OiOj - g[i].r - g[neighbours[j]].r;
        if (dn >= 0) {
          fji.f1 = 0;
          fji.f2 = 0;
          fji.f3 = 0;
        } else {
          // relative normal velocity
          vxOiOj = g[i].v1 - g[neighbours[j]].v1;
          vyOiOj = g[i].v2 - g[neighbours[j]].v2;
          xn = xOiOj / OiOj;
          yn = yOiOj / OiOj;
          vn = vxOiOj * xn + vyOiOj * yn;
          vt = -vxOiOj * yn + vyOiOj * xn - g[i].v3 * g[i].r -
               g[neighbours[j]].v3 * g[neighbours[j]].r;
          c1[i][j].fn = -kg * dn - nug * vn;
          if (c1[i][j].fn < 0) c1[i][j].fn = 0.0;
          c1[i][j].ft = c1[i][j].ft - kt * vt * dt;
          ftest = mu * c1[i][j].ft;
          if (fabs(c1[i][j].ft) > ftest) {
            if (c1[i][j].ft > 0.0)
              c1[i][j].ft = ftest;
            else
              c1[i][j].ft = -ftest;
          }
          //calculate the normal force
          fji.f1 = c1[i][j].fn * xn - c1[i][j].ft * yn;
          fji.f2 = c1[i][j].fn * yn + c1[i][j].ft * xn;
          fji.f3 = -c1[i][j].ft * g[i].r * murf;

          g[i].p += c1[i][j].fn;
          g[neighbours[j]].p += c1[i][j].fn;
          g[i].s += c1[i][j].ft;
          g[neighbours[j]].s += c1[i][j].ft;
          g[i].slip += fabs(c1[i][j].ft) *
                       (fabs(vt * dt) + (fabs(c1[i][j].ft - pft)) / kt);
          g[neighbours[j]].slip += fabs(c1[i][j].ft) *
                               (fabs(vt * dt) + (fabs(c1[i][j].ft - pft)) / kt);
          g[i].rw +=
              fabs(fji.f3) * (fabs(g[i].v3 * dt) + (fabs(fji.f3 - pff)) / kt);
          g[neighbours[j]].rw +=
              fabs(fji.f3) * (fabs(g[i].v3 * dt) + (fabs(fji.f3 - pff)) / kt);
          g[i].z += 1;
          pff = fji.f3;
          pft = c1[i][j].ft;
          // Stress computations
          g[i].M11 += fji.f1 * xOiOj;
          g[i].M12 += fji.f1 * yOiOj;
          g[i].M21 += fji.f2 * xOiOj;
          g[i].M22 += fji.f2 * yOiOj;
        }
        // end force_grains
        g[i].a1 = g[i].a1 + fji.f1;
        g[i].a2 = g[i].a2 + fji.f2;
        g[i].a3 = g[i].a3 + fji.f3;
        g[neighbours[j]].a1 = g[neighbours[j]].a1 - fji.f1;
        g[neighbours[j]].a2 = g[neighbours[j]].a2 - fji.f2;
        g[neighbours[j]].a3 = g[neighbours[j]].a3 + fji.f3;
      }
    }
  } else {  // Calculate normal
    for (i = 0; i <= nbgrains - 1; i++) {
      g[i].a1 = fhf1[i];
      g[i].a2 = fhf2[i];
      g[i].a3 = fhf3[i];
    }
    // summation of forces between the grains
    for (i = 0; i <= nbgrains - 1; i++) {
      //  printf("cumul(%d)= %d\n",i,cumul[i]);
      if (i == 0)
        jdep = 0;
      else
        jdep = cumul[i - 1];
      for (j = jdep; j < cumul[i]; j++) {
        // printf("grain(%d), neighbours(%d)= %d\n",i,i,neighbours[j]);
        fji = force_grains(i, neighbours[j]);
        g[i].a1 = g[i].a1 + fji.f1;
        g[i].a2 = g[i].a2 + fji.f2;
        g[i].a3 = g[i].a3 + fji.f3;
        g[neighbours[j]].a1 = g[neighbours[j]].a1 - fji.f1;
        g[neighbours[j]].a2 = g[neighbours[j]].a2 - fji.f2;
        g[neighbours[j]].a3 = g[neighbours[j]].a3 + fji.f3;
      }
    }
  }

  // Forces on the botton wall

  for (i = 0; i < nNeighWallb; i++) {
    dn = g[neighbourWallB[i]].x2 - g[neighbourWallB[i]].r - Mby;
    if (dn < 0) {
      fji = force_WallB(neighbourWallB[i], dn);
      g[neighbourWallB[i]].a1 = g[neighbourWallB[i]].a1 + fji.f1;
      g[neighbourWallB[i]].a2 = g[neighbourWallB[i]].a2 + fji.f2;
      g[neighbourWallB[i]].a3 = g[neighbourWallB[i]].a3 + fji.f3;
      g[neighbourWallB[i]].fr +=
          fabs(fji.f1) * (fabs(dt * g[i].v1) + fabs(dt2 * g[i].a1) +
                          (fabs(fji.f1 - pf)) /
                              kt);  // Friction work at the wall dt2*g[i].a1/2
      pf = fji.f1;  // Previous force fji.f1
    }
  }

  // Forces on le Top Wall

  for (i = 0; i < nNeighWallt; i++) {
    dn = -g[neighbourWallT[i]].x2 - g[neighbourWallT[i]].r + Mhy;
    if (dn < 0) {
      fji = force_WallT(neighbourWallT[i], dn);
      g[neighbourWallT[i]].a1 = g[neighbourWallT[i]].a1 + fji.f1;
      g[neighbourWallT[i]].a2 = g[neighbourWallT[i]].a2 + fji.f2;
      g[neighbourWallT[i]].a3 = g[neighbourWallT[i]].a3 + fji.f3;
    }
  }
  // Forces on le Wall left

  for (i = 0; i < nNeighWallL; i++) {
    dn = g[neighbourWallL[i]].x1 - g[neighbourWallL[i]].r - Mgx;
    if (dn < 0) {
      fji = force_WallL(neighbourWallL[i], dn);
      g[neighbourWallL[i]].a1 = g[neighbourWallL[i]].a1 + fji.f1;
      g[neighbourWallL[i]].a2 = g[neighbourWallL[i]].a2 + fji.f2;
      g[neighbourWallL[i]].a3 = g[neighbourWallL[i]].a3 + fji.f3;
      g[neighbourWallL[i]].fr +=
          fabs(fji.f2) *
          (fabs(dt * g[i].v1) + fabs(dt2 * g[i].a1) +
           (fabs(fji.f2 - pf)) / kt);  // Friction work at the wall
      pf = fji.f2;  // Previous force fji.f1
    }
  }

  // Forces on the right Wall

  for (i = 0; i < nNeighWallR; i++) {
    dn = -g[neighbourWallR[i]].x1 - g[neighbourWallR[i]].r + Mdx;
    if (dn < 0) {
      fji = force_WallR(neighbourWallR[i], dn);
      g[neighbourWallR[i]].a1 = g[neighbourWallR[i]].a1 + fji.f1;
      g[neighbourWallR[i]].a2 = g[neighbourWallR[i]].a2 + fji.f2;
      g[neighbourWallR[i]].a3 = g[neighbourWallR[i]].a3 + fji.f3;
    }
  }

  // calcul acc
  for (i = 0; i <= nbgrains - 1; i++) {
    // printf("moment (%d)= %f\n",i,g[i].a3,i,g[i].a3/g[i].It);
    g[i].a1 = g[i].a1 / g[i].m + ((g[i].m - g[i].mw) / g[i].m) * xG;
    g[i].a2 = (g[i].a2 / g[i].m) + ((g[i].m - g[i].mw) / g[i].m) * yG;
    g[i].a3 = g[i].a3 / g[i].It;
    // printf("acc(%d)= %f\n",i,g[i].a3);
  }
  // printf("acceleration grain (nbgrains-1) %f\n",g[nbgrains-1].a3);
}
//**********************************************************************

//-------------------------------------------------
void initVerlet() {
  int i, j;
  int jneighbours;
  double distx, disty;
  //	double distVerlet=.1e-7;
  jneighbours = 0;
  for (i = 0; i < nbgrains; i++) {
    for (j = i + 1; j < nbgrains; j++) {
      distx = g[i].x1 - g[j].x1;
      disty = g[i].x2 - g[j].x2;
      if (((fabs(distx) - g[i].r - g[j].r) <= distVerlet) &&
          ((fabs(disty) - g[i].r - g[j].r) <= distVerlet)) {
        if ((sqrt(distx * distx + disty * disty) - g[i].r - g[j].r) <=
            distVerlet) {
          neighbours[jneighbours] = j;
          jneighbours++;
          if (jneighbours == (nbgrains * 6 - 1))
            printf("error! size of vector verlet neighbors  is outdated");
        }
      }
      cumul[i] = jneighbours;
    }
    // printf("cumul(%d)= %d\n",i,cumul[i]);
  }
}

void VerletWall() {
  int i;
  nNeighWallb = 0;
  nNeighWallL = 0;
  nNeighWallt = 0;
  nNeighWallR = 0;
  double dn;
  // double distVerlet=.1e-7;
  // Verlet WallB

  if (nbsteps * dt < dtt) {
    Mdx = 1.e-3 * lx / 10;
    Mhy = (1.e-3 * ly / 10);
  } else {
    Mdx = 1.e-3 * lx;
    Mhy = 1.e-3 * ly;
  }

  for (i = 0; i < nbgrains; ++i) {
    dn = g[i].x2 - g[i].r - Mby;
    if (dn < distVerlet) {
      neighbourWallB[nNeighWallb] = i;
      ++nNeighWallb;
    }
  }
  // Verlet WallT
  for (i = 0; i < nbgrains; ++i) {
    dn = -g[i].x2 - g[i].r + Mhy;
    if (dn < distVerlet) {
      neighbourWallT[nNeighWallt] = i;
      ++nNeighWallt;
    }
  }
  // Verlet WallL
  for (i = 0; i < nbgrains; ++i) {
    dn = g[i].x1 - g[i].r - Mgx;
    if (dn < distVerlet) {
      neighbourWallL[nNeighWallL] = i;
      ++nNeighWallL;
    }
  }
  // Verlet WallR
  for (i = 0; i < nbgrains; ++i) {
    dn = -g[i].x1 - g[i].r + Mdx;
    if (dn < distVerlet) {
      neighbourWallR[nNeighWallR] = i;
      ++nNeighWallR;
    }
  }
}

//---------------------------------------------------
void WallPosition() {

  Mgx = 0.;
  Mdx = 1.e-3 * lx / 10;
  Mhy = 1.e-3 * ly / 10;
  Mby = 0.;

  //
  aff_transx = (Mdx - Mgx) * 0.25;  //*aff_zoom;
  aff_transy = (Mhy - Mby) * 0.5;   //*aff_zoom;

}

// *******************************************************************************************
// *   writing obstacle arrays 
// *
// *******************************************************************************************
void obst_writing() {
  int x, y, i;
  char filename1[] = "obst_LB.dat";
  FILE* outfile1;
  outfile1 = fopen(filename1, "w");
  for (y = 0; y < ly; y++) {
    for (x = 0; x < lx; x++) {

      fprintf(outfile1, "%d ", obst[x][y]);
    }
    fprintf(outfile1, "\n");
  }
  fclose(outfile1);

  char filename2[] = "active_nodes.dat";
  FILE* outfile2;
  outfile2 = fopen(filename2, "w");

  for (y = 0; y < ly; y++) {
    for (x = 0; x < lx; x++) {
      fprintf(outfile2, "%d ", act[x][y]);
    }
    fprintf(outfile2, "\n");
  }
  fclose(outfile2);

  char filename3[] = "links.dat";
  FILE* outfile3;
  outfile3 = fopen(filename3, "w");

  for (y = 0; y < ly; y++) {
    for (x = 0; x < lx; x++) {
      for (i = 1; i < Q; i++) {
        if (delta[x][y][i] != 0) {
          fprintf(outfile3, "%d  %d  %d  %f\n", x, y, i, delta[x][y][i]);
        }
      }
    }
  }
  fclose(outfile3);
}

// ****************************************************************************
// *   Output of results to velocity files                                    *
// *   Distribution of verlocity x - y                        *
// ****************************************************************************
void velocity_profile() {
  int x, y, i;
  double u_y, d_loc;
  double u_y1;

  char filename1[] = "yvel_vs_x.dat";

  FILE* outfile1;
  outfile1 = fopen(filename1, "w");
  fprintf(outfile1, "# vitesse u_x      ordonnée \n");

  y = (int)((g[0].x2 - Mby) / dx);
  for (x = 0; x < lx; x++) {
    if (obst[x][y] != -1 && obst[x][y] != nbgrains)
      u_y1 = g[obst[x][y]].v2 / c;
    else {
      u_y = 0;
      d_loc = 0.;
      for (i = 0; i < Q; i++) {
        d_loc = d_loc + f[x][y][i];
      }
      for (i = 0; i < Q; i++) {
        u_y = u_y + f[x][y][i] * ey[i];
      }
      u_y1 = u_y / d_loc;
    }
    fprintf(outfile1, "%d    %.10lf\n", x, u_y1);
  }
  fclose(outfile1);
}

//************************************************************
//   Calculate pressure                                   *
//************************************************************
void pressures() {
  int x, y;
  for (x = 0; x < lx; x++) {
    for (y = 0; y < ly; y++) {
      if (obst[x][y] == -1) {
        press[x][y] =
            (f[x][y][0] + f[x][y][1] + f[x][y][2] + f[x][y][3] + f[x][y][4] +
             f[x][y][5] + f[x][y][6] + f[x][y][7] + f[x][y][8] - rho_moy) *
            c_squ;
      } else
        press[x][y] = 0.;
    }
  }
}
//----------------------------------------------------------

void renderScene(void) {
  long i;

  if (start == 1) {
    if (vib == 1) {
      t = t + dt;
      //  Mby=Mby+0.1*amp*sin(freq*t);
      Mgx = Mgx + amp * sin(freq * t);
      Mdx = Mdx + amp * sin(freq * t);
    }

// FORCES FLUIDES !!!
#ifdef _FLUIDE_
    if (nbsteps % npDEM == 0) forces_fluid();
#endif

    if (nbsteps % UpdateVerlet == 0) initVerlet();
    if (nbsteps % UpdateVerlet == 0) VerletWall();

    /*
     for(i=0; i<=nbgrains_bas-1; i++) {
       g[i].a1=0.;
       g[i].a2=0.;
       g[i].a3=0.;
     }
    */
    for (i = nbgrains_bas - 1; i <= nbgrains - 1; i++) {
      g[i].p = 0;  // reset pressure
      g[i].s = 0.;
      g[i].ifm = 0;
      g[i].f1 = 0.;
      g[i].f2 = 0.;
      g[i].ice = 0;
      g[i].fr = 0.;
      g[i].slip = 0;
      g[i].rw = 0.;
      ic = 0.;
      g[i].M11 = g[i].M12 = g[i].M21 = g[i].M22 = 0.;  // Moments
      g[i].z = 0;   // reset coordination numbers
      g[i].zz = 0;  

      g[i].x1 = g[i].x1 + dt * g[i].v1 + dt2 * g[i].a1 / 2.;
      g[i].x2 = g[i].x2 + dt * g[i].v2 + dt2 * g[i].a2 / 2.;
      g[i].x3 = g[i].x3 + dt * g[i].v3 + dt2 * g[i].a3 / 2.;
      g[i].v1 = g[i].v1 + dt * g[i].a1 / 2.;
      g[i].v2 = g[i].v2 + dt * g[i].a2 / 2.;
      g[i].v3 = g[i].v3 + dt * g[i].a3 / 2.;
    }

    acceleration_grains();

#pragma omp parallel for
    for (i = 0; i <= nbgrains - 1; i++) {
      //	g[i].p=g[i].p/(2.*M_PI*g[i].r); // pressure on grains
      g[i].v1 = g[i].v1 + dt * g[i].a1 / 2.;
      g[i].v2 = g[i].v2 + dt * g[i].a2 / 2.;
      g[i].v3 = g[i].v3 + dt * g[i].a3 / 2.;
    }
    nbsteps++;
  }

  if (nbsteps % stepFilm == 0 && start == 1) {
    write_densities();
    write_grains();
    nFile++;
  }
  if (nbsteps % stepStrob == 0 && start == 1) write_DEM();
  if (nbsteps % (1 * stepStrob) == 0 && start == 1) write_forces();
}

//////////////////////////////////////////////////////////////////////////////
////////////////////////////////  MAIN   /////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {
  time_t time_raw_format;
  struct tm* ptr_time;

  printf("2D LBM-DEM code\n");
  int i;
  double dtmax;

  if (argc != 2) {    
    printf("Error! Usage should be: ./lbmdem filename\n");
    exit(EXIT_FAILURE);
  }

  printf("Opening file : %s\n", argv[1]);
  strcpy(filename_sample, argv[1]);

  c_squ = 1. / 3.;
  sample();
  stat();
  init_density();
  WallPosition();
  xG = -G * sin(angleG);
  yG = -G * cos(angleG);

  dx = (Mdx - Mgx) / (lx - 1);
  printf("no space %le\n", dx);

  // Compute the time step for DEM
  dtmax = (1 / iterDEM) * pi * sqrt(pi * rMin * rMin * rhoS / kg);
  dtLB = dx * dx * (tau - 0.5) / (3 * nu);
  npDEM = (dtLB / dtmax + 1);
  c = dx / dtLB;
  dt = dtLB / npDEM;
  dt2 = dt * dt;

  // write_strob(_INIT_); //JYD

  printf("dtLB=%le,  dtmax=%le,   dt=%le,   npDEM=%d,   c=%lf\n", dtLB, dtmax,
         dt, npDEM, c);
  for (i = 0; i <= nbgrains - 1; i++) {
    rLB[i] = reductionR * g[i].r / dx;
  }
  init_obst();

  //	VerletWall();
  time(&time_raw_format);
  ptr_time = localtime(&time_raw_format);
  printf("Current local time and date: %s", asctime(ptr_time));
  char filename_stats[] = "stats.data";
  s_stats = fopen(filename_stats, "w");

  fprintf(s_stats,
          "#1_t 2_xfront 3_xgrainmax 4_height 5_zmean 6_energie_x 7_energie_y "
          "8_energie_teta 9_energie_cin 10_N0 11_N1 12_N2 13_N3 14_N4 15_N5 "
          "16_energy_Potential 17_Strain_Energy 18_Frictional_Work "
          "19_Internal_Friction 20_Inelastic_Collision 21_Slip "
          "22_Rotational_Work\n");

  fclose(s_stats);

  start = 1;
  do {
    renderScene();
    time(&time_raw_format);
    ptr_time = localtime(&time_raw_format);
    if (nbsteps % UpdateVerlet == 0)
      printf(
          "steps %li steps %le KE %le PE %le SE %le WF %le INCE %le SLIP %le "
          "RW %le Time %s \n",
          nbsteps, nbsteps * dt, energie_cin, energy_p, SE, WF, INCE, TSLIP,
          TRW, asctime(ptr_time));  // IFR TSE TBW
  } while (nbsteps * dt <= duration);

  time(&time_raw_format);
  ptr_time = localtime(&time_raw_format);
  printf("End local time and date: %s", asctime(ptr_time));

  return 0;
}
