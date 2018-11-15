/*  Facultad de Física de la Pontificia Universidad Católica de Chile
    Física Computacional (FIZ1431) - Segundo Semestre de 2018

    Alumnos:  Nicolás Van Sint Jan Campos (nicovsj@uc.cl)
              Juan Manuel González Brantes (jmgonzalez4@uc.cl)

    Profesor: Edgardo Dörner Yaksic (endorner@fis.uc.cl)

    --- TAREA 3: . ---

    Enunciado: 

    Objetivos propuestos: 

    BIBLIOGRAFÍA: Libro guía del curso
    Stickler, Benjamin and Schachinger, Ewald: "Basic Concepts in Computational
    Physics", Editorial Springer Verlag, Ebook de 2014. */

// ----------------------------  INICIO TAREA  ------------------------------ //


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Constantes del material
const double L = 10.0;                      // Largo del cilíndro
const double kappa = 1.0;                   // Constante

// Condiciones de borde
const double T0 = 0.0;                      // Temperatura T(0, t)
const double Tf  = 2.0;                     // Temperatura T(N, t)

// Constantes del sumidero
const double Theta = -0.4;
const double l = 1.0;

// Elección de discretizaciones para resolver 
// numéricamente
const int nnodes = 60;
                    // Número de pasos
const double dx = L / (nnodes-1);                      // Ancho espacial (debe ser >= 1.0)
const double dt = dx * dx / (kappa * 2.0) + 0.1;                      // Ancho temporal
const int nsteps = 1200; 

double Gamma(double number);
void F(double* Yin, double* Yout, int dim);
void tridiag(double* a, double* b, double* c, double* r, double* u, int n);

int main(int argc, const char * argv[]) {



    printf("**********************************************\n");
    printf("Resolución de la ecuacion de calor u_t = k*u_xx - f \n");
    printf("mediante el método de Runge-Kutta explícito\n");
    printf("**********************************************\n");



    double u[nsteps][nnodes];

    // Condición inicial (t=0)
    for (int i = 0; i < nnodes; ++i) {
      u[0][i] = 0.0;
    }

    // Condicion de borde (x=0; x=N)
    for (int i = 0; i < nsteps; ++i) {
          u[i][0] = T0;
          u[i][nnodes-1] = Tf;
    }



    // Start allocating a, b, c and r vectors
    double *a = (double*)malloc(nnodes*sizeof(double));
    double *b = (double*)malloc(nnodes*sizeof(double));
    double *c = (double*)malloc(nnodes*sizeof(double));
    double *r = (double*)malloc(nnodes*sizeof(double));


    for (int i = 0; i < nnodes; ++i) {
      a[i] = - kappa / (2*dx*dx);
      c[i] = a[i];
      b[i] = 1/dt - (a[i] + c[i]);
    }

    // Arreglar para condiciones de borde
    b[0] = 1; b[nnodes-1] = 1;
    c[0] = 0; a[nnodes-1] = 0;


    // IMPLICIT CRANK - NICOLSON

    for (int t = 0; t < nsteps-1; ++t) {
      r[0] = T0;
      for (int x = 1; x < nnodes-1; ++x) {
        r[x] = - a[x]*u[t][x-1] + 
                (1.0/dt + c[x] + a[x])*u[t][x] 
               - c[x]*u[t][x+1] 
               - Gamma(x*dx);
      }
      r[nnodes-1] = Tf;



      tridiag(a, b, c, r, u[t+1], nnodes);
    }

    // Guardamos los resultados en un archivo para su analisis posterior.
    FILE *fp = fopen("res_implicit.csv", "w");
    
    fprintf(fp, "t,");
    for (int i = 0; i < nnodes-1; ++i) {
        fprintf(fp, "x=%2.3f,", i*dx);
    }
    fprintf(fp, "x=%2.3f\n", dx*(nnodes-1));

    
    for (int t=0; t<nsteps; t++) {

        fprintf(fp, "%3d,", t);

        for (int x = 0; x < nnodes-1; ++x) {
            fprintf(fp, "%10.5f,", u[t][x]);
        }
        fprintf(fp, "%10.5f\n", u[t][nnodes-1]);
    }

   // Imprimimos los valores finales como muestra.
    printf("Simulación terminada exitosamente\n");

    free(fp);
    free(a);
    free(b);
    free(c);
    free(r);
    
    return 0;
}

  

double Gamma(double number) {
  /* Sumidero/fuente en L/2 para la ecuación de calor inhomogénea */
  return Theta / l * exp(-1.0 * pow((number - L/2)/l, 2));
}

void F(double* Yin, double* Yout, int dim) {
    Yout[0] = 0.0; Yout[dim-1] = 0.0;  // Truco para que las condiciones de borde funcionen
    for (int i = 1; i < dim-1; ++i) {
        Yout[i] = kappa/(dx*dx) * (Yin[i-1] - 2*Yin[i] + Yin[i+1]) - Gamma(i*dx);
    }
}


void tridiag(double* a, double* b, double* c, double* r, double* u, int n)
/* Resuelve para un vector u[1..n] el sistema de equaciones tridiagonal
   dado por

   |b_1 c_1  0  ...                   |   | u_1 |   | r_1 |
   |a_2 b_2 c_2 ...                   |   | u_2 |   | r_2 |
   |            ...                   | . | ... | = | ... |
   |            ... a_n-1  b_n-1 c_n-1|   |u_n-1|   |r_n-1|
   |            ...   0     a_n   b_n |   | u_n |   | r_n |

   a[1..n], b[1..n], c[1..n] y r[1..n] son vectores de input y no son
   modificados. */
{
  int j;
  double bet, *gam;



  gam = (double*)malloc(n*sizeof(double));
  u[0] = r[0] / (bet=b[1]);
  for (j=1;j<=n-1; ++j)
  {
    gam[j] = c[j-1]/bet;
    bet = b[j]-a[j]*gam[j];

    u[j] = (r[j]- a[j]*u[j-1])/bet;
  }

  for (j=(n-2);j>=0;j--)
    u[j] -= gam[j+1]*u[j+1];

  free(gam);
}
