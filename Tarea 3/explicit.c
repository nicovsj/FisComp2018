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
const double dt = dx * dx / (kappa * 2.0);                      // Ancho temporal
const int nsteps = 1200; 

double Gamma(double number);
void F(double* Yin, double* Yout, int dim);
void eRK4(double* Yin, double* Yout, double dt, int dim);

int main(int argc, const char * argv[]) {

    printf("**********************************************\n");
    printf("Resolución de la ecuacion de calor u_t = k*u_xx - f \n");
    printf("mediante el método de Runge-Kutta explícito\n");
    printf("**********************************************\n");

    double u[nsteps][nnodes];

    // Condición inicial (t = 0)
    for (int i=0; i<nnodes; i++) {
        u[0][i] = 0.0;
    }

    // Condición de borde (x = 0; x = N)
    for (int i = 0; i < nsteps; ++i) {
        u[i][0] = T0;
        u[i][nnodes-1] = Tf;
    }

    // EXPLICIT RUNGE-KUTTA

    for (int t=0; t<nsteps-1; t++) {
      eRK4(u[t], u[t+1], dt, nnodes);
    }

    // Guardamos los resultados en un archivo para su analisis posterior.
    FILE *fp = fopen("res_explicit.csv", "w");

    fprintf(fp, "t,");
    for (int i = 0; i < nnodes-1; ++i) {
        fprintf(fp, "x=%2.3f,", dx*i);
    }
    fprintf(fp, "x=%2.3f\n", dx*(nnodes-1));

    
    for (int t=0; t<nsteps; t++) {

        fprintf(fp, "%3d,", t);

        for (int x = 0; x < nnodes-1; ++x) {
            fprintf(fp, "%10.5f,", u[t][x]);
        }
        fprintf(fp, "%10.5f\n", u[t][nnodes-1]);
    }
    fclose(fp);

    // Imprimimos los valores finales como muestra.
    printf("Simulación terminada exitosamente\n");
    
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

void eRK4(double* Yin, double* Yout, double dt, int dim) {
    /* Implementación del Método Runge-Kutta a orden 4

    Args:
        - double* Yin: Arreglo input a realizar (Y_n)
        - double* Yout: Arreglo output esperado luego de eRK4 (Y_(n+1))
        - double dt: Paso temporal
        - double dim: Número de entradas que posee Y_n 
    */

  double Y1[dim], Y2[dim], Y3[dim], Y4[dim];
  double _Y1[dim], _Y2[dim], _Y3[dim], _Y4[dim];

  for (size_t i = 0; i < dim; i++) {  // Construcción de Y1[i] por eRK-4
    Y1[i] = Yin[i];
  }

  F(Y1, _Y1, dim);  // Guardamos F(Y1) en un arreglo auxiliar _Y1
  for (size_t i = 0; i < dim; i++) {
    Y2[i] = dt/2 * _Y1[i] + Yin[i];
  }

  F(Y2, _Y2, dim); // Guardamos F(Y2) en un arreglo auxiliar _Y2
  for (size_t i = 0; i < dim; i++) {
    Y3[i] = dt/2 * _Y2[i] + Yin[i];
  }

  F(Y3, _Y3, dim); // Guardamos F(Y3) en un arreglo auxiliar _Y3
  for (size_t i = 0; i < dim; i++) {
    Y4[i] = dt*_Y3[i] + Yin[i];
  }

  F(Y4, _Y4, dim); // Guardamos F(Y4) en un arreglo auxiliar _Y4
  for (size_t i = 0; i < dim; i++) {
    Yout[i] = Yin[i] + dt/6 * (_Y1[i] + 2*_Y2[i] + 2*_Y3[i] + _Y4[i]);
  }
}
