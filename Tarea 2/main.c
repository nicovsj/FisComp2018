/*  Facultad de Física de la Pontificia Universidad Católica de Chile
    Física Computacional (FIZ1431) - Segundo Semestre de 2018

    Alumnos:  Nicolás Van Sint Jan Campos (nicovsj@uc.cl)
              Juan Manuel González Brantes (jmgonzalez4@uc.cl)

    Profesor: Edgardo Dörner Yaksic (endorner@fis.uc.cl)

    --- TAREA 2: ECUACIONES DIFERENCIALES CON VALORES INICIALES. ---

    Enunciado: El Péndulo doble. Resolver, numéricamente, el problema de péndulo
    doble, es decir, resolver las ecuaciones de Hamilton para el movimiento del
    sistema mediante el método de Runge-Kutta explícito de cuarto orden (eRK-4).
    Utilizar las condiciones iniciales empleadas en la sección 6.2 del libro
    guía, y comprobar resultados con los ejemplos allí mostrados.
    Objetivos: Entregar los valores de phi_i (en el programa qi con i = 1, 2) y
    p_i (en el programa pi con i = 1, 2) con respecto al tiempo. Entregar un
    segundo archivo de salida que entregue las posiciones de cada masa en coor-
    denadas cartesianas.

    Objetivos propuestos: Entregar animaciones y gráficos interactivos que ayuden
    a comprender la dinámica del péndulo doble, problema característico de la
    mecánica intermedia, analítica, o de Hamilton.

    BIBLIOGRAFÍA: Libro guía del curso
    Stickler, Benjamin and Schachinger, Ewald: "Basic Concepts in Computational
    Physics", Editorial Springer Verlag, Ebook de 2014. */

// ----------------------------  INICIO TAREA  ------------------------------ //

// Importación de los headers necesarios para poder trabajar en la tarea.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Definición de constantes para el problema físico. En particular, se han
   definido de acuerdo a valores sugeridos por el texto guía del curso.
   Observación: Se asumirá conocida la Mecánica de Hamilton y Lagrange, de modo
   que no se definirán los ángulos ni los momentos generalizados. */


const double m = 1.0;     // masa de los péndulos.
const double l = 1.0;     // largo de los péndulos.
const double g = 9.8067;  // magnitud de la aceleración de gravedad terrestre.
const double q1_0 = 0.0;  // q1 es phi_1, obertura angular del primer péndulo.
const double q2_0 = 0.0;  // q2 es phi_2, obertura angular del segundo péndulo.
const double p1_0 = 4.0;  // p1 es p_1, momento generalizado para la masa 1.
const double p2_0 = 2.0;  // p2 es p_2, momento generalizado para la masa 2.

void F(double* Yin, double* Yout);
double f1(double Yin[4]);   // Creación de la función que albergará la expresión para q1.
double f2(double Yin[4]);   // Creación de la función que albergará la expresión para q2.
double f3(double Yin[4]);   // Creación de la función que albergará la expresión para p1.
double f4(double Yin[4]);   // Creación de la función que albergará la expresión para p2.
void eRK4(double* Yin, double* Yout, double dt);

int main(int argc, const char * argv[]) {

    // Definición de parámetros de tiempo para el sistema.
    const double t0 = 0.0;    // instante inicial, medido en segundos.
    const double tf  = 60.0;  // instante final, medido en segundos.
    const double dt = .01;   // paso temporal, división de la grilla de tiempo.

    // Definición del número de pasos temporales.
    const int nsteps = 1 + (tf-t0)/dt;    // Número de pasos temporales.

    // Mensaje de interacción con el usuario.
    printf("**************************************************** \n");
    printf("*     RESOLUCIÓN DEL PROBLEMA DE PÉNDULO DOBLE     * \n");
    printf("*      MEDIANTE MÉTODO DE RUNGE-KUTTA4 (eRK-4)     * \n");
    printf("**************************************************** \n");

    // Discretización del tiempo: Construcción de la grilla de tiempo.
    double *t = (double*) malloc(nsteps*sizeof(double));
    for (int i=0; i<nsteps; i++) {
        t[i] = i*dt;
    }

    // Definición de arreglos para guardar variables calculadas.
    double *q1 = (double*) malloc(sizeof(double)*nsteps);
    double *q2 = (double*) malloc(sizeof(double)*nsteps);
    double *p1 = (double*) malloc(sizeof(double)*nsteps);
    double *p2 = (double*) malloc(sizeof(double)*nsteps);

    // Definición de valores iniciales.
    q1[0] = q1_0;
    q2[0] = q2_0;
    p1[0] = p1_0;
    p2[0] = p2_0;

    // Arreglos auxiliares para inicializar el programa.
    // En esta parte del código, el programa comienza a calcular recursivamente.
    double yn[4] = {0.0, 0.0, 0.0, 0.0};
    double ynp[4] = {0.0, 0.0, 0.0, 0.0};

    for (int i=0; i<nsteps-1; i++) {
        yn[0] = q1[i];
        yn[1] = q2[i];
        yn[2] = p1[i];
        yn[3] = p2[i];

        eRK4(yn, ynp, dt);

        q1[i+1] = ynp[0];
        q2[i+1] = ynp[1];
        p1[i+1] = ynp[2];
        p2[i+1] = ynp[3];
    }

    // Anexión de los resultados en un archivo .csv para su posterior análisis.
    FILE *fp = fopen("res_erk4.csv", "w");
    fprintf(fp, "t,phi1,phi2,p1,p2,x1,z1,x2,z2\n");
    double x1, z1, x2, z2;
    for (int i=0; i<nsteps; i++) {

        x1 = l* sin(q1[i]);
        z1 = 2*l - l*cos(q1[i]);
        x2 = l*(sin(q1[i]) + sin(q2[i]));
        z2 = 2*l - l*(cos(q1[i]) + cos(q2[i]));

        fprintf(fp, "%2.5f, %2.5f, %2.5f, %2.5f, %2.5f, %2.5f, %2.5f, %2.5f, %2.5f\n",
                t[i], q1[i], q2[i], p1[i], p2[i], x1, z1, x2, z2);
    }
    fclose(fp);

    // Impresión de los valores finales calculados por el método para mostrar
    // la funcionalidad del código.
    printf("Valores finales Runge-Kuta 4\n");
    printf("tf = %2.5f, q1 = %2.5f, q2 = %2.5f, p1 = %2.5f, p2 = %2.5f\n",
           t[nsteps-1], q1[nsteps-1], q2[nsteps-1], p1[nsteps-1], p2[nsteps-1]);

    // Limpieza de variables.
    free(t);
    free(p1);
    free(p2);
    free(q1);
    free(q2);

    return 0;
}

// Utilización de Y_in = {q1, q2, p1, p2} para definir los Y_out[i], a partir
// de la inicialización de más arriba.
void F(double Yin[4], double Yout[4]) {
  Yout[0] = f1(Yin);
  Yout[1] = f2(Yin);
  Yout[2] = f3(Yin);
  Yout[3] = f4(Yin);
}

// Implementación del Método Runge-Kutta a orden 4
void eRK4(double* Yin, double* Yout, double dt) {
  double Y1[4], Y2[4], Y3[4], Y4[4];
  double _Y1[4], _Y2[4], _Y3[4], _Y4[4];

  for (size_t i = 0; i < 4; i++) {  // Construcción de Y1[i] por eRK-4
    Y1[i] = Yin[i];
  }

  F(Y1, _Y1);  // Guardamos F(Y1) en un arreglo auxiliar _Y1
  for (size_t i = 0; i < 4; i++) {
    Y2[i] = dt/2 * _Y1[i] + Yin[i];
  }

  F(Y2, _Y2);
  for (size_t i = 0; i < 4; i++) {
    Y3[i] = dt/2 * _Y2[i] + Yin[i];
  }

  F(Y3, _Y3);
  for (size_t i = 0; i < 4; i++) {
    Y4[i] = dt*_Y3[i] + Yin[i];
  }

  F(Y4, _Y4);
  for (size_t i = 0; i < 4; i++) {
    Yout[i] = Yin[i] + dt/6 * (_Y1[i] + 2*_Y2[i] + 2*_Y3[i] + _Y4[i]);
  }
}

// Implementación del primer ángulo de obertura phi_1 = q1.
double f1(double Yin[4]) {
  double p1, p2, q1, q2, ret;
  q1 = Yin[0];
  q2 = Yin[1];
  p1 = Yin[2];
  p2 = Yin[3];
  ret = p1 - p2*cos(q1-q2);
  ret /= m*pow(l,2)*(1 + pow(sin(q1-q2), 2));
  return ret;
}

// Implementación del primer ángulo de obertura phi_2 = q2.
double f2(double Yin[4]) {
  double p1, p2, q1, q2, ret;
  q1 = Yin[0];
  q2 = Yin[1];
  p1 = Yin[2];
  p2 = Yin[3];
  ret = 2*p2 - p1*cos(q1-q2);
  ret /= m*pow(l,2)*(1 + pow(sin(q1-q2), 2));
  return ret;
}

// Implementación de la derivada del momento generalizado p1.
double f3(double Yin[4]) {
  double p1, p2, q1, q2, ret;
  q1 = Yin[0];
  q2 = Yin[1];
  p1 = Yin[2];
  p2 = Yin[3];
  ret = pow(p1, 2) + 2*pow(p2, 2) - 2*p1*p2*cos(q1-q2);
  ret /= 1 + pow(sin(q1-q2), 2);
  ret *= sin(q1-q2) * cos(q1 -q2);
  ret -= p1*p2*sin(q1-q2);
  ret /= m*pow(l,2)*(1 + pow(sin(q1-q2), 2));
  ret -= 2*m*g*l*sin(q1);
  return ret;
}

// Implementación de la derivada del momento generalizado p2.
double f4(double Yin[4]) {
  double p1, p2, q1, q2, ret;
  q1 = Yin[0];
  q2 = Yin[1];
  p1 = Yin[2];
  p2 = Yin[3];
  ret = pow(p1, 2) + 2*pow(p2, 2) - 2*p1*p2*cos(q1-q2);
  ret /= 1 + pow(sin(q1-q2), 2);
  ret *= sin(q1-q2) * cos(q1 -q2) *-1;
  ret += p1*p2*sin(q1-q2);
  ret /= m*pow(l,2)*(1 + pow(sin(q1-q2), 2));
  ret -= m*g*l*sin(q2);
  return ret;
}
