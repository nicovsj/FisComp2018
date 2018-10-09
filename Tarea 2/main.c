#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Parametros del problema

const double m = 1.0;
const double l = 1.0;
const double g = 9.8067;
const double q1_0 = 0.0;
const double q2_0 = 0.0;
const double p1_0 = 0.0;
const double p2_0 = 6.5;

void F(double* Yin, double* Yout);
double f1(double Yin[4]);
double f2(double Yin[4]);
double f3(double Yin[4]);
double f4(double Yin[4]);
void eRK4(double* Yin, double* Yout, double dt);

int main(int argc, const char * argv[]) {
    // Parametros del sistema.
    const double t0 = 0.0;              // instante inicial (s)
    const double tf  = 60.0;            // instante final (s)
    const double dt = .001;              // ancho temporal (s)

    // A partir de los parametros del sistema calculamos el numero de pasos
    // necesario.
    const int nsteps = 1 + (tf-t0)/dt;      // numero de pasos temporales

    printf("******************************************\n");
    printf(" Resolución del pendulo doble \n");
    printf(" mediante Runge-Kutta         \n");
    printf("******************************************\n");

    // Discretizacion del eje temporal.
    double *t = (double*) malloc(nsteps*sizeof(double));
    for (int i=0; i<nsteps; i++) {
        t[i] = i*dt;
    }

    // Arreglos para guardar variables.
    double *q1 = (double*) malloc(sizeof(double)*nsteps);
    double *q2 = (double*) malloc(sizeof(double)*nsteps);
    double *p1 = (double*) malloc(sizeof(double)*nsteps);
    double *p2 = (double*) malloc(sizeof(double)*nsteps);

    // Valores iniciales
    q1[0] = q1_0;
    q2[0] = q2_0;
    p1[0] = p1_0;
    p2[0] = p2_0;

    // Arreglos auxiliares
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

    // Guardamos los resultados en un archivo para su analisis posterior.
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

    // Imprimimos los valores finales como muestra.
    printf("Valores finales Runge-Kuta 4\n");
    printf("tf = %2.5f, q1 = %2.5f, q2 = %2.5f, p1 = %2.5f, p2 = %2.5f\n",
           t[nsteps-1], q1[nsteps-1], q2[nsteps-1], p1[nsteps-1], p2[nsteps-1]);

    // Limpieza
    free(t);
    free(p1);
    free(p2);
    free(q1);
    free(q2);

    return 0;
}

// Y_in = {q1, q2, p1, p2}
void F(double Yin[4], double Yout[4]) {
  Yout[0] = f1(Yin);
  Yout[1] = f2(Yin);
  Yout[2] = f3(Yin);
  Yout[3] = f4(Yin);
}

void eRK4(double* Yin, double* Yout, double dt) {
  double Y1[4], Y2[4], Y3[4], Y4[4];
  double _Y1[4], _Y2[4], _Y3[4], _Y4[4];
  for (size_t i = 0; i < 4; i++) {
    Y1[i] = Yin[i];
  }
  // printf("F(Y1)\n");
  F(Y1, _Y1);
  for (size_t i = 0; i < 4; i++) {
    Y2[i] = dt/2 * _Y1[i] + Yin[i];
  }
  // printf("F(Y2)\n");
  F(Y2, _Y2);
  for (size_t i = 0; i < 4; i++) {
    Y3[i] = dt/2 * _Y2[i] + Yin[i];
  }
  // printf("F(Y3)\n");
  F(Y3, _Y3);
  for (size_t i = 0; i < 4; i++) {
    Y4[i] = dt*_Y3[i] + Yin[i];
  }
  // printf("F(Y4)\n");
  F(Y4, _Y4);
  for (size_t i = 0; i < 4; i++) {
    Yout[i] = Yin[i] + dt/6 * (_Y1[i] + 2*_Y2[i] + 2*_Y3[i] + _Y4[i]);
  }
}

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
