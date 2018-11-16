/*  FACULTAD DE FÍSICA DE LA PONTIFICIA UNIVERSIDAD CATÓLICA DE CHILE
    FÍSICA COMPUTACIONAL (FIZ1431) - SEGUNDO SEMESTRE DE 2018

    ALUMNOS:  NICOLÁS ANDRÉ VAN SINT JAN CAMPOS (nicovsj@uc.cl)
              JUAN MANUEL GONZÁLEZ BRANTES (jmgonzalez4@uc.cl)

    PROFESOR: EDGARDO ANDRÉS DÖRNER YAKSIC (endorner@fis.uc.cl)

    FECHA LÍMITE DE ENTREGA: JUEVES 15 DE NOVIEMBRE DE 2018 A LAS 23:59 HORAS

-------------------------------------------------------------------------------
              TAREA 3 - LA ECUACIÓN DEL CALOR DEPENDIENTE DEL TIEMPO

    Calcule la evolución temporal de la distribución de temperatura T(x,t) a
    través de una barra cilíndrica, que está expuesta a un sumidero de calor
    alrededor de x = L/2. Este sumidero de calor está descrito por una función
    Gamma(x) de forma gaussiana. La barra cilíndrica está expuesta a temperaturas
    constantes T_0 Y T_N en sus extremos. Utilice los siguientes parámetros:

    L = 10; k = 1; Theta (mayúscula) = -0.4; l = 1; T_0 = 0; T_N = 2; dt = 0.5

    Adicionalmente, para la condición inicial considere T(x,0) = 0, x \in [0, L].
    Puede obtener gráficos para diferentes pasos temporales, o bien, una simula-
    ción continua (animación) de hasta 300 pasos.

    Para resolver el problema, debe utilizar un método explícito y un método
    implícito de resolución numérica, según lo visto en clases. Justifique ade-
    cuadamente su elección. Por supuesto, no puede emplear ni el método de Euler
    implícito como el explícito. Una resolución implícita implica resolver un
    sistema de ecuaciones lineales. Para esto, debe buscar una librería o código
    externo que le permita resolver el problema. Una buena referencia es:

    Press, William., et. al., "Numerical Recipes in C" (Cambridge University Press,
    First Edition, 1988)

    Sea cual sea el método escogido, debe referenciarlo correctamente, indicando
    las ventajas y desventajas en función del problema numérico que está intentan
    do resolver. No puede utilizar el código visto en clases para este propósito.

--------------------------------------------------------------------------------

    OBJETIVOS PROPUESTOS:
    -Resolver, tanto mediante condiciones de borde definidas como parámetros fí-
    sicos ya indicados, la ecuación del calor dependiente del tiempo para una
    barra cilíndrica, utilizando para estos efectos la física computacional. En
    particular, con métodos numéricos de programación implícito y explícito.
    -Utilizar, de manera óptima, métodos numéricos de programación implícitos y
    explícitos para el propósito anterior; la idea será no reciclar el código
    exhibido en clases, sino, la implementación de un código propio que incluya
    esta teoría.
    -Evidenciar, por medio de animación computacional, la evolución temporal que
    sufre la barra cilíndrica en cuestión, tal que permita establecer una visua-
    lización sobre la fenomenología del problema físico.
    -Referenciar e indicar, mediante la bibliografía apropiada, la implementación
    del código en virtud de sus ventajas y desventajas.

--------------------------------------------------------------------------------

    REFERENCIAS:
    [1] Stickler, Benjamin; Schachinger, Ewald: "Basic Concepts in Computational
    Physics". Primera Edición, Editorial Springer-Verlag, 2014.
    [2] Press, William H., et. al.; “Numerical recipes in C.”. Primera Edición,
    Cambridge University Press, 1988.
    [3] Recktenwald, Gerald: Crank-Nicholson Solution to the Heat Equation:
    ME 448/548 Notes. Portland State University; Mechanical Engineering Department.
    http://web.cecs.pdx.edu/~gerry/class/ME448/notes/pdf/CN_slides.pdf (visitado
    el 15 de noviembre de 2018).
    [4] Dörner, Edgardo. "Clase 11: ODE Valores Iniciales". Diapositivas del curso
    Física Computacional, Segundo Semestre de 2018, Pontificia Universidad Católica
    de Chile.
    [5] Dörner, Edgardo. "Clase 18: Ecuación de Calor 1D". Diapositivas del curso
    Física Computacional, Segundo Semestre de 2018, Pontificia Universidad Católica
    de Chile.

--------------------------------------------------------------------------------

    A continuación, se presenta el código que alberga un método explícito e implícito 
    de resolución numérica para el problema de ecuación de calor en una barra cilín-
    drica. 

    Para el caso del método explícito, no se requirió de búsqueda de librerías o de
    una decisión más pensada para establecer el código, pues ya se contaba con el 
    método de Runge-Kutta explícito para poder resolver el problema.

    Por otro lado, respecto a la elección e implementación del método implícito.
    En primer lugar, se buscó en [1], el libro guía del curso, alguna indicación 
    que pudiera servir, terminando este libro por referenciar a [2]. En [2], se
    encontró un algoritmo para resolver sistemas de ecuaciones con matriz tridiag-
    onal en el capítulo 3 de sistemas de ecuaciones lineales para, posteriormente, 
    encontrar [3] y decidir finalmente, que el método Crank-Nicolson era el indicado
     para resolver implícitamente la ecuación del calor.

    Esta elección tuvo un asidero teórico visto en clases, pues, en [4], se obser-
    vó que es un método implícito, para el que, en [5], existía un argumento muy
    contundente que permitió finalmente decidir usarlo: si bien es un método implí
    cito, es incondicionalmente estable, queriendo decir que, sobre la base de ser
    predictor-corrector, no es alterado por condiciones externas lo suficientemente
    desestabilizadoras.

    No obstante, este método cuenta con la desventaja de que, como ya se señaló,
    es implícito, es decir, su solución no es hallada de forma directa: de ahí,
    la dificultad de implementar el código (de más abajo) que define la matriz
    tridiagonal junto a sus vectores, el vector solución, las condiciones de borde,
    y los índices que definen las posiciones de los inputs.

*/
// ----------------------------  INICIO TAREA  ------------------------------ //

// Importación de las librerías necesarias para la compilación correcta de los
// módulos que contienen y que se utilizarán.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Constantes del material que constituye la barra cilíndrica.
const double L = 10.0;                      // Longitud del cilindro
const double kappa = 1.0;                   // Constante de conductividad térmica

// Establecimiento de condiciones de borde para el problema.
const double T0 = 0.0;                      // Temperatura en x = 0; T(0,t)
const double Tf  = 2.0;                     // Temperatura en x = N; T(N,t)

// Constantes del sumidero que afecta a la barra cilíndrica.
const double Theta = -0.4;
const double l = 1.0;

// Discretización para poder resolver numéricamente.
const int nnodes = 60;                        // Número de nodos espaciales
const double dx = L / (nnodes-1);             // Ancho espacial (debe ser >= 1.0)
const double dt = dx * dx / (kappa * 2.0);    // Ancho temporal (ELEGIDO A PARTIR DE dx 
                                              // PARA ASEGURAR ESTABILDAD EN MÉTODO EXPLÍCITO)
const int nsteps = 300;                      // Número de pasos temporales

// Definición computacional del sumidero de calor.
double **initUMatrix();
void freeUMatrix(double ** u);
void printSolution(double **u, FILE* fp);
double Gamma(double number);
void F(double* Yin, double* Yout, int dim);
void eRK4(double* Yin, double* Yout, double dt, int dim);
void tridiag(double* a, double* b, double* c, double* r, double* u, int n);


int main(int argc, const char * argv[]) {

    printf("**********************************************\n");
    printf("Resolución de la ecuacion de calor u_t = k*u_xx - f \n");
    printf("mediante el método de Runge-Kutta explícito\n");
    printf("**********************************************\n");

    double **u_exp = initUMatrix();
    double **u_imp = initUMatrix();


    // INICIALIZACIÓN PARA MÉTODO IMPLÍCITO
    double *a = (double*)malloc(nnodes*sizeof(double));
    double *b = (double*)malloc(nnodes*sizeof(double));
    double *c = (double*)malloc(nnodes*sizeof(double));
    double *r = (double*)malloc(nnodes*sizeof(double));

    // Seteo de la matriz tridiagonal
    for (int i = 0; i < nnodes; ++i) {
      a[i] = - kappa / (2*dx*dx);
      c[i] = a[i];
      b[i] = 1/dt - (a[i] + c[i]);
    }

    // Fijación de las condiciones de borde para garantizar matriz cuadrada de n x n
    b[0] = 1; b[nnodes-1] = 1;
    c[0] = 0; a[nnodes-1] = 0;

    // IMPLEMENTACIÓN DE AMBOS MÉTODOS
    for (int t=0; t<nsteps-1; t++) {
      // EXPLÍCITO RUNGE-KUTTA 
      eRK4(u_exp[t], u_exp[t+1], dt, nnodes);
      // IMPLÍCITO CRANK-NICOLSON
      r[0] = T0;
      for (int x = 1; x < nnodes-1; ++x) {
        r[x] = - a[x]*u_imp[t][x-1] + (1.0/dt + c[x] + a[x])*u_imp[t][x]
               - c[x]*u_imp[t][x+1] - Gamma(x*dx);
      }
      r[nnodes-1] = Tf;


      tridiag(a, b, c, r, u_imp[t+1], nnodes);
    }

    // Registro y retención de los datos obtenidos mediante el método numérico
    // para el posterior análisis y uso gráfico animado.
    FILE *fp_exp = fopen("res_explicit.csv", "w");
    FILE *fp_imp = fopen("res_implicit.csv", "w");

    printSolution(u_exp, fp_exp);
    printSolution(u_imp, fp_imp);

    // Limpieza

    fclose(fp_exp);
    fclose(fp_imp);

    freeUMatrix(u_exp);
    freeUMatrix(u_imp);

    free(a);
    free(b);
    free(c);
    free(r);

    // Exhibición de los valores finales como muestra de la efectividad del código
    printf("Simulación terminada exitosamente\n");

    return 0;
}

// Función sumidero-fuente en L/2 para la barra del calor descrita por la ecuación
// del calor no homogénea.
double Gamma(double number) {
  return Theta / l * exp(-1.0 * pow((number - L/2)/l, 2));
}

void freeUMatrix(double ** u) {
/* Función que realiza la limpieza de la memoria utlizada por la matriz u */

  for (int i = 0; i < nsteps; ++i)
    {
      free(u[i]);
    }
    free(u);
}

double ** initUMatrix() {
/* Inicialización de la matriz u que guarda las soluciones encontradas a través del tiempo */
  double **u = (double **)malloc(nsteps * sizeof(double*));
  for (int i = 0; i < nsteps; ++i) {
    u[i] = (double *)malloc(nnodes * sizeof(double));
  }

  // Condición inicial (t = 0)
  for (int i=0; i<nnodes; i++) {
      u[0][i] = 0.0;
  }

  // Condiciones de borde (x = 0; x = N)
  for (int i = 0; i < nsteps; ++i) {
      u[i][0] = T0;
      u[i][nnodes-1] = Tf;
  }

  return u;
}

void printSolution(double **u, FILE* fp) {
/* Dada la matriz u RESUELTA imprime los resultados obtenidos en el file pointer fp */

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
}

void F(double* Yin, double* Yout, int dim) {
    Yout[0] = 0.0; Yout[dim-1] = 0.0;     // Truco para que las condiciones de borde funcionen
    for (int i = 1; i < dim-1; ++i) {
        Yout[i] = kappa/(dx*dx) * (Yin[i-1] - 2*Yin[i] + Yin[i+1]) - Gamma(i*dx);
    }
}

// Inicio - Método eRK-4
void eRK4(double* Yin, double* Yout, double dt, int dim) {

    /*
    Argumentos necesarios:
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

void tridiag(double* a, double* b, double* c, double* r, double* u, int n)
/* ALGORITMO OBTENIDO DESDE [2].
   Resuelve para un vector u[1..n] el sistema de equaciones tridiagonal
   dado por

   |b_1 c_1  0  ...                   |   | u_1 |   | r_1 |
   |a_2 b_2 c_2 ...                   |   | u_2 |   | r_2 |
   |            ...                   | . | ... | = | ... |
   |            ... a_n-1  b_n-1 c_n-1|   |u_n-1|   |r_n-1|
   |            ...   0     a_n   b_n |   | u_n |   | r_n |

   a[1..n], b[1..n], c[1..n] y r[1..n] son vectores de input y no son
   modificados.
*/

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
