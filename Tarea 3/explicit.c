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

    A continuación, se presenta el código que alberga un método explícito de re-
    solución numérica para el problema de ecuación de calor en una barra cilíndrica.
    En este caso, es decir, en la resolución del problema mediante el método explí-
    cito, no se requirió de búsqueda de librerías o de una decisión más pensada para
    establecer el código, pues ya se contaba con el método de Runge-Kutta explícito
    para poder resolver el problema.

    Nuevamente, esta elección tuvo un contacto directo con lo visto en clases, pues
    de éstas, tanto por el sustento teórico como por la práctica misma de ejercicios
    en clases, ya era sabido que el método explícito garantizaba efectividad para el
    problema presente acá.

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
const int nnodes = 60;                        // Número de pasos temporales
const double dx = L / (nnodes-1);             // Ancho espacial (debe ser >= 1.0)
const double dt = dx * dx / (kappa * 2.0);    // Ancho temporal
const int nsteps = 1200;                      // Número de pasos espaciales

// Definición computacional del sumidero de calor.
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

    // Condiciones de borde (x = 0; x = N)
    for (int i = 0; i < nsteps; ++i) {
        u[i][0] = T0;
        u[i][nnodes-1] = Tf;
    }

    // IMPLEMENTACIÓN DEL MÉTODO RUNGE-KUTTA EXPLÍCITO

    for (int t=0; t<nsteps-1; t++) {
      eRK4(u[t], u[t+1], dt, nnodes);
    }

    // Registro y retención de los datos obtenidos mediante el método numérico
    // para el posterior análisis y uso gráfico animado.
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

    // Exhibición de los valores finales como muestra de la efectividad del código
    printf("Simulación terminada exitosamente\n");

    return 0;
}

// Función sumidero-fuente en L/2 para la barra del calor descrita por la ecuación
// del calor no homogénea.
double Gamma(double number) {
  return Theta / l * exp(-1.0 * pow((number - L/2)/l, 2));
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
