/*
Pontificia Universidad Católica de Chile - Segundo Semestre de 2018.
Instituto de Física - Curso: Física Computacional - FIZ1431.
Nombre de los alumnos autores:
        -Nicolás Van Sint Jan Campos.
        -Juan Manuel González Brantes.
Nombre del profesor:
        -Edgardo Dörner Yaksic.
Fecha de entrega: 13 Septiembre 2018, 23:59 horas.

------------------------------------ TAREA 1 -----------------------------------
                 RESOLVIENDO EL PROBLEMA DE DOS CUERPOS DE KEPLER
--------------------------------------------------------------------------------

ENUNCIADO:
  Resolver de forma numérica el problema de Kepler, es decir, resolver las
  ecuaciones diferenciales ordinarias para las coordenadas phi y rho.
  Utilizar métodos explícitos para resolverlas.

CONSIDERACIONES:
  (1) Recordar que tanto la energía E como el momento angular l del problema
  son cantidades conservadas.
  (2) Reescalamiento del problema: las distancias Tierra - Sol serán expresadas
  en unidades astronómicas (UA) y el tiempo será considerado en años.

RESULTADO ESPERADO:
  El programa será capaz de manipular los datos, resolver los modelos
  matemáticos, y entregar un archivo del tipo .txt que contendrá los valores de
  phi y rho para todo tiempo t. Éstos serán útiles para luego estudiar
  y analizar este problema.

--------------------------------------------------------------------------------
----------------------------- INICIO DEL CÓDIGO --------------------------------
--------------------------------------------------------------------------------
*/
// Se importan los módulos stdio y stdlib, pues contienen funciones predefinidas
// que pueden utilizarse.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double const SEC_YEAR = 3.154e7;        // Segundos en un año
double const MET_AU = 1.496e11;         // Metros en una unidad astronómica (AU)
double const KG_SUN = 1.988e30;         // Masa del sol
double const KG_EARTH = 5.972e24;       // Masa tierra
double const G = 6.67408e-11;           // Constante de gravitación universal
                                        // en m^3/(kg*s^2)

double U_eff(double rho, double l, double mu, double alpha);

// Se ingresan los parámetros del sistema.
int main(int argc, const char * argv[]) {
    // Parametros del sistema.
    const double t0 = 0.0;                // Tiempo inicial en años
    const double tf  = 300.0;             // Instante final, en años
    const double dt = 0.001;               // Ancho temporal, en años

    const double rho0 = 1.0;              // Distancia inicial, en UA
    const double rhop0 = 5.5;             // Velocidad radial inicial, en UA/año
    const double phip0 = M_PI*2;          // Velocidad angular, en rad/año
    const double phi0 = 0.0;              // Ángulo inicial, en radianes.
    const double m1 = KG_SUN/KG_EARTH;    // Masa 1 (sol), en masas terrestres.
    const double m2 = 1.0;                // Masa 2 (tierra), en masas terrestres.

    const double mu = m1 * m2 / (m1 + m2);  // Masa reducida
    const double alpha = (G*pow(SEC_YEAR, 2)*KG_EARTH/pow(MET_AU, 3))*m1*m2;// Alpha del potencial U(rho) = -alpha/rho

    const double l = phip0 * mu * pow(rho0, 2); // Momento angular
    const double E = pow(rhop0, 2) * mu * 0.5 + U_eff(rho0,l,mu,alpha);  // Energía

    // A partir de los parametros del sistema calculamos el numero de pasos necesarios.
    // es decir, acá se decide la cantidad de subintervalos que utilizará el programa.
    const int nsteps = 1 + (tf-t0)/dt;      // número de pasos temporales.

    printf("|-------------------------------------------------|\n");
    printf("|Resolución del problema de Kepler                |\n");
    printf("|mediante el método numérico de Euler (explícito) |\n");
    printf("|-------------------------------------------------|\n");

    // Discretizacion del eje temporal.
    double *t = (double*) malloc(nsteps*sizeof(double));

    for (int i=0; i<nsteps; i++) {
        t[i] = i*dt;
    }

    printf("Parámetros del problema:\n");
    printf("\t t0 \t\t = \t%10.2f (años)\n", t0);
    printf("\t tf \t\t = \t%10.2f (años)\n", t[nsteps-1]);
    printf("\t nsteps \t = \t%10d\n", nsteps);
    printf("\t dt \t\t = \t%10.2f (años)\n", dt);

    // Resolucion de la ecuacion diferencial.
    double *rho = (double*) malloc(nsteps*sizeof(double));
    double *phi = (double*) malloc(nsteps*sizeof(double));
    double *ueff = (double*) malloc(nsteps*sizeof(double));

    rho[0] = rho0;
    phi[0] = phi0;

    double de, sgn;
    if (rhop0 == 0) sgn = 1;
    else sgn = rhop0 / fabs(rhop0);
    for (int i=0; i<nsteps-1; i++) {
        ueff[i] = U_eff(rho[i], l, mu, alpha);
        // Obtenemos rho para i+1
        de = E - U_eff(rho[i],l,mu,alpha);
        if (de <= 0 && i != 0) {
          sgn *= -1;
        }
        rho[i+1] = sgn*sqrt(2/mu *fabs(E - U_eff(rho[i],l,mu,alpha))) * dt + rho[i];

        phi[i+1] = l /(mu * pow(rho[i+1], 2)) * dt + phi[i];

    }
    ueff[nsteps-1] = U_eff(rho[nsteps-1],l,mu,alpha);
    // Guardamos los resultados en un archivo para su analisis posterior.
    FILE *fp = fopen("res_explicit.txt", "w");

    fprintf(fp, "t,rho,phi,ueff,rhox,ueffteo,E\n");
    double x;
    for (int i=0; i<nsteps; i++) {
        x = 0.1*(i+1);
        fprintf(fp, "%10.5f, %10.5f, %10.5f, %10.5f, %10.5f, %10.5f, %10.5f\n",
        t[i], rho[i], phi[i], ueff[i], x, U_eff(x,l,mu,alpha), E);
    }

    // Imprimimos los valores finales como muestra.
    printf("|----------------------------------------------|\n");
    printf("|Simulación terminada de forma exitosa:        |\n");
    printf("|tf = %2.5f, rhof = %2.5f, phif = %2.5f |\n",
          t[nsteps-1], rho[nsteps-1], phi[nsteps-1]);
    printf("|----------------------------------------------| \n ");

    // Limpieza
    free(t);
    free(rho);
    free(phi);
    fclose(fp);

    return 0;
}

double U_eff(double rho, double l, double mu, double alpha) {
  return -alpha / rho + pow(l,2) / (2 * mu * pow(rho, 2));
}
