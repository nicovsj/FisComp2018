//
//  main.c
//  rng_dispersion
//
//  Modelo simple del transporte de particulas mediante el metodo Monte Carlo.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

struct Rng {
    /* Parametros del pseudo generador de numeros aleatorios */
    long int a; // mucho cuidado con la eleccion de 'long int'!
    long int c;
    long int m;
    
    /* Numero pseudo aleatorio (estado del generador) */
    long int x;
};

void initRng(long int seed, long int a, long int c, long int m,
             struct Rng *rng);
long int setRng(struct Rng *rng);

double getRng(struct Rng *rng);
double getR(double x, double y, double z);

// Parametros del problema
const int nreg = 3; // numero de regiones (sin considerar region anterior y
                    // posterior a la placa)

// Los siguientes parametros se entregan para cada region.
const double d[] = {0.5, 0.1, 0.1, 0.1};    // espesor de las regiones del
                                        // atenuador (cm)
const double sigma_t[] = {0.05, 0.1, 10.0, 100.0}; // seccion eficaz total (cm^-1)
const double sigma_a[] = {0.005, 0.01, 0.1, 10.0};//seccion eficaz de absorcion (cm^-1)
const int nhist = 10000; // numero de historias
const int nbatch = 10;  // grupos estadisticos
const int ptracks = 1;  // guardar posicion particulas

int main(int argc, const char * argv[]) {
    // Inicio del programa.
    clock_t start = clock();
    
    // Archivo para guardar las trayectorias
    FILE *fp;
    
    if (ptracks) {
        fp = fopen("ptracks.lst","w");
        
        // Guardamos el numero de historias que seran registradas.
        fprintf(fp, "%d\n", nhist);
    }
    
    // Primero declaramos la estructura del generador y lo inicializamos con
    // una semilla concreta.
    struct Rng rng;
    
    // Parametros sugeridos de Numerical Recipes in C.
    initRng(0, 1664525, 1013904223, 4294967296, &rng);
    
    // Fronteras de las regiones del atenuador. Si hay tres regiones de
    // atenuador significa que tenemos cuatro bordes en la geometria.
    double *zbounds = (double*) malloc((nreg+1)*sizeof(double));
    
    // Registro de las particulas absorbidasmybib en cada region. +2 por la region
    // anterior y posterior a la placa.
    double *pdep = (double*) malloc((nreg+2)*sizeof(double));
    double *apdep = (double*) malloc((nreg+2)*sizeof(double));
    double *apdep2 = (double*) malloc((nreg+2)*sizeof(double));
    memset(pdep, 0.0, (nreg+2)*sizeof(double));
    memset(apdep, 0.0, (nreg+2)*sizeof(double));
    memset(apdep2, 0.0, (nreg+2)*sizeof(double));

    // Inicializar geometria y arreglo de scoring
    zbounds[0] = d[0];
    for (int i=0; i<nreg; i++) {
        zbounds[i+1] = zbounds[i] + d[i+1];
    }
    
    // Inicio de la simulacion
    
    double x, y, z;             // posicion de la particula
    double u, v, w;             // direccion de la particula
    double theta0, phi0;        // angulos de dispersion
    double xnew, ynew, znew;    // posiciones después de aplicar el camino libre medio
    double unew, vnew, wnew;    // direccion despues de la dispersion
    double b, c;             // coeficientes de una ecuación cuadrática
    double s;                   // sqrt(1-w^2)
    
    int ir;  // region actual de la particula.
                // 0: region anterior a la placa
                // 1: primera region del atenuador
                // 2: segunda region del atenuador
                // ...
                // nreg: ultima region del atenuador
                // nreg+1: region posterior al atenuador
    
    double tstep; // distancia a la siguiente interaccion
    double rnno;
    
    int ptrans;     // indica si la particula debe seguir transportandose antes
                    // de interactuar.
    
    int idisc;      // descarte de la particula si escapa del atenuador
    int irnew;      // indice de la nueva region de la particula
    
    int nperbatch = nhist/nbatch;
    
    for (int ibatch=0; ibatch<nbatch; ibatch++) {
        
        for (int ihist=0; ihist<nperbatch; ihist++) {

            // Radio del núcleo
            double R = zbounds[0];

            // Elegimos el punto en el núcleo aleatoriamente (coord esféricas)
            double r = getRng(&rng) * R;
            double theta = getRng(&rng) * M_PI;
            double phi = getRng(&rng) * M_PI * 2;

            
            // Pasamos a cartesianas
            x = r*sin(theta)*cos(phi); y = r*sin(theta)*sin(phi); z = r*cos(theta);  
            // Elegimos la dirección (componente entre -1..1)
            double u0 = getRng(&rng)*2 - 1; 
            double v0 = getRng(&rng)*2 - 1; 
            double w0 = getRng(&rng)*2 - 1; 

            s = getR(u0, v0, w0);
            u = u0/s; v = v0/s; w = w0/s;

            ir = 0;                               // en el núcleo
            idisc = 0; ptrans=1;
            
            if (ptracks) {
                // Marcamos el inicio de cada historia con un caracter especial
                fprintf(fp, "%s\n", "&");
            }
            
            while (1) {
                
                // Comenzamos el transporte del foton.
                do {
                    // Primero debemos muestrear el camino libre.
                    rnno = getRng(&rng);
                    tstep = -(1.0/sigma_t[ir-1])*log(1.0 - rnno);
                    
                    // Guardamos la region actual de la particula.
                    irnew = ir;

                    // Calculamos a donde debe parar la partícula si se mueve la distancia
                    s = getR(u, v, w);
                    xnew = x + (u/s)*tstep; ynew = y + (v/s)*tstep; znew = z + (w/s)*tstep;

                    // Hacemos que R sea cero de tal manera de realizar un control de flujo
                    // modificando este valor después si es que se cambia de región
                    R = 0;

                    // Dado que tenemos a donde debería parar la partícula, podemos revisar
                    // si es que se sale del cascarón revisando la distancia del nuevo punto al
                    // centro de los cascarones
                    if (ir != nreg+1) { // Si estamos en los cascarones o el núcleo

                        if (getR(xnew, ynew, znew) >= zbounds[ir]) { // Revisa si hay que cambiar al cascarón siguiente
                            R = zbounds[ir]; 
                            irnew += 1;
                        }

                        else if (ir != 0) { // Si estamos en un cascarón
                            if (getR(xnew, ynew, znew) <= zbounds[ir-1]) { // Revisa si hay que cambiar al cascarón anterior
                                R = zbounds[ir-1];
                                irnew -= 1;
                            }
                        }
                        
                    }

                    else { // Si estamos en la región de afuera 
                        if (getR(xnew, ynew, znew) >= zbounds[ir-1]) { // Nos estamos saliendo de 
                            idisc = 1; // Nos salimos de la geometría
                        }
                    }

                    if (R) { // Si nos cambiamos de región tenemos que calcular la distancia hasta la frontera

                        /* Acá tenemos que resolver una ecuación cuadrática que se consigue al intersectar la línea

                                (x, y, z) = (x0, y0, z0) + tstep * (u, v, w)    con tstep = [0, inf)

                            Con la esfera

                                x^2 + y^2 + z^2 = R^2

                            Obtenemos la siguiente ecuación cuadrática para tstep

                                tstep^2 *(u^2 + v^2 + y^2) + tstep*2*(u*x0 + v*y0+ w*z0) + (x0^2 + y0^2 + z0^2 - R^2) = 0

                                a * tstep^2 + b * tstep + c = 0

                            La solución positiva será (dado que (u^2 + v^2 + y^2) = 1)

                                tstep = 1/2 * (-b + sqrt(b^2 - 4*c ))

                            Que es lo que hacemos justamente*/
                        // a = 1.0 por lo que no contribuye en la resolución de la ecuación cuadrática

                        b = 2*(u*x+v*y+w*z);
                        c = x*x+y*y+z*z - R*R;
                        tstep = (-b + sqrt(b*b-4*c)) / 2;
                    }
                    
 
                    if (idisc == 1) {
                        // La particula abandona la geometria.
                        break;
                    }
                    
                    if (ptracks) {
                        // Guardamos la posicion de la particula antes del
                        // transporte.
                        fprintf(fp, "%10.5f,%10.5f,%10.5f,", x, y, z);
                    }
                    
                    // Transportamos la particula.
                    x += u*tstep;
                    y += v*tstep;
                    z += w*tstep;

                    
                    if (ptracks) {
                        // Guardamos la posicion de la particula despues del
                        // transporte.
                        fprintf(fp, "%10.5f,%10.5f,%10.5f\n", x, y, z);
                    }
                    
                    // Si llegamos a este punto y la particula no ha cambiado de
                    // region, significa que debe interactuar.
                    if (ir == irnew) {
                        ptrans = 0;
                    }
                    else {
                        // Actualizamos la region de la particula y continuamos
                        // el transporte
                        ir = irnew;
                    }
                    
                } while (ptrans);
                
                if (idisc == 1) {
                    // La particula fue descartada. La depositamos y finalizamos
                    // su rastreo.
                    ir = irnew;
                    pdep[ir] += 1;
                    break;
                }
                
                // Debemos seleccionar una interaccion
                rnno = getRng(&rng);
                if (rnno <= sigma_a[ir-1]/sigma_t[ir-1]) {
                    // La particula es absorbida
                    pdep[ir] += 1.0;
                    break;
                } else {
                    // La particula es dispersada, debemos muestrear el angulo
                    // de dispersion. Muestreo isotropico.
                    rnno = getRng(&rng);
                    theta0 = acos(2.0*rnno - 1.0);
                    rnno = getRng(&rng);
                    phi0 = 2.0*M_PI*rnno;
                    
                    // Calculamos la nueva direccion de la particula. Si el
                    // cambio del angulo polar es muy pequeño, w ~ 1 y s ~ 0,
                    // por lo tanto agregamos un caso especial para manejar
                    // esa situacion.
                    if (pow(u, 2.0) + pow(v, 2.0) < 1.0E-20) {
                        unew = sin(theta0)*cos(phi0);
                        vnew = sin(theta0)*sin(phi0);
                        wnew = w*cos(theta0);
                    } else {
                        s = sqrt(1 - pow(w, 2.0));
                        unew = (w*u/s)*sin(theta0)*cos(phi0) - (v/s)*sin(theta0)*sin(phi0) + u*cos(theta0);
                        vnew = (w*v/s)*sin(theta0)*cos(phi0) + (u/s)*sin(theta0)*sin(phi0) + v*cos(theta0);
                        wnew = -s*(sin(theta0)*cos(phi0)) + w*cos(theta0);
                    }
                    
                    // Actualizamos la direccion de la particula.
                    u = unew;
                    v = vnew;
                    w = wnew;
                    
                    // La particula debe ser nuevamente transportada
                    ptrans = 1;
                }
            }
        }
        // AL finalizar el grupo estadistico acumulamos los resultados en los
        // arreglos respectivos.
        for (int i=0; i<nreg+2; i++) {
            apdep[i] += pdep[i];
            apdep2[i] += pow(pdep[i], 2.0);
        }
        
        // Anulamos pdep y comenzamos el nuevo grupo estadistico.
        memset(pdep, 0.0, (nreg+2)*sizeof(double));
    }
    
    // Analisis de resultados. Obtenemos el promedio y error asociado a las
    // fracciones de reflexion, absorcion y transmision.
    double *updep = (double*) malloc((nreg+2)*sizeof(double));
    for (int i=0; i<nreg+2; i++) {
        pdep[i] = apdep[i]/nbatch;
        
        updep[i] = (apdep2[i] - nbatch*pow(pdep[i], 2.0));
        updep[i] /= (nbatch*(nbatch-1));
        updep[i] = sqrt(updep[i]);
    }
    
    double pdep_abs = 0.0;
    double updep_abs = 0.0;
    for (int i=0; i<=nreg; i++) {
        pdep_abs += pdep[i];
        updep_abs += pow(updep[i], 2.0);
    }
    updep_abs = sqrt(updep_abs);
    

    printf("Probabilidad de absorcion : %6.5f +/- %6.5f \n",
           pdep_abs/nperbatch, updep_abs/nperbatch);
    printf("Probabilidad de transmision : %6.5f +/- %6.5f \n",
           pdep[nreg+1]/nperbatch, updep[nreg+1]/nperbatch);
    
    // Limpieza
    free(pdep);
    free(zbounds);
    if (ptracks) {
        fclose(fp);
    }
    
    // Calculamos el tiempo transcurrido y reportamos el resultado.
    printf("\nTiempo transcurrido = %10.2f (s)\n",
           (double)(clock()-start)/(double)CLOCKS_PER_SEC);
    
    return 0;
}

void initRng(long int seed, long int a, long int c, long int m,
             struct Rng *rng) {
    /* Inicializacion del generador */
    rng->a = a;
    rng->c = c;
    rng->m = m;
    
    rng->x = seed;
    
    printf("Generador pseudo-aleatorio LCG con parametros:\n");
    printf("\t a = %ld\n", rng->a);
    printf("\t c = %ld\n", rng->c);
    printf("\t m = %ld\n", rng->m);
    printf("\t seed = %ld\n", rng->x);
    
    return;
}

long int setRng(struct Rng *rng) {
    /* Esta funcion entrega un numero aleatorio en [0, m) */
    long int rnno = 0;
    
    rnno = (rng->a*rng->x + rng->c) % rng->m;
    
    /* Actualizacion del estado del generador */
    rng->x = rnno;
    
    return rnno;
}

double getRng(struct Rng *rng) {
    /* Esta funcion entrega un numero pseudo aleatorio en [0, 1) */
    double rnno = 0.0;
    
    rnno = (double)setRng(rng)/rng->m;
    
    return rnno;
}

double getR(double x, double y, double z) {
    return sqrt(x*x+y*y+z*z);
}
