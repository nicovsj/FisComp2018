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

// Parametros del problema
const int nreg = 3; // numero de regiones (sin considerar region anterior y
                    // posterior a la placa)

// Los siguientes parametros se entregan para cada region.
const double d[] = {0.10, 0.10, 0.10};    // espesor de las regiones del
                                        // atenuador (cm)
const double sigma_t[] = {0.1, 10.0, 100.0}; // seccion eficaz total (cm^-1)
const double sigma_a[] = {0.01, 0.1, 10.0};//seccion eficaz de absorcion (cm^-1)
const int nhist = 1000000; // numero de historias
const int nbatch = 10;  // grupos estadisticos
const int ptracks = 0;  // guardar posicion particulas

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
    zbounds[0] = 0.0;
    for (int i=0; i<nreg; i++) {
        zbounds[i+1] = zbounds[i] + d[i];
    }
    
    // Inicio de la simulacion
    
    double x, y, z;             // posicion de la particula
    double u, v, w;             // direccion de la particula
    double theta0, phi0;        // angulos de dispersion
    double unew, vnew, wnew;    // direccion despues de la dispersion
    double s;                   // sqrt(1-w^2)
    
    int ir;  // region actual de la particula.
                // 0: region anterior a la placa
                // 1: primera region del atenuador
                // 2: segunda region del atenuador
                // ...
                // nreg: ultima region del atenuador
                // nreg+1: region posterior al atenuador
    
    double tstep; // distancia a la siguiente interaccion
    double dist;  // distancia hacia uno de los bordes
    double rnno;
    
    int ptrans;     // indica si la particula debe seguir transportandose antes
                    // de interactuar.
    
    int idisc;      // descarte de la particula si escapa del atenuador
    int irnew;      // indice de la nueva region de la particula
    
    int nperbatch = nhist/nbatch;
    
    for (int ibatch=0; ibatch<nbatch; ibatch++) {
        
        for (int ihist=0; ihist<nperbatch; ihist++) {
            
            // Inicializacion de la particula
            x = 0.0; y = 0.0; z = 0.0;       // i.e., frente al atenuador
            u = 0.0; v= 0.0; w = 1.0;             // perpendicular a la placa
            ir = 1;                               // en la primera region de la placa
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
                    
                    // Ahora debemos calcular la distancia a la siguiente
                    // frontera en la direccion de propagacion de la particula.
                    // Notese que al ser una placa 1D nos interesa solamente
                    // analizar la direccion z.
                    
                    if(w < 0.0) {
                        // La particula viaja hacia el borde anterior de la
                        // placa.
                        dist = -(zbounds[ir-1] - z)/w;
                        
                        // Ahora revisamos si la particula cambia de region.
                        if (dist < tstep) {
                            tstep = dist;
                            if (irnew != 0) {
                                irnew -= 1;
                            } else {
                                // Indicamos que la particula abandona la
                                // geometria.
                                idisc = 1;
                            }
                        }
                    }
                    else if(w > 0.0) 
{                        // La particula viaja hacia el borde posterior de la
                        // placa.
                        dist = (zbounds[ir] - z)/w;
    
                        // Ahora revisamos si la particula cambia de region.
                        if (dist < tstep) {
                            tstep = dist;
                            if (irnew != nreg+1) {
                                irnew += 1;
                            } else {
                                // Indicamos que la particula abandona la
                                // geometria.
                                idisc = 1;
                            }
                        }
                    }
                    
                    if (idisc == 1) {
                        // La particula abandona la geometria.
                        break;
                    }
                    
                    if (ptracks) {
                        // Guardamos la posicion de la particula antes del
                        // transporte.
                        fprintf(fp, "%10.5f %10.5f %10.5f ", x, y, z);
                    }
                    
                    // Transportamos la particula.
                    x += u*tstep;
                    y += v*tstep;
                    z += w*tstep;
                    
                    if (ptracks) {
                        // Guardamos la posicion de la particula despues del
                        // transporte.
                        fprintf(fp, "%10.5f %10.5f %10.5f\n", x, y, z);
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
    for (int i=1; i<=nreg; i++) {
        pdep_abs += pdep[i];
        updep_abs += pow(updep[i], 2.0);
    }
    updep_abs = sqrt(updep_abs);
    
    printf("Probabilidad de reflexion : %6.5f +/- %6.5f \n",
           pdep[0]/nperbatch, updep[0]/nperbatch);
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
