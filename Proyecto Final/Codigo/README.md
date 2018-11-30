# __Proyecto Final - Modelo Monte Carlo Del Transporte De Partículas__

Se compila y ejecuta el archivo `main.c` para obtener los resultados por consola. Si se coloca la opción `ptracks = 1` dentro de `main.c` se creará un archivo `ptracks.lst` que listará las historias simuladas, recomendamos hacerlo para no más de 10 000 historias simuladas. Con este archivo se puede hacer uso del script `check.py` que eligirá una historia al azar presente en `ptracks.lst` y la graficará en un gráfico 3D en conjunto a la geometría del problema.

Dentro de `main.c` se definen las constantes relevantes al problema (condiciones iniciales, constantes inherentes al sistema). El código es altamente inspirado en el programa `rnd_exp.c` entregado por el prof. Dörner durante la última clase del curso.