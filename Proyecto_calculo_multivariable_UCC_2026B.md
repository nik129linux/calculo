                    Proyecto de Cálculo Multivariable
Regresión lineal y optimización mediante funciones
                de varias variables
1.    Planteamiento del problema
En diferentes áreas de la ciencia y la ingeniería es necesario encontrar relaciones matemáticas
(modelos matemáticos) que permitan describir un conjunto de datos experimentales o simulados.
Una de las aproximaciones más sencillas consiste en buscar una relación lineal de la forma

                                               y = ax + b,                                                 (1)
donde a representa la pendiente y b el intercepto. Sin embargo, determinar los valores de a y
b que producen el mejor ajuste puede formularse como un problema de optimización de una
función de dos variables.
Dado un conjunto de datos

                                  (x1 , y1 ), (x2 , y2 ), . . . , (xn , yn ),                              (2)
se propone utilizar el modelo

                                               y = ax + b                                                  (3)
y definir la función de error cuadrático
                                                 n
                                                ∑︂
                                 S(a, b) =            (yi − axi − b)2 .                                    (4)
                                                i=1

El objetivo será determinar los valores a∗ y b∗ que minimizan S(a, b), utilizando inicialmente
herramientas de Cálculo Multivariable y, posteriormente, herramientas computacionales para
verificar y analizar los resultados.


2.    Pregunta central del proyecto
El trabajo deberá responder, a partir del análisis matemático y computacional, la siguiente
pregunta:

       ¿Cómo se puede utilizar el Cálculo Multivariable para determinar y
       analizar el modelo lineal que mejor representa un conjunto de datos?

La respuesta deberá integrar tres perspectivas:

       análisis analítico   +   interpretación geométrica                 +     análisis computacional .   (5)
El computador no deberá sustituir el desarrollo matemático. Su función será verificar los
resultados obtenidos analíticamente, visualizar las funciones de varias variables y permitir un
análisis más profundo de los datos.
3.    Objetivo general
Aplicar conceptos de Cálculo Multivariable para formular, resolver y analizar un problema
de regresión lineal mediante la construcción y minimización de una función de error de dos
variables, comparando los resultados obtenidos analíticamente de forma manual con los resultados
obtenidos mediante herramientas computacionales.


4.    Objetivos específicos
  1. Analizar un conjunto de datos asociado a un problema real de su elección e identificar una
     posible relación lineal entre las variables involucradas (Aclaración: aquí hay que buscar un
     problema real y recolectar los datos cuyo comportamiento se ajuste a un modelo lineal)

  2. Formular la función de error cuadrático
                                                   n
                                                  ∑︂
                                        S(a, b) =    (yi − axi − b)2 .                        (6)
                                                  i=1


  3. Calcular manualmente las derivadas parciales de S(a, b).

  4. Utilizar el gradiente para determinar las condiciones necesarias para encontrar el punto
     crítico.

  5. Resolver manualmente el sistema de ecuaciones obtenido para determinar a∗ y b∗ .

  6. Calcular y analizar el Hessiano de S(a, b).

  7. Determinar mediante el Hessiano si el punto crítico corresponde a un mínimo.

  8. Implementar computacionalmente el procedimiento de ajuste.

  9. Comparar los resultados analíticos y computacionales.

 10. Construir representaciones gráficas que permitan interpretar geométricamente el proceso
     de minimización.

 11. Analizar la calidad del ajuste y las limitaciones del modelo lineal.


5.    Exploración de los datos
Buscar un problema real y recolectar los datos cuyo comportamiento se ajuste a un modelo
lineal. Inicialmente deberán Analizar gráficamente la relación entre las variables mediante un
diagrama de dispersión, con el propósito de determinar si existe una tendencia aproximadamente
lineal.

Actividades
Los estudiantes deberán:

     organizar los datos en una tabla;

     identificar las variables x y y;
       determinar las unidades de las variables, cuando corresponda;

       construir un gráfico de dispersión de los datos;

       analizar visualmente si existe una posible tendencia aproximadamente lineal.

El gráfico de dispersión deberá permitir discutir si un modelo de la forma

                                                y = ax + b                                     (7)
es razonable para representar los datos.


6.     Construcción del modelo
Se propone inicialmente el modelo lineal

                                               y = ax + b.                                     (8)
Los datos medidos u observados se representan mediante pares (xi , yi ), donde xi es el valor
medido de la variable independiente y yi es el valor observado de la variable dependiente.
Para cada xi , el modelo proporciona un valor predicho

                                          ypred,i = axi + b.                                   (9)
La diferencia entre el valor observado y el valor predicho se denomina error o residuo y se define
como

                                  ei = yi − ypred,i = yi − (axi + b).                         (10)
Como algunos errores pueden ser positivos y otros negativos, se utilizará el error cuadrático
total
                                          n
                                         ∑︂             n
                                                       ∑︂
                                                  2
                             S(a, b) =         (ei ) =    (yi − axi − b)2 .                   (11)
                                         i=1             i=1



7.     Procedimiento analítico
Esta etapa deberá realizarse manualmente antes de utilizar herramientas computacionales (solo
para realizar alguna grafica puede usar computadora).

7.1.     Derivadas parciales
A partir de la función de error
                                                  n
                                                 ∑︂
                                    S(a, b) =          (yi − axi − b)2 ,                      (12)
                                                 i=1

los estudiantes deberán calcular las derivadas parciales respecto de a y b, mostrando el procedi-
miento empleado.
7.2.     Condición de punto crítico
Utilizando las derivadas parciales obtenidas, los estudiantes deberán establecer la condición de
punto crítico

                                          ∇S(a, b) = 0                                        (13)
y obtener el sistema de ecuaciones correspondiente.

7.3.     Determinación de los parámetros
Los estudiantes deberán resolver algebraicamente el sistema obtenido y determinar los valores
de a∗ y b∗ que minimizan la función de error.
El procedimiento deberá realizarse manualmente y presentarse de manera clara en su defensa
del trabajo.

7.4.     Análisis mediante el Hessiano
       Una vez determinado el punto crítico, los estudiantes deberán calcular el Hessiano de la
       función de error.

       Deberán utilizar el criterio del Hessiano para analizar la naturaleza del punto crítico y
       justificar matemáticamente si corresponde a un mínimo.

       En la defensa, deberá presentar el procedimiento y los argumentos utilizados para sustentar
       la conclusión.

7.5.     Interpretación geométrica
S(a, b) puede representarse como una superficie de tres dimensiones, z = S(a, b). Los estudiantes
deberán construir computacionalmente esta superficie y analizar su comportamiento en función
de los parámetros a y b.

7.6.     Superficie de error
Deberán representar gráficamente la superficie de error para un intervalo adecuado de valores
de a y b, e identificar la región en la que se encuentra el mínimo.

7.7.     Trazas y secciones transversales
       Deberán obtener y analizar algunas trazas de la superficie.

       Deberán seleccionar algunas secciones transversales y representarlas gráficamente, expli-
       cando qué información proporcionan sobre el comportamiento de la función de error.

7.8.     Curvas de nivel
       Deberán construir las curvas de nivel de la función S(a, b), y analizar cómo se distribuyen
       alrededor del punto mínimo.

       Deberán representar gráficamente algunas curvas de nivel, señalar sobre ellas el punto
       (a∗ , b∗ ) obtenido y analizar su correspondencia con la región de mínimo error.
A partir de la superficie, las trazas, las secciones transversales y las curvas de nivel, deberán
explicar geométricamente cómo se comporta la función de error y cómo puede identificarse el
punto correspondiente a los parámetros que producen el mejor ajuste.


8.    Verificación computacional
Después de completar el procedimiento analítico, implementar el problema mediante una
herramienta computacional.
El programa deberá:

  1. cargar los datos;

  2. calcular a∗ y b∗ ;

  3. calcular S(a∗ , b∗ );

  4. Construir un gráfico superpuesto que muestre los datos observados y la recta correspon-
     diente al modelo ajustado. Representar también rectas obtenidas con otros valores (a, b),
     y comparar su ajuste con el modelo que minimiza la función de error.

  5. construir la superficie S(a, b); y representar las curvas de nivel, trazas y secciones trans-
     versales y el minimo (a*, b*);

  6. calcular el Hessiano;

  7. verificar numéricamente las propiedades del Hessiano.

Aclaracion: en éste trabajo, no deberán limitarse a utilizar una función incorporada de
regresión lineal. Primero deberán implementar el procedimiento que dedujeron analíticamente.
Posteriormente podrán utilizar una función de regresión incorporada como método adicional de
verificación y comparacion (es decir, deben programar el metodo de ajuste, despues comparan
con alguno que haya en las librerias del lenguaje usado).


9.    Comparación entre el procedimiento analítico y compu-
      tacional
Los resultados obtenidos mediante el procedimiento analítico realizado manualmente deberán
compararse con los resultados obtenidos mediante la implementación computacional desarrollada
por los estudiantes y con los resultados proporcionados por las funciones disponibles en las
bibliotecas del lenguje utilizado. Deberán comparar los valores de a∗ , b∗ y Smı́n obtenidos por
cada método, y analizar las posibles diferencias.
Se deberá construir una tabla similar a:

     Parámetro      Analítico   Código propio   Biblioteca de Python    Diferencia máxima
          a           aman          acod                alib                    ∆a
          b           bman          bcod                blib                    ∆b
         Smı́n        Sman          Scod                Slib                    ∆S
10.     Análisis de la recta de ajuste
Una vez obtenidos los parámetros óptimos se construirá la recta

                                           y = a∗ x + b∗ .                                  (14)
La recta deberá superponerse al conjunto de datos y los estudiantes deberán analizar:

      qué tan bien representa la recta la tendencia general de los datos;

      cuáles puntos presentan mayor desviación respecto a la recta;

      si existen puntos alejados de la tendencia general que puedan influir en el ajuste.


11.     Análisis de los residuos
Para cada dato deberán calcularse los residuos

                                       ei = yi − (a∗ xi + b∗ ).                             (15)
Deberán construirse:

      un gráfico de los residuos ei en función de xi ;

      un histograma de los residuos.

A partir de estas representaciones, analizar:

      si los residuos están aproximadamente centrados alrededor de cero;

      si presentan alguna tendencia o patrón;

      si aumenta o disminuye su dispersión con xi ;

      si existen valores atípicos;

      si existe alguna región donde el modelo lineal presenta desviaciones sistemáticas.

Aclaración: Se considerará que existe una desviación sistemática cuando los residuos presenten
algún patrón o tendencia en función de xi , en lugar de distribuirse aproximadamente alrededor
de cero sin una estructura evidente. Por ejemplo, puede ocurrir que los residuos sean predo-
minantemente negativos en una región de los datos y positivos en otra, o que presenten un
patrón curvo. Este tipo de comportamiento puede indicar que el modelo lineal no representa
adecuadamente alguna región del conjunto de datos.
El análisis deberá ser interpretativo y permitir determinar si el comportamiento de los residuos
es compatible con un modelo lineal adecuado para los datos.
12.     Evaluación de la calidad del modelo
Calcular al menos una medida cuantitativa de la calidad del ajuste. Por ejemplo, la raíz del
error cuadrático medio:
                                           ⌜
                                           ⃓ n
                                           ⃓ 1 ∑︂
                                 RMSE = ⎷           e2 .                                (16)
                                              n i=1 i

También podrán calcular el coeficiente de determinación:

                                                 (yi − ŷi )2
                                             ∑︁
                                   R2 = 1 − ∑︁i             2
                                                              .                            (17)
                                                i (yi − ȳ)
Explicar qué representa cada medida para el conjunto de datos estudiado.
No será suficiente afirmar que un determinado valor de R2 es “bueno” o “malo”. Deberán
relacionarlo con la distribución observada de los datos y con el propósito del modelo.


13.     Análisis adicional
Cada grupo deberá realizar al menos uno de los siguientes análisis.

13.1.     Opción A: Sensibilidad de los parámetros
Modificar ligeramente algunos de los datos y estudiar cómo cambian

                                          a∗       y    b∗ .                               (18)
Discutir qué tan sensible es el modelo a pequeñas modificaciones de los datos.

13.2.     Opción B: Influencia de un dato
Eliminar temporalmente uno de los puntos y volver a calcular la regresión.
Deberán analizar cómo cambia la recta de ajuste y determinar si existen observaciones con una
influencia particularmente importante.

13.3.     Opción C: Comparación con otro modelo
Comparar el modelo lineal

                                               y = ax + b                                  (19)
con otro modelo sencillo apropiado para los datos, por ejemplo,

                                         y = ax2 + bx + c.                                 (20)
Deberán comparar los errores y discutir si el modelo más complejo realmente proporciona una
representación significativamente diferente de los datos.

13.4.     Opción D: Efecto del número de datos
Realizar el ajuste utilizando diferentes cantidades de datos y analizar cómo cambian los paráme-
tros obtenidos.
