# Cálculo Numérico para Ingeniería Química y de Alimentos

[![Pruebas](https://github.com/jhabriel/calculo_numerico_iq_ia/actions/workflows/tests.yml/badge.svg)](https://github.com/jhabriel/calculo_numerico_iq_ia/actions/workflows/tests.yml)

> [!IMPORTANT]
> **Antes de usar las funciones o los ejemplos, ejecuta el archivo `iniciar.m`.**
>
> MATLAB y Octave olvidan el path cada vez que se cierran, así que al abrir el programa:
>
> 1. Descarga este repositorio como un `.zip` y descomprímelo.
> 2. Ingresa a OCTAVE/MATLAB y sitúate en la ventana de comandos.
> 3. Sitúate en la carpeta descomprimida (tienes que poder ver que contiene `iniciar.m`). Para verificar esto, escribe `ls` en la ventana de comandos y fíjate si aparece `iniciar.m`.
> 4. En la ventana de comandos, escribe y ejecuta:
>
>    ```matlab
>    iniciar
>    ```
>
> A partir de ahí, todas las funciones de este repositorio podrán ser usadas como 
> funciones estándar, sin importar de donde se ejecuten tus scripts.
>
> Si no lo haces, verás errores como `'biseccion' undefined`: no es un error de tu
> código, solo falta ejecutar `iniciar`.

Implementaciones didácticas de métodos de cálculo numérico en **MATLAB/Octave**,
pensadas para estudiantes de Ingeniería Química e Ingeniería de Alimentos de la
Facultad de Ciencias Químicas, Universidad Nacional de Asunción. El código prioriza
la claridad: cada método de `src/` está escrito "a mano" (sin llamar a `\`, `lu`,
`interp1`, etc.); los ejemplos solo usan funciones nativas como referencia, para
graficar o para evaluar (`polyval`). Los métodos se validan con pruebas automáticas
y los ejemplos con la solución exacta, además de una prueba de humo que los ejecuta.

**Requisitos:** probado con GNU Octave 9.2. Está escrito para funcionar también en
MATLAB, pero esa combinación todavía no se ha probado.

## Estructura

| Carpeta | Contenido |
|---|---|
| `src/01_ecuaciones_no_lineales` | `biseccion`, `newton` |
| `src/02_sistemas_lineales` | `sustitucion_adelante`, `sustitucion_atras`, `gauss_naive`, `gauss_pivoteo`, `lu_factorizacion`, `resolver_lu`, `jacobi`, `gauss_seidel` |
| `src/03_optimizacion` | `seccion_aurea` (1D, sin restricciones) |
| `src/04_interpolacion_regresion` | `lagrange` |
| `tests/` | Pruebas automáticas |
| `ejemplos/` | Ejemplos con problemas de ingeniería química y alimentos |

## Uso

Desde la raíz del repositorio, después de ejecutar `iniciar`:

```matlab
iniciar                                   % una sola vez por sesión (ver aviso al inicio)

f  = @(x) x.^2 - 2;
x  = biseccion(f, 0, 2, 1e-10);
x2 = newton(f, @(x) 2*x, 1);

A = [4 -2 1; -2 4 -2; 1 -2 4];  b = [11; -16; 17];
x = resolver_lu(A, b);
```

Cada función documenta sus argumentos con `help nombre_funcion`.

### Ejemplos

Un script por sección en `ejemplos/`, listos para ejecutar **después de `iniciar`** (imprimen
resultados y, en varios casos, grafican):

| Script | Tema |
|---|---|
| `ej_01_ecuacion_no_lineal` | Bisección y Newton: x³ − 2x − 5 = 0 |
| `ej_02_ax_b_gauss_y_lu_simple` | Gauss y LU para una A, b fijas |
| `ej_03_ax_b_gauss_y_lu_pivoteo_parcial` | Pivote nulo y pivote diminuto |
| `ej_04_ax_b_jacobi_y_gauss_seidel` | Iterativos en sistema diagonalmente dominante |
| `ej_05_sistema_no_lineal_newton` | Newton para sistemas (armado en el script), sistema 3×3 clásico de Burden & Faires |
| `ej_06_optimizacion_seccion_dorada_y_newton` | Máximo de 2 sin x − x²/10 |
| `ej_07_interpolacion_lagrange_y_vandermonde` | Presión de vapor del agua: Vandermonde (armado en el script) y Lagrange |
| `ej_08_ajuste_de_curvas` | Regresión simple, polinómica y múltiple con ecuaciones normales |
| `ej_09_diferenciacion_numerica` | 1ª y 2ª derivada con diferencias finitas, orden de convergencia (fórmulas en el script) |
| `ej_10_integracion_numerica` | Trapecio, Simpson 1/3, 3/8 y Boole (fórmulas de una aplicación) con error relativo verdadero |
| `ej_11_edo_lineal` | EDO lineal (reacción de 1er orden): Euler explícito, implícito y Heun |
| `ej_12_edo_no_lineal` | EDO no lineal (reacción de 2º orden): Euler implícito con Newton |
| `ej_13_edo_segundo_orden_reduccion` | Oscilador amortiguado: reducción de orden (u1, u2) |
| `ej_14_sistema_edo_lineal` | Sistema lineal: modelo compartimental |

## Pruebas (ignorar esto si eres alumno)

```matlab
run('tests/run_all_tests.m')
```

o desde la terminal con Octave:

```bash
octave --no-gui --eval "run('tests/run_all_tests.m')"
```

## Convenciones

- Nombres de funciones en español, `snake_case`, un método por archivo.
- Los scripts no llevan `addpath`: se asume que se ejecutó `iniciar` al abrir MATLAB/Octave.
- Compatible con MATLAB y Octave (sin sintaxis exclusiva de uno u otro).
- Diferenciación, integración y EDOs no tienen funciones en `src/`: los métodos se escriben directamente en los scripts (`ej_09` a `ej_14`). En EDOs la historia se guarda en una matriz de (n+1) × m: una fila por instante y una columna por variable.

## Agregar un método nuevo (Ignorar esto si eres alumno)

1. Crear `src/<tema>/nombre.m` con encabezado de ayuda (uso, argumentos, salidas).
2. Agregar pruebas en `tests/test_<tema>.m` (comparar contra solución exacta o contra el operador nativo).
3. Listarlo en esta tabla y ejecutar `run_all_tests`.

## Licencia

Copyright © 2026 Equipo Docente - Cálculo Numérico - FCQ-UNA

Este proyecto se distribuye bajo la licencia **GNU General Public License v3.0**
(GPL-3.0). Puedes usar, estudiar, modificar y redistribuir el código, siempre que
las obras derivadas se distribuyan bajo la misma licencia. Ver el archivo
[`LICENSE`](LICENSE) para el texto completo.
