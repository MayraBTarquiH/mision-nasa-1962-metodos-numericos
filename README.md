# 🚀 Misión espacial NASA 1962 — Simulación del despegue de un cohete

## Descripción del proyecto

Este proyecto consiste en una aplicación web interactiva desarrollada con Python y Streamlit para simular el movimiento vertical de un cohete considerando el empuje, la gravedad y la resistencia lineal del aire.

Se implementan tres métodos numéricos: Euler, Heun y Runge-Kutta de cuarto orden (RK4). Sus resultados se comparan con la solución analítica para evaluar la precisión de cada método.

**Aplicación web:** https://mision-nasa-1962-metodos-numericos-nujrhzunjrefwq3xd89gpn.streamlit.app/

## Modelo matemático

El movimiento del cohete se representa mediante la segunda ley de Newton:

\[
m\frac{dv}{dt}=T-mg-kv
\]

El sistema de ecuaciones diferenciales es:

\[
\frac{dv}{dt}=a-bv
\]

\[
\frac{dy}{dt}=v
\]

Donde:

- \(m\): masa del cohete (kg).
- \(T\): fuerza de empuje (N).
- \(g\): aceleración de la gravedad (m/s²).
- \(k\): coeficiente de resistencia lineal (kg/s).
- \(v\): velocidad vertical (m/s).
- \(y\): altura (m).
- \(a=T/m-g\).
- \(b=k/m\).

## Métodos numéricos utilizados

La aplicación implementa:

1. **Método de Euler:** utiliza la pendiente al inicio de cada intervalo.
2. **Método de Heun:** utiliza el promedio de dos pendientes para mejorar la aproximación.
3. **Método Runge-Kutta de cuarto orden (RK4):** utiliza cuatro evaluaciones de la pendiente para obtener una mayor precisión.

Los resultados numéricos se comparan con la solución analítica del modelo.

## Parámetros iniciales

Los valores predeterminados de la simulación son:

| Parámetro | Valor |
|---|---:|
| Masa | 1000 kg |
| Empuje | 15000 N |
| Gravedad | 9.81 m/s² |
| Resistencia lineal | 100 kg/s |
| Velocidad inicial | 0 m/s |
| Altura inicial | 0 m |
| Tiempo final | 20 s |
| Paso de integración | 2 s |

La aplicación permite modificar estos parámetros y las tolerancias de error desde la barra lateral.

## Resultados obtenidos

Con los parámetros iniciales y un paso de integración de 2 segundos, se obtuvieron los siguientes resultados:

| Método | Velocidad final (m/s) | Altura final (m) |
|---|---:|---:|
| Euler | 46.3273 | 574.7272 |
| Heun | 44.7664 | 590.3355 |
| RK4 | 44.8759 | 589.2412 |
| Solución analítica | 44.8761 | 589.2390 |

Al reducir el paso de integración a 0.5 segundos, se observa una disminución de los errores y una mejora en la precisión de los tres métodos.

## Funcionalidades de la aplicación

La aplicación cuenta con cuatro pestañas:

- **Resultados:** presenta los valores finales de velocidad y altura, los errores y el cumplimiento de las tolerancias.
- **Tabla de iteraciones:** muestra los resultados de cada instante y permite descargar la información en formato Excel.
- **Gráficas:** compara las aproximaciones numéricas con la solución analítica y muestra los errores absolutos.
- **Efecto del paso h:** permite analizar cómo influye el tamaño del paso de integración en la precisión de los métodos.

## Tecnologías utilizadas

- Python
- Streamlit
- NumPy
- Pandas
- Matplotlib
- OpenPyXL

## Instalación y ejecución

Para ejecutar el proyecto localmente:

1. Clonar el repositorio:

```bash
git clone https://github.com/MayraBTarquiH/mision-nasa-1962-metodos-numericos.git
```

2. Ingresar a la carpeta:

```bash
cd mision-nasa-1962-metodos-numericos
```

3. Instalar las dependencias:

```bash
pip install -r requirements.txt
```

4. Ejecutar la aplicación:

```bash
streamlit run app.py
```

## Conclusión

La aplicación permite observar cómo los métodos numéricos aproximan la solución de un sistema de ecuaciones diferenciales. Los resultados muestran que, al reducir el paso de integración, disminuyen los errores de aproximación. Asimismo, RK4 proporciona una precisión superior a Euler y Heun para los pasos analizados.
