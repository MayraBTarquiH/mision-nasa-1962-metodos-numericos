# 🚀 Misión espacial NASA 1962 — Simulación del despegue de un cohete

**Autora:** [Tu nombre completo]
**Materia:** Métodos Numéricos
**Carrera:** Ingeniería de Sistemas
**Universidad:** [Nombre de tu universidad]
**Docente:** [Nombre de tu docente]
**Fecha:** [Fecha de entrega]

## Descripción del proyecto

Este proyecto consiste en una aplicación web interactiva desarrollada con Python y Streamlit para simular el movimiento vertical de un cohete considerando el empuje, la gravedad y la resistencia lineal del aire.

Se implementan tres métodos numéricos: Euler, Heun y Runge-Kutta de cuarto orden (RK4). Sus resultados se comparan con la solución analítica para evaluar la precisión de cada método.

**Aplicación web:** https://mision-nasa-1962-metodos-numericos-nujrhzunjrefwq3xd89gpn.streamlit.app/

## Modelo matemático

El movimiento del cohete se representa mediante la segunda ley de Newton:

$$
m\frac{dv}{dt}=T-mg-kv
$$

El sistema de ecuaciones diferenciales es:

$$
\frac{dv}{dt}=a-bv
$$

$$
\frac{dy}{dt}=v
$$

Donde:

- $m$: masa del cohete (kg).
- $T$: fuerza de empuje (N).
- $g$: aceleración de la gravedad (m/s²).
- $k$: coeficiente de resistencia lineal (kg/s).
- $v$: velocidad vertical (m/s).
- $y$: altura (m).
- $a=T/m-g$.
- $b=k/m$.

### Solución analítica

Con $v(0)=v_0$ y $y(0)=y_0$:

$$
v(t)=\frac{a}{b}\left(1-e^{-bt}\right)+v_0e^{-bt}
$$

$$
y(t)=y_0+\frac{a}{b}t+\left(v_0-\frac{a}{b}\right)\frac{1-e^{-bt}}{b}
$$

Con los parámetros iniciales y $v_0=y_0=0$ se obtiene $v(t)=51.9\left(1-e^{-0.1t}\right)$ y $y(t)=51.9\,t-519\left(1-e^{-0.1t}\right)$.

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

Con los parámetros iniciales y un paso de integración de 2 segundos, se obtuvieron los siguientes resultados a los 20 s:

| Método | Velocidad final (m/s) | Altura final (m) | Error en velocidad (m/s) | Error en altura (m) |
|---|---:|---:|---:|---:|
| Euler | 46.3273 | 574.7272 | 1.4512 | 14.5118 |
| Heun | 44.7664 | 590.3355 | 0.1097 | 1.0965 |
| RK4 | 44.8759 | 589.2412 | 0.0002 | 0.0022 |
| Solución analítica | 44.8761 | 589.2390 | — | — |

### Cumplimiento de la precisión requerida

El desafío exige un error menor a 0.5 m/s en velocidad y menor a 5 m en altura:

- Con **h = 2 s**, Heun y RK4 cumplen la precisión requerida, mientras que Euler no.
- Euler cumple recién con **h = 0.5 s**, donde los errores son 0.3541 m/s en velocidad y 3.5412 m en altura.

Al reducir el paso de integración se observa una disminución de los errores en los tres métodos. En Euler, al dividir h entre 2 el error baja aproximadamente a la mitad, lo que corresponde a un método de primer orden.

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

La aplicación permite observar cómo los métodos numéricos aproximan la solución de un sistema de ecuaciones diferenciales. Los resultados muestran que, al reducir el paso de integración, disminuyen los errores de aproximación. Asimismo, RK4 proporciona una precisión superior a Euler y Heun para los pasos analizados, y con h = 2 s solo Heun y RK4 satisfacen los requisitos de precisión del desafío.
