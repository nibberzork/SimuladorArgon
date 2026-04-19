# Simulador de dinámica molecular para argón

Este paquete implementa un simulador de dinámica molecular para un sistema de átomos de argón en unidades reducidas de Lennard-Jones. El núcleo numérico está desarrollado en C++ y se expone a Python mediante `pybind11`. Sobre dicha base se proporciona una interfaz de alto nivel que devuelve las magnitudes termodinámicas muestreadas como un `DataFrame` de `pandas` y, de forma opcional, un array de `NumPy` con los módulos de velocidad.

## Características

- Integración temporal mediante el algoritmo velocity-Verlet.
- Interacciones de Lennard-Jones truncadas a un radio de corte fijo.
- Cálculo de fuerzas optimizado mediante listas de celdas (`cell lists`) cuando la discretización espacial lo permite.
- Condiciones de contorno periódicas e imagen mínima.
- Corrección opcional del momento lineal total para evitar la deriva del centro de masas.
- Reescalado opcional de velocidades durante la fase de equilibrado.
- Muestreo de temperatura, presión y energías.
- Exportación opcional de resultados termodinámicos en formato CSV y de módulos de velocidad en formato NPY.
- Funciones de análisis para energía total, resumen termodinámico e histogramas de velocidades.

## Estructura del proyecto

```text
SimuladorArgon/
├── simulador_dm/
│   ├── __init__.py
│   ├── _wrapper_simulador.py
│   ├── analisis.py
│   ├── _simulador.hpp
│   ├── _simulador.cpp
│   └── bindings.cpp
├── pyproject.toml
├── setup.py
└── CMakeLists.txt
```

## Requisitos

- Python 3.7 o superior
- Un compilador compatible con C++17
- `pip`
- `git`

Dependencias Python del proyecto:

- `pybind11`
- `numpy`
- `pandas`

## Obtención del código fuente

Para trabajar con el proyecto, el primer paso consiste en obtener una copia local del repositorio. La forma recomendada de hacerlo es mediante `git`.

1. Abra una terminal en la carpeta donde desee guardar el proyecto.
2. Ejecute el siguiente comando para clonar el repositorio:

```bash
git clone https://github.com/nibberzork/SimuladorArgon.git
```

3. Acceda a la carpeta descargada:

```bash
cd SimuladorArgon
```

Si el proyecto ya se encuentra descargado en su equipo, este paso puede omitirse.

## Instalación

Una vez situado en la raíz del proyecto, la vía recomendada de instalación consiste en ejecutar:

```bash
pip install .
```

Si se desea trabajar en modo desarrollo, de forma que los cambios realizados en el código fuente se reflejen sin reinstalar el paquete, puede utilizarse:

```bash
pip install -e .
```

En sistemas Windows puede ser necesario disponer de las herramientas de compilación de Visual Studio para C++.

## Módulos públicos

El paquete `simulador_dm` exporta directamente:

- `Simulador`: alias público de la clase `WraperSimulador`.
- `graficar_energia`
- `graficar_resumen_termodinamico`
- `graficar_histograma_velocidades`

Importación típica:

```python
from simulador_dm import (
    Simulador,
    graficar_energia,
    graficar_resumen_termodinamico,
    graficar_histograma_velocidades,
)
```

## Uso básico

```python
from simulador_dm import Simulador

sim = Simulador(
    particulas_por_lado=8,
    densidad_reducida=0.84,
    paso_tiempo=0.005,
    temp_objetivo=1.002,
)

df, velocidades = sim.ejecutar(
    num_pasos=25000,
    pasos_equilibrado=1000,
    frecuencia_muestreo=50,
    frecuencia_muestreo_velocidades=100,
    muestrear_velocidades=True,
    csv="resultados.csv",
    npy_velocidades="velocidades.npy",
)

print(df.head())
print(velocidades.shape)
```

Si no se requiere el análisis de velocidades, puede desactivarse dicho muestreo. En ese caso, el segundo valor devuelto será un array vacío:

```python
df, velocidades = sim.ejecutar(
    num_pasos=10000,
    pasos_equilibrado=2000,
    muestrear_velocidades=False,
)
```

## API principal

### `Simulador(...)`

Constructor de alto nivel en Python. Internamente crea una instancia de `ArgonSimulator`, implementada en C++.

Parámetros:

- `particulas_por_lado`: número de partículas por lado de la red inicial. El total es `N = n^3`.
- `densidad_reducida`: densidad reducida del sistema.
- `paso_tiempo`: paso temporal de integración.
- `temp_objetivo`: temperatura objetivo usada durante el equilibrado.
- `seed`: semilla del generador aleatorio. Si vale `0`, se usa entropía del sistema.
- `corregir_cm`: activa la corrección de deriva del centro de masas.
- `correccion_presion_cola`: activa la corrección de cola en la presión.
- `reescalar_velocidades`: activa el termostato por reescalado durante el equilibrado.

### `Simulador.ejecutar(...)`

Ejecuta la simulación y devuelve una tupla `(df, velocidades)`.

Parámetros:

- `num_pasos`: número total de pasos de simulación.
- `pasos_equilibrado`: pasos iniciales en los que se reescalan velocidades para estabilizar la temperatura.
- `frecuencia_muestreo`: cada cuántos pasos se guardan observables.
- `frecuencia_muestreo_velocidades`: cada cuántos pasos se muestrean módulos de velocidad.
- `muestrear_velocidades`: activa el muestreo de módulos de velocidad para análisis posteriores.
- `csv`: ruta opcional para guardar los resultados en disco.
- `npy_velocidades`: ruta opcional para guardar los módulos de velocidad con `np.save`.
- `forzar_calculo`: si es `False` y el archivo CSV ya existe, reutiliza los resultados almacenados en lugar de relanzar la simulación.

Comportamiento adicional:

- Si `csv` ya existe y `forzar_calculo=False`, los resultados se leen desde disco y no se ejecuta una nueva simulación.
- Si `npy_velocidades` apunta a un archivo existente, también se cargan los módulos de velocidad previamente almacenados.
- Si el núcleo C++ detecta una inestabilidad numérica, el contenedor Python emite un `RuntimeWarning` e intenta devolver los datos parciales ya muestreados.

## Salida

El primer elemento de retorno, `df`, es un `DataFrame` con las columnas:

- `paso`
- `tiempo`
- `temperatura`
- `presion`
- `energia_potencial`
- `energia_cinetica`
- `energia_total`

El segundo elemento, `velocidades`, es un `numpy.ndarray` con los módulos de velocidad muestreados. Si `muestrear_velocidades=False`, se devuelve un array vacío.

Si se proporciona el parámetro `csv`, se genera además un archivo con estas mismas magnitudes. Si se proporciona `npy_velocidades` y hubo muestreo de velocidades, se guarda también el array correspondiente en formato `.npy`.

## Funciones públicas de análisis

Además de `Simulador`, el paquete exporta estas utilidades:

- `graficar_energia(df, cutoff=None, *, ax=None, figsize=(9, 4.5))`: representa la energía total en función del tiempo o del paso de integración.
- `graficar_resumen_termodinamico(df, cutoff=None, *, axes=None, figsize=(11, 8))`: genera una rejilla `2x2` con temperatura, presión, energía cinética y energía potencial.
- `graficar_histograma_velocidades(velocidades, *, temperatura=1.002, bins=60, figsize=(12, 5), filepath=None)`: dibuja el histograma normalizado de velocidades y la distribución teórica de Maxwell-Boltzmann.

Las funciones de análisis esperan un `DataFrame` con las columnas termodinámicas necesarias y emplean `tiempo` como eje X cuando dicha columna está disponible; en caso contrario, utilizan `paso`.

Ejemplo de uso:

```python
fig1, axes = graficar_resumen_termodinamico(df, cutoff=2000)
ax = graficar_energia(df, cutoff=2000)

if velocidades.size > 0:
    fig2 = graficar_histograma_velocidades(
        velocidades,
        temperatura=df["temperatura"].iloc[-1],
        filepath="hist_velocidades.png",
    )
```

## Detalles físicos y numéricos

- Las magnitudes se manejan en unidades reducidas de Lennard-Jones.
- Las partículas se inicializan sobre una red cúbica simple.
- Las velocidades iniciales se generan aleatoriamente en distribución uniforme y se reescalan para arrancar cerca de la temperatura objetivo.
- El momento lineal total se corrige para evitar deriva del centro de masas.
- La presión incluye una corrección de cola por truncamiento del potencial.
- El cálculo de fuerzas utiliza un barrido directo de pares cuando el número de celdas por lado es pequeño y recurre a `cell lists` cuando la discretización espacial es suficiente.
- El reescalado de velocidades se aplica únicamente durante la fase de equilibrado, si esta opción se encuentra activada.

## Alcance físico del modelo

El simulador está orientado al estudio de un sistema monoatómico de argón descrito mediante el potencial de Lennard-Jones en unidades reducidas. La dinámica integra la evolución microscópica de las partículas bajo condiciones de contorno periódicas. Durante la fase de equilibrado puede aplicarse un control sencillo de temperatura mediante reescalado de velocidades; fuera de dicha fase, la evolución prosigue con las ecuaciones de movimiento sin la aplicación de ese ajuste, salvo que el usuario modifique el código fuente.

## Nota sobre CMake

El repositorio incluye un `CMakeLists.txt`, pero actualmente los nombres de archivo definidos ahí no coinciden con los del módulo real dentro de `simulador_dm/`. Por eso, para compilar e instalar el proyecto se recomienda usar `pip install .` hasta que ese archivo se actualice.

