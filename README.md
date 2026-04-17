# Simulador de dinámica molecular para Argón

Este proyecto implementa un simulador de dinámica molecular para un sistema de átomos de argón en unidades reducidas de Lennard-Jones. La parte numérica está escrita en C++ y se expone a Python mediante `pybind11`, con una interfaz de alto nivel que devuelve los resultados termodinámicos como un `DataFrame` de `pandas` y, opcionalmente, un array de NumPy con módulos de velocidad.

## Características

- Simulación de Argón en un ensamble `NVT`.
- Integración temporal mediante Verlet de velocidad.
- Interacciones de Lennard-Jones con radio de corte.
- Optimización del cálculo de fuerzas con `cell lists`.
- Condiciones de contorno periódicas.
- Muestreo de temperatura, presión y energías.
- Exportación opcional de resultados a CSV y de módulos de velocidad a NPY.

## Estructura del proyecto

```text
SimuladorArgon/
├── simulador_dm/
│   ├── __init__.py
│   ├── _wrapper_simulator.py
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

Dependencias Python del proyecto:

- `pybind11`
- `numpy`
- `pandas`

## Instalación

La vía más consistente con el estado actual del repositorio es instalar el paquete con `pip` desde la raíz del proyecto:

```bash
pip install .
```

Si quieres trabajar en modo desarrollo:

```bash
pip install -e .
```

En Windows puede ser necesario tener instaladas las herramientas de compilación de Visual Studio para C++.

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

## API principal

### `Simulador(...)`

Constructor de alto nivel en Python.

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
- `forzar_calculo`: si `False` y el CSV ya existe, reutiliza el archivo existente en lugar de recalcular.

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

Si `csv` ya existe y `forzar_calculo=False`, el wrapper carga ese archivo y devuelve sus datos sin lanzar una nueva simulación. Si además existe el archivo indicado por `npy_velocidades`, también se carga.

En configuraciones numéricamente inestables, el wrapper intenta devolver los datos parciales ya muestreados y emite un `RuntimeWarning`.

## Funciones públicas de análisis

Además de `Simulador`, el paquete exporta estas utilidades:

- `graficar_energia(df, cutoff=None, *, ax=None, figsize=(9, 4.5))`: representa la energía total frente a tiempo o paso.
- `graficar_resumen_termodinamico(df, cutoff=None, *, axes=None, figsize=(11, 8))`: genera una rejilla 2x2 con temperatura, presión, energía cinética y energía potencial.
- `graficar_histograma_velocidades(velocidades, *, temperatura=1.002, bins=60, figsize=(12, 5), filepath=None)`: dibuja el histograma normalizado de velocidades y la distribución teórica de Maxwell-Boltzmann.

Las funciones de análisis esperan un `DataFrame` con las columnas termodinámicas necesarias y usan `tiempo` como eje X cuando está disponible; si no, usan `paso`.

## Detalles físicos y numéricos

- Las magnitudes se manejan en unidades reducidas de Lennard-Jones.
- Las partículas se inicializan sobre una red cúbica simple.
- Las velocidades iniciales siguen una distribución de Maxwell-Boltzmann.
- El momento lineal total se corrige para evitar deriva del centro de masas.
- La presión incluye una corrección de cola por truncamiento del potencial.

## Nota sobre CMake

El repositorio incluye un `CMakeLists.txt`, pero actualmente los nombres de archivo definidos ahí no coinciden con los del módulo real dentro de `simulador_dm/`. Por eso, para compilar e instalar el proyecto se recomienda usar `pip install .` hasta que ese archivo se actualice.

## Posibles mejoras

- Añadir ejemplos o notebooks de análisis de resultados.
- Incorporar pruebas automáticas para validar energía, temperatura y presión.
- Corregir y sincronizar la configuración de compilación con CMake.
- Documentar casos de uso y parámetros recomendados para diferentes regímenes físicos.
