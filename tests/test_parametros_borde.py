from pathlib import Path

import numpy as np
import pandas as pd
import pytest

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")

from simulador_dm import Simulador
from simulador_dm.analisis import (
    graficar_energia,
    graficar_histograma_velocidades,
    graficar_resumen_termodinamico,
)


OUTPUT_DIR = Path("tests/outputs")
COLUMNAS_ESPERADAS = [
    "paso",
    "tiempo",
    "temperatura",
    "presion",
    "energia_potencial",
    "energia_cinetica",
    "energia_total",
]


def test_ejecutar_con_cero_pasos_devuelve_estructuras_vacias():
    sim = Simulador(particulas_por_lado=3, seed=11)

    df, velocidades = sim.ejecutar(num_pasos=0, pasos_equilibrado=0)

    assert isinstance(df, pd.DataFrame)
    assert list(df.columns) == COLUMNAS_ESPERADAS
    assert df.empty
    assert isinstance(velocidades, np.ndarray)
    assert velocidades.size == 0


def test_ejecutar_devuelve_velocidades_y_guarda_npy_sin_extension():
    OUTPUT_DIR.mkdir(exist_ok=True)
    npy_base = OUTPUT_DIR / "velocidades_borde"

    sim = Simulador(particulas_por_lado=3, seed=7)
    _, velocidades = sim.ejecutar(
        num_pasos=40,
        pasos_equilibrado=10,
        frecuencia_muestreo=5,
        muestrear_velocidades=True,
        frecuencia_muestreo_velocidades=10,
        npy_velocidades=str(npy_base),
        forzar_calculo=True,
    )

    ruta_npy = npy_base.with_suffix(".npy")
    assert velocidades.size > 0
    assert ruta_npy.exists()

    velocidades_guardadas = np.load(ruta_npy)
    assert np.array_equal(velocidades, velocidades_guardadas)


def test_reutiliza_csv_y_npy_existentes_cuando_no_fuerza_calculo():
    OUTPUT_DIR.mkdir(exist_ok=True)
    csv_path = OUTPUT_DIR / "cache_parametros.csv"
    npy_path = OUTPUT_DIR / "cache_parametros_vel.npy"

    sim = Simulador(particulas_por_lado=3, seed=5)
    df_generado, vel_generadas = sim.ejecutar(
        num_pasos=40,
        pasos_equilibrado=10,
        frecuencia_muestreo=5,
        muestrear_velocidades=True,
        frecuencia_muestreo_velocidades=10,
        csv=str(csv_path),
        npy_velocidades=str(npy_path),
        forzar_calculo=True,
    )

    df_cache, vel_cache = sim.ejecutar(
        num_pasos=999,
        pasos_equilibrado=0,
        frecuencia_muestreo=1,
        muestrear_velocidades=True,
        frecuencia_muestreo_velocidades=1,
        csv=str(csv_path),
        npy_velocidades=str(npy_path),
        forzar_calculo=False,
    )

    assert csv_path.exists()
    assert npy_path.exists()
    assert df_cache.equals(pd.read_csv(csv_path))
    assert len(df_cache) == len(df_generado)
    assert np.array_equal(vel_cache, vel_generadas)


def test_graficar_energia_falla_si_falta_energia_total():
    df = pd.DataFrame({"paso": [0, 1], "temperatura": [1.0, 1.1]})

    with pytest.raises(ValueError, match="energia_total"):
        graficar_energia(df)


def test_graficar_resumen_falla_si_faltan_columnas_requeridas():
    df = pd.DataFrame({"tiempo": [0.0, 0.1], "temperatura": [1.0, 1.1]})

    with pytest.raises(ValueError, match="presion"):
        graficar_resumen_termodinamico(df)


def test_graficar_energia_usa_paso_si_no_hay_tiempo():
    df = pd.DataFrame(
        {
            "paso": [0, 10, 20],
            "energia_total": [-1.0, -0.9, -0.95],
        }
    )

    ax = graficar_energia(df)

    assert ax.get_xlabel() == "Paso"
    matplotlib.pyplot.close(ax.figure)


def test_histograma_velocidades_guarda_figura():
    OUTPUT_DIR.mkdir(exist_ok=True)
    ruta_figura = OUTPUT_DIR / "hist_velocidades.png"
    velocidades = np.array([0.8, 1.0, 1.2, 1.4, 1.6], dtype=np.float64)

    fig = graficar_histograma_velocidades(
        velocidades,
        temperatura=1.0,
        filepath=ruta_figura,
    )

    assert ruta_figura.exists()
    assert len(fig.axes) == 2
    matplotlib.pyplot.close(fig)
