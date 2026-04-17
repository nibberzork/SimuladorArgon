# tests/test_inestabilidad.py
import warnings
import pytest
import numpy as np
from simulador_dm import Simulador  # ajusta al nombre de tu paquete

def test_explosion_numerica():
    """
    Con la política actual, el fallback solo se activa ante NaN/Inf.
    Si la simulación produce valores muy grandes pero finitos, no debe
    emitir warning ni cortar los datos.
    """
    sim = Simulador(
        particulas_por_lado=4,   # 64 partículas, rápido
        densidad_reducida=1.2,   # densidad alta → más colisiones
        paso_tiempo=0.08,        # dt* agresivo, pero no necesariamente NaN/Inf
        temp_objetivo=1.0,
    )

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        df, velocidades = sim.ejecutar(
            num_pasos=5000,
            pasos_equilibrado=100,
            frecuencia_muestreo=10,
            muestrear_velocidades=False,
        )

    runtime_warnings = [x for x in w if issubclass(x.category, RuntimeWarning)]
    assert len(runtime_warnings) == 0
    assert df is not None
    assert velocidades.size == 0
    assert len(df) == 5000 // 10

    columnas_esperadas = [
        "paso", "tiempo", "temperatura",
        "presion", "energia_potencial",
        "energia_cinetica", "energia_total"
    ]
    for col in columnas_esperadas:
        assert col in df.columns

    assert not df.isnull().any().any()


def test_simulacion_ordinaria():
    """Configuración estable no debe emitir warnings de inestabilidad."""
    sim = Simulador(
        particulas_por_lado=4,
        densidad_reducida=0.84,
        paso_tiempo=0.005,       # dt* conservador
        temp_objetivo=1.002,
    )

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        df, _ = sim.ejecutar(num_pasos=2000, pasos_equilibrado=500)

    runtime_warnings = [x for x in w if issubclass(x.category, RuntimeWarning)]
    assert len(runtime_warnings) == 0
    assert len(df) > 0
