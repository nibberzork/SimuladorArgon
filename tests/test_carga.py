import pytest
import numpy as np
from simulador_dm import Simulador

def test_heavy_load_profiling():
    """
    Test de alta carga diseñado para el perfilado con Scalene.
    Usa 1331 partículas (11^3) para asegurar que el cálculo de fuerzas
    sea el cuello de botella predominante.
    """
    # Configuramos una carga lo suficientemente pesada
    sim = Simulador(
        particulas_por_lado=11,  # 1331 partículas
        densidad_reducida=0.8,
        paso_tiempo=0.005,
        temp_objetivo=1.0,
        seed=42
    )

    # Ejecutamos 5000 pasos SIN escribir en disco (csv=None)
    # Esto elimina el tiempo 'System' (amarillo) y resalta el 'Native' (morado)
    df, _ = sim.ejecutar(
        num_pasos=5000,
        pasos_equilibrado=500,
        frecuencia_muestreo=100, # Muestreo espaciado para no saturar Python
        csv=None,                # CRÍTICO: No escribir CSV durante el profiling
        muestrear_velocidades=False
    )
    
    assert len(df) > 0