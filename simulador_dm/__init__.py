"""API pública del paquete `simulador_dm`.

Este módulo expone la interfaz de alto nivel del simulador y las
funciones de análisis asociadas a los resultados termodinámicos.
"""

from ._wrapper_simulador import WraperSimulador as Simulador
from .analisis import (
    graficar_energia,
    graficar_resumen_termodinamico,
    graficar_histograma_velocidades,
)

__all__ = [
    "Simulador",
    "graficar_energia",
    "graficar_resumen_termodinamico",
    "graficar_histograma_velocidades",
]
