# _wrapper_simulator.py
import pandas as pd
from .simulador import ArgonSimulator, ConfiguracionSimulacion

# Nota: ArgonSimulator y ConfiguracionSimulacion se importan del módulo privado _simulador
# compilado como simulador.cp312-win_amd64.pyd

class Simulador:
    def __init__(self,
                 particulas_por_lado: int = 8,
                 densidad_reducida: float = 0.84,
                 paso_tiempo: float = 0.005,
                 temp_objetivo: float = 1.002,
                 seed: int = 0):
        self._sim = ArgonSimulator(particulas_por_lado, densidad_reducida, paso_tiempo, temp_objetivo, seed)

    def ejecutar(self,
                 num_pasos: int = 25000,
                 pasos_equilibrado: int = 1000,
                 frecuencia_muestreo: int = 50,
                 csv: str = None) -> pd.DataFrame:

        config = ConfiguracionSimulacion()
        config.num_pasos = num_pasos
        config.pasos_equilibrado = pasos_equilibrado
        config.frecuencia_muestreo = frecuencia_muestreo

        resultados = self._sim.ejecutar(config, csv)

        return pd.DataFrame({
            'paso':               resultados.pasos,
            'tiempo':             resultados.tiempos,
            'temperatura':        resultados.temperaturas,
            'presion':            resultados.presiones,
            'energia_potencial':  resultados.energias_potenciales,
            'energia_cinetica':   resultados.energias_cineticas,
            'energia_total':      resultados.energias_totales,
        })
