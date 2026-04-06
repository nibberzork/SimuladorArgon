# _wrapper_simulator.py
import pandas as pd
from .simulador import ArgonSimulator, ConfiguracionSimulacion

# Nota: ArgonSimulator y ConfiguracionSimulacion se importan del módulo privado _simulador
# compilado como simulador.cp312-win_amd64.pyd

class Simulador:
    """
    Simulador de dinamica molecular para Argon en unidades reducidas de Lennard-Jones.

    Envuelve el simulador C++ (ArgonSimulator) exponiendo una interfaz Python
    con resultados como DataFrame de pandas.

    Unidades reducidas:
        - Longitud   : sigma (diametro atomico)
        - Energia    : epsilon (profundidad del pozo LJ)
        - Masa       : m (masa atomica)
        - Tiempo     : sigma * sqrt(m / epsilon)
        - Temperatura: epsilon / k_B

    Parameters
    ----------
    particulas_por_lado : int, optional
        Particulas por dimension. El sistema tendra n^3 particulas en total.
        Por defecto 8.
    densidad_reducida : float, optional
        Densidad reducida rho* = rho * sigma^3. Tipicamente entre 0.8 y 1.2.
        Por defecto 0.84.
    paso_tiempo : float, optional
        Paso de integracion dt*. Tipicamente entre 0.001 y 0.01.
        Por defecto 0.005.
    temp_objetivo : float, optional
        Temperatura objetivo T* para el termostato. Por defecto 1.002.
    seed : int, optional
        Semilla del generador aleatorio. 0 usa entropia del hardware.
        Por defecto 0.
    corregir_cm : bool, optional
        Corrige la deriva del centro de masa en cada paso. Por defecto True.
    correccion_presion_cola : bool, optional
        Añade correccion de largo alcance a la presion. Por defecto True.
    reescalar_velocidades : bool, optional
        Activa termostato por reescalado durante el equilibrado. Por defecto True.

    Examples
    --------
    >>> sim = Simulador(particulas_por_lado=6, densidad_reducida=0.84, temp_objetivo=1.0)
    >>> df = sim.ejecutar(num_pasos=10000, pasos_equilibrado=2000)
    >>> df['temperatura'].mean()
    """
    def __init__(self,
                 particulas_por_lado: int = 8,
                 densidad_reducida: float = 0.84,
                 paso_tiempo: float = 0.005,
                 temp_objetivo: float = 1.002,
                 seed: int = 0,
                 corregir_cm: bool = True,
                 correccion_presion_cola: bool = True,
                 reescalar_velocidades: bool = True):
        self._sim = ArgonSimulator(
            particulas_por_lado,
            densidad_reducida,
            paso_tiempo,
            temp_objetivo,
            seed,
            corregir_cm,
            correccion_presion_cola,
            reescalar_velocidades,
        )

    def ejecutar(self,
                 num_pasos: int = 25000,
                 pasos_equilibrado: int = 1000,
                 frecuencia_muestreo: int = 50,
                 csv: str = None,
                 desde_csv: bool = True) -> pd.DataFrame:
        """
        Ejecuta la simulacion y devuelve los resultados muestreados.

        Realiza la inicializacion del sistema, el equilibrado con termostato
        activo (si reescalar_velocidades=True) y la produccion, guardando
        una muestra cada frecuencia_muestreo pasos.

        Parameters
        ----------
        num_pasos : int, optional
            Numero total de pasos de integracion. Por defecto 25000.
        pasos_equilibrado : int, optional
            Pasos iniciales con termostato activo. Debe ser menor que
            num_pasos. Por defecto 1000.
        frecuencia_muestreo : int, optional
            Intervalo de pasos entre muestras guardadas. Por defecto 50.
        csv : str, optional
            Ruta del archivo CSV de salida. None para no guardar en disco.
            Por defecto None.
        desde_csv : bool, optional
            Si True y csv apunta a un archivo existente, carga los datos
            sin simular. Por defecto True.

        Returns
        -------
        pd.DataFrame
            DataFrame con una fila por muestra (cada frecuencia_muestreo
            pasos) y las columnas:

            - paso              : int,   numero de paso de integracion
            - tiempo            : float, tiempo reducido t*
            - temperatura       : float, temperatura instantanea T*
            - presion           : float, presion instantanea P*
            - energia_potencial : float, energia potencial U*
            - energia_cinetica  : float, energia cinetica K*
            - energia_total     : float, energia total E* = U* + K*

        Raises
        ------
        RuntimeError
            Si csv no es None y no se puede escribir el archivo indicado.
        """
        if desde_csv and csv is not None:
            from os.path import exists
            if exists(csv):
                return pd.read_csv(csv)
        

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
