import pandas as pd
import numpy as np
import warnings
from . import simulador as _simulador_mod

ArgonSimulator = _simulador_mod.ArgonSimulator
ConfiguracionSimulacion = _simulador_mod.ConfiguracionSimulacion
# Compatibilidad con binarios antiguos que aún no exponen la excepción específica.
ErrorInestabilidadNumerica = getattr(_simulador_mod, "ErrorInestabilidadNumerica", RuntimeError)


# Nota: ArgonSimulator y ConfiguracionSimulacion se importan del módulo privado _simulador
# compilado como simulador.cp312-win_amd64.pyd

class WraperSimulador:
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
                 frecuencia_muestreo_velocidades: int = 100,
                 muestrear_velocidades: bool = False,
                 csv: str = None,
                 npy_velocidades: str = None,
                 forzar_calculo: bool = False) -> tuple[pd.DataFrame, np.ndarray]:
        """
        Ejecuta la simulacion y devuelve los resultados muestreados.

        Realiza la inicializacion del sistema, el equilibrado con termostato
        activo (si reescalar_velocidades=True) y la produccion, guardando
        una muestra cada frecuencia_muestreo pasos.

        Parameters
        ----------
        num_pasos : int, optional
            Número total de pasos de integración. Por defecto 25000.
        pasos_equilibrado : int, optional
            Pasos iniciales con termostato activo. Debe ser menor que
            num_pasos. Por defecto 1000.
        frecuencia_muestreo : int, optional
            Intervalo de pasos entre muestras guardadas (termodinámicas). Por defecto 50.
        frecuencia_muestreo_velocidades : int, optional
            Intervalo de pasos entre muestras de módulos de velocidad. Por defecto 100.
        muestrear_velocidades : bool, optional
            Si activar muestreo de módulos de velocidad. Por defecto False.
        csv : str, optional
            Ruta del archivo CSV de salida para termodinámicas. None para no guardar en disco.
            Por defecto None.
        npy_velocidades : str, optional
            Ruta del archivo .npy para guardar módulos de velocidad con np.save. None para no guardar.
            Por defecto None.
        forzar_calculo : bool, optional
            Si True, fuerza la simulación incluso si el archivo CSV existe.
            Por defecto False.

        Returns
        -------
        tuple[pd.DataFrame, np.ndarray]
            Tupla con el DataFrame de termodinámicas y un array NumPy con
            los módulos de velocidad. El array estará vacío si no se muestrearon.

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
        if not forzar_calculo and csv is not None:
            from os.path import exists
            if exists(csv):
                df = pd.read_csv(csv)
                velocidades_array = np.empty((0,), dtype=np.float64)
                if npy_velocidades is not None:
                    ruta_vel = npy_velocidades
                    if not ruta_vel.lower().endswith('.npy'):
                        ruta_vel += '.npy'
                    if exists(ruta_vel):
                        velocidades_array = np.load(ruta_vel)
                return df, velocidades_array

        config = ConfiguracionSimulacion()
        config.num_pasos = num_pasos
        config.pasos_equilibrado = pasos_equilibrado
        config.frecuencia_muestreo = frecuencia_muestreo
        config.frecuencia_velocidades = frecuencia_muestreo_velocidades
        config.muestrear_velocidades = muestrear_velocidades

        # TODO: Capturar el error y extraer los datos parciales relanzando un warning para que no corte la ejecución
        try:
            print("[DEBUG] Iniciando simulación C++...")
            resultados = self._sim.ejecutar(config, csv)
            print("[DEBUG] Simulación completada sin errores.")
        except ErrorInestabilidadNumerica as e:
            print(f"[DEBUG] ErrorInestabilidadNumerica capturado: {e}")
            warnings.warn(
                f"Simulación terminada prematuramente: {e}. "
                "Se devuelven los datos parciales hasta el momento de la explosión numérica.",
                RuntimeWarning,
                stacklevel=2,
            )
            resultados = e.resultados_parciales
        except Exception as e:
            print(f"[DEBUG] Excepción inesperada tipo {type(e).__name__}: {e}")
            raise
        
        # Crear DataFrame con termodinámicas
        df = pd.DataFrame({
            'paso':               resultados.pasos,
            'tiempo':             resultados.tiempos,
            'temperatura':        resultados.temperaturas,
            'presion':            resultados.presiones,
            'energia_potencial':  resultados.energias_potenciales,
            'energia_cinetica':   resultados.energias_cineticas,
            'energia_total':      resultados.energias_totales,
        })
        
        # Exponer módulos de velocidad como array NumPy
        if muestrear_velocidades:
            velocidades_array = np.array(resultados.modulos_velocidades, dtype=np.float64)
        else:
            velocidades_array = np.empty((0,), dtype=np.float64)
        
        
        # Guardar DataFrame en CSV si se especifica
        if csv is not None:
            df.to_csv(csv, index=False)
            print(f"Datos termodinámicos guardados en {csv}")

        # Guardar velocidades con np.save si se especifica y se muestrearon
        if npy_velocidades is not None and velocidades_array.size > 0:
            ruta_vel = npy_velocidades
            if not ruta_vel.lower().endswith('.npy'):
                ruta_vel += '.npy'
            np.save(ruta_vel, velocidades_array)
            print(f"Módulos de velocidad guardados en {ruta_vel}")
        
        return df, velocidades_array
