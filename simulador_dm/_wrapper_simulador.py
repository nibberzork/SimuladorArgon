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
    Simulador de dinámica molecular para argón en unidades reducidas de Lennard-Jones.

    Esta clase envuelve el simulador implementado en C++ (`ArgonSimulator`)
    y proporciona una interfaz de alto nivel en Python para ejecutar
    simulaciones y recuperar sus resultados en estructuras de datos
    adecuadas para análisis posterior.

    Unidades reducidas:
        - Longitud   : sigma (diámetro atómico)
        - Energia    : epsilon (profundidad del pozo LJ)
        - Masa       : m (masa atómica)
        - Tiempo     : sigma * sqrt(m / epsilon)
        - Temperatura: epsilon / k_B

    Parameters
    ----------
    particulas_por_lado : int, optional
        Número de partículas por dimensión. El sistema tendrá n^3
        partículas en total.
        Por defecto 8.
    densidad_reducida : float, optional
        Densidad reducida rho* = rho * sigma^3. Típicamente entre 0.8 y 1.2.
        Por defecto 0.84.
    paso_tiempo : float, optional
        Paso de integración dt*. Típicamente entre 0.001 y 0.01.
        Por defecto 0.005.
    temp_objetivo : float, optional
        Temperatura objetivo T* para el termostato. Por defecto 1.002.
    seed : int, optional
        Semilla del generador aleatorio. Si vale 0, se utiliza entropía
        del hardware.
        Por defecto 0.
    corregir_cm : bool, optional
        Corrige la deriva del centro de masas en cada paso. Por defecto True.
    correccion_presion_cola : bool, optional
        Añade la corrección de largo alcance a la presión. Por defecto True.
    reescalar_velocidades : bool, optional
        Activa el termostato por reescalado durante el equilibrado.
        Por defecto True.

    Examples
    --------
    >>> sim = Simulador(particulas_por_lado=6, densidad_reducida=0.84, temp_objetivo=1.0)
    >>> df, velocidades = sim.ejecutar(num_pasos=10000, pasos_equilibrado=2000)
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
        Ejecuta la simulación y devuelve los resultados muestreados.

        El método realiza la inicialización del sistema, la fase de
        equilibrado con termostato activo (si `reescalar_velocidades=True`)
        y la fase de producción, almacenando una muestra cada
        `frecuencia_muestreo` pasos.

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
            Indica si debe activarse el muestreo de módulos de velocidad.
            Por defecto False.
        csv : str, optional
            Ruta del archivo CSV de salida para magnitudes termodinámicas.
            Si es None, no se guarda salida en disco.
            Por defecto None.
        npy_velocidades : str, optional
            Ruta del archivo `.npy` para guardar módulos de velocidad con
            `np.save`. Si es None, no se guardan.
            Por defecto None.
        forzar_calculo : bool, optional
            Si es True, fuerza la simulación incluso si el archivo CSV ya existe.
            Por defecto False.

        Returns
        -------
        tuple[pd.DataFrame, np.ndarray]
            Tupla con el `DataFrame` de magnitudes termodinámicas y un array
            de `NumPy` con los módulos de velocidad. El array estará vacío
            si no se activó dicho muestreo.

            - paso              : int,   número de paso de integración
            - tiempo            : float, tiempo reducido t*
            - temperatura       : float, temperatura instantánea T*
            - presion           : float, presión instantánea P*
            - energia_potencial : float, energía potencial U*
            - energia_cinetica  : float, energía cinética K*
            - energia_total     : float, energía total E* = U* + K*

        Raises
        ------
        RuntimeError
            Si `csv` no es None y no se puede escribir el archivo indicado.
        
        Notes
        -----
        Si `csv` existe y `forzar_calculo` es False, los resultados se cargan
        directamente desde disco y no se lanza una nueva simulación.

        Si el núcleo C++ detecta una inestabilidad numérica, este método emite
        un `RuntimeWarning` y devuelve los datos parciales disponibles.
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
            resultados = self._sim.ejecutar(config, csv)
        except ErrorInestabilidadNumerica as e:
            warnings.warn(
                f"Simulación terminada prematuramente: {e}. "
                "Se devuelven los datos parciales hasta el momento de la explosión numérica.",
                RuntimeWarning,
                stacklevel=2,
            )
            resultados = e.resultados_parciales
        except Exception as e:
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
