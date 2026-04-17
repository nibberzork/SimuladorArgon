import numpy as np
import pandas as pd
from pathlib import Path
from simulador_dm import Simulador

OUTPUT_DIR = Path("tests/outputs")

def test_simulacion_exportar_csv(tmp_path):
    csv_path = tmp_path / "ejemplo_ejecucion.csv"

    sim = Simulador(particulas_por_lado=4)  # 64 partículas, rápido
    df, velocidades = sim.ejecutar(
        num_pasos=10000,
        pasos_equilibrado=1000,
        frecuencia_muestreo=10,
        csv=str(csv_path),
        forzar_calculo=True,
    )

    # ¿devuelve un DataFrame con las columnas esperadas?
    columnas = ['paso', 'tiempo', 'temperatura', 'presion',
                'energia_potencial', 'energia_cinetica', 'energia_total']
    assert isinstance(df, pd.DataFrame)
    assert isinstance(velocidades, np.ndarray)
    assert velocidades.size == 0
    assert all(c in df.columns for c in columnas)
    assert len(df) > 0

    # ¿se guardó el CSV?
    assert csv_path.exists()
    df_csv = pd.read_csv(csv_path)
    assert len(df_csv) > 0

    # Checks físicos mínimos
    assert (df['temperatura'] > 0).all()
    # Energía total debe ser aproximadamente constante tras equilibrado
    df_prod = df[df['paso'] > 100]
    variacion_relativa = df_prod['energia_total'].std() / abs(df_prod['energia_total'].mean())
    assert variacion_relativa < 0.05, f"Energía total varía demasiado: {variacion_relativa:.3f}"

    print(df.describe())  # visible con pytest -s


def test_ejecutar_forzar_simulacion_si_csv_existe():
    OUTPUT_DIR.mkdir(exist_ok=True)
    csv_path = str(OUTPUT_DIR / "ejemplo_ejecucion_forzar.csv")

    sim = Simulador(particulas_por_lado=4)
    df1, vel1 = sim.ejecutar(
        num_pasos=1000,
        pasos_equilibrado=100,
        frecuencia_muestreo=10,
        csv=csv_path,
        muestrear_velocidades=False,
    )

    assert Path(csv_path).exists()
    assert isinstance(df1, pd.DataFrame)
    assert isinstance(vel1, np.ndarray)
    assert vel1.size == 0

    df2, vel2 = sim.ejecutar(
        num_pasos=1000,
        pasos_equilibrado=100,
        frecuencia_muestreo=10,
        csv=csv_path,
        forzar_calculo=True,
        muestrear_velocidades=False,
    )

    assert isinstance(df2, pd.DataFrame)
    assert isinstance(vel2, np.ndarray)
    assert vel2.size == 0
