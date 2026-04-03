# tests/test_simulacion_basica.py
import pytest
import pandas as pd
from pathlib import Path
from simulador_dm import Simulador

OUTPUT_DIR = Path("tests/outputs")

def test_simulacion_exportar_csv():
    OUTPUT_DIR.mkdir(exist_ok=True)
    csv_path = str(OUTPUT_DIR / "ejemplo_ejecucion.csv")

    sim = Simulador(particulas_por_lado=4)  # 64 partículas, rápido
    df = sim.ejecutar(
        num_pasos=10000,
        pasos_equilibrado=1000,
        frecuencia_muestreo=10,
        csv=csv_path
    )

    # ¿devuelve un DataFrame con las columnas esperadas?
    columnas = ['paso', 'tiempo', 'temperatura', 'presion',
                'energia_potencial', 'energia_cinetica', 'energia_total']
    assert isinstance(df, pd.DataFrame)
    assert all(c in df.columns for c in columnas)
    assert len(df) > 0

    # ¿se guardó el CSV?
    assert Path(csv_path).exists()
    df_csv = pd.read_csv(csv_path)
    assert len(df_csv) > 0

    # Checks físicos mínimos
    assert (df['temperatura'] > 0).all()
    # Energía total debe ser aproximadamente constante tras equilibrado
    df_prod = df[df['paso'] > 100]
    variacion_relativa = df_prod['energia_total'].std() / abs(df_prod['energia_total'].mean())
    assert variacion_relativa < 0.05, f"Energía total varía demasiado: {variacion_relativa:.3f}"

    print(df.describe())  # visible con pytest -s