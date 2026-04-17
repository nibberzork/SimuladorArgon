import pytest

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")

from simulador_dm import Simulador
from simulador_dm.analisis import graficar_energia, graficar_resumen_termodinamico


PARAMETROS_ESTANDAR = {
    "particulas_por_lado": 10,
    "densidad_reducida": 0.84,
    "paso_tiempo": 0.005,
    "temp_objetivo": 1.002,
    "seed": 12345,
}


def _simulacion_estandar():
    sim = Simulador(**PARAMETROS_ESTANDAR)
    return sim.ejecutar(
        num_pasos=500,
        pasos_equilibrado=100,
        frecuencia_muestreo=2,
        csv=None,
    )


def test_graficas_analisis_con_parametros_estandar():
    df, _ = _simulacion_estandar()

    ax_energia = graficar_energia(df, cutoff=100)
    fig_resumen, axes_resumen = graficar_resumen_termodinamico(df, cutoff=100)

    assert ax_energia.get_title() == "Energía total"
    assert len(ax_energia.lines) == 5

    assert fig_resumen is not None
    assert axes_resumen.shape == (2, 2)

    for ax in axes_resumen.flat:
        assert len(ax.lines) == 5

    matplotlib.pyplot.close(ax_energia.figure)
    matplotlib.pyplot.close(fig_resumen)
