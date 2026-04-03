"""Funciones de analisis y representacion para los resultados del simulador."""

from __future__ import annotations

import matplotlib.pyplot as plt


def _validar_columnas(df, columnas):
    faltantes = [columna for columna in columnas if columna not in df.columns]
    if faltantes:
        columnas_txt = ", ".join(faltantes)
        raise ValueError(
            f"Faltan columnas necesarias en el DataFrame: {columnas_txt}"
        )


def _obtener_eje_x(df):
    if "tiempo" in df.columns:
        return "tiempo", "Tiempo reducido"
    if "paso" in df.columns:
        return "paso", "Paso"
    raise ValueError("El DataFrame debe contener la columna 'tiempo' o 'paso'.")


def _obtener_posicion_cutoff(df, cutoff, columna_x):
    if cutoff is None:
        return None

    if "paso" not in df.columns:
        return cutoff if columna_x == "tiempo" else cutoff

    filas_post_cutoff = df[df["paso"] >= cutoff]
    if filas_post_cutoff.empty:
        return None

    return filas_post_cutoff.iloc[0][columna_x]


def _dibujar_linea_cutoff(ax, df, cutoff, columna_x):
    posicion_cutoff = _obtener_posicion_cutoff(df, cutoff, columna_x)
    if posicion_cutoff is None:
        return

    ax.axvline(
        posicion_cutoff,
        color="crimson",
        linestyle="--",
        linewidth=1.5,
        label="Fin equilibrado",
    )


def graficar_energia(df, ax=None, cutoff=1000):
    """
    Representa la energia total del sistema respecto al tiempo o al paso.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame devuelto por el simulador.
    ax : matplotlib.axes.Axes, optional
        Eje sobre el que dibujar. Si no se pasa, se crea uno nuevo.
    cutoff : int, optional
        Paso a partir del cual se considera que termina el equilibrado.

    Returns
    -------
    matplotlib.axes.Axes
        Eje con la grafica generada.
    """
    _validar_columnas(df, ["energia_total"])
    columna_x, etiqueta_x = _obtener_eje_x(df)

    if ax is None:
        _, ax = plt.subplots(figsize=(9, 4.5))

    ax.plot(df[columna_x], df["energia_total"], color="steelblue", linewidth=1.6)
    _dibujar_linea_cutoff(ax, df, cutoff, columna_x)

    ax.set_title("Energia total")
    ax.set_xlabel(etiqueta_x)
    ax.set_ylabel("Energia reducida")
    ax.grid(True, alpha=0.3)

    if ax.get_legend_handles_labels()[0]:
        ax.legend()

    return ax


def graficar_resumen_termodinamico(df, cutoff=1000, figsize=(11, 8)):
    """
    Dibuja una rejilla 2x2 con temperatura, presion, energia cinetica
    y energia potencial.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame devuelto por el simulador.
    cutoff : int, optional
        Paso a partir del cual se considera terminado el equilibrado.
    figsize : tuple, optional
        Tamano de la figura creada.

    Returns
    -------
    tuple[matplotlib.figure.Figure, numpy.ndarray]
        Figura y ejes de la rejilla 2x2.
    """
    _validar_columnas(
        df,
        [
            "temperatura",
            "presion",
            "energia_cinetica",
            "energia_potencial",
        ],
    )
    columna_x, etiqueta_x = _obtener_eje_x(df)

    fig, axes = plt.subplots(2, 2, figsize=figsize, sharex=True)

    series = [
        ("temperatura", "Temperatura", "firebrick"),
        ("presion", "Presion", "darkgreen"),
        ("energia_cinetica", "Energia cinetica", "darkorange"),
        ("energia_potencial", "Energia potencial", "slateblue"),
    ]

    for ax, (columna, titulo, color) in zip(axes.flat, series):
        ax.plot(df[columna_x], df[columna], color=color, linewidth=1.4)
        _dibujar_linea_cutoff(ax, df, cutoff, columna_x)
        ax.set_title(titulo)
        ax.set_ylabel(titulo)
        ax.grid(True, alpha=0.3)

        if ax.get_legend_handles_labels()[0]:
            ax.legend()

    axes[1, 0].set_xlabel(etiqueta_x)
    axes[1, 1].set_xlabel(etiqueta_x)
    fig.suptitle("Resumen termodinamico", fontsize=14)
    fig.tight_layout()

    return fig, axes
