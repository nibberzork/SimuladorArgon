"""Funciones de analisis y representacion para los resultados del simulador."""

from __future__ import annotations

import pandas as pd
from typing import Sequence, Optional
import matplotlib.pyplot as plt

_UNIDADES = {
    "temperatura":       "T*",
    "presion":           "P*",
    "energia_cinetica":  "K*",
    "energia_potencial": "U*",
    "energia_total":     "E*",
}

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


def _dibujar_estadisticas(ax, df, cutoff,  columna):
    """Dibuja media y banda min/max en la region estable de cualquier magnitud."""
    filas_estables = df[df["paso"] >= cutoff] if cutoff and "paso" in df.columns else df

    if filas_estables.empty:
        return None

    media = filas_estables[columna].mean()
    minimo = filas_estables[columna].min()
    maximo = filas_estables[columna].max()

    unidad = _UNIDADES.get(columna, columna)

    ax.axhline(media,  color="black",   linestyle="-.",  linewidth=1.4,
               label=f"{unidad} media = {media:.4f}")
    ax.axhline(minimo, color="dimgray", linestyle=":",   linewidth=1.2,
               label=f"{unidad} min = {minimo:.4f}")
    ax.axhline(maximo, color="dimgray", linestyle=":",   linewidth=1.2,
               label=f"{unidad} max = {maximo:.4f}")

    return media

def graficar_energia(
    df : pd.DataFrame,
    cutoff=None,
    *,
    ax=None,
    figsize=(9, 4.5),
) -> plt.Axes:
    """
    Representa la energía total del sistema respecto al tiempo o al paso.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame devuelto por el simulador.
    cutoff : int, optional
        Paso a partir del cual se considera que termina el equilibrado.
        Si es None, no se marca región de equilibrio.
    ax : matplotlib.axes.Axes, optional
        Eje sobre el que dibujar. Si no se pasa, se crea uno nuevo.
    figsize : tuple, optional
        Tamaño de la figura si se crea un eje nuevo.

    Returns
    -------
    matplotlib.axes.Axes
        Eje con la gráfica generada.
    """
    _validar_columnas(df, ["energia_total"])
    columna_x, etiqueta_x = _obtener_eje_x(df)

    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    ax.plot(
        df[columna_x],
        df["energia_total"],
        color="steelblue",
        linewidth=1.6,
    )

    _dibujar_linea_cutoff(ax, df, cutoff, columna_x)

    _dibujar_estadisticas(ax, df, cutoff, "energia_total")

    ax.set_title("Energía total")
    ax.set_xlabel(etiqueta_x)
    ax.set_ylabel(_UNIDADES.get("energia_total", "E*"))
    ax.grid(True, alpha=0.3)

    if ax.get_legend_handles_labels()[0]:
        ax.legend()

    return ax


def graficar_resumen_termodinamico(
    df : pd.DataFrame,
    cutoff=None,
    *,
    axes=None,
    figsize=(11, 8),
) -> tuple[plt.Figure, Sequence[Sequence[plt.Axes]]]:
    """
    Dibuja una rejilla 2x2 con temperatura, presión, energía cinética
    y energía potencial.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame devuelto por el simulador.
    cutoff : int, optional
        Paso a partir del cual se considera terminado el equilibrado.
        Si es None, no se marca región de equilibrio.
    axes : array-like of matplotlib.axes.Axes, optional
        Ejes sobre los que dibujar. Si no se pasan, se crean nuevos.
    figsize : tuple, optional
        Tamaño de la figura si se crean nuevos ejes.

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

    if axes is None:
        fig, axes = plt.subplots(2, 2, figsize=figsize, sharex=True)
    else:
        fig = axes.flat[0].figure  # inferir figura a partir de los ejes

    series = [
        ("temperatura", "Temperatura", "firebrick"),
        ("presion", "Presión", "darkgreen"),
        ("energia_cinetica", "Energía cinética", "darkorange"),
        ("energia_potencial", "Energía potencial", "slateblue"),
    ]

    for ax, (columna, titulo, color) in zip(axes.flat, series):
        ax.plot(
            df[columna_x],
            df[columna],
            color=color,
            linewidth=1.4,
        )

        _dibujar_linea_cutoff(ax, df, cutoff, columna_x)

        _dibujar_estadisticas(ax, df, cutoff, columna)
        
        ax.set_title(titulo)
        ax.set_ylabel(_UNIDADES.get(columna, columna))
        ax.grid(True, alpha=0.3)

        # Escalar los ejes para dar espacio a la leyenda
        y_min, y_max = ax.get_ylim()
        margen = (y_max - y_min) * 0.15  # 15% de margen arriba
        ax.set_ylim(y_min, y_max + margen)

        if ax.get_legend_handles_labels()[0]:
            ax.legend(loc='upper right')

    axes[1, 0].set_xlabel(etiqueta_x)
    axes[1, 1].set_xlabel(etiqueta_x)

    fig.suptitle("Resumen termodinámico", fontsize=14)
    fig.tight_layout()

    return fig, axes

