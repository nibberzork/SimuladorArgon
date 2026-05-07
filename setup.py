from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension, build_ext
import os, sys


DEBUG = os.environ.get("DEBUG", "0") == "1"

# Configuración de flags según el sistema operativo y modo
if sys.platform == "win32":
    # Flags para Windows (MSVC)
    cpp_std = '/std:c++17'
    base_compile_args = ['/EHsc', '/bigobj']
    if DEBUG:
        rebuild_flags = ['/Zi', '/Od'] if DEBUG else ['/O2'] # /Zi: Debug symbols, /Od: No optimizar
        link_args = ['/DEBUG', '/OPT:REF', '/OPT:ICF']
    else:
        rebuild_flags = base_compile_args + ['/O2', '/Zi']

        # /DEBUG: Genera el archivo .pdb que lee Scalene
        # /OPT:REF y /OPT:ICF: Son obligatorios cuando usas /DEBUG en Release, 
        # de lo contrario, el enlazador desactiva las optimizaciones.
        link_args = ['/DEBUG', '/OPT:REF', '/OPT:ICF']

ext_modules = [
    Pybind11Extension(
        'simulador_dm.simulador',
        ['simulador_dm/bindings.cpp', 'simulador_dm/_simulador.cpp'],
        extra_compile_args=[cpp_std] + rebuild_flags,
        extra_link_args=link_args,
    ),
]

setup(
    name='simulador-dm',
    version='1.0.1',
    description='Simulador de dinámica molecular para argón en unidades reducidas de Lennard-Jones',
    packages=['simulador_dm'],
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
    install_requires=[
        'pybind11>=2.6.0',
        'pandas>=1.0',
        'numpy>=1.20',
    ],
    python_requires='>=3.7',
)