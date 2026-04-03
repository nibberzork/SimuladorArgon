from setuptools import setup, Extension
import pybind11

# Extensión para el módulo simulador compilado con pybind11
ext_modules = [
    Extension(
        'simulador_dm.simulador',
        ['simulador_dm/bindings.cpp', 'simulador_dm/_simulador.cpp'],
        include_dirs=[pybind11.get_include(), 'simulador_dm'],
        language='c++',
        extra_compile_args=['/std:c++17', '/O2'],
    ),
]

setup(
    name='simulador-dm',
    version='1.0.0',
    description='Simulador de dinámica molecular para Argón con interfaz Python',
    author='Usuario',
    url='https://example.com',
    packages=['simulador_dm'],
    ext_modules=ext_modules,
    install_requires=[
        'pybind11>=2.6.0',
        'pandas>=1.0',
        'numpy>=1.20',
    ],
    python_requires='>=3.7',
)
