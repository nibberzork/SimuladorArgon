from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension, build_ext

ext_modules = [
    Pybind11Extension(
        'simulador_dm.simulador',
        ['simulador_dm/bindings.cpp', 'simulador_dm/_simulador.cpp'],
        include_dirs=['simulador_dm'],
        extra_compile_args=['/std:c++17', '/O3'],
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
    cmdclass={"build_ext": build_ext},
    install_requires=[
        'pybind11>=2.6.0',
        'pandas>=1.0',
        'numpy>=1.20',
    ],
    python_requires='>=3.7',
)