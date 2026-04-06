#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "_simulador.hpp"

namespace py = pybind11;

PYBIND11_MODULE(simulador, m) {
    m.doc() = "Módulo de simulación molecular para Argón";

    py::class_<ConfiguracionSimulacion>(m, "ConfiguracionSimulacion")
        .def(py::init<>())
        .def_readwrite("num_pasos", &ConfiguracionSimulacion::num_pasos)
        .def_readwrite("pasos_equilibrado", &ConfiguracionSimulacion::pasos_equilibrado)
        .def_readwrite("frecuencia_muestreo", &ConfiguracionSimulacion::frecuencia_muestreo);

    py::class_<ResultadosSimulacion>(m, "ResultadosSimulacion")
        .def(py::init<>())
        .def_readonly("pasos", &ResultadosSimulacion::pasos)
        .def_readwrite("tiempos", &ResultadosSimulacion::tiempos)
        .def_readwrite("temperaturas", &ResultadosSimulacion::temperaturas)
        .def_readwrite("presiones", &ResultadosSimulacion::presiones)
        .def_readwrite("energias_potenciales", &ResultadosSimulacion::energias_potenciales)
        .def_readwrite("energias_cineticas", &ResultadosSimulacion::energias_cineticas)
        .def_readwrite("energias_totales", &ResultadosSimulacion::energias_totales);

    py::class_<ArgonSimulator>(m, "ArgonSimulator")
        .def(py::init<int, double, double, double, unsigned int, bool, bool, bool>(),
             py::arg("particulas_por_lado") = 8,
             py::arg("densidad_reducida") = 0.84,
             py::arg("paso_tiempo") = 0.005,
             py::arg("temp_objetivo") = 1.002,
             py::arg("seed") = 0,
             py::arg("corregir_cm") = true,
             py::arg("correccion_presion_cola") = true,
             py::arg("reescalar_velocidades") = true)
        .def("ejecutar", &ArgonSimulator::ejecutar);
}

