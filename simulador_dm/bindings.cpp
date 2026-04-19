#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <iostream>
#include "_simulador.hpp"

#ifndef _WIN32
#include <windows.h>
#endif


namespace py = pybind11;

PYBIND11_MODULE(simulador, m) {

    #ifndef _WIN32
        SetConsoleOutputCP(CP_UTF8);
        SetConsoleCP(CP_UTF8);
    #endif


    m.doc() = "Módulo de simulación molecular para Argón";

    py::class_<ConfiguracionSimulacion>(m, "ConfiguracionSimulacion")
        .def(py::init<>())
        .def_readwrite("num_pasos", &ConfiguracionSimulacion::num_pasos)
        .def_readwrite("pasos_equilibrado", &ConfiguracionSimulacion::pasos_equilibrado)
        .def_readwrite("frecuencia_muestreo", &ConfiguracionSimulacion::frecuencia_muestreo)
        .def_readwrite("frecuencia_velocidades", &ConfiguracionSimulacion::frecuencia_velocidades)
        .def_readwrite("muestrear_velocidades", &ConfiguracionSimulacion::muestrear_velocidades);

    py::class_<ResultadosSimulacion>(m, "ResultadosSimulacion")
        .def(py::init<>())
        .def_readonly("pasos", &ResultadosSimulacion::pasos)
        .def_readwrite("tiempos", &ResultadosSimulacion::tiempos)
        .def_readwrite("temperaturas", &ResultadosSimulacion::temperaturas)
        .def_readwrite("presiones", &ResultadosSimulacion::presiones)
        .def_readwrite("energias_potenciales", &ResultadosSimulacion::energias_potenciales)
        .def_readwrite("energias_cineticas", &ResultadosSimulacion::energias_cineticas)
        .def_readwrite("energias_totales", &ResultadosSimulacion::energias_totales)
        .def_readonly("modulos_velocidades", &ResultadosSimulacion::modulos_velocidades);

    py::register_exception<ErrorInestabilidadNumerica>(
        m,
        "ErrorInestabilidadNumerica"
    );
    
    py::class_<ArgonSimulator>(m, "ArgonSimulator")
        .def(py::init<
             int,
             double,
             double,
             double,
             unsigned int,
             bool,
             bool,
             bool>(),
             py::arg("particulas_por_lado"),
             py::arg("densidad_reducida"),
             py::arg("paso_tiempo"),
             py::arg("temp_objetivo"),
             py::arg("semilla") = 0,
             py::arg("corregir_cm") = true,
             py::arg("correccion_presion_cola") = true,
             py::arg("reescalar_velocidades") = true)
        .def("ejecutar", 
            [](ArgonSimulator& self,
               const ConfiguracionSimulacion& cfg,
               const std::optional<std::string>& archivo) {
                try {
                    return self.ejecutar(cfg, archivo);
                } catch (ErrorInestabilidadNumerica& e) {
                    std::cout << "\n[BIND] Excepción capturada en bindings, what(): " 
                            << e.what() << std::flush;
                    std::cout << "\n[BIND] Accediendo a resultados_parciales..." << std::flush;
                    auto rp = e.resultados_parciales;  // copia explícita antes de tocar Python
                    std::cout << "\n[BIND] resultados_parciales.pasos.size() = " 
                            << rp.pasos.size() << std::flush;
                    py::object exc = py::module_::import("simulador_dm.simulador")
                        .attr("ErrorInestabilidadNumerica");
                    py::object exc_inst = exc(e.what());
                    exc_inst.attr("resultados_parciales") = rp;
                    PyErr_SetObject(exc.ptr(), exc_inst.ptr());
                    throw py::error_already_set();
                }
            },
            py::arg("config"),
            py::arg("nombre_archivo") = py::none()
        );
}
