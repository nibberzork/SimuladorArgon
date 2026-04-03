#include "_simulador.hpp"

#include <cmath>
#include <random>
#include <iostream>
#include <fstream>
#include <iomanip>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/**
 * @brief Constructor del simulador de dinámica molecular de Argón
 * 
 * Inicializa el simulador con los parámetros dados en unidades reducidas de Lennard-Jones.
 * Calcula automáticamente el tamaño de la caja, reserva memoria para todas las partículas
 * y precalcula los parámetros de cell lists para optimización O(n).
 * 
 * @param particulas_por_lado Número de partículas por dimensión del cubo (n³ = N total)
 * @param densidad_reducida Densidad reducida ρ* = ρ σ³ (típicamente 0.8-1.2)
 * @param paso_tiempo Paso de tiempo reducido Δt* (típicamente 0.001-0.01)
 * @param temp_objetivo Temperatura objetivo reducida T* (típicamente 1.0-2.0)
 * 
 * @note Las unidades están en términos reducidos de Lennard-Jones donde:
 *       - Longitud: σ (diámetro atómico)
 *       - Energía: ε (profundidad del pozo)
 *       - Masa: m (masa atómica)
 *       - Tiempo: σ√(m/ε)
 */
ArgonSimulator::ArgonSimulator(
    int particulas_por_lado = 8,
    double densidad_reducida = 0.84, 
    double paso_tiempo = 0.005, 
    double temp_objetivo = 1.002,
    unsigned int semilla
) {
    this->dt = paso_tiempo;
    this->temp_referencia = temp_objetivo;
    this->semilla = semilla;
    
    
    sistema.num_particulas = particulas_por_lado * particulas_por_lado * particulas_por_lado; // N = n^3
    sistema.longitud_caja = std::pow(static_cast<double>(sistema.num_particulas) / densidad_reducida, 1.0/3.0); // L* = (N / rho*)^(1/3)

    // Precalcular parámetros de las listas enlazadas
    sistema.nc = static_cast<int>(std::floor(sistema.longitud_caja / R_CORTE));
    sistema.celda_tam = sistema.longitud_caja / sistema.nc;

    const int N = sistema.num_particulas; // Número de partículas
    const int num_celdas = sistema.nc * sistema.nc * sistema.nc; // Número de celdas
    

    // Reservar memoria para los vectores (SoA)
    sistema.rx.resize(N); sistema.ry.resize(N); sistema.rz.resize(N);
    sistema.vx.resize(N); sistema.vy.resize(N); sistema.vz.resize(N);
    sistema.ax.resize(N); sistema.ay.resize(N); sistema.az.resize(N);

    // Calcular corrección de cola una sola vez
    calcular_presion_cola();

    // Asignar espacio para la estructura de listas enlazadas
    sistema.cabeza.resize(num_celdas);
    sistema.lista.resize(N);
}



/**
 * @brief Inicializa posiciones y velocidades del sistema
 * 
 * Coloca las partículas en una red cúbica simple centrada en el origen
 * y asigna velocidades aleatorias con distribución Maxwell-Boltzmann.
 * Corrige el momento lineal total para evitar deriva del centro de masa.
 * 
 * @note Las posiciones se generan en una malla cúbica uniforme
 * @note Las velocidades siguen distribución gaussiana con varianza T*
 * @warning Debe llamarse después del constructor antes de cualquier simulación
 */
void ArgonSimulator::inicializar_sistema() {
    const int N = sistema.num_particulas; // Número de partículas
    const double L = sistema.longitud_caja; // Longitud de la caja
    
    // Inicializar posiciones en una red cúbica simple
    int n = std::cbrt(sistema.num_particulas); // Número de partículas por lado
    double paso = sistema.longitud_caja / n; // Paso entre partículas
    
    int index = 0;
    for (int i = 0; i < n; ++i) { // Bucle sobre x
        for (int j = 0; j < n; ++j) { // Bucle sobre y
            for (int k = 0; k < n; ++k) { // Bucle sobre z
                if (index < sistema.num_particulas) {
                    sistema.rx[index] = i * paso;
                    sistema.ry[index] = j * paso;
                    sistema.rz[index] = k * paso;
                    ++index;
                }
            }
        }
    }
    

    std::random_device rd; // Semilla aleatoria en base a la entropía del hardware
    std::mt19937 gen(this->semilla != 0 ? this->semilla : rd());
    std::uniform_real_distribution<double> distribucion(0.0, 1.0);

    // Inicializar velocidades aleatorias con distribución de Maxwell-Boltzmann
    for (int i = 0; i < sistema.num_particulas; ++i) {
        sistema.vx[i] = distribucion(gen);
        sistema.vy[i] = distribucion(gen);
        sistema.vz[i] = distribucion(gen);
    }

    // Corregir momento lineal para evitar deriva del centro de masa
    corregir_momento_lineal();

    // Reescalar velocidades para compensar la pérdida de temperatura por corrección del momento
    // La corrección del momento reduce la temperatura, así que la reescalamos
    double temp_actual = 0.0;
    for (int i = 0; i < N; ++i) {
        temp_actual += sistema.vx[i]*sistema.vx[i] +
                       sistema.vy[i]*sistema.vy[i] + 
                       sistema.vz[i]*sistema.vz[i];
    }
    temp_actual *= 0.5 / N; // (1/2)<v^2> = (3/2)T

    if (temp_actual > 0.0) {
        double factor = std::sqrt(temp_referencia / temp_actual);
        for (int i = 0; i < N; ++i) {
            sistema.vx[i] *= factor;
            sistema.vy[i] *= factor;
            sistema.vz[i] *= factor;
        }
        corregir_momento_lineal();
    }

    // Inicializar las cell list
    std::fill(sistema.cabeza.begin(), sistema.cabeza.end(), -1);
    std::fill(sistema.lista.begin(), sistema.lista.end(), -1);
}

/**
 * @brief Corrige el momento lineal total del sistema
 * 
 * Calcula la velocidad del centro de masa y la resta a todas las partículas
 * para asegurar que el momento lineal total sea cero. Esto evita la deriva
 * no física del centro de masa durante la simulación.
 * 
 * @note Asume que todas las partículas tienen la misma masa (masa reducida = 1)
 * @note Se ejecuta automáticamente en inicializar_sistema()
 */
void ArgonSimulator::corregir_momento_lineal() {
    double vx_cm = 0.0, vy_cm = 0.0, vz_cm = 0.0;
    
    // Calcular velocidade del centro de masa (Asumiendo masa igual para todas las partículas)
    for (int i = 0; i < sistema.num_particulas; ++i) {
        vx_cm += sistema.vx[i];
        vy_cm += sistema.vy[i];
        vz_cm += sistema.vz[i];
    }

    vx_cm /= sistema.num_particulas;
    vy_cm /= sistema.num_particulas;
    vz_cm /= sistema.num_particulas;

    for (int i = 0; i < sistema.num_particulas; ++i) {
        sistema.vx[i] -= vx_cm;
        sistema.vy[i] -= vy_cm;
        sistema.vz[i] -= vz_cm;
    }
}


/**
 * @brief Calcula fuerzas entre todas las partículas usando cell lists O(n)
 * 
 * Implementa el algoritmo de cell lists 3D optimizado para calcular las fuerzas
 * de Lennard-Jones. Divide el espacio en celdas de tamaño ~rc y solo evalúa
 * interacciones entre partículas en celdas adyacentes. Utiliza la tercera ley
 * de Newton para calcular cada par solo una vez.
 * 
 * @note Complejidad: O(n) en lugar de O(n²) gracias a cell lists
 * @note Fuerza de Lennard-Jones: F = 48 * [(1/r)^14 - 0.5*(1/r)^8] * r_vec
 * @note Energía potencial: U = 4 * [(1/r)^12 - (1/r)^6]
 * @note Calcula también el virial para la presión: W = Σ r_ij · f_ij
 * 
 * @warning Requiere posiciones inicializadas
 * @warning Modifica ax, ay, az, energia_potencial y virial del sistema
 */
void ArgonSimulator::calcular_fuerzas() {

    const int N = sistema.num_particulas;
    const double L = sistema.longitud_caja;


    // Reiniciar aceleraciones y energía potencial
    std::fill(sistema.ax.begin(), sistema.ax.end(), 0.0);
    std::fill(sistema.ay.begin(), sistema.ay.end(), 0.0);
    std::fill(sistema.az.begin(), sistema.az.end(), 0.0);
    sistema.energia_potencial = 0.0;
    sistema.virial = 0.0;   //  r_ij * f^ij, necesario para la presión


    // Función para calcular fuerzas entre un par de partículas
    auto calcular_par_fuerzas = [&](int i, int j) {

        double dx = sistema.rx[i] - sistema.rx[j];
        double dy = sistema.ry[i] - sistema.ry[j];
        double dz = sistema.rz[i] - sistema.rz[j];

        // Imagen mínima usando redondeo para manejar desplazamientos
        // mayores que una longitud de caja de forma robusta.
        dx -= L * std::round(dx / L);
        dy -= L * std::round(dy / L);
        dz -= L * std::round(dz / L);


        const double r2 = dx*dx + dy*dy + dz*dz;
        if (r2 >= R_CORTE_SQ || r2 < 1e-10) return;

        const double r2inv  = 1.0 / r2;
        const double r6inv  = r2inv * r2inv * r2inv;   // 1/r^6
        const double r12inv = r6inv * r6inv;            // 1/r^12

        const double factor = COEF_FUERZA * r2inv * (r12inv - 0.5 * r6inv);
                
        // Limitar fuerzas extremadamente grandes para estabilidad
        double fax = factor * dx;
        double fay = factor * dy;
        double faz = factor * dz;

        // Segun la 3ra ley de newton, toda acción tiene
        sistema.ax[i] += fax;  sistema.ay[i] += fay;  sistema.az[i] += faz;
        sistema.ax[j] -= fax;  sistema.ay[j] -= fay;  sistema.az[j] -= faz;

        // Energía potencial: U* = 4*(r^{-12} - r^{-6})
        sistema.energia_potencial += COEF_ENERGIA * (r12inv - r6inv);

        // Virial: r_ij * f_ij = factor * r^2       (f_ij = factor * r_ij)
        sistema.virial += factor * r2;
    };

    const int nc = sistema.nc;

    if (nc < 3) {
        // Si los pares estan muy juntos y no hay espacio se calculan todos los pares
        for (int i = 0; i < N; ++i)
        for (int j = i + 1; j < N; ++j)
            calcular_par_fuerzas(i, j);
        return;
    } else {
        // Listas de celdas
        // Aprovechando que las interacciones solo ocurren a < rc, se discretiza el espacio
        // en celdas de esa longitud de forma que solo se miran las vecinas
        const double celda_tam = sistema.celda_tam;

        // Resetear cell lists (inicializar a -1)
        std::fill(sistema.cabeza.begin(), sistema.cabeza.end(), -1);
        std::fill(sistema.lista.begin(), sistema.lista.end(), -1);



        // Construcción de la estructura de listas enlazadas, se toma la última particula
        // como cabecera y se va enlazando con la partícula anterior
        for (int i = 0; i < N; ++i) {
            int cx = static_cast<int>(sistema.rx[i] / celda_tam); // Indice de celda respecto a x
            int cy = static_cast<int>(sistema.ry[i] / celda_tam); // Indice de celda respecto a y
            int cz = static_cast<int>(sistema.rz[i] / celda_tam); // Indice de celda respecto a z

            // Si hay particulas en el borde exacto, se debe restar 1 para evitar desbordamiento en lista[]
            cx = std::min(cx, nc - 1);
            cy = std::min(cy, nc - 1);
            cz = std::min(cz, nc - 1);

            const int idx = cx + cy * nc + cz * nc * nc;
            sistema.lista[i]    = sistema.cabeza[idx];   // el nuevo nodo apunta al anterior primero
            sistema.cabeza[idx] = i;                     // ahora i es la cabeza de la celda
        }


        // Semi-vecindad: 13 desplazamientos + celda propia
        static const int SEMI_NB[13][3] = {
            // dz = +1 (9 celdas)
            {-1, -1, +1}, { 0, -1, +1}, {+1, -1, +1},
            {-1,  0, +1}, { 0,  0, +1}, {+1,  0, +1},
            {-1, +1, +1}, { 0, +1, +1}, {+1, +1, +1},
            // dz = 0, dy = +1 (3 celdas)
            {-1, +1,  0}, { 0, +1,  0}, {+1, +1,  0},
            // dz = 0, dy = 0, dx = +1 (1 celda)
            {+1,  0,  0}
        };

        // Bucle principal sobre celdas
        for (int cz = 0; cz < nc; ++cz) {
        for (int cy = 0; cy < nc; ++cy) {
        for (int cx = 0; cx < nc; ++cx) {

            const int idx_self = cx + cy * nc + cz * nc * nc;

            // Pares dentro de la misma celda (j > i)
            for (int i = sistema.cabeza[idx_self]; i != -1; i = sistema.lista[i]) {
            for (int j = sistema.lista[i];         j != -1; j = sistema.lista[j]) {
                calcular_par_fuerzas(i, j);
            }}

            // Pares con las 13 celdas de la semi-vecindad
            for (int nb = 0; nb < 13; ++nb) {
                const int nx = (cx + SEMI_NB[nb][0] + nc) % nc;
                const int ny = (cy + SEMI_NB[nb][1] + nc) % nc;
                const int nz = (cz + SEMI_NB[nb][2] + nc) % nc;
                const int idx_nb = nx + ny * nc + nz * nc * nc;

                for (int i = sistema.cabeza[idx_self]; i != -1; i = sistema.lista[i]) {
                for (int j = sistema.cabeza[idx_nb];   j != -1; j = sistema.lista[j]) {
                    calcular_par_fuerzas(i, j);
                }}
            }

        }}} // fin bucle celdas cx, cy, cz
    }

    

    // Presión instantánea (sin corrección de cola)
    const double V = L * L * L;
    sistema.presion_inst = sistema.virial / (3.0 * V);
    // Nota: el término N*T*/V* se suma en propiedades_termodinamicas()
}

/**
 * @brief Calcula la corrección de presión de largo alcance
 * 
 * Corrige la presión por el truncamiento del potencial en rc. Esta corrección
 * es necesaria porque el potencial de Lennard-Jones tiene cola de largo alcance.
 * Se calcula una sola vez ya que solo depende de ρ* y rc*.
 * 
 * @note Fórmula: P_cola* = (16 * pi/3) * rho*^2 * [(2/3)*(1/rc*)^9 - (1/rc*)^3]
 * @note Se ejecuta automáticamente en el constructor
 * @note Corrige el efecto de truncar el potencial en r_corte
 */
void ArgonSimulator::calcular_presion_cola() {
    const double rho    = sistema.num_particulas
                        / (sistema.longitud_caja * sistema.longitud_caja * sistema.longitud_caja);
    const double rc3    = R_CORTE * R_CORTE * R_CORTE;   // rc^3
    const double rc9    = rc3 * rc3 * rc3;               // rc^9

    sistema.presion_correccion_cola =
        (16.0 * M_PI / 3.0) * rho * rho
        * ((2.0 / 3.0) / rc9 - 1.0 / rc3);
}

/**
 * @brief Aplica condiciones de contorno periódicas
 * 
 * Asegura que todas las partículas permanezcan dentro de la caja
 * aplicando la convención de imagen mínima.
 */
void ArgonSimulator::aplicar_condiciones_contorno() {
    const double L = sistema.longitud_caja;
    
    for (int i = 0; i < sistema.num_particulas; ++i) {
        // Mantener posiciones en [0, L).
        sistema.rx[i] -= L * std::floor(sistema.rx[i] / L);
        sistema.ry[i] -= L * std::floor(sistema.ry[i] / L);
        sistema.rz[i] -= L * std::floor(sistema.rz[i] / L);
    }
}

/**
 * @brief Integra las ecuaciones de movimiento usando Verlet
 * 
 * Actualiza posiciones y velocidades usando el algoritmo de Verlet:
 * r(t+dt) = 2r(t) - r(t-Δt) + a(t)*dt^2
 * v(t+dt) = [r(t+dt) - r(t-dt)] / (2*dt)
 */
void ArgonSimulator::integracion_verlet() {
    // TODO: cambiar el algoritmo por el de verlet que usa v(dt/2)
    const double half_dt = 0.5 * dt;
    const double L = sistema.longitud_caja;
    
    // Medio paso de la velocidad
    // v(t + dt/2) = v(t) + (1/2)*dt*a(t)
    for (int i = 0; i < sistema.num_particulas; ++i) {
        sistema.vx[i] += sistema.ax[i] * half_dt;
        sistema.vy[i] += sistema.ay[i] * half_dt;
        sistema.vz[i] += sistema.az[i] * half_dt;
    };

    // Paso completo para la posición con velocidad intermedia
    // r(t + dt) = r(t) + v(t + dt/2)*dt
    for (int i = 0; i < sistema.num_particulas; ++i) {
        sistema.rx[i] += sistema.vx[i] * dt;
        sistema.ry[i] += sistema.vy[i] * dt;
        sistema.rz[i] += sistema.vz[i] * dt;
        
        // Aplicar PBC en la caja para mantener en [0, L).
        sistema.rx[i] -= L * std::floor(sistema.rx[i] / L);
        sistema.ry[i] -= L * std::floor(sistema.ry[i] / L);
        sistema.rz[i] -= L * std::floor(sistema.rz[i] / L);
    };

    // Calcular nuevas fuerzas usando r(t + dt)
    calcular_fuerzas();

    // Integrar la velocidad usando las nuevas aceleraciones
    // v(t + dt) = v(t + dt/2) + (1 / 2)*dt*a(t + dt)
    for (int i = 0; i < sistema.num_particulas; i++) {
        sistema.vx[i] += sistema.ax[i] * half_dt;
        sistema.vy[i] += sistema.ay[i] * half_dt;
        sistema.vz[i] += sistema.az[i] * half_dt;
    };

    corregir_momento_lineal(); // Para corregir los errores numéricos
}

/**
 * @brief Calcula propiedades termodinámicas del sistema
 * 
 * Computa temperatura, energía cinética y presión usando las
 * definiciones estadísticas en el ensemble NVT.
 */
void ArgonSimulator::propiedades_termodinamicas() {
    const int N = sistema.num_particulas;
    const double V = sistema.longitud_caja * sistema.longitud_caja * sistema.longitud_caja;
    
    // Energía cinética: K = (1/2) Σ m v² = (3/2) N T (para m=1)
    double kin_energy = 0.0;
    for (int i = 0; i < N; ++i) {
        kin_energy += sistema.vx[i]*sistema.vx[i] + 
                     sistema.vy[i]*sistema.vy[i] + 
                     sistema.vz[i]*sistema.vz[i];
    }
    kin_energy *= 0.5; // (1/2) sum(v_i^2)
    sistema.energia_cinetica = kin_energy;
    
    // Temperatura: T = (2/3) K / N
    sistema.temperatura_inst = (2.0 / 3.0) * kin_energy / N;
    
    // Presión: P = ρ T + (1/(3V)) Σ r_ij · f_ij
    // El término virial ya se calculó en calcular_fuerzas()
    const double rho = N / V;
    sistema.presion_inst = rho * sistema.temperatura_inst + sistema.virial / (3.0 * V);
    sistema.presion_inst += sistema.presion_correccion_cola; // Corrección de largo alcance
}

/**
 * @brief Escala velocidades para mantener temperatura constante
 * 
 * Implementa un termostato simple que reescala suavemente las velocidades
 * hacia la temperatura objetivo con un factor de amortiguamiento.
 */
void ArgonSimulator::escalar_velocidades() {
    // Si por lo que sea hay un bug sobre la temperatura se salta el escalado
    if (sistema.temperatura_inst > 0.0) {
        // Forzar la temperatura según el factor
        const double factor = std::sqrt(temp_referencia / sistema.temperatura_inst);
        
        // Escalado a las particulas
        for (int i = 0; i < sistema.num_particulas; ++i) {
            sistema.vx[i] *= factor;
            sistema.vy[i] *= factor;
            sistema.vz[i] *= factor;
        }
    }
}

/**
 * @brief Ejecuta la simulación completa de dinámica molecular
 * 
 * Realiza una simulación completa con los siguientes pasos:
 * 1. Inicialización del sistema
 * 2. Bucle principal: fuerzas -> integración -> propiedades
 * 3. Guardado de trayectorias cada cierto número de pasos
 * 
 * @param num_pasos Número total de pasos de tiempo a simular
 * @param nombre_archivo Archivo donde guardar las trayectorias
 */
ResultadosSimulacion ArgonSimulator::ejecutar(const ConfiguracionSimulacion& config,
                                              const std::optional<std::string>& nombre_archivo) {
    // Inicializar sistema
    inicializar_sistema();
    
    ResultadosSimulacion resultados;  // Objeto para almacenar resultados
    
    // Abrir archivo de salida solo si se proporciona nombre_archivo
    std::optional<std::ofstream> archivo_salida;
    if (nombre_archivo.has_value()) {
        archivo_salida.emplace(*nombre_archivo);
        if (!archivo_salida->is_open()) {
            throw std::runtime_error("Error: No se pudo abrir el archivo " + *nombre_archivo);
        }
        *archivo_salida << "paso,tiempo,temperatura,presion,energia_pot,energia_cin,energia_tot\n";
    }
    
    // Bucle principal de simulación
    for (int paso = 0; paso < config.num_pasos; ++paso) {
        integracion_verlet();  // Avanzar dt
        
        propiedades_termodinamicas(); // Necesario para poder escalar durante el equilibrado
        // Control de temperatura
        if (paso < config.pasos_equilibrado) {
            escalar_velocidades();
            propiedades_termodinamicas();  // Actualizar tras el escalado
        }

        if (paso % config.frecuencia_muestreo == 0) {
            const double tiempo = paso * dt;
            const double energia_total = sistema.energia_potencial + sistema.energia_cinetica;
            
            // Guardar los datos en la estructura
            resultados.pasos.push_back(paso);
            resultados.tiempos.push_back(tiempo);
            resultados.temperaturas.push_back(sistema.temperatura_inst);
            resultados.presiones.push_back(sistema.presion_inst);
            resultados.energias_potenciales.push_back(sistema.energia_potencial);
            resultados.energias_cineticas.push_back(sistema.energia_cinetica);
            resultados.energias_totales.push_back(energia_total);
            
            // Escribir al archivo si está abierto
            if (archivo_salida.has_value()) {
                *archivo_salida << std::fixed << std::setprecision(6)
                                << paso << ","
                                << tiempo << ","
                                << sistema.temperatura_inst << ","
                                << sistema.presion_inst << ","
                                << sistema.energia_potencial << ","
                                << sistema.energia_cinetica << ","
                                << energia_total << "\n";
            }
            
            // Mostrar progreso
            std::cout << "Paso " << paso << "/" << config.num_pasos 
                      << " - T=" << sistema.temperatura_inst 
                      << " - P=" << sistema.presion_inst << std::endl;
        }
    }
    
    if (archivo_salida.has_value()) {
        archivo_salida->close();
        std::cout << "Simulación completada. Datos guardados en " << *nombre_archivo << std::endl;
    } else {
        std::cout << "Simulación completada. Resultados en memoria." << std::endl;
    }
    
    return resultados;  // Devolver los resultados estructurados
}


