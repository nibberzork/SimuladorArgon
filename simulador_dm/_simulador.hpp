#ifndef SIMULADOR_HPP   // Si no se ha definido SIMULADOR_HPP
#define SIMULADOR_HPP

#include <vector>
#include <string>
#include <optional>

/**
 * @brief Estructura que contiene todas las variables del sistema de partículas
 * 
 * Almacena posiciones, velocidades, aceleraciones y propiedades termodinámicas
 * en unidades reducidas de Lennard-Jones.
 */
struct SistemaParticulas {
    int num_particulas; ///< Número total de partículas (N)
    double longitud_caja; ///< Longitud del lado de la caja cúbica (L*)

    // Arrays SoA para posiciones, velocidades y aceleraciones
    std::vector<double> rx, ry, rz; ///< Posiciones x,y,z de cada partícula
    std::vector<double> vx, vy, vz; ///< Velocidades x,y,z de cada partícula  
    std::vector<double> ax, ay, az; ///< Aceleraciones x,y,z de cada partícula


    // Variables termodinámicas
    double energia_potencial = 0.0;     ///< Energía potencial total (U*)
    double energia_cinetica = 0.0;      ///< Energía cinética total (K*)
    double temperatura_inst = 0.0;      ///< Temperatura instantánea (T*)
    double presion_inst = 0.0;          ///< Presión instantánea (P*)
    double presion_correccion_cola = 0.0; ///< Corrección de presión de largo alcance
    double virial = 0.0;                ///< Virial para cálculo de presión (W*)
    
    // Parámetros precalculados de cell lists
    int nc = 0;         ///< Número de celdas por lado (nc × nc × nc celdas totales)
    double celda_tam = 0.0; ///< Tamaño de cada celda (L*/nc)

    std::vector<int> cabeza; ///< Lista que apunta a la primera particula de una celda
    std::vector<int> lista; ///< Lista que apunta a la siguiente particula de una celda, -1 si no hay más
};

struct ConfiguracionSimulacion {
    int num_pasos = 25000; ///< Número total de pasos de tiempo
    int pasos_equilibrado = 1000; ///< Número de pasos de tiempo para equilibrar
    int frecuencia_muestreo = 50; ///< Frecuencia de muestreo de variables
};

struct ResultadosSimulacion {
    std::vector <int>       pasos;
    std::vector <double>    tiempos;
    std::vector <double>    temperaturas;
    std::vector <double>    presiones;
    std::vector <double>    energias_potenciales;
    std::vector <double>    energias_cineticas;
    std::vector <double>    energias_totales;
};

/**
 * @brief Simulador de dinámica molecular para sistemas de Argón
 * 
 * Implementa integración de Verlet con optimización de cell lists O(n)
 * para simulaciones de dinámica molecular en el ensemble NVT.
 * Utiliza unidades reducidas de Lennard-Jones.
 */
class ArgonSimulator {
public:
    /**
     * @brief Constructor del simulador
     * @param particulas_por_lado Número de partículas por dimensión (n³ = N)
     * @param densidad_reducida Densidad reducida ρ* = ρ σ³
     * @param paso_tiempo Paso de tiempo dt*
     * @param temp_objetivo Temperatura objetivo T*
     * @param semilla Semilla opcional para reproducibilidad; si vale 0 se usa una aleatoria
     * @param corregir_cm Activa la corrección del momento lineal total
     * @param correccion_presion_cola Activa la corrección de cola en la presión
     * @param reescalar_velocidades Activa el termostato por reescalado durante el equilibrado
     */
    ArgonSimulator(
        int particulas_por_lado,
        double densidad_reducida,
        double paso_tiempo,
        double temp_objetivo,
        unsigned int semilla = 0,
        bool corregir_cm = true,
        bool correccion_presion_cola = true,
        bool reescalar_velocidades = true
    );

    /**
     * @brief Ejecuta la simulación completa
     * @param config Configuración de la simulación
     * @param nombre_archivo Archivo de salida opcional para guardar trayectorias
     * @return ResultadosSimulacion con los datos muestreados
     */
    ResultadosSimulacion ejecutar(const ConfiguracionSimulacion& config, const std::optional<std::string>& nombre_archivo = std::nullopt);

    /**
     * @brief Obtiene la temperatura instantánea actual
     * @return Temperatura reducida T*
     */
    double get_current_temperature() const { return sistema.temperatura_inst; };
    
    /**
     * @brief Obtiene la presión instantánea actual
     * @return Presión reducida P*
     */
    double get_current_pressure() const { return sistema.presion_inst; };

private:
    SistemaParticulas sistema;     ///< Estado completo del sistema
    double dt;                     ///< Paso de tiempo reducido
    double temp_referencia;        ///< Temperatura objetivo para termostato
    unsigned int semilla;          ///< Semilla opcional para reproducibilidad
    bool corregir_cm;              ///< Activa o desactiva la corrección del centro de masa
    bool correccion_presion_cola;  ///< Activa o desactiva la corrección de cola en presión
    bool reescalar_velocidades;    ///< Activa o desactiva el reescalado de velocidades en equilibrado
    
    // Constantes físicas de Lennard-Jones
    static constexpr double R_CORTE = 2.5;                  ///< Radio de corte rc* 
    static constexpr double R_CORTE_SQ = R_CORTE * R_CORTE; ///< rc*^2 (para optimizar cálculos)
    static constexpr double COEF_FUERZA = 48.0;             ///< Coeficiente para fuerza LJ
    static constexpr double COEF_ENERGIA = 4.0;             ///< Coeficiente para energía potencial LJ
    
    // Métodos privados para las diferentes etapas de la simulación
    /**
     * @brief Inicializa posiciones y velocidades del sistema
     * 
     * Coloca las partículas en una red cúbica simple y asigna velocidades
     * aleatorias. Puede corregir el momento lineal total y reescalar las
     * velocidades iniciales para arrancar cerca de la temperatura objetivo.
     * 
     * @note Se ejecuta automáticamente al iniciar la simulación
     * @warning Requiere que el constructor haya sido llamado primero
     */
    void inicializar_sistema();
    
    /**
     * @brief Corrige el momento lineal total del sistema
     * 
     * Calcula la velocidad del centro de masa y la resta a todas las
     * partículas para asegurar que el momento lineal total sea cero.
     * Evita la deriva no física del centro de masa cuando esta corrección
     * está habilitada.
     * 
     * @note Asume que todas las partículas tienen la misma masa (m = 1)
     * @note Se ejecuta automáticamente en inicializar_sistema()
     */
    void corregir_momento_lineal();
    
    /**
     * @brief Calcula fuerzas entre todas las partículas usando cell lists
     * 
     * Implementa el cálculo de fuerzas de Lennard-Jones. Usa cell lists
     * cuando la discretización espacial es segura y recurre a un barrido
     * exhaustivo de pares cuando el número de celdas por lado es pequeño.
     * Actualiza las aceleraciones, energía potencial y virial del sistema.
     * 
     * @note Fuerza LJ: F = 48 * [(1/r)^14 - 0.5*(1/r)^8] * r_vec
     * @note Energía: U = 4 * [(1/r)^12 - (1/r)^6]
     * @note Calcula también el virial: W = Σ r_ij · f_ij
     * @warning Requiere posiciones inicializadas
     * @warning Se ejecuta antes y después de la integración Verlet
     */
    void calcular_fuerzas();
    
    
    /**
     * @brief Integra las ecuaciones de movimiento usando el algoritmo de Verlet
     * 
     * Actualiza posiciones usando: r(t+dt) = 2r(t) - r(t-Δt) + a(t)*dt²
     * Y velocidades usando: v(t+dt) = [r(t+dt) - r(t-dt)] / (2*dt)
     * 
     * @note Requiere posiciones actuales, anteriores y aceleraciones
     * @note Se ejecuta dos veces por ciclo de simulación
     * @warning Debe alternar con cálculo de fuerzas para convergencia correcta
     */
    void integracion_verlet();
    
    /**
     * @brief Calcula propiedades termodinámicas del sistema
     * 
     * Computa temperatura, energía cinética y presión utilizando las
     * definiciones estadísticas en el ensemble NVT. Realiza promedios
     * sobre todas las partículas e incluye la corrección de cola solo si
     * está activada.
     * 
     * @note K = (1/2) Σ m*v² = (1/2) Σ v² (para m = 1)
     * @note T = (2/3) K / N
     * @note P = ρ*T + (1/3V) Σ r_ij·f_ij + P_cola
     * @note Se ejecuta al final de cada paso de integración
     */
    void propiedades_termodinamicas();
    
    /**
     * @brief Escala velocidades para mantener temperatura constante
     * 
     * Implementa un termostato simple que reescala suavemente las
     * velocidades de todas las partículas hacia la temperatura objetivo.
     * Utilizado durante la fase de equilibración cuando el reescalado está
     * habilitado para estabilizar el sistema.
     * 
     * @note factor = √(T_ref / T_actual)
     * @note Se aplica solo durante los primeros pasos (fase de equilibrado)
     * @warning No debe usarse después del equilibrado para evitar sesgos
     */
    void escalar_velocidades();
    
    /**
     * @brief Calcula la corrección de presión por truncamiento del potencial
     * 
     * Corrige el efecto de truncar el potencial de Lennard-Jones en el
     * radio de corte. Si la funcionalidad está desactivada, deja la
     * corrección almacenada a cero como salvaguarda.
     * 
     * @note Fórmula: P_cola = (16π/3) * ρ² * [(2/3)/rc⁹ - 1/rc³]
     * @note Se ejecuta una sola vez en el constructor
     * @note Depende únicamente de ρ* y rc*, no cambia durante la simulación
     */
    void calcular_presion_cola();
};

#endif // SIMULADOR_HPP
