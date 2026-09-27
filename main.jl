using Gridap
using GridapGmsh
using LinearAlgebra

# Bemholtz ===================================================
#-
# @brief Simulación de interacción EM para ADCS en satelite cubesat
# utilizando elementos finitos (FEM) en Julia usando Gridap.
#
# Pertinencia para el grupo:
#  - Control: Software in The Loop para pruebas de Control
#  - Sensores: Distribución de campo esperada para validación con Magnetómetros
#
# @author: Sofia Vera y Andrés Morales
# @project: Kuntur
#
#
# --- Parámetros físicos -------------------------------------

const μ₀ = 4π * 1e-7
const I₀ = 1.0           # corriente por bobina [A]
const R  = 0.100          # radio de las bobinas [m]
const dz = R/2            # posición Helmholtz: z = ±R/2
const ν  = 1.0 / μ₀
const ε  = ν * 1e-6      # regularización de gauge
const σ  = 0.006          # grosor del anillo gaussiano [m]

const amp = I₀ / (2π * σ^2) # Normalizamos
println("amp = ", round(amp, sigdigits=4), " A/m²")

# --- Planteamiento ------------------------------------------
include(joinpath(@__DIR__, "src", "stationary", "helmholtz.jl"))

# --- Solución -----------------------------------------------
#
println("Ensamblando y resolviendo...")
op = AffineFEOperator(a, l, U, V)
Ah = solve(op)

# --- Resultados ---------------------------------------------
#
Bh = ∇ × Ah
Hh = (1/μ₀) * Bh
normB(B) = sqrt(B⋅B)
Bmag = normB ∘ Bh

W = 0.5 * sum(∫( Bh⋅Hh )dΩ)
println("Energía magnética     : W = $(round(W, sigdigits=4)) J")

B_centro = μ₀ * (4/5)^(3/2) * I₀ / R
println("B centro (teórico)    : $(round(B_centro*1e6, sigdigits=4)) µT")

mkpath("results")
writevtk(Ω, "results/magnetostatic",
    cellfields = ["A"=>Ah, "B"=>Bh, "H"=>Hh, "normB"=>Bmag, "J"=>J_cf])

println("Simulación finalizada. Revisar en PawaraView...");
