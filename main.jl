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
const I₀ = 100.0         # corriente por bobina de Helmholtz [A]
const R  = 0.60           # radio [m] (60 cm)
const dz = R/2            # bobinas en z = ±R/2; separación entre planos = R (60 cm)
const ν  = 1.0 / μ₀
const ε  = ν * 1e-6      # regularización de gauge
const σ  = 0.030          # grosor del anillo [m]; ancho para que la malla lo resuelva

const amp = I₀ / (2π * σ^2) # Normalizamos
println("amp = ", round(amp, sigdigits=4), " A/m²")

# Magnetorquers. Diámetro, largo y holgura coinciden con gmsh/magnetorquers.geo
const mt_diameter = 0.010   # [m]  varilla, 10 mm
const mt_length   = 0.140   # [m]  140 mm
const mt_gap      = 0.003   # [m]
const N_mt  = 100.0         # espiras por solenoide
const I_mtx = 100.0         # corriente DC, eje X [A]
const I_mtz = 100.0         # corriente DC, eje Z [A]

# --- Planteamiento ------------------------------------------
include(joinpath(@__DIR__, "src", "stationary", "helmholtz_magnetorquers.jl"))

# Flecha de longitud fija en la dirección del torque.
# El vector guardado es el torque real, en N·m.
function write_torque_arrow(path, τ)
    mag = sqrt(τ ⋅ τ)
    u = mag > 0 ? τ / mag : VectorValue(0.0, 1.0, 0.0)
    side = u × VectorValue(0.0, 0.0, 1.0)
    if sqrt(side ⋅ side) < 1e-8
        side = u × VectorValue(1.0, 0.0, 0.0)
    end
    side = side / sqrt(side ⋅ side)
    L = 0.40
    tip = L * u
    base = tip - 0.10 * u
    left = base + 0.05 * side
    right = base - 0.05 * side
    pts = (VectorValue(0.0, 0.0, 0.0), tip, left, right)
    open(path, "w") do io
        println(io, "# vtk DataFile Version 3.0")
        println(io, "Torque resultante")
        println(io, "ASCII")
        println(io, "DATASET POLYDATA")
        println(io, "POINTS 4 float")
        for p in pts
            println(io, p[1], " ", p[2], " ", p[3])
        end
        println(io, "LINES 3 9")
        println(io, "2 0 1")
        println(io, "2 1 2")
        println(io, "2 1 3")
        println(io, "POINT_DATA 4")
        println(io, "VECTORS torque float")
        for _ in 1:4
            println(io, τ[1], " ", τ[2], " ", τ[3])
        end
        println(io, "SCALARS torque_Nm float 1")
        println(io, "LOOKUP_TABLE default")
        for _ in 1:4
            println(io, mag)
        end
    end
end

# --- Solución -----------------------------------------------
println("Resolviendo campo de Helmholtz...")
Ah_H = solve(AffineFEOperator(a, l_H, U, V))
Bh_H = ∇ × Ah_H

println("Resolviendo Helmholtz + magnetorquers...")
Ah = solve(AffineFEOperator(a, l, U, V))

# --- Resultados ---------------------------------------------
Bh = ∇ × Ah
Hh = (1/μ₀) * Bh
normB(B) = sqrt(B⋅B)
logB(B) = log10(sqrt(B⋅B) + 1e-8)
Bmag = normB ∘ Bh
Bmag_H = normB ∘ Bh_H
Blog = logB ∘ Bh

r_cf = CellField(x -> VectorValue(x[1], x[2], x[3]), Ω)
J_mt = CellField(x -> J_mtx(x) + J_mtz(x), Ω)
tau_density = r_cf × (J_mt × Bh_H)

W = 0.5 * sum(∫( Bh⋅Hh )dΩ)
println("Energía magnética     : W = $(round(W, sigdigits=4)) J")

B_centro = μ₀ * (4/5)^(3/2) * I₀ / R
println("B centro (teórico)    : $(round(B_centro*1e6, sigdigits=4)) µT")

# Torque de interacción: la J de cada magnetorquer con el campo del par
τx = torque_on(J_mtx, Bh_H, Ωx)
τz = torque_on(J_mtz, Bh_H, Ωz)
τ  = τx + τz
B_H = VectorValue(0.0, 0.0, B_centro)
τ_m = (m_x + m_z) × B_H

fmt(v) = "($(round(v[1], sigdigits=4)), $(round(v[2], sigdigits=4)), $(round(v[3], sigdigits=4))) N·m"
println("Torque magnetorquer X : ", fmt(τx))
println("Torque magnetorquer Z : ", fmt(τz))
println("Torque resultante     : ", fmt(τ))
println("Torque m × B          : ", fmt(τ_m))

mkpath("results")
writevtk(Ω, "results/helmholtz_magnetorquers",
    cellfields = ["A"=>Ah, "B"=>Bh, "B_helmholtz"=>Bh_H,
                  "H"=>Hh, "normB"=>Bmag, "normB_helmholtz"=>Bmag_H,
                  "logB"=>Blog, "J"=>J_cf, "torque_density"=>tau_density])
write_torque_arrow("results/torque.vtk", τ)

println("Simulación finalizada.")
println("ParaView: normB_helmholtz muestra el par; logB muestra los dos campos.")
println("          torque.vtk es la flecha del torque. torque_density, en las varillas.")
