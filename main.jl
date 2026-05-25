using Gridap
using GridapGmsh
using LinearAlgebra

# Hola sofi

# ============================================================
# Bobina de Helmholtz — Formulación magnetostática
# ∇×(1/μ ∇×A) = J   con gauge de Coulomb ∇·A = 0
# ============================================================

# --- Malla --------------------------------------------------
model = GmshDiscreteModel("meshes/torus_coil.msh")

# --- Espacio de funciones -----------------------------------
order = 1
reffe = ReferenceFE(nedelec, Float64, order)
V = TestFESpace(model, reffe;
                conformity     = :HCurl,
                dirichlet_tags = ["boundary"])
U = TrialFESpace(V, VectorValue(0.0, 0.0, 0.0))

println("Grados de libertad libres : ", num_free_dofs(V))

# --- Triangulaciones ----------------------------------------
Ω      = Triangulation(model)
dΩ     = Measure(Ω, 2*order + 1)

Ωcoil  = Triangulation(model; tags = ["coil1", "coil2"])
dΩcoil = Measure(Ωcoil, 2*order + 1)

Ωair   = Triangulation(model; tags = ["air"])
dΩair  = Measure(Ωair, 2*order + 1)

# --- Parámetros físicos -------------------------------------
const μ₀  = 4π * 1e-7        # permeabilidad vacío     [H/m]
const I₀  = 1.0               # corriente por espira    [A]
const R   = 0.100             # radio de la bobina      [m]
const r_w = 0.008             # radio sección del hilo  [m]
const J₀  = I₀ / (π * r_w^2) # densidad de corriente   [A/m²]
const α   = 1.0 / μ₀         # penalización gauge       [m/H]

const μ_air  = μ₀
const μ_coil = μ₀   # cobre: μ_r ≈ 1

# --- Fuente de corriente ------------------------------------
# J azimutal: ê_φ = (-y/r, x/r, 0) — dirección real en una espira
function J(x)
    r = sqrt(x[1]^2 + x[2]^2)
    r < 1e-10 && return VectorValue(0.0, 0.0, 0.0)
    return J₀ * VectorValue(-x[2]/r, x[1]/r, 0.0)
end

# --- Formas variacionales -----------------------------------
a(u,v) = ∫( (1/μ_air)  * (∇×u)⋅(∇×v)
           + α          * (∇⋅u) * (∇⋅v) )dΩ

l(v)   = ∫( v⋅J )dΩcoil

# --- Ensamblaje y solución ----------------------------------
println("Ensamblando sistema...")
op = AffineFEOperator(a, l, U, V)

println("Resolviendo...")
Ah = solve(op)
println("✓ Solución obtenida")

# --- Postproceso físico -------------------------------------
Bh = ∇ × Ah                  # densidad de flujo  B = ∇×A    [T]
Hh = (1/μ₀) * Bh             # intensidad de campo H = B/μ   [A/m]

# Energía magnética almacenada: W = (1/2)∫ B·H dΩ
W = 0.5 * sum(∫( Bh⋅Hh )dΩ)
println("Energía magnética almacenada : W = $(round(W, sigdigits=4)) J")

# Campo B teórico en el centro (fórmula analítica Helmholtz)
B_teorico = μ₀ * (4/5)^(3/2) * I₀ / R
println("B en centro (teórico)        : $(round(B_teorico*1e6, sigdigits=4)) μT")
println("Nota: comparar con |B| en (0,0,0) en ParaView")

# --- Exportar a ParaView ------------------------------------
mkpath("results")
writevtk(Ω, "results/magnetostatic",
    cellfields = [
        "A"  => Ah,    # potencial vector         [Wb/m]
        "B"  => Bh,    # densidad de flujo        [T]
        "H"  => Hh,    # intensidad de campo      [A/m]
    ])

println("✓ Exportado : results/magnetostatic.vtu")
println("  Campos disponibles en ParaView: A, B, H")
println("  Para |B|: Filters → Calculator → mag(B)")
println("  Para isosuperficies: Filters → Contour → selecciona mag(B)")