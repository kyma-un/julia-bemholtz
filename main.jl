using Gridap
using GridapGmsh
using LinearAlgebra

# ============================================================
# Bobina de Helmholtz — fuente volumétrica gaussiana (robusta)
# ∇×(ν ∇×A) + ε A = J,  J = anillo gaussiano azimutal
# ============================================================

# --- Malla --------------------------------------------------
model = GmshDiscreteModel("gmsh/torus_coil.msh")

order = 1
reffe = ReferenceFE(nedelec, Float64, order)
V = TestFESpace(model, reffe; conformity = :HCurl, dirichlet_tags = ["boundary"])
U = TrialFESpace(V, VectorValue(0.0, 0.0, 0.0))

Ω  = Triangulation(model)
dΩ = Measure(Ω, 2*order + 1)
println("Celdas: ", num_cells(Ω), " | DOFs libres: ", num_free_dofs(V))

# --- Parámetros físicos -------------------------------------
const μ₀ = 4π * 1e-7
const I₀ = 1.0           # corriente por bobina [A]
const R  = 0.100          # radio de las bobinas [m]
const dz = R/2            # posición Helmholtz: z = ±R/2
const ν  = 1.0 / μ₀
const ε  = ν * 1e-6      # regularización de gauge
const σ  = 0.006          # grosor del anillo gaussiano [m]

# Normalización: ∫∫ amp·exp(-ρ²/2σ²) dA = I₀  sobre la sección
# La integral de una gaussiana 2D es 2πσ², por tanto amp = I₀/(2πσ²)
const amp = I₀ / (2π * σ^2)

# --- Fuente: dos anillos azimutales gaussianos --------------
function Jsrc(x)
    r = sqrt(x[1]^2 + x[2]^2)
    r < 1e-12 && return VectorValue(0.0, 0.0, 0.0)
    φ̂ = VectorValue(-x[2]/r, x[1]/r, 0.0)   # dirección azimutal
    # distancia² al centro de cada anillo en el plano (r, z)
    d1 = (r - R)^2 + (x[3] - dz)^2          # anillo superior
    d2 = (r - R)^2 + (x[3] + dz)^2          # anillo inferior
    g  = exp(-d1/(2σ^2)) + exp(-d2/(2σ^2))
    return amp * g * φ̂
end
J_cf = CellField(Jsrc, Ω)

# Verifica la corriente total inyectada (debería ≈ I₀ por bobina)
# Flujo de J a través del plano y=0, x>0 (corta ambos anillos una vez cada uno)
println("amp = ", round(amp, sigdigits=4), " A/m²")

# --- Formas variacionales -----------------------------------
a(u,v) = ∫( ν*((∇×u)⋅(∇×v)) + ε*(u⋅v) )dΩ
l(v)   = ∫( v⋅J_cf )dΩ

# --- Solución -----------------------------------------------
println("Ensamblando y resolviendo...")
op = AffineFEOperator(a, l, U, V)
Ah = solve(op)
println("✓ Resuelto")

# --- Postproceso --------------------------------------------
Bh = ∇ × Ah
Hh = (1/μ₀) * Bh
normB(B) = sqrt(B⋅B)
Bmag = normB ∘ Bh

W = 0.5 * sum(∫( Bh⋅Hh )dΩ)
println("Energía magnética     : W = $(round(W, sigdigits=4)) J")

B_centro = μ₀ * (4/5)^(3/2) * I₀ / R
println("B centro (teórico)    : $(round(B_centro*1e6, sigdigits=4)) µT")

# --- Exportar -----------------------------------------------
mkpath("results")
writevtk(Ω, "results/magnetostatic",
    cellfields = ["A"=>Ah, "B"=>Bh, "H"=>Hh, "normB"=>Bmag, "J"=>J_cf])
println("✓ Exportado results/magnetostatic.vtu")
println("  Verificar: Plot Over Line eje z (0,0,-0.15)→(0,0,0.15), componente B_Z")
println("  Debe verse meseta plana en el centro ≈ $(round(B_centro*1e6,sigdigits=3)) µT")