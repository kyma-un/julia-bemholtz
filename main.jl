using Gridap
using GridapGmsh
using LinearAlgebra

model = GmshDiscreteModel("meshes/torus_coil.msh")

order = 1
reffe = ReferenceFE(nedelec, Float64, order)
V = TestFESpace(model, reffe;
                conformity    = :HCurl,
                dirichlet_tags = ["boundary"])
U = TrialFESpace(V, VectorValue(0.0, 0.0, 0.0))

# Dominio completo — para la forma bilineal
Ω  = Triangulation(model)
dΩ = Measure(Ω, 2*order + 1)

# Dominio solo de las bobinas — para la fuente J
Ωcoil  = Triangulation(model; tags = ["coil1", "coil2"])
dΩcoil = Measure(Ωcoil, 2*order + 1)

# Parámetros físicos
const μ₀ = 4π * 1e-7   # permeabilidad del vacío [H/m]
const I₀ = 1.0          # corriente total [A]
const r_w = 0.008       # radio sección del hilo [m]  ← igual que en el .geo
const J₀  = I₀ / (π * r_w^2)   # densidad de corriente [A/m²]

# Dirección azimutal: ê_φ = (-y/r, x/r, 0)
# Es la dirección real de la corriente en una bobina circular
function J(x)
    r = sqrt(x[1]^2 + x[2]^2)
    r < 1e-10 && return VectorValue(0.0, 0.0, 0.0)   # evitar singularidad en r=0
    return J₀ * VectorValue(-x[2]/r, x[1]/r, 0.0)
end

# Forma bilineal — sobre TODO el dominio (aire + bobinas)
a(u,v) = ∫( (1/μ₀) * (∇×u)⋅(∇×v) )dΩ

# Forma lineal — SOLO sobre las bobinas
l(v) = ∫( v⋅J )dΩcoil

# Resolver
op = AffineFEOperator(a, l, U, V)
Ah = solve(op)

# Postproceso
Bh = ∇ × Ah
Hh = (1/μ₀) * Bh

mkpath("results")
writevtk(Ω, "results/magnetostatic",
    cellfields = [
        "A" => Ah,
        "B" => Bh,
        "H" => Hh
    ])

println("✓ Exportado: results/magnetostatic.vtu")

