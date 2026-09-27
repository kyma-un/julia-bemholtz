# Planteamiento de sistema estacionario...
# --- Malla y espacios ---------------------------------------
model = GmshDiscreteModel("meshes/torus_coil.msh")

order = 1
reffe = ReferenceFE(nedelec, Float64, order)
V = TestFESpace(model, reffe; conformity = :HCurl, dirichlet_tags = ["boundary"])
U = TrialFESpace(V, VectorValue(0.0, 0.0, 0.0))

Ω  = Triangulation(model)
dΩ = Measure(Ω, 2*order + 1)
println("Celdas: ", num_cells(Ω), " | DOFs libres: ", num_free_dofs(V))

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

# --- Formas variacionales -----------------------------------
a(u,v) = ∫( ν*((∇×u)⋅(∇×v)) + ε*(u⋅v) )dΩ
l(v)   = ∫( v⋅J_cf )dΩ
