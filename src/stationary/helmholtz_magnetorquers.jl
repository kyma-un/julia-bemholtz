using Gridap.CellData: change_domain

# ============================================================
# Estacionario — Helmholtz + magnetorquers en DC
# ∇×(ν ∇×A) + ε A = J
#   J = anillos de Helmholtz + corriente azimutal de cada solenoide
#
# Lo incluye main.jl, que define antes:
#   ν, ε, R, dz, σ, amp
#   mt_diameter, mt_length, mt_gap, N_mt, I_mtx, I_mtz
# ============================================================

model = GmshDiscreteModel("gmsh/helmholtz_magnetorquers.msh")

order = 1
reffe = ReferenceFE(nedelec, Float64, order)
V = TestFESpace(model, reffe; conformity = :HCurl, dirichlet_tags = ["boundary"])
U = TrialFESpace(V, VectorValue(0.0, 0.0, 0.0))

Ω  = Triangulation(model)
dΩ = Measure(Ω, 2*order + 1)
println("Celdas: ", num_cells(Ω), " | DOFs libres: ", num_free_dofs(V))

labels = get_face_labeling(model)
Ωx = Triangulation(model, labels, tags="magnetorquer_x")
Ωz = Triangulation(model, labels, tags="magnetorquer_z")

# Geometría de los solenoides, la misma que gmsh/magnetorquers.geo
const mt_radius = mt_diameter / 2
const z_join = mt_length/2 + mt_radius + mt_gap
const z_mid  = z_join / 2
const c_x = VectorValue(0.0, 0.0,  z_mid)
const c_z = VectorValue(0.0, 0.0, -z_mid)
const ax_x = VectorValue(1.0, 0.0, 0.0)
const ax_z = VectorValue(0.0, 0.0, 1.0)
const half_L = mt_length / 2

# Jφ uniforme en el cilindro, normalizada a N·I espiras-amperio
const Jφ_x = N_mt * I_mtx / (mt_length * mt_radius)
const Jφ_z = N_mt * I_mtz / (mt_length * mt_radius)

function J_helmholtz(x)
    r = sqrt(x[1]^2 + x[2]^2)
    r < 1e-12 && return VectorValue(0.0, 0.0, 0.0)
    φ̂ = VectorValue(-x[2]/r, x[1]/r, 0.0)
    d1 = (r - R)^2 + (x[3] - dz)^2
    d2 = (r - R)^2 + (x[3] + dz)^2
    g  = exp(-d1/(2σ^2)) + exp(-d2/(2σ^2))
    return amp * g * φ̂
end

function J_solenoid(x, center, axis, Jφ)
    rel = VectorValue(x[1], x[2], x[3]) - center
    s = rel ⋅ axis
    abs(s) > half_L && return VectorValue(0.0, 0.0, 0.0)
    radial = rel - s * axis
    ρ2 = radial ⋅ radial
    (ρ2 > mt_radius^2 || ρ2 < 1e-16) && return VectorValue(0.0, 0.0, 0.0)
    return Jφ * ((axis × radial) / sqrt(ρ2))
end

J_mtx(x) = J_solenoid(x, c_x, ax_x, Jφ_x)
J_mtz(x) = J_solenoid(x, c_z, ax_z, Jφ_z)
Jsrc(x)  = J_helmholtz(x) + J_mtx(x) + J_mtz(x)

J_H  = CellField(J_helmholtz, Ω)
J_cf = CellField(Jsrc, Ω)

a(u,v) = ∫( ν*((∇×u)⋅(∇×v)) + ε*(u⋅v) )dΩ
l_H(v) = ∫( v⋅J_H )dΩ
l(v)   = ∫( v⋅J_cf )dΩ

# Momento del solenoide con esta J: m = N I π R² / 3 , a lo largo del eje
const m_x = (N_mt * I_mtx * π * mt_radius^2 / 3) * ax_x
const m_z = (N_mt * I_mtz * π * mt_radius^2 / 3) * ax_z

# τ = ∫ r × (J × B) dV sobre el volumen del magnetorquer
function torque_on(Jfun, Bfield, tri)
    dV = Measure(tri, 2*order + 1)
    r = CellField(x -> VectorValue(x[1], x[2], x[3]), tri)
    J = CellField(Jfun, tri)
    B = change_domain(Bfield, tri, PhysicalDomain())
    t = r × (J × B)
    ex = VectorValue(1.0, 0.0, 0.0)
    ey = VectorValue(0.0, 1.0, 0.0)
    ez = VectorValue(0.0, 0.0, 1.0)
    VectorValue(sum(∫(t⋅ex)dV), sum(∫(t⋅ey)dV), sum(∫(t⋅ez)dV))
end
