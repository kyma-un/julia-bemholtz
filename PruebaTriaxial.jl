using Gridap
using Gridap.Algebra
using Gridap.FESpaces
using GridapGmsh
using Gmsh: gmsh
using LinearAlgebra
using SparseArrays
using IterativeSolvers   # Pkg.add("IterativeSolvers") si no lo tienes

# ============================================================
# Bobina de Helmholtz TRIAXIAL — fuente volumétrica gaussiana
# ∇×(ν ∇×A) + ε A = J,  J = suma de 3 pares de anillos (X,Y,Z)
# Barrido de R generando la malla vía Gmsh.jl (sin binario externo)
# ============================================================

const μ₀ = 4π * 1e-7   # constante física real, no depende de R

# --- Rutas (edita SOLO aquí si cambias de carpeta) -----------
const GEO_PATH = "gmsh/triaxial_coil.geo"   # ruta al .geo triaxial
const MESH_DIR = "gmsh"                     # dónde se guardan los .msh
const OUT_DIR  = "results"                  # dónde se guardan los .vtu

# --- Factores de anidado (deben coincidir con el .geo) -------
const FX = 1.00
const FY = 1.08
const FZ = 1.16

# --- Direcciones de eje (vectores unitarios) ------------------
const EX = VectorValue(1.0, 0.0, 0.0)
const EY = VectorValue(0.0, 1.0, 0.0)
const EZ = VectorValue(0.0, 0.0, 1.0)

"""
    ring_pair_field(x, ê, Reje, d, σ, amp)

Campo fuente (gaussiano azimutal) de UN par Helmholtz cuyo eje de simetría
es `ê` (vector unitario), con bobinas centradas en ±d a lo largo de `ê`.
Generaliza tu `Jsrc` original (que asumía ê = ẑ) a cualquier eje.
"""
function ring_pair_field(x::VectorValue, ê::VectorValue, Reje::Float64,
                          d::Float64, σ::Float64, amp::Float64)
    xax = (x ⋅ ê)                 # coordenada a lo largo del eje
    ρvec = x - xax * ê             # componente perpendicular al eje
    r = sqrt(ρvec ⋅ ρvec)
    r < 1e-12 && return VectorValue(0.0, 0.0, 0.0)
    ρ̂ = ρvec / r
    # dirección azimutal = ê × ρ̂  (right-hand rule respecto al eje)
    φ̂ = VectorValue(
        ê[2]*ρ̂[3] - ê[3]*ρ̂[2],
        ê[3]*ρ̂[1] - ê[1]*ρ̂[3],
        ê[1]*ρ̂[2] - ê[2]*ρ̂[1],
    )
    d1 = (r - Reje)^2 + (xax - d)^2   # anillo "+"
    d2 = (r - Reje)^2 + (xax + d)^2   # anillo "-"
    g  = exp(-d1/(2σ^2)) + exp(-d2/(2σ^2))
    return amp * g * φ̂
end

"""
    generar_malla(R, geo_path, msh_path)

Genera la malla triaxial para un radio base `R` usando Gmsh.jl.
IMPORTANTE: `R` se fija con `onelab.setNumber("Parameters/R", ...)` ANTES
de abrir el .geo (que lo declara con `R = DefineNumber[...]`).
"""
function generar_malla(R::Float64, geo_path::String, msh_path::String)
    # gmsh.open NO da error si el archivo no existe: deja el modelo vacío.
    isfile(geo_path) || error("No encuentro el .geo en: $(abspath(geo_path))")

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)   # 0 para silenciar

    gmsh.onelab.setNumber("Parameters/R", [R])     # 1) fijar R
    gmsh.open(geo_path)                            # 2) abrir el .geo con ese R

    length(gmsh.model.getEntities(3)) == 0 &&
        error("El .geo no generó ningún volumen (revisa errores de sintaxis arriba)")

    gmsh.model.mesh.generate(3)
    gmsh.write(msh_path)
    gmsh.finalize()
end

"""
    run_helmholtz_triaxial(R; Ix, Iy, Iz, σ, ...)

Corre el pipeline completo (malla -> FEM -> B -> export) para un radio
base `R` y corrientes `Ix,Iy,Iz` en cada eje. Para simular un solo eje
(como tu script uniaxial original) deja las otras dos corrientes en 0.0.
"""
function run_helmholtz_triaxial(R::Float64;
        Ix::Float64 = 1.0, Iy::Float64 = 1.0, Iz::Float64 = 1.0,
        σ::Float64 = 0.006,
        geo_path::String = GEO_PATH, meshdir::String = MESH_DIR, outdir::String = OUT_DIR)

    mkpath(meshdir); mkpath(outdir)

    tag       = "R$(round(R*1000, digits=1))mm"
    msh_path  = joinpath(meshdir, "triaxial_coil_$tag.msh")

    # --- (a) Malla ---------------------------------------------
    generar_malla(R, geo_path, msh_path)
    model = GmshDiscreteModel(msh_path)

    order = 1
    reffe = ReferenceFE(nedelec, Float64, order)
    V = TestFESpace(model, reffe; conformity = :HCurl, dirichlet_tags = ["boundary"])
    U = TrialFESpace(V, VectorValue(0.0, 0.0, 0.0))

    Ω  = Triangulation(model)
    dΩ = Measure(Ω, 2*order + 1)
    println("R = $R m | Celdas: ", num_cells(Ω), " | DOFs libres: ", num_free_dofs(V))

    # --- (b) Geometría de cada par (debe coincidir con el .geo) --
    Rx, Ry, Rz = R*FX, R*FY, R*FZ
    dx, dy, dz = Rx/2, Ry/2, Rz/2

    ν = 1.0 / μ₀
    ε = ν * 1e-6

    amp_x = Ix / (2π * σ^2)
    amp_y = Iy / (2π * σ^2)
    amp_z = Iz / (2π * σ^2)

    # --- (c) Fuente total: superposición de los 3 pares ----------
    function Jsrc(x)
        Jx = ring_pair_field(x, EX, Rx, dx, σ, amp_x)
        Jy = ring_pair_field(x, EY, Ry, dy, σ, amp_y)
        Jz = ring_pair_field(x, EZ, Rz, dz, σ, amp_z)
        return Jx + Jy + Jz
    end
    J_cf = CellField(Jsrc, Ω)

    # --- (d) Forma variacional -----------------------------------
    a(u,v) = ∫( ν*((∇×u)⋅(∇×v)) + ε*(u⋅v) )dΩ
    l(v)   = ∫( v⋅J_cf )dΩ

    println("Ensamblando...")
    op = AffineFEOperator(a, l, U, V)
    A = get_matrix(op)
    b = get_vector(op)
    println("  Sistema: $(size(A,1)) DOFs, $(nnz(A)) no-ceros ",
            "(~$(round(nnz(A)*12/1e6, digits=1)) MB en CSC)")

    # --- Solver iterativo (CG) en vez de LU directo ---------------
    # El sistema es simétrico definido positivo (curl-curl + término
    # de masa ε>0), así que CG es aplicable. Esto evita el relleno
    # (fill-in) del solver directo, que en 3D agota la RAM fácilmente
    # incluso con unos cientos de miles de DOFs.
    println("Resolviendo con CG (Jacobi precond.)...")
    Pl = Diagonal(A)   # precondicionador de Jacobi; si converge muy
                        # lento, considera AlgebraicMultigrid.jl (AMG)
    x, hist = IterativeSolvers.cg(A, b; Pl=Pl, reltol=1e-8,
                                   maxiter=5000, log=true, verbose=true)
    println("  CG: convergió=$(hist.isconverged) en $(hist.iters) iteraciones")
    Ah = FEFunction(U, x)
    println("✓ Resuelto")

    # --- (e) Postproceso ------------------------------------------
    Bh = ∇ × Ah
    Hh = (1/μ₀) * Bh
    normB(B) = sqrt(B⋅B)
    Bmag = normB ∘ Bh

    W = 0.5 * sum(∫( Bh⋅Hh )dΩ)

    # B teórico en el centro para cada eje, asumiendo pares
    # desacoplados ideales (aproximación: el anidado introduce un
    # pequeño acoplamiento geométrico que este valor no captura)
    Bc_x = μ₀ * (4/5)^(3/2) * Ix / Rx
    Bc_y = μ₀ * (4/5)^(3/2) * Iy / Ry
    Bc_z = μ₀ * (4/5)^(3/2) * Iz / Rz

    println("  W = $(round(W, sigdigits=4)) J")
    println("  B centro teórico: X=$(round(Bc_x*1e6,sigdigits=4)) µT | ",
            "Y=$(round(Bc_y*1e6,sigdigits=4)) µT | ",
            "Z=$(round(Bc_z*1e6,sigdigits=4)) µT")

    writevtk(Ω, joinpath(outdir, "magnetostatic_$tag"),
        cellfields = ["A"=>Ah, "B"=>Bh, "H"=>Hh, "normB"=>Bmag, "J"=>J_cf])
    println("✓ Exportado $(joinpath(outdir, "magnetostatic_$tag.vtu"))")

    return (R=R, W=W, Bc_x=Bc_x, Bc_y=Bc_y, Bc_z=Bc_z)
end

# ============================================================
# --- Barrido de R --------------------------------------------
# ============================================================
R_values = [0.05, 0.075, 0.100, 0.125, 0.150]

resultados = [run_helmholtz_triaxial(R; Ix=1.0, Iy=1.0, Iz=1.0) for R in R_values]

println("\n--- Resumen del barrido ---")
for r in resultados
    println("R=$(r.R) m -> W=$(round(r.W, sigdigits=4)) J | ",
            "Bx=$(round(r.Bc_x*1e6,sigdigits=4)) µT, ",
            "By=$(round(r.Bc_y*1e6,sigdigits=4)) µT, ",
            "Bz=$(round(r.Bc_z*1e6,sigdigits=4)) µT")
end

# ============================================================
# Para simular UN SOLO EJE (equivalente a tu script uniaxial
# original), pon las otras dos corrientes en 0, por ejemplo:
#
#   run_helmholtz_triaxial(0.100; Ix=0.0, Iy=0.0, Iz=1.0)
#
# ============================================================
