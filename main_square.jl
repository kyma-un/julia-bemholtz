using Gridap
using Gridap.Algebra
using Gridap.FESpaces
using GridapGmsh
using Gmsh: gmsh
using LinearAlgebra
using SparseArrays
using IterativeSolvers    # Pkg.add("IterativeSolvers") si no lo tienes
using AlgebraicMultigrid  # Pkg.add("AlgebraicMultigrid") si no lo tienes

# ============================================================
# Bobina de Helmholtz TRIAXIAL CUADRADA — fuente gaussiana
# ∇×(ν ∇×A) + ε A = J,  J = suma de 3 pares de lazos CUADRADOS
# ============================================================

const μ₀ = 4π * 1e-7   # constante física real, no depende de L

# --- Factores de anidado (deben coincidir EXACTO con el .geo) ---
const FX = 1.00
const FY = 1.08
const FZ = 1.16

# --- Separación Helmholtz óptima para bobina CUADRADA ------------
# (NO es L/2 como en la circular; es ≈0.5445*L, ver .geo)
const KOPT = 0.5445

# --- Direcciones de eje y de los dos lados locales (u,v) ---------
# Deben coincidir con cómo se orientó cada cuadrado en el .geo:
#   Par Z: cuadrado en plano XY  -> (û,v̂) = (X,Y)
#   Par X: cuadrado en plano YZ  -> (û,v̂) = (Y,Z)
#   Par Y: cuadrado en plano XZ  -> (û,v̂) = (Z,X)
const EX = VectorValue(1.0, 0.0, 0.0)
const EY = VectorValue(0.0, 1.0, 0.0)
const EZ = VectorValue(0.0, 0.0, 1.0)

"""
    nearest_square_boundary(u, v, s)

Distancia mínima desde el punto (u,v) al contorno de un cuadrado de
semi-lado `s` centrado en el origen, y la dirección tangente (sentido
antihorario, regla de la mano derecha respecto a ê) en el punto más
cercano del contorno. Generaliza la dirección azimutal φ̂ que se usaba
para el anillo circular.
"""
function nearest_square_boundary(u::Float64, v::Float64, s::Float64)
    vc = clamp(v, -s, s)
    uc = clamp(u, -s, s)
    d_right  = sqrt((u-s)^2  + (v-vc)^2)
    d_left   = sqrt((u+s)^2  + (v-vc)^2)
    d_top    = sqrt((u-uc)^2 + (v-s)^2)
    d_bottom = sqrt((u-uc)^2 + (v+s)^2)
    dists = (d_right, d_left, d_top, d_bottom)
    tans  = ((0.0,1.0), (0.0,-1.0), (-1.0,0.0), (1.0,0.0))
    imin = argmin(dists)
    return dists[imin], tans[imin][1], tans[imin][2]
end

"""
    square_pair_field(x, ê, û, v̂, Lside, d, σ, amp)

Campo fuente (gaussiano) de UN par Helmholtz CUADRADO cuyo eje de
simetría es `ê`, con el lazo definido en el plano (û,v̂), lado `Lside`,
y las dos bobinas centradas en ±d a lo largo de `ê`. Generaliza
`ring_pair_field` (circular) reemplazando el anillo por un cuadrado.
"""
function square_pair_field(x::VectorValue, ê::VectorValue, û::VectorValue, v̂::VectorValue,
                            Lside::Float64, d::Float64, σ::Float64, amp::Float64)
    xax = x ⋅ ê
    u = x ⋅ û
    v = x ⋅ v̂
    s = Lside / 2
    dist2d, tu, tv = nearest_square_boundary(u, v, s)

    d1 = dist2d^2 + (xax - d)^2   # lazo "+"
    d2 = dist2d^2 + (xax + d)^2   # lazo "-"
    g  = exp(-d1/(2σ^2)) + exp(-d2/(2σ^2))

    tangent = tu*û + tv*v̂
    return amp * g * tangent
end

"""
    B_theo_square(x, L, D, N, I)

Ecuación (8): campo axial teórico de un par Helmholtz CUADRADO de
lado `L`, separación total `D` (= 2*d, distancia entre las dos
bobinas), `N` vueltas y corriente `I`, evaluado en la posición `x`
a lo largo del eje del par.
"""
function B_theo_square(x::Float64, L::Float64, D::Float64, N::Float64, I::Float64)
    halfL = L/2
    τp = (x + D/2) / halfL
    τm = (x - D/2) / halfL
    denom_p = (τp^2 + 1) * sqrt(τp^2 + 2)
    denom_m = (τm^2 + 1) * sqrt(τm^2 + 2)
    return (4*μ₀*N*I/(π*L)) * (1/denom_p + 1/denom_m)
end

"""
    generar_malla(L, geo_path, msh_path)

Genera la malla triaxial cuadrada para un lado base `L` usando Gmsh.jl.
IMPORTANTE: `L` se fija con `onelab.setNumber("Parameters/L", ...)` ANTES
de abrir el .geo (que lo declara con `L = DefineNumber[...]`).
"""
function generar_malla(L::Float64, geo_path::String, msh_path::String)
    isfile(geo_path) || error("No encuentro el .geo en: $(abspath(geo_path))")

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)   # 0 para silenciar

    gmsh.onelab.setNumber("Parameters/L", [L])     # 1) fijar L
    gmsh.open(geo_path)                            # 2) abrir el .geo con esa L

    length(gmsh.model.getEntities(3)) == 0 &&
        error("El .geo no generó ningún volumen (revisa errores de sintaxis arriba)")

    gmsh.model.mesh.generate(3)
    gmsh.write(msh_path)
    gmsh.finalize()
end

# --- Rutas (edita SOLO aquí si cambias de carpeta) ---------------
const GEO_PATH = "gmsh/triaxial_square.geo"
const MESH_DIR = "gmsh"
const OUT_DIR  = "results"

"""
    run_helmholtz_square(L; Ix, Iy, Iz, Nx, Ny, Nz, σ_frac, ...)

Corre el pipeline completo (malla -> FEM -> B -> export) para un lado
base `L` (bobina Z), corrientes `Ix,Iy,Iz` y número de vueltas
`Nx,Ny,Nz` en cada eje. Para un solo eje, deja las otras dos
corrientes en 0.0.

`σ_frac` es la fracción de L usada como ancho físico de la fuente
gaussiana (ver .geo: sigma_src = L/σ_frac). DEBE coincidir con el
valor usado en el .geo para que la malla resuelva bien la fuente.
"""
function run_helmholtz_square(L::Float64;
        Ix::Float64 = 1.0, Iy::Float64 = 1.0, Iz::Float64 = 1.0,
        Nx::Float64 = 1.0, Ny::Float64 = 1.0, Nz::Float64 = 1.0,
        σ_frac::Float64 = 40.0,
        geo_path::String = GEO_PATH, meshdir::String = MESH_DIR, outdir::String = OUT_DIR)

    mkpath(meshdir); mkpath(outdir)

    tag      = "L$(round(L*1000, digits=1))mm"
    msh_path = joinpath(meshdir, "triaxial_square_$tag.msh")

    # --- (a) Malla ---------------------------------------------
    generar_malla(L, geo_path, msh_path)
    model = GmshDiscreteModel(msh_path)

    order = 1
    reffe = ReferenceFE(nedelec, Float64, order)
    V = TestFESpace(model, reffe; conformity = :HCurl, dirichlet_tags = ["boundary"])
    U = TrialFESpace(V, VectorValue(0.0, 0.0, 0.0))

    Ω  = Triangulation(model)
    dΩ = Measure(Ω, 2*order + 1)
    println("L = $L m | Celdas: ", num_cells(Ω), " | DOFs libres: ", num_free_dofs(V))

    # --- (b) Geometría de cada par (debe coincidir con el .geo) --
    Lx, Ly, Lz = L*FX, L*FY, L*FZ
    dx, dy, dz = KOPT*Lx, KOPT*Ly, KOPT*Lz   # separación (centro a cada bobina)
    σ = L / σ_frac

    ν = 1.0 / μ₀
    ε = ν * 1e-6

    amp_x = Ix / (2π * σ^2)
    amp_y = Iy / (2π * σ^2)
    amp_z = Iz / (2π * σ^2)

    # --- (c) Fuente total: superposición de los 3 pares cuadrados --
    function Jsrc(x)
        Jx = square_pair_field(x, EX, EY, EZ, Lx, dx, σ, amp_x)   # plano YZ
        Jy = square_pair_field(x, EY, EZ, EX, Ly, dy, σ, amp_y)   # plano ZX
        Jz = square_pair_field(x, EZ, EX, EY, Lz, dz, σ, amp_z)   # plano XY
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

    # --- Solver iterativo (CG + AMG) ------------------------------
    println("Construyendo precondicionador AMG...")
    ml = ruge_stuben(A)
    Pl = aspreconditioner(ml)

    println("Resolviendo con CG (AMG precond.)...")
    x, hist = IterativeSolvers.cg(A, b; Pl=Pl, reltol=1e-8,
                                   maxiter=2000, log=true, verbose=true)
    println("  CG: convergió=$(hist.isconverged) en $(hist.iters) iteraciones")

    if !hist.isconverged
        error("El solver NO convergió: los resultados de este L no son válidos. " *
              "Sube maxiter, relaja reltol, o revisa la calidad de la malla.")
    end

    Ah = FEFunction(U, x)
    println("✓ Resuelto")

    # --- (e) Postproceso ------------------------------------------
    Bh = ∇ × Ah
    Hh = (1/μ₀) * Bh
    normB(B) = sqrt(B⋅B)
    Bmag = normB ∘ Bh

    W = 0.5 * sum(∫( Bh⋅Hh )dΩ)

    # --- Validación: B simulado en el centro vs. teórico (ec. 8) --
    pt = Point(0.0, 0.0, 0.0)
    Bsim = Bh(pt)

    Bc_x = B_theo_square(0.0, Lx, 2*dx, Nx, Ix)
    Bc_y = B_theo_square(0.0, Ly, 2*dy, Ny, Iy)
    Bc_z = B_theo_square(0.0, Lz, 2*dz, Nz, Iz)

    err_x = 100*abs(Bsim[1]-Bc_x)/Bc_x
    err_y = 100*abs(Bsim[2]-Bc_y)/Bc_y
    err_z = 100*abs(Bsim[3]-Bc_z)/Bc_z

    println("  W = $(round(W, sigdigits=4)) J")
    println("  B centro TEÓRICO  : X=$(round(Bc_x*1e6,sigdigits=4)) µT | ",
            "Y=$(round(Bc_y*1e6,sigdigits=4)) µT | ",
            "Z=$(round(Bc_z*1e6,sigdigits=4)) µT")
    println("  B centro SIMULADO : X=$(round(Bsim[1]*1e6,sigdigits=4)) µT | ",
            "Y=$(round(Bsim[2]*1e6,sigdigits=4)) µT | ",
            "Z=$(round(Bsim[3]*1e6,sigdigits=4)) µT")
    println("  Error relativo    : X=$(round(err_x,sigdigits=3))% | ",
            "Y=$(round(err_y,sigdigits=3))% | ",
            "Z=$(round(err_z,sigdigits=3))%")
    if max(err_x,err_y,err_z) > 15
        @warn "Error > 15% en L=$L: sospecha de malla insuficiente cerca de las bobinas"
    end

    writevtk(Ω, joinpath(outdir, "magnetostatic_square_$tag"),
        cellfields = ["A"=>Ah, "B"=>Bh, "H"=>Hh, "normB"=>Bmag, "J"=>J_cf])
    println("✓ Exportado $(joinpath(outdir, "magnetostatic_square_$tag.vtu"))")

    return (L=L, W=W,
            Bc_x=Bc_x, Bc_y=Bc_y, Bc_z=Bc_z,
            Bsim_x=Bsim[1], Bsim_y=Bsim[2], Bsim_z=Bsim[3],
            err_x=err_x, err_y=err_y, err_z=err_z,
            cg_iters=hist.iters, cg_converged=hist.isconverged)
end

# ============================================================
# --- Corrida principal: bobina de L=1 m ------------------------
# ============================================================
resultado = run_helmholtz_square(1.0; Ix=1.0, Iy=1.0, Iz=1.0)

println("\n--- Resultado (SIMULADO vs teórico) ---")
println("L=$(resultado.L) m | CG: conv=$(resultado.cg_converged) iters=$(resultado.cg_iters) | W=$(round(resultado.W, sigdigits=4)) J")
println("  Bz: sim=$(round(resultado.Bsim_z*1e6,sigdigits=4)) µT  teo=$(round(resultado.Bc_z*1e6,sigdigits=4)) µT  err=$(round(resultado.err_z,sigdigits=3))%")
println("  Bx: sim=$(round(resultado.Bsim_x*1e6,sigdigits=4)) µT  teo=$(round(resultado.Bc_x*1e6,sigdigits=4)) µT  err=$(round(resultado.err_x,sigdigits=3))%")
println("  By: sim=$(round(resultado.Bsim_y*1e6,sigdigits=4)) µT  teo=$(round(resultado.Bc_y*1e6,sigdigits=4)) µT  err=$(round(resultado.err_y,sigdigits=3))%")

# ============================================================
# Para barrer varios tamaños de bobina, descomenta:
#
# L_values = [0.5, 0.75, 1.0, 1.25, 1.5]
# resultados = [run_helmholtz_square(L) for L in L_values]
#
# Para simular UN SOLO EJE, pon las otras dos corrientes en 0:
#   run_helmholtz_square(1.0; Ix=0.0, Iy=0.0, Iz=1.0)
# ============================================================
