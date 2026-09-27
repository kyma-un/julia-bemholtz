# Kuntur

Simulación magnetostática por elementos finitos de una bobina de Helmholtz, pensada como banco de campo magnético uniforme para el ADCS de un CubeSat.

`main.jl` concentra los parámetros físicos, la solución y la exportación del campo. El planteamiento de elementos finitos de cada caso vive en `src`. Hoy solo está el estacionario, en `src/stationary/helmholtz.jl`. `src/transient` y `src/closed_loop` siguen reservados.

## Formulación

Se resuelve la forma débil regularizada

$$\nabla \times (\nu \nabla \times \mathbf{A}) + \varepsilon \mathbf{A} = \mathbf{J}$$

con elementos de Nédélec de orden 1 en \(H(\mathrm{curl})\) y \(\mathbf{A} = \mathbf{0}\) en la frontera del dominio de aire. La corriente de cada bobina es un anillo azimutal de perfil gaussiano, normalizado para que la integral de la sección valga \(I_0\).

| Símbolo | Valor | Significado |
| --- | --- | --- |
| \(R\) | 0.100 m | Radio de cada bobina |
| \(d\) | \(R/2\) | Separación Helmholtz |
| \(I_0\) | 1 A | Corriente por bobina |
| \(\sigma\) | 6 mm | Grosor del anillo gaussiano |
| \(\nu\) | \(1/\mu_0\) | Reluctividad del aire |
| \(\varepsilon\) | \(\nu \times 10^{-6}\) | Regularización de gauge |

El campo en el centro tiene solución analítica

$$B_z(0) = \mu_0 \left(\frac{4}{5}\right)^{3/2} \frac{I_0}{R}$$

y el script la imprime junto con la energía magnética \(\tfrac12 \int \mathbf{B}\cdot\mathbf{H}\,\mathrm{d}V\).

## Estructura

```
.
├── main.jl                  # parámetros, solve y exportación VTK
├── Project.toml
├── Manifest.toml            # entorno fijado en Julia 1.11.7
├── gmsh/torus_coil.geo      # dominio esférico y curvas de las bobinas
├── src/
│   ├── stationary/
│   │   └── helmholtz.jl     # malla, espacios, fuente y formas débiles
│   ├── transient/           # reservado
│   └── closed_loop/         # reservado
├── config/                  # reservado
├── cad/                     # reservado
└── results/                 # salidas VTK
```

`src/stationary/helmholtz.jl` arma la malla, los espacios de Nédélec y la fuente de dos anillos gaussianos, y deja listas las formas `a` y `l`. `main.jl` las incluye, resuelve y escribe el VTK. La geometría de Gmsh describe las bobinas como curvas filamentarias embebidas en una esfera de aire de radio \(3R\), con refinamiento junto a las bobinas y en el centro del par.

## Requisitos

- [Julia](https://julialang.org/) 1.11
- [Gmsh](https://gmsh.info/) 4.x, solo si hay que regenerar la malla
- [ParaView](https://www.paraview.org/), para revisar el campo

Dependencias del entorno Julia, declaradas en `Project.toml`:

- [Gridap](https://github.com/gridap/Gridap.jl) y [GridapGmsh](https://github.com/gridap/GridapGmsh.jl) para el problema de elementos finitos
- CairoMakie, disponible para postproceso; el caso actual exporta VTK y no lo usa

## Ejecución

Desde la raíz del repositorio:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'

mkdir -p meshes
cp gmsh/torus_coil.msh meshes/torus_coil.msh

julia --project=. main.jl
```

El solver lee `meshes/torus_coil.msh`. El archivo versionado está en `gmsh/`; la copia anterior lo deja donde `main.jl` lo espera.

Para regenerar la malla:

```bash
mkdir -p meshes
gmsh gmsh/torus_coil.geo -3 -o meshes/torus_coil.msh
```

La corrida imprime el número de celdas y grados de libertad, la amplitud de la fuente, la energía magnética y \(B_z\) teórico en el centro. Escribe `results/magnetostatic.vtu` con \(\mathbf{A}\), \(\mathbf{B}\), \(\mathbf{H}\), \(|\mathbf{B}|\) y \(\mathbf{J}\).

## Comprobación

Abrir el VTK en ParaView y trazar **Plot Over Line** sobre el eje \(z\), de \((0,0,-0.15)\) a \((0,0,0.15)\), componente \(B_z\). En el centro debe aparecer una meseta cercana al valor teórico que imprime el script (del orden de unos µT para \(I_0 = 1\,\mathrm{A}\) y \(R = 0.1\,\mathrm{m}\)).

## Licencia

MIT. Copyright (c) 2026 KYMA UN. Ver [LICENSE](LICENSE).
