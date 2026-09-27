# Kuntur

Magnetostática por elementos finitos de una bobina de Helmholtz, como banco de campo uniforme para el ADCS de un CubeSat.

## Resultados

<img src="assets/field-h.png" alt="Magnitud de H en el dominio de aire" width="100%">

*Magnitud de H. El par concentra el campo entre las bobinas.*

<img src="assets/mesh-oblique.png" alt="Malla, vista oblicua" width="49%"> <img src="assets/mesh-axial.png" alt="Malla, vista axial" width="49%">

*Malla del dominio de aire. Vista oblicua y vista axial, con refinamiento sobre las bobinas y en el centro del par.*

## Ejecución

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. main.jl
```

Los parámetros, la solución y el VTK están en `main.jl`. El planteamiento estacionario está en `src/stationary/helmholtz.jl`.

MIT · KYMA UN
