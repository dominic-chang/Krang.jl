# cameras/IntensityCamera.jl
# Reactant specialisations mirroring the corresponding file in `src/`.

function _reactant_intensity_pixel(met::Krang.Kerr, α, β, θo)
    tempη = Krang.η(met, α, β, θo)
    tempλ = Krang.λ(met, α, θo)
    roots = get_radial_roots(met, tempη, tempλ)

    I0_inf = Krang.Ir_inf(met, roots)
    τ_total = Krang.total_mino_time(met, roots)
    absGθo_Gθhat = Krang._absGθo_Gθhat(met, θo, tempη, tempλ)

    return Krang.IntensityPixel(
        met,
        (α, β),
        roots,
        I0_inf,
        τ_total,
        absGθo_Gθhat,
        θo,
        tempη,
        tempλ,
    )
end
