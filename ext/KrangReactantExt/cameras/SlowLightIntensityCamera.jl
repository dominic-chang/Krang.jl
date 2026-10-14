# cameras/SlowLightIntensityCamera.jl
# Reactant specialisations mirroring the corresponding file in `src/`.

function _reactant_slow_light_intensity_pixel(met::Krang.Kerr, α, β, θo)
    tempη = Krang.η(met, α, β, θo)
    tempλ = Krang.λ(met, α, θo)
    roots = Krang.get_radial_roots(met, tempη, tempλ)

    I1, I2, Ip, Im = Krang.radial_inf_integrals(met, roots)
    I0_inf = Krang.Ir_inf(met, roots)
    τ_total = Krang.total_mino_time(met, roots)
    Iϕ_inf_temp = Krang.Iϕ_inf(met, roots, tempλ)
    It_inf_temp = Krang.It_inf(met, roots, tempλ)
    absGθo_Gθhat = Krang._absGθo_Gθhat(met, θo, tempη, tempλ)
    absGϕo_Gϕhat = Krang._absGϕo_Gϕhat(met, θo, tempη, tempλ)
    absGto_Gthat = Krang._absGto_Gthat(met, θo, tempη, tempλ)

    return Krang.SlowLightIntensityPixel(
        met,
        (α, β),
        roots,
        I0_inf,
        τ_total,
        Iϕ_inf_temp,
        It_inf_temp,
        I1,
        I2,
        Ip,
        Im,
        absGθo_Gθhat,
        absGϕo_Gϕhat,
        absGto_Gthat,
        θo,
        tempη,
        tempλ,
    )
end

Reactant.@reactant_overlay function Krang.SlowLightIntensityPixel(met::Krang.Kerr, α, β, θo)
    return _reactant_slow_light_intensity_pixel(met, α, β, θo)
end
