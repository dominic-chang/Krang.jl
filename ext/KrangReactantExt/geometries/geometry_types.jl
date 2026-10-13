# geometries/geometry_types.jl
# Reactant specialisations mirroring the corresponding file in `src/`.

Reactant.@reactant_overlay function Krang._raytrace(
    observation,
    pix::Krang.SlowLightIntensityPixel,
    mesh::Krang.Mesh{<:Krang.ConeGeometry,<:Krang.AbstractMaterial};
    res,
)
    geometry = mesh.geometry
    material = mesh.material
    θs = geometry.opening_angle
    subimgs = material.subimgs

    isindir = false
    for _ = 1:2
        isindir ⊻= true
        for n in subimgs
            rs = zero(typeof(metric(pix).spin))
            issuccess = false
            rs, issuccess = _reactant_emission_radius_theta(pix, θs, isindir, n)
            θs_tr = zero(rs) + θs
            intersection = Intersection(zero(rs), rs, θs_tr, zero(rs), true, true)

            cond = issuccess & (horizon(metric(pix)) < rs) & (rs < T(Inf))
            contribution = material(pix, intersection)
            observation += Base.ifelse(cond, contribution, zero(contribution))
        end
    end

    return observation
end
