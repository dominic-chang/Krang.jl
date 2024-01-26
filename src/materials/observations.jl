function observe(geometry::AbstractGeometry, material::AbstractMaterial, observer::AbstractCamera) 
   error("Observations by type $(typeof(observer)) on type $(typeof(material)) is not defined for geometry type $(typeof(geometry)).") 
end

#function observe(geometry::ConeGeometry, material::PowerLawPolarization{T}, observer::AbstractCamera; α=one(T), subimgs=[0,], profile=x->one(T)) where T
#    observation = Array{Union{Missing, Vector{T}}}(missing, size(observer.screen.pixels))
#    θs = geometry.opening_angle
#    θo = inclination(observer.screen.pixels[1])
#    # is θ̇s increasing or decreasing?
#
#    for (I, pix) in enumerate(observer.screen.pixels)
#        observation[I] = zeros(4)
#
#        for n in subimgs
#            for isindir in [true, false]
#                νθ = cos(θs) < abs(cos(θo)) ? (θo > θs) ⊻ (n % 2 == 1) : !isindir
#                rs, νr, _ =  emission_radius(pix, geometry.opening_angle, isindir, n)
#                eα, eβ, redshift, lp = material(pix, rs, θs, νr, νθ) 
#            
#                prof = profile(rs)*max(redshift , eps(T))^(3+α)
#                q = T(-(eα^2 - eβ^2)*lp*prof + eps(T))
#                u = T(-2*eα*eβ*lp*prof + eps(T))
#                i = profile(rs)#T(hypot(q, u))
#                nan2zero = x -> isnan(x) ? zero(T) : x
#                observation[I] += [nan2zero(i), nan2zero(q), nan2zero(u), zero(T)]
#            end
#        end
#    end
#    return observation
#end
#