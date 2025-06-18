@testset "Finite Observer" begin
    metric = Krang.Kerr(0.01); # Kerr metric with a spin of 0.99
    θo = 45 * π / 180; # inclination angle of the observer. θo ∈ (0, π)
    θs = π/2
    ro  = 1e8
    sze = 40; # resolution of the screen is sze x sze
    rmin = Krang.horizon(metric); # minimum radius to be ray traced
    rmax = 6.0; # maximum radius to be ray traced
    ρmax = 10; # horizontal and vertical limits of the screen

    infinite_coordinates = [zeros(sze, sze) for _ = 1:3]
    finite_coordinates = [zeros(sze, sze) for _ = 1:3]
    infinite_camera = Krang.SlowLightIntensityCamera(metric, θo, -ρmax, ρmax, -ρmax, ρmax, sze);
    finite_camera = Krang.SlowLightIntensityCamera(metric, θo, ro, -ρmax/ro, ρmax/ro, -ρmax/ro, ρmax/ro, sze);

    @test all([i.screen_coordinate[1] .* ro for i in finite_camera.screen.pixels]  .≈ [i.screen_coordinate[1]  for i in infinite_camera.screen.pixels])

    @test maximum([abs.((finite_camera.screen.pixels[i].η / infinite_camera.screen.pixels[i].η))-1 for i in range(1, length(finite_camera.screen.pixels))]) ≈ 0.0 atol = 1e-5
    @test maximum([abs.((finite_camera.screen.pixels[i].λ / infinite_camera.screen.pixels[i].λ))-1 for i in range(1, length(finite_camera.screen.pixels))]) ≈ 0.0 atol = 1e-5

    function coordinate_point(
        pix::Krang.AbstractPixel,
        geometry::Krang.ConeGeometry{T,A},
    ) where {T,A}
        n, rmin, rmax = geometry.attributes
        θs = geometry.opening_angle

        coords = ntuple(j -> zero(T), Val(4))

        isindir = false
        for _ = 1:2 # Looping over isindir this way is needed to get Metal.jl to work
            isindir ⊻= true
            ts, rs, ϕs = emission_coordinates(pix, geometry.opening_angle, isindir, n)
            if rs ≤ rmin || rs ≥ rmax
                continue
            end
            coords = isnan(rs) ? observation : (ts, rs, θs, ϕs)
        end
        return coords
    end


    function infinite_draw!(camera, coordinates, rmin, rmax, θs, n)
        times, radii, azimuths = coordinates

        geometry = Krang.ConeGeometry(θs, (n, rmin, rmax))

        rendered_scene = coordinate_point.(camera.screen.pixels, Ref(geometry))
        for I in CartesianIndices(rendered_scene)
            temp = rendered_scene[I][1]
            times[I] =  temp 
            radii[I] = rendered_scene[I][2]
            azimuths[I] = mod2pi(rendered_scene[I][4]) # azimuths are in radians
        end
        coordinates[1] .= times
        coordinates[2] .= radii
        coordinates[3] .= azimuths

    end

    function finite_draw!(camera, coordinates, rmin, rmax, θs, n)
        times, radii, azimuths = coordinates

        geometry = Krang.ConeGeometry(θs, (n, rmin, rmax))
        rendered_scene = coordinate_point.(camera.screen.pixels, Ref(geometry))
        for I in CartesianIndices(rendered_scene)
            temp = rendered_scene[I][1]
            times[I] =  temp  == 0 ? 0 : temp - 2log(ro) -ro
            radii[I] = rendered_scene[I][2]
            azimuths[I] = mod2pi(rendered_scene[I][4]) # azimuths are in radians
        end
        coordinates[1] .= times
        coordinates[2] .= radii
        coordinates[3] .= azimuths

    end

    @testset "n=$n" for n in 0:2
        finite_draw!(finite_camera, finite_coordinates, rmin, rmax, θs, n)
        infinite_draw!(infinite_camera, infinite_coordinates, rmin, rmax, θs, n)

        @testset "$coord" for (i,coord) in enumerate(["radius", "azimuth"])
            @test maximum(abs.(finite_coordinates[i+1] .- infinite_coordinates[i+1] )) .≈ 0.0 atol = 1e-4
        end
    end

end