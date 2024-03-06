using Krang
using GLMakie 
using Colors
Makie.inline!(true)

curr_theme = Theme(
    Axis = (
        xticksvisible = false, 
        xticklabelsvisible = false,
        xgridvisible = false,
        yticksvisible = false,
        yticklabelsvisible = false,
        ygridvisible = false,
        gridvisible = false,
        topspinevisible = false,
        bottomspinevisible = false,
        leftspinevisible = false,
        rightspinevisible = false,
        ),
    Text = (
        align = (:left,:center),
        ),
)
set_theme!(merge!(curr_theme, theme_latexfonts()))

sze = 200;
ρmax = 18e0;
# Create the material

fig = Figure(resolution = (1000, 1000));
ax = Axis(fig[1, 1:2], aspect=1)
textTag(loc, text) = begin
    ax = Axis(loc, height = 30, width=160)
    xlims!(ax,0,1)
    WGLMakie.text!(ax, 0,0; text=text)
end

ivals = Matrix{Float64}(undef, sze, sze)
θo = 85e0/180*π
θs = 40/180*π
a = 1.0

struct Annulus <: Krang.AbstractMaterial end

function (m::Annulus)(pix::Krang.AbstractPixel, geometry::Krang.ConeGeometry)
    subimgs,mul, rmin, rmax = geometry.attributes
    
    observation = false
    isindir = false
    for _ in 1:2 # Looping over isindir this way is needed to get Metal to work
        isindir ⊻= true
        for n in subimgs
            rs, _, _ = emission_radius(pix, geometry.opening_angle, isindir, n)
            if !(rmin < rs ≤ rmax) 
                continue
            end
            observation |= true
        end
    end

    return observation .* mul
end
function (m::Annulus)(pix::Krang.AbstractPixel, geometry::Krang.UnionGeometry) 
    return max(m(pix, geometry.geometry1),m(pix, geometry.geometry2))
end

function ishorizon(pix)  
    #roots = filter(Krang._isreal2, Krang.roots(pix))
    rad = Krang.emission_radius(pix, 0.01Krang.mino_time(pix, Krang.horizon(Krang.metric(pix)), false))[1]
    if !isnan(rad)
        return true
    else
        return false
    end
end

metric = Krang.Kerr(a);
camera = Krang.IntensityCamera(metric, θo, -ρmax, ρmax, -ρmax, ρmax, sze);
material = Annulus()
rmin = 18
rmax = 20

geometry = Krang.ConeGeometry(π/2, ((0, 1, 2),2, 5, 15))

mesh = Krang.Mesh(geometry, material)
ivals .= ishorizon.(camera.screen.pixels)
ivals .= max.(ivals, Array(mesh.material.(camera.screen.pixels, Ref(mesh.geometry))))

heatmap!(ax, ivals, colormap=[:black, colorant"#CB3C33",colorant"rgba(1,1,1,1)"])#, colorant"#389826",colorant"#9558B2",colorant"#CB3C33"])#], show_axis=false)
cl = append!([i for i in range(start=colorant"#389826", stop=colorant"#9558B2", length=10)], range(start=colorant"#9558B2", stop=colorant"#CB3C33", length=10)[2:end])
cl2 = append!([i for i in range(start=colorant"#389826", stop=colorant"#9558B2", length=10)], range(start=colorant"#9558B2", stop=colorant"#CB3C33", length=10)[2:end], range(start=colorant"#CB3C33", stop=colorant"#389826", length=10)[2:end-1])
text!(ax, 200/1000*sze*20/ρmax,180/1000*sze*20/ρmax,text="K.R.A.N.G", fontsize=104, color=:white)#cgrad(cl)[LinRange(0, 1, 9)])#colorant"#4063D8")
#lines!(ax, 425/1000*sze*20/ρmax .* sin.([i for i in 0:0.05:(2π+0.05)]) .+ sze/2, 425/1000*sze*20/ρmax .* cos.([i for i in 0:0.05:(2π+0.05)]) .+ sze/2, color=cgrad(cl2)[LinRange(0, 1, 127)], linewidth=10)


save("Krang_logo.png", fig)
fig
