using Krang
using CairoMakie 
using Colors

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

sze = 4000;
ρmax = 17e0;


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
a = 1-1e-8

struct Annulus <: Krang.AbstractMaterial
end

function (m::Annulus)(pix::Krang.AbstractPixel, geometry::Krang.ConeGeometry)
    subimgs,mul, rmin, rmax = geometry.attributes
    
    observation = false
    isindir = false
    for _ in 1:2 # Looping over isindir this way is needed to get Metal to work
        isindir ⊻= true
        for n in subimgs
            
            rs, _, ϕs,_ = emission_coordinates_fast_light(pix, geometry.opening_angle, isindir, n)
            #if (n == 1 && isindir) || (n == 1 && (1.5π/2 > ϕs >0 ||  2π>ϕs > 1.2π) )
            #    continue
            #end
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
camera = Krang.SlowLightIntensityCamera(metric, θo, -ρmax, ρmax, -ρmax, ρmax, sze);
material = Annulus()
rmin = 18
rmax = 20

cl_1 = colorant"#586BA6"
cl_2 = colorant"#4590BFBF"
#cl_2 = colorant"#4590BF"
cl_3 = colorant"#F2CB05"
cl_4 = colorant"#F27F3D"
cl_5 = colorant"#D95A4E"
r_cl = colorant"#CB3C33"
p_cl = colorant"#9558B2"
g_cl = colorant"#389826"
cl = append!([i for i in range(start=g_cl, stop=p_cl, length=10)], range(start=p_cl, stop=r_cl, length=10))
cl2 = append!([i for i in range(start=r_cl, stop=p_cl, length=10)], range(start=p_cl, stop=g_cl, length=10)[2:end-1], range(start=g_cl, stop=r_cl, length=10))

geometry = Krang.ConeGeometry(π/2, ((0,),2, 2, 15))

mesh = Krang.Mesh(geometry, material)
horizon_vals = ishorizon.(camera.screen.pixels)
heatmap!(ax, horizon_vals,colormap=[:white, cl_5])
ivals = map(x->x[1] >0 &&  x[2]>0 ? 1.0 : 0.0, zip(horizon_vals, Array(mesh.material.(camera.screen.pixels, Ref(mesh.geometry)))))
#heatmap!(ax, ivals, colormap=[colorant"rgba(256,256,256,0)",  colorant"rgba(0,0,0,0.25)"])#, colorant"#389826",colorant"#9558B2",colorant"#CB3C33"])#], show_axis=false)
fig

##for i in (2,1,0)
#    #geometry5 = Krang.ConeGeometry(π/2, ((i), 5, 15))
#    geometry5 = Krang.ConeGeometry(π/2, ((0, 1, ),2, 5, 15))
#
#
#    #mesh = Krang.Mesh(geometry5, Krang.CoordinatePoint())
#    mesh = Krang.Mesh(geometry5, Annulus())
#    #ivals .= map(x->begin t=x[4];iszero(t) ? NaN : t end,Array(mesh.material.(camera.screen.pixels, Ref(mesh.geometry))))
#    ivals .= mesh.material.(camera.screen.pixels, Ref(mesh.geometry))
#    heatmap!(ax, ivals, colormap=[colorant"rgba(0,0,0,0)",cl_2])#[colorant"rgba(256,245,256,0)",:gray])#, colorant"#389826",colorant"#9558B2",colorant"#CB3C33"])#], show_axis=false)
##end

g_cr = begin
    t = cgrad(:grays)[LinRange(0, 1, 9)]
    append!(t,reverse(t)[begin+1:end])
end


for i in (1, 0)
    geometry5 = Krang.ConeGeometry(π/2, ((i), 6, 15))

    mesh = Krang.Mesh(geometry5, Krang.CoordinatePoint())
    #mesh = Krang.Mesh(geometry5, Annulus())
    ivals .= map(x->begin t=x[4];iszero(t) ? NaN : t end,Array(mesh.material.(camera.screen.pixels, Ref(mesh.geometry))))
    #ivals = mesh.material.(camera.screen.pixels, Ref(mesh.geometry))
    heatmap!(ax, ivals, colormap=g_cr)#[colorant"rgba(256,245,256,0)",:gray])#, colorant"#389826",colorant"#9558B2",colorant"#CB3C33"])#], show_axis=false)
end

#for i in (2,1,0)
    #geometry5 = Krang.ConeGeometry(π/2, ((i), 5, 15))
    geometry5 = Krang.ConeGeometry(π/2, ((0, 1, ),2, 6, 15))
    #mesh = Krang.Mesh(geometry5, Krang.CoordinatePoint())
    mesh = Krang.Mesh(geometry5, Annulus())
    #ivals .= map(x->begin t=x[4];iszero(t) ? NaN : t end,Array(mesh.material.(camera.screen.pixels, Ref(mesh.geometry))))
    ivals .= mesh.material.(camera.screen.pixels, Ref(mesh.geometry))
    heatmap!(ax, ivals, colormap=[colorant"rgba(0,0,0,0)",cl_2])#[colorant"rgba(256,245,256,0)",:gray])#, colorant"#389826",colorant"#9558B2",colorant"#CB3C33"])#], show_axis=false)
#end




fig

#text!(ax, 220/1000*sze*20/ρmax,180/1000*sze*20/ρmax,text="K.R.A.N.G", fontsize=104, color=cgrad(cl)[LinRange(0, 1, 9)])#colorant"#4063D8")


save("Krang_logo.png", fig)
fig
