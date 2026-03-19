using Reactant 
using Krang
using VLBISkyModels
to_rarray = Reactant.to_rarray

function _isreal2(num) 
	ren, imn = reim(num)
	ren2 = ren^2
	imn2 = imn^2
	return Base.:&(1, (imn2 / (imn2 + ren2))  < eps(real(num)))
end

function test(a, b)
	numreals = sum(_isreal2, b)
	ans = zero(real(eltype(b)))
	@trace if numreals == 4 
		ans = sum(Base.real, b)
		nothing
	else	
		ans = sum(Base.real, b)* a
		nothing	
	end
	return ans
end

a = Reactant.to_rarray(rand(2))
b = Reactant.to_rarray([rand(4) .+ rand(4)im, [-5.0 + 0.0im, 1.0 + 2.0im, 1.0 - 2.0im, -5.0 + 0.0im] ])

@compile test.(a,b)

a =rand(2)
b = [rand(4) .+ rand(4)im, [-5.0 + 0.0im, 1.0 + 2.0im, 1.0 - 2.0im, -5.0 + 0.0im] ]
test.(a,b)


function Base._rf_findmax(a::Tuple{Union{<:Reactant.TracedRNumber{T1}, T1}, A}, b::Tuple{Union{<:Reactant.TracedRNumber{T2}, T2}, B}) where {T1, T2, A, B}
	fm, im = a
	fx, ix = b
	Base.ifelse(fm < fx, (fx, ix), (fm, im))
end

function Base._rf_findmin(a::Tuple{Union{<:Reactant.TracedRNumber{T1}, T1}, A}, b::Tuple{Union{<:Reactant.TracedRNumber{T2}, T2}, B}) where {T1, T2, A, B}
	fm, im = a
	fx, ix = b
	Base.ifelse(fm > fx, (fx, ix), (fm, im))
end

struct JuKeBOX{T, F} <: VLBISkyModels.ComradeBase.AbstractModel
	spin::T
	θo::T
	scene::F
end

function JuKeBOX(θ::NamedTuple)
	(;
		spin,
		θo,
		θs,
		rpeak,
		p1,
		p2,
		χ,
		ι,
		βv,
		spec,
		η,
	) = θ
	T = typeof(θo)
	magfield1 = [sin(ι) * cos(η), sin(ι) * sin(η), cos(ι)]
	vel = [βv, T(π / 2), χ]

	subimgs = (1,)

	magfield1 = [sin(ι) * cos(η), sin(ι) * sin(η), cos(ι)]
	material1 = Krang.ElectronSynchrotronPowerLawIntensity(magfield1..., vel..., spec, rpeak, p1, p2, subimgs)
	geometry1 = Krang.ConeGeometry(θs*π/180, nothing)
	mesh1 = Krang.Mesh(geometry1, material1)

	magfield2 = [-sin(ι) * cos(η), -sin(ι) * sin(η), cos(ι)]
	material2 = Krang.ElectronSynchrotronPowerLawIntensity(magfield2..., vel..., spec, rpeak, p1, p2, subimgs)
	geometry2 = Krang.ConeGeometry(π-θs*π/180, nothing)
	mesh2 = Krang.Mesh(geometry2, material2)

	scene = Krang.Scene([mesh1, mesh2])

	return JuKeBOX(
		spin,
		θo,
		scene,
	)
end

@inline function VLBISkyModels.ComradeBase.intensity_point(m::JuKeBOX, p) 
	(; X, Y) = p
	(; spin, scene, θo) = m

	pix = Krang.IntensityPixel(Krang.Kerr(spin), -X, Y, θo*π/180)
	ans = Krang.render(pix, scene)
	return ans #+ one(T)
end


θo = 17 * π / 180;
ρmax = 10.0;
spin = 0.94
#camera = Krang.IntensityCamera(spin, θo, -ρmax, ρmax, -ρmax, ρmax, 100);

αval = 5.0
βval = 5.0
pixel = Krang.IntensityPixel(Kerr(spin), αval, βval, θo)
rpixel = Reactant.to_rarray(pixel)


θs = π/4
p1 = 1.5
p2 = 1.5
χ = π/2
ι = π/2
βv = 0.5
rpeak = 5.0
spec = 1.0
η = π/4

model = JuKeBOX((; spin,
	θo,
	θs,
	rpeak,
	p1,
	p2,
	χ,
	ι,
	βv,
	spec,
	η)
)
rmodel = Reactant.to_rarray(model; track_numbers=Number)

fovy = fovx = (20.0)
npix = 10
grid = imagepixels(fovx, fovy, npix, npix)
rgrid = Reactant.to_rarray(grid)

#@jit ComradeBase.allocate_imgmap(model, grid)
intensitymap(model, grid)  
@code_hlo intensitymap(rmodel, rgrid)  
@code_hlo intensitymap(rmodel, rgrid)  
rp = Reactant.to_rarray((X = 10.0, Y =10.0))
@code_hlo VLBISkyModels.ComradeBase.intensity_point(rmodel, rp)



heatmap(intmap)

function condinf0(roots)
	return 4.0* +(real,roots)
end

function condinf(result, cond, numreals, roots)
	cond2 = cond & (numreals ==4)
	@trace if cond2
		result[] = condinf0(numreals, roots)
		nothing
	end
	nothing
end

function foo(rs, roots) 
	numreals = reduce(+, ntuple(n->isreal(roots[n]), Val(4)))
	ans = 0.0
	result = Ref(ans)
	cond = isinf(rs)

	condinf(result, cond, numreals, roots)
    return result[]
end

rs = ConcreteRNumber(1.0)
roots = (ConcreteRNumber(2.0-2im), ConcreteRNumber(2.0+2im), ConcreteRNumber(3.0), ConcreteRNumber(4.0))

@jit foo(ConcreteRNumber(true), roots)