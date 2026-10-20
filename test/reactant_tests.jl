push!(LOAD_PATH, dirname(@__DIR__))
using Krang
using Test
using Reactant
import Reactant: to_rarray


struct ReactantRaytraceMaterial{N} <: Krang.AbstractMaterial
    subimgs::NTuple{N,Int}
end
function (mat::ReactantRaytraceMaterial)(pix, intersection)
    rs = intersection.rs
    return rs / (one(rs) + rs)
end

const MET = Krang.Kerr(0.94)
const ΘO = 45 / 180 * π
const ΘS = 80 / 180 * π

"Screen points grouped by how many of the four radial roots are real."
function screen_grid(met, θo; npts = 3)
    pts = Dict(
        4 => Tuple{Float64,Float64}[],
        2 => Tuple{Float64,Float64}[],
        0 => Tuple{Float64,Float64}[],
    )
    for α in range(-9.0, 9.0, 40), β in range(-9.0, 9.0, 40)
        η = Krang.η(met, α, β, θo)
        λ = Krang.λ(met, α, θo)
        n = sum(Krang._isreal2, Krang.get_radial_roots(met, η, λ))
        length(pts[n]) < npts && push!(pts[n], (α, β))
    end
    return pts
end

const GRID = screen_grid(MET, ΘO)

# Shared helpers. `rootsof` recomputes the roots inside each scalar closure instead of
# broadcasting a tuple-returning function, because a broadcast whose element type is a
# `Tuple` cannot be materialised by Reactant.
ηof(m, α, β) = Krang.η(m, α, β, ΘO)
λof(m, α) = Krang.λ(m, α, ΘO)
rootsof(m, α, β) = Krang.get_radial_roots(m, ηof(m, α, β), λof(m, α))
pixof(m, α, β) = Krang.SlowLightIntensityPixel(m, α, β, ΘO)

"A radius between the horizon and the outermost root, so the source lies on the ray."
function rsof(m, α, β)
    rh = Krang.horizon(m)
    rmax = maximum(real, rootsof(m, α, β))
    return rh + (max(rmax, rh + 2 * one(rh)) - rh) / 2
end

"Half the ray's own total Mino time, so τ names a point the geodesic reaches."
τof(m, α, β) = Krang.total_mino_time(m, rootsof(m, α, β)) / 2

"""
Broadcast slot `i` of `f` over the screen points, for each `i`, as a tuple of arrays.
`ntuple(..., Val(N))` keeps each index a compile-time constant, which matters because the
return tuples are heterogeneous -- `emission_radius` returns `(Float64, Bool, Bool, Int64,
Bool)` -- and a runtime index would make the closure uninferrable.
"""
slots(f, ::Val{N}, m, as, bs) where {N} =
    ntuple(i -> ((mm, a, b) -> f(mm, a, b)[i]).(Ref(m), as, bs), Val(N))

"Single-value variant for the functions that return a bare scalar."
only_value(f, m, as, bs) = (((mm, a, b) -> f(mm, a, b)).(Ref(m), as, bs),)

"""
    compare(name, labels, make, m, cpu_pts, traced_pts)

`make(m, as, bs)` returns a tuple of arrays. Run it on CPU inputs and on traced inputs and
compare slot by slot. A throw during tracing is reported as a failure carrying the
exception message rather than escaping, so one broken overlay does not suppress the
verdict on every overlay after it.
"""
function compare(name, labels, make, m, cpu_pts, traced_pts; atol = 1e-8, rtol = 1e-6)
    @testset "$name" begin
        base = make(m, cpu_pts...)
        local traced
        # Each bundle is a separate XLA compilation taking tens of seconds, and a bundle
        # that never finishes compiling is otherwise indistinguishable from a slow one, so
        # announce which one is being traced before entering it.
        t0 = time()
        print(stderr, "  tracing $name ... ")
        flush(stderr)
        try
            traced = @jit make(m, traced_pts...)
            println(stderr, round(time() - t0, digits = 1), "s")
        catch err
            println(stderr, "threw after ", round(time() - t0, digits = 1), "s")
            @test sprint(showerror, err) == "<traced without error>"
            return
        end
        for (i, label) in enumerate(labels)
            got = Array(traced[i])
            want = base[i]
            @testset "$label" begin
                # An overlay that lets an integer count collapse to `Bool` is wrong even
                # when the values happen to compare equal, so the element type is checked.
                @test eltype(got) === eltype(want)
                if eltype(want) <: Union{Integer,Bool}
                    @test got == want
                else
                    @test all(isapprox.(got, want; atol = atol, rtol = rtol))
                end
            end
        end
    end
end

@testset "Reactant overlays" begin

    @testset "input coverage" begin
        # Without this the rest of the file silently tests less than it claims to.
        @test length(GRID[4]) > 0
        @test length(GRID[2]) > 0
        @test length(GRID[0]) > 0
    end

    for nreal in (4, 2, 0)
        isempty(GRID[nreal]) && continue
        αs = first.(GRID[nreal])
        βs = last.(GRID[nreal])
        cpu_pts = (αs, βs)
        traced_pts = (to_rarray(αs), to_rarray(βs))
        cmp(name, labels, make) = compare(name, labels, make, MET, cpu_pts, traced_pts)

        @testset "$nreal real radial roots" begin

            # --- conserved quantities and the roots themselves ---
            cmp(
                "η and λ",
                ("η", "λ"),
                (m, as, bs) ->
                    (ηof.(Ref(m), as, bs), ((mm, a, b) -> λof(mm, a)).(Ref(m), as, bs)),
            )

            # Only the real parts are compared: the sqrt branch cut can return a conjugate
            # pair in either order, which is a labelling difference, not a wrong root.
            cmp(
                "get_radial_roots",
                ("r1", "r2", "r3", "r4", "numreals"),
                (m, as, bs) -> (
                    ntuple(
                        i ->
                            ((mm, a, b) -> real(rootsof(mm, a, b)[i])).(Ref(m), as, bs),
                        Val(4),
                    )...,
                    # The `ifelse` is the form the overlays themselves use, and it is
                    # load-bearing rather than stylistic: `sum(_isreal2, roots)` reduces
                    # with `Bool + Bool`, which Reactant does not promote to `Int` (see
                    # the regression guard at the bottom of this file). Counting the way
                    # the overlays count is what makes this a test of Krang.
                    (
                        (mm, a, b) -> sum(
                            x -> Base.ifelse(Krang._isreal2(x), 1, 0),
                            rootsof(mm, a, b),
                        )
                    ).(
                        Ref(m),
                        as,
                        bs,
                    ),
                ),
            )

            # --- radial integrals evaluated at infinity ---
            cmp(
                "Ir_inf",
                ("Ir_inf",),
                (m, as, bs) -> only_value(
                    (mm, a, b) -> Krang.Ir_inf(mm, rootsof(mm, a, b)),
                    m,
                    as,
                    bs,
                ),
            )
            cmp(
                "total_mino_time",
                ("τ_total",),
                (m, as, bs) -> only_value(
                    (mm, a, b) -> Krang.total_mino_time(mm, rootsof(mm, a, b)),
                    m,
                    as,
                    bs,
                ),
            )
            # `Iϕ_inf` and `It_inf` are the two kernels that lean hardest on the elliptic
            # Π integral, and tracing either of them produces an XLA module that was still
            # compiling after twenty-five minutes while every other bundle here finishes in
            # under twenty seconds. They are opt-in rather than deleted, so the gap stays
            # visible: set KRANG_REACTANT_SLOW_TESTS=true to attempt them.
            if get(ENV, "KRANG_REACTANT_SLOW_TESTS", "false") == "true"
                cmp(
                    "Iϕ_inf",
                    ("Iϕ_inf",),
                    (m, as, bs) -> only_value(
                        (mm, a, b) -> Krang.Iϕ_inf(mm, rootsof(mm, a, b), λof(mm, a)),
                        m,
                        as,
                        bs,
                    ),
                )
                cmp(
                    "It_inf",
                    ("It_inf",),
                    (m, as, bs) -> only_value(
                        (mm, a, b) -> Krang.It_inf(mm, rootsof(mm, a, b), λof(mm, a)),
                        m,
                        as,
                        bs,
                    ),
                )
            end
            cmp(
                "radial_inf_integrals",
                ("I1", "I2", "Ip", "Im"),
                (m, as, bs) -> slots(
                    (mm, a, b) -> Krang.radial_inf_integrals(mm, rootsof(mm, a, b)),
                    Val(4),
                    m,
                    as,
                    bs,
                ),
            )

            # --- angular antiderivative pairs cached on every pixel ---
            cmp(
                "_absGθo_Gθhat",
                ("Gθo", "Gθhat"),
                (m, as, bs) -> slots(
                    (mm, a, b) ->
                        Krang._absGθo_Gθhat(mm, ΘO, ηof(mm, a, b), λof(mm, a)),
                    Val(2),
                    m,
                    as,
                    bs,
                ),
            )
            cmp(
                "_absGϕo_Gϕhat",
                ("Gϕo", "Gϕhat"),
                (m, as, bs) -> slots(
                    (mm, a, b) ->
                        Krang._absGϕo_Gϕhat(mm, ΘO, ηof(mm, a, b), λof(mm, a)),
                    Val(2),
                    m,
                    as,
                    bs,
                ),
            )
            cmp(
                "_absGto_Gthat",
                ("Gto", "Gthat"),
                (m, as, bs) -> slots(
                    (mm, a, b) ->
                        Krang._absGto_Gthat(mm, ΘO, ηof(mm, a, b), λof(mm, a)),
                    Val(2),
                    m,
                    as,
                    bs,
                ),
            )

            # --- integrals at a finite emission radius ---
            cmp(
                "Ir_s",
                ("Ir_s",),
                (m, as, bs) -> only_value(
                    (mm, a, b) ->
                        Krang.Ir_s(mm, rsof(mm, a, b), rootsof(mm, a, b), true),
                    m,
                    as,
                    bs,
                ),
            )
            # Same story as `Iϕ_inf`: these two carry the finite-radius elliptic integrals
            # and do not finish compiling. Opt in with KRANG_REACTANT_SLOW_TESTS=true.
            if get(ENV, "KRANG_REACTANT_SLOW_TESTS", "false") == "true"
                cmp(
                    "radial_w_I0_terms_integrals",
                    ("I1", "I2", "Ip", "Im"),
                    (m, as, bs) -> slots(
                        (mm, a, b) -> Krang.radial_w_I0_terms_integrals(
                            mm,
                            rsof(mm, a, b),
                            rootsof(mm, a, b),
                            τof(mm, a, b),
                            true,
                        ),
                        Val(4),
                        m,
                        as,
                        bs,
                    ),
                )
                cmp(
                    "radial_integrals",
                    ("I0", "I1", "I2", "Ip", "Im"),
                    (m, as, bs) -> slots(
                        (mm, a, b) -> Krang.radial_integrals(
                            pixof(mm, a, b),
                            rsof(mm, a, b),
                            τof(mm, a, b),
                            true,
                        ),
                        Val(5),
                        m,
                        as,
                        bs,
                    ),
                )

            end

            # --- pixel-level overlays ---
            # Each extra slot in a bundle is another full copy of the traced program, and
            # a `SlowLightIntensityPixel` traces its whole cache (roots, I0_inf, Mino time
            # and all three angular pairs). Bundling six slots of `Gt` that way produced an
            # XLA module that had not finished compiling after ten minutes, so the
            # pixel-level overlays are spot-checked a slot or two at a time rather than
            # swept exhaustively. `emission_radius` is the one driven across all three root
            # cases, since it is what selects `_rs_case1_and_2` / `_rs_case3` / `_rs_case4`.
            cmp(
                "emission_radius rs",
                ("rs",),
                (m, as, bs) -> only_value(
                    (mm, a, b) ->
                        Krang.emission_radius(pixof(mm, a, b), ΘS, false, 0)[1],
                    m,
                    as,
                    bs,
                ),
            )
            cmp(
                "emission_radius issuccess",
                ("issuccess",),
                (m, as, bs) -> only_value(
                    (mm, a, b) ->
                        Krang.emission_radius(pixof(mm, a, b), ΘS, false, 0)[5],
                    m,
                    as,
                    bs,
                ),
            )
        end
    end

    # The remaining pixel-level overlays are checked once, on four real roots, for the
    # cost reason above.
    # Tracing a `SlowLightIntensityPixel` compiles its entire cache -- roots, I0_inf,
    # Mino time and all three angular pairs -- and these bundles produced XLA modules that
    # had not finished compiling after ten minutes. They are opt-in so the rest of the file
    # stays runnable; set KRANG_REACTANT_PIXEL_TESTS=true to include them.
    @testset "pixel-level overlays (4 real roots)" begin
        # `return` must not be used to skip here: a `@testset` body is not a function, so a
        # `return` propagates out of the *enclosing* testset and silently drops every test
        # after it (it previously swallowed the two testsets below).
        if get(ENV, "KRANG_REACTANT_PIXEL_TESTS", "false") != "true"
            @test_skip "set KRANG_REACTANT_PIXEL_TESTS=true to run"
        else
            αs = first.(GRID[4])
            βs = last.(GRID[4])
            cpu_pts = (αs, βs)
            traced_pts = (to_rarray(αs), to_rarray(βs))
            cmp(name, labels, make) = compare(name, labels, make, MET, cpu_pts, traced_pts)

            cmp(
                "SlowLightIntensityPixel I0_inf",
                ("I0_inf",),
                (m, as, bs) ->
                    only_value((mm, a, b) -> Krang.I0_inf(pixof(mm, a, b)), m, as, bs),
            )
            cmp(
                "emission_coordinates ts",
                ("ts",),
                (m, as, bs) -> only_value(
                    (mm, a, b) ->
                        Krang.emission_coordinates(pixof(mm, a, b), ΘS, false, 0)[1],
                    m,
                    as,
                    bs,
                ),
            )
            cmp(
                "emission_inclination θs",
                ("θs",),
                (m, as, bs) -> only_value(
                    (mm, a, b) -> Krang.emission_inclination(
                        pixof(mm, a, b),
                        rsof(mm, a, b),
                        true,
                    )[1],
                    m,
                    as,
                    bs,
                ),
            )
        end
    end

    # Not a Krang test: this pins the upstream Reactant behaviour the `ifelse` above
    # exists to dodge, so that if Reactant ever promotes `Bool + Bool` to `Int` again these
    # flip to "Unexpectedly Pass" and the workaround can be dropped.
    #
    #   julia> f(x) = (x > 0) + (x > 0)
    #   julia> f(1.0)                                 # 2    :: Int64
    #   julia> Array(@jit f.(to_rarray([1.0])))[1]    # true :: Bool
    #
    # Base promotes `Bool + Bool` to `Int`; Reactant keeps it `Bool`, so every reduction
    # built on `+` answers "did any match" instead of "how many matched". Verified against
    # Reactant v0.2.295 (#main).
    @testset "Reactant Bool-sum regression guard" begin
        xs = [1.0, -2.0, 3.0]
        xr = to_rarray(xs)
        pos(x) = x > 0
        four(x) = (x, -x, 2x, -2x) # two of the four are positive

        traced_eltype(f) = eltype(Array(@jit f.(xr)))

        @test eltype((x -> pos(x) + pos(x)).(xs)) === Int64       # Base promotes
        @test_broken traced_eltype(x -> pos(x) + pos(x)) === Int64
        @test_broken traced_eltype(x -> sum(pos, four(x))) === Int64
        @test_broken traced_eltype(x -> reduce(+, pos.(four(x)))) === Int64

        # Mixed-type promotion is fine, so this is specific to `Bool + Bool`.
        @test traced_eltype(x -> pos(x) + 1) === Int64

        # The workaround the overlays use: the summands are already `Int`, so `+` is never
        # asked to promote.
        @test traced_eltype(x -> sum(y -> Base.ifelse(pos(y), 1, 0), four(x))) === Int64

        # Two neighbouring gaps: `count` throws rather than collapsing, and the obvious
        # manual conversion has no method.
        @test_throws Exception (@jit (x -> count(pos, four(x))).(xr))
        @test_throws Exception (@jit (x -> Int(pos(x))).(xr))
    end

    # --- the cone render path, which composes every overlay above ---
    @testset "cone _raytrace" begin
        αs = [-6.0, -2.0, 3.0, 7.0]
        βs = [-5.0, 1.0, 4.0, 6.0]
        mesh = Krang.Mesh(Krang.ConeGeometry(ΘS), ReactantRaytraceMaterial((0, 1)))
        compare(
            "raytrace cone",
            ("observation",),
            (m, as, bs) ->
                only_value((mm, a, b) -> Krang.raytrace(pixof(mm, a, b), mesh), m, as, bs),
            MET,
            (αs, βs),
            (to_rarray(αs), to_rarray(βs)),
        )
    end
end
