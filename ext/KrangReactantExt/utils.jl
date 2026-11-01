# Shared helpers
# Reactant specialisations mirroring the corresponding file in `src/`.

function _regularize_zero(x, ::Type{T}) where {T}
    return Base.ifelse(x == zero(T), eps(T), x)
end

@inline function _argmax_real3(x1, x2, x3)
    return Base.ifelse(
        real(x1) >= real(x2),
        Base.ifelse(real(x1) >= real(x3), x1, x3),
        Base.ifelse(real(x2) >= real(x3), x2, x3),
    )
end
