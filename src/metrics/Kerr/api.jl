abstract type AbstractSpaceTimeEvent{T} end
function get_coordinates(event::AbstractSpaceTimeEvent{T}) where T
    throw("get_coordinates has not been defined for even of type: $(typeof(event))")
end
