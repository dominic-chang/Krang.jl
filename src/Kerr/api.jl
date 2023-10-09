abstract type AbstractSpaceTimeEvent{T} end
function get_coordinates(event::AbstractSpaceTimeEvent{T}) where T
    throw("get_coordinates has not been defined for even of type: $(typeof(event))")
end

export AssympototicObserver, get_coordinates

abstract type AbstractSpaceTimeEvent{T} end
struct AssymptoticObserver{T} <: AbstractSpaceTimeEvent{T}
    azimuth::T
    inclination::T
end
function get_coordinates(observer::AssymptoticObserver{T}) where T
    return (azimuth=observer.azimuth, inclination=observer.inclination)
end