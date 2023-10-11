"""
    $TYPEDEF

Abstract Metric Type
"""
abstract type AbstractMetric end

"""
Returns the metric in some representation (usually as an nxn matrix).
"""
function metric_dd(metric::AbstractMetric, args...)
    throw(MethodError(metric_dd, (metric, args...)))
end

"""
Returns the inverse metric in some representation (usually as an nxn matrix).
"""
function metric_uu(metric::AbstractMetric, args...)
    try return inv(metric_dd(metric, args...)) catch err throw(err) end 
end

