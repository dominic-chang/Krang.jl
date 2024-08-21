export render

function render!(::RayTrace, camera::AbstractCamera, scene::Scene)
    mapreduce(mesh -> mesh.material.(camera.screen.pixels, Ref(mesh.geometry)), +, scene)
end
