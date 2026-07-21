## Meshes
Meshes are the simplest objects in Krang that can be rendered to an image. 
Each mesh is made from material that defines emission mechanics, and a geometry that defines where the emission mechanics occur in the spacetime.

### Geometries
Geometries define the regions in spacetime where the emission physics originates from.
`Krang` provides two concrete geometries and an interface for defining level-set geometries.

* `ConeGeometry` : A spin axis centered cone with its apex placed at the coordinate origin.

* `MeshGeometry` : A geometry made from a [triangle mesh](https://en.wikipedia.org/wiki/Triangle_mesh). The mesh is embedded by placing vertices at points in the Cartesian Kerr-Schild coordinate system. Check out the [Raytracing a Triangular Mesh](@ref)
There are convenience functions defines to `translate`, `rotate` and `scale` these geometries.

* `AbstractLevelSetGeometry` : Subtype this interface to define a geometry whose zero set is expressed in the Cartesian Kerr-Schild coordinate system. See the [level-set example](examples/level-set-example.md).

### Materials
Materials define the local emission physics necessary to render geometries.
Materials may sometimes need additional information, which can be stored in geometry fields. `ConeGeometry`, for example, accepts attributes as its second positional argument.
