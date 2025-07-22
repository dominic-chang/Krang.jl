## Meshes
Meshes are the simplest objects in Krang that can be rendered to an image. 
Each mesh is made from material that defines emission mechanics, and a geometry that defines where the emission mechanics occur in the spacetime.

### Geometries
Geometries define the regions in spacetime where the emission physics originates from.
There are three basic geometries currently available in `Krang`. 

* `ConeGeometry` : A spin axis centered cone with its apex placed at the coordinate origin.

* `MeshGeometry` : A geometry made from a [triangle mesh](https://en.wikipedia.org/wiki/Triangle_mesh). The mesh is embedded by placing vertices at points in the Cartesian Kerr-Schild coordinate system. Check out a mesh example [here](examples/raytracing-mesh-example)
There are convenience functions defines to `translate`, `rotate` and `scale` these geometries.

* `LevelSetGeometry` : A geometry defined by a level set function expressed in the Cartesian Kerr-Schild coordinate system.

### Materials
Materials define the local emission physics necessary to render geometries.
Materials may sometimes need additional information that can be stored in geometries by passing `attributes` to the geometry constructor.

