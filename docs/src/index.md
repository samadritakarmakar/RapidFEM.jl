# RapidFEM.jl Documentation

## Basic Usage

First we shall discuss basic usage of the library. By basic usage we mean Types, variables and methods that may be used as a end-user. One of the best methods to uderstand the library is to look at the examples in the directory "src/Examples". Start with "poisson.jl" then go to "LinearElastic.jl" and keep on moving ahead from then on.  

#### Mesh Reading

```@docs
    RapidFEM.Mesh
```

```@docs
    RapidFEM.readMesh
```
#### Initializing a Finite Element Space

```@docs
    RapidFEM.createFeSpace
```

#### Assembling the Finite Element Matrix

```@docs
    RapidFEM.assembleMatrix
```
#### Assembling the Finite Element Vector

```@docs
    RapidFEM.assembleVector
```
#### Calculating the Parameters for the model

```@docs
    RapidFEM.assembleScalar
```

#### Applying Dirichlet Boundary Condition

```@docs
    RapidFEM.RapidFEM.applyDirichletBC!
```
#### Post-Processing Solution Data

```@docs
    RapidFEM.InvDistInterpolation
```

```@docs
    RapidFEM.voigtToTensor
```
## Advanced Usage
In this part we try to delve into more details about the library so that the advanced user may find it useful to add new mathematical models.  

#### Mesh Related Data

```@docs
    RapidFEM.getNoOfElements
```

```@docs
    RapidFEM.getCoordArray
```
#### Types of Elements

```@docs
    RapidFEM.AbstractElement
```

```@docs
    RapidFEM.LineElement
```

```@docs
    RapidFEM.TriElement
```

```@docs
    RapidFEM.QuadElement
```

```@docs
    RapidFEM.TetElement
```
```@docs
    RapidFEM.HexElement
```
