# ```BaseModel```

## Variables

```name```
: The name of the particular `Model` instance. If a specific name is not present then this will be a random generated string. This variable
will be set based on the `exportname` argumnet given to the constructor.

### Paths and files

`dir_current`
: This is the path where the subclass of this class is being instantiated. 

`dir_resources`
: This directory contains any additional files that is needed for the  model creation.
For example: different parts of the geometry in different *.dxf or *.svg files or B-H curves.

`dir_snapshots`
: This is the main directory where all the models are going to be exported.

`dir_media`
: If any figure is being created then you can save it under this directory.

`dir_data`
: If you have to save any data during or after the calculation, use this directory.

`dir_export`
: All the solver realted files and results will be saved under this directory.

`file_solver_script`
: This is the name of the file that the particular FEM solver will execute.

`file_solution`
: The will export the computation results into this csv file.


### Model realted

`geom`
: This is an empty Geometry instance. It's not necessary to use this variable but highly recommended. If a model has complicated geometry, then
  building it with different pieces is a preferable way. After transforming the different pieces you can just simply merge them into this variable.
  After the geometry is done you can add it to the `snapshot`.

`snapshot`
: This is the variable that the `BaseModel` class wraps. It holds the material/boundary definitions, geometry, solver setups, post-processing steps.
  You have to interact with this variable, usually at the end of the methods.

### Label Queues

Label queues are simple python lists. They hold the the label positions and label names for the different regions of the models. You can reach
these variables from any method. This is useful because some labels are easier to specify during the initialization of the model while others during the
geometry generation. To add a new label use ```.append(newlabel)``` format on any of the label queues. The `newlabel` is a 3 element tuple with the
structure: 
```py
newlabel= (label_position_x, label_position_y, label_name)

# Example:
self.label_queue.append((0.5, 1.5, "air"))
```

`label_queue`
: This list contains the labels for the different materials. 

`boundary_queue`
: This list holds the points for the boundary conditions.

`boundary_arc_queue`
: This list holds the points for the circle arc boundary conditions.

!!! note "Specifying boundary conditions"
    To specify a boundary condition you have to provide a point and a name. The closest line to the point will be selected and the boundary condition
    will be assigned to it. There are 2 types of boundary label queues because the distance calculation is different for lines and circle arcs. For lines 
    the perpendicular distance is used. For circle arcs the distance is going to be the minimum distance from the arcs characteristic points, namely the 
    start point, end point and the apex point (halfway point on the arc).


##  Class Reference

::: adze_modeler.model.BaseModel
    handler: python
    rendering:
      heading_level: 3

