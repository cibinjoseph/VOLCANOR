OpenDatabase("localhost:/home/cibin/WorkInProgress/VOLCANOR/Results/r01force.txt", 0)
AddPlot("Scatter", "var00", 1, 1)
ScatterAtts = ScatterAttributes()
ScatterAtts.var1 = "var00"
ScatterAtts.var2 = "var01"
ScatterAtts.pointSize = 0.05
ScatterAtts.pointSizePixels = 5
ScatterAtts.pointType = ScatterAtts.Point  # Box, Axis, Icosahedron, Octahedron, Tetrahedron, SphereGeometry, Point, Sphere
ScatterAtts.scaleCube = 0
ScatterAtts.legendFlag = 1
SetPlotOptions(ScatterAtts)

SetPlotDescription(0, "Thrust Coefficient")
AnnotationAtts = AnnotationAttributes()
AnnotationAtts.axes2D.xAxis.title.userTitle = 1
AnnotationAtts.axes2D.xAxis.title.title = "Iterations"
AnnotationAtts.axes2D.yAxis.title.userTitle = 1
AnnotationAtts.axes2D.yAxis.title.title = "Thrust Coefficient"
SetAnnotationAttributes(AnnotationAtts)
DrawPlots()
