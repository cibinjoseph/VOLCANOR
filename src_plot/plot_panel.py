#Plot wing PC
OpenDatabase("Results/wingPC.tec")
AddPlot("Mesh", "mesh", 1, 1)

#Plot wing PC
OpenDatabase("Results/wingVR.tec")
AddPlot("Mesh", "mesh", 1, 1)
MeshAtts = MeshAttributes()
MeshAtts.meshColor = (255, 0, 0, 255)  # Red color
MeshAtts.meshColorSource = MeshAtts.MeshCustom  # Foreground, MeshCustom
MeshAtts.lineStyle = MeshAtts.DOTDASH  # SOLID, DASH, DOT, DOTDASH
SetPlotOptions(MeshAtts)

#Plot wing PC
OpenDatabase("Results/wingCP.tec")
ScatterAtts = ScatterAttributes()
ScatterAtts.var1Role = ScatterAtts.Coordinate0  # Coordinate0, Coordinate1, Coordinate2, Color, None
ScatterAtts.var1 = "X"
ScatterAtts.var2Role = ScatterAtts.Coordinate1  # Coordinate0, Coordinate1, Coordinate2, Color, None
ScatterAtts.var2 = "Y"
ScatterAtts.var3Role = ScatterAtts.Coordinate2  # Coordinate0, Coordinate1, Coordinate2, Color, None
ScatterAtts.var3 = "Z"
ScatterAtts.pointSizePixels = 8
ScatterAtts.pointType = ScatterAtts.Sphere  # Box, Axis, Icosahedron, Octahedron, Tetrahedron, SphereGeometry, Point, Sphere
ScatterAtts.scaleCube = 0
ScatterAtts.colorType = ScatterAtts.ColorBySingleColor  # ColorByForegroundColor, ColorBySingleColor, ColorByColorTable
ScatterAtts.singleColor = (0, 0, 255, 255)  # Blue color
SetDefaultPlotOptions(ScatterAtts)
AddPlot("Scatter", "X", 1, 1)

DrawPlots()
