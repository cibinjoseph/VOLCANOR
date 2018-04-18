OpenDatabase("Results/gam*.curve database")
AddPlot("Curve", "# Gamma", 1, 1)
CurveAtts = CurveAttributes()
CurveAtts.showPoints = 1
CurveAtts.symbol = CurveAtts.Point  # Point, TriangleUp, TriangleDown, Square, Circle, Plus, X
CurveAtts.pointSize = 10
CurveAtts.curveColor = (255, 0, 0, 255)
CurveAtts.showLabels = 0
SetPlotOptions(CurveAtts)
DrawPlots()
