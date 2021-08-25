var documenterSearchIndex = {"docs":
[{"location":"library/#Public-API","page":"Library","title":"Public API","text":"","category":"section"},{"location":"library/#Airfoil-geometries","page":"Library","title":"Airfoil geometries","text":"","category":"section"},{"location":"library/","page":"Library","title":"Library","text":"naca00\ngaw1\ncircle\nellipse","category":"page"},{"location":"library/#VSFlow.naca00","page":"Library","title":"VSFlow.naca00","text":"naca00(x; kwargs...)\n\nReturn the coordinates of a NACA00ZZ airfoil.\n\nKeyword arguments:\n\nZZ: max thickness (in % of the chord) of the airfoil.\n\n\n\n\n\n","category":"function"},{"location":"library/#VSFlow.gaw1","page":"Library","title":"VSFlow.gaw1","text":"gaw1(x)\n\nReturn the coordinates of a GAW1 airfoil.\n\n\n\n\n\n","category":"function"},{"location":"library/#VSFlow.circle","page":"Library","title":"VSFlow.circle","text":"circle(x; kwargs...)\n\nReturn the coordinates of a circle.\n\nKeyword arguments:\n\nR=.5: radius of the circle.\n\n\n\n\n\n","category":"function"},{"location":"library/#VSFlow.ellipse","page":"Library","title":"VSFlow.ellipse","text":"ellipse(x; kwargs...)\n\nReturn the coordinates of an ellipse.\n\nKeyword arguments:\n\nB=.5: half height of the ellipse.\n\n\n\n\n\n","category":"function"},{"location":"library/#Airfoil-motion","page":"Library","title":"Airfoil motion","text":"","category":"section"},{"location":"library/","page":"Library","title":"Library","text":"heavepitch\ncircularmotion","category":"page"},{"location":"library/#VSFlow.heavepitch","page":"Library","title":"VSFlow.heavepitch","text":"heavepitch(t, h0, αmax, ψ1, ψ2, strouhal)\n\nReturn the vector [x, ẋ, ẍ] corresponding to a heaving & pitching motion where x = [X, Y, α].\n\n\n\n\n\n","category":"function"},{"location":"library/#VSFlow.circularmotion","page":"Library","title":"VSFlow.circularmotion","text":"circularmotion(t, R, strouhal)\n\nReturn the vector [x, ẋ, ẍ] corresponding to a circular motion of radius R where x = [X, Y, α].\n\n\n\n\n\n","category":"function"},{"location":"library/#Running-the-simulation","page":"Library","title":"Running the simulation","text":"","category":"section"},{"location":"library/","page":"Library","title":"Library","text":"Profile\nsetinitvelocity\nprofilerun","category":"page"},{"location":"library/#VSFlow.Profile","page":"Library","title":"VSFlow.Profile","text":"PROFILE object\n\n\n\n\n\n","category":"type"},{"location":"library/#VSFlow.setinitvelocity","page":"Library","title":"VSFlow.setinitvelocity","text":"setinitvelocity(p::Profile, u, v, α̇)\n\n\n\n\n\n","category":"function"},{"location":"library/#VSFlow.profilerun","page":"Library","title":"VSFlow.profilerun","text":"\n\n\n\n","category":"function"},{"location":"library/#Post-Processing","page":"Library","title":"Post Processing","text":"","category":"section"},{"location":"library/","page":"Library","title":"Library","text":"","category":"page"},{"location":"library/#Private-API","page":"Library","title":"Private API","text":"","category":"section"},{"location":"library/#Vortex-Elements","page":"Library","title":"Vortex Elements","text":"","category":"section"},{"location":"library/","page":"Library","title":"Library","text":"","category":"page"},{"location":"library/#Airfoil-geometries-2","page":"Library","title":"Airfoil geometries","text":"","category":"section"},{"location":"library/","page":"Library","title":"Library","text":"VSFlow.initiator\nVSFlow.carlingfish\nVSFlow.splinefish\nVSFlow.swimmerthickness\nVSFlow.getbodymass\nVSFlow.getinertia\nVSFlow.getcenterofmass\nVSFlow.getcenterofmassvelocity\nVSFlow.correctvelocity","category":"page"},{"location":"library/#VSFlow.initiator","page":"Library","title":"VSFlow.initiator","text":"initiator(ξ, T)\n\nSmooth increase from 0 to 1 in half a period T and its time derivative.\n\n\n\n\n\n","category":"function"},{"location":"library/#VSFlow.carlingfish","page":"Library","title":"VSFlow.carlingfish","text":"carlingfish(s, t, T)\n\nReturn the curvature of a curve of unit length associated with the motion of a carling fish of period T. Its time derivative is computed too.\n\n\n\n\n\n","category":"function"},{"location":"library/#VSFlow.splinefish","page":"Library","title":"VSFlow.splinefish","text":"splinefish(s, t, T, X, K)\n\nReturn the curvature of a curve of unit length associated with the motion of a fish following a cubic spline.\n\n\n\n\n\n","category":"function"},{"location":"library/#VSFlow.swimmerthickness","page":"Library","title":"VSFlow.swimmerthickness","text":"swimmerthickness(s)\n\nReturn the thickness of the swimming body along its coordinate s.\n\n\n\n\n\n","category":"function"},{"location":"library/#VSFlow.getbodymass","page":"Library","title":"VSFlow.getbodymass","text":"getbodymass(XYm, XYt, XYb)\n\nIntegrate the constant and uniform unit density on the body. Integral computed with 25-point Gauss quadrature.\n\n\n\n\n\n","category":"function"},{"location":"library/#VSFlow.getinertia","page":"Library","title":"VSFlow.getinertia","text":"getinertia(XYm, XYt, XYb)\n\nReturn the moment of inertia of the body through 25-point Gauss quatradure on triangles.\n\n\n\n\n\n","category":"function"},{"location":"library/#VSFlow.getcenterofmass","page":"Library","title":"VSFlow.getcenterofmass","text":"getcenterofmass(xm, ym, xb, yb, xt, yt)\n\nCompute the center of mass of the body. Integrals computed with 25-point Gauss quadrature.\n\n\n\n\n\n","category":"function"},{"location":"library/#VSFlow.getcenterofmassvelocity","page":"Library","title":"VSFlow.getcenterofmassvelocity","text":"getcenterofmassvelocity(xm, ym, xb, yb, xt, yt, ub, vb, ut, vt)\n\nCompute the velocity of the center of mass of the body. Integrals computed with bi-linear assumption on quadrilaterals.\n\n\n\n\n\n","category":"function"},{"location":"library/#VSFlow.correctvelocity","page":"Library","title":"VSFlow.correctvelocity","text":"correctvelocity(xm, ym, xb, yb, xt, yt, ub, vb, ut, vt)\n\nCompute the angular velocity (ratio between angular impulse and scalar inertia moment). Integrals computed with bi-linear assumption on quadrilaterals, and 25-point Gauss quadrature.\n\n\n\n\n\n","category":"function"},{"location":"library/#Maths","page":"Library","title":"Maths","text":"","category":"section"},{"location":"library/","page":"Library","title":"Library","text":"VSFlow.rotate\nVSFlow.simps\nVSFlow.trapz\nVSFlow.backwarddifference\nVSFlow.centereddifference\nVSFlow.checkangle\nVSFlow.gausslegendretriangle\nVSFlow.get∫fdV\nVSFlow.cubicspline","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = VSFlow","category":"page"},{"location":"#VSFlow","page":"Home","title":"VSFlow","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for VSFlow.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"}]
}
