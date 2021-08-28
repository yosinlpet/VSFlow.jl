var documenterSearchIndex = {"docs":
[{"location":"library/","page":"Library","title":"Library","text":"Pages = [\"library.md\"]\nDepth = 2","category":"page"},{"location":"library/#Public-API","page":"Library","title":"Public API","text":"","category":"section"},{"location":"library/#Airfoil-geometries","page":"Library","title":"Airfoil geometries","text":"","category":"section"},{"location":"library/","page":"Library","title":"Library","text":"circle\nellipse\ngaw1\nnaca00","category":"page"},{"location":"library/#VSFlow.circle","page":"Library","title":"VSFlow.circle","text":"circle(R)(x)\n\nReturn the y coordinates of a circle of radius R.\n\n\n\n\n\n","category":"function"},{"location":"library/#VSFlow.ellipse","page":"Library","title":"VSFlow.ellipse","text":"ellipse(B)(x)\n\nReturn the y coordinates of an ellipse with half height B.\n\n\n\n\n\n","category":"function"},{"location":"library/#VSFlow.gaw1","page":"Library","title":"VSFlow.gaw1","text":"gaw1(x)\n\nReturn the y coordinates of a GAW1 airfoil.\n\n\n\n\n\n","category":"function"},{"location":"library/#VSFlow.naca00","page":"Library","title":"VSFlow.naca00","text":"naca00(ZZ)(x)\n\nReturn the y coordinates of a NACA00ZZ airfoil.\n\n\n\n\n\n","category":"function"},{"location":"library/#Airfoil-motion","page":"Library","title":"Airfoil motion","text":"","category":"section"},{"location":"library/","page":"Library","title":"Library","text":"uniform\nheavepitch\ncircularmotion","category":"page"},{"location":"library/#VSFlow.uniform","page":"Library","title":"VSFlow.uniform","text":"uniform(vx, vy, α̇)(t)\n\nReturn the vector [x, ẋ, ẍ] corresponding to a uniform straight motion in the [-X, Y, -θ] direction  with velocity [vx, vy, α̇].\n\n\n\n\n\n","category":"function"},{"location":"library/#VSFlow.heavepitch","page":"Library","title":"VSFlow.heavepitch","text":"heavepitch(h0, αmax, ψ1, ψ2, strouhal)(t)\n\nReturn the vector [x, ẋ, ẍ] corresponding to a heaving & pitching motion where x = [-X, Y, α].\n\n\n\n\n\n","category":"function"},{"location":"library/#VSFlow.circularmotion","page":"Library","title":"VSFlow.circularmotion","text":"circularmotion(R, strouhal)(t)\n\nReturn the vector [x, ẋ, ẍ] corresponding to a circular motion of radius R where x = [X, Y, α].\n\n\n\n\n\n","category":"function"},{"location":"library/#Running-the-simulation","page":"Library","title":"Running the simulation","text":"","category":"section"},{"location":"library/","page":"Library","title":"Library","text":"Profile(;id, profile::Function, N, position, dt, T, δ, ϵ, lumpargs)\nprofilerun","category":"page"},{"location":"library/#VSFlow.Profile-Tuple{}","page":"Library","title":"VSFlow.Profile","text":"Profile(; id, profileshape::Function, x0, ẋ0, N, dt, T, δ = 1e-2, ϵ = 1e-2, (η, Tmin, Smin) = zeros(3))\n\nConstructs a Profile object. If eta, Tmin, Smin are all 0, then no lumping operation occurs.\n\nKeyword Arguments\n\nid: String used to identify the profile.\nprofileshape: function determining the profile shape.\nN: number of panels on the surface of the profile.\nx0: initial position [X0, Y0, α0] of the quarter chord of the profile.\nẋ0: initial velocity [Ẋ0, Ẏ0, α̇0] of the quarter chord of the profile.\ndt: timestep of the simulation.\nT: horizon time.\nδ = 1e-2: kernel cut-off width.\nϵ = 1e-2: percentage of the average panel length to chop the trailing edge.\nη = 0: maximum error due to lumping.\nTmin = 0: minimum time between two successive active vortices.\nSmin = 0: minimum vortex sheet length in the wake.\n\n\n\n\n\n","category":"method"},{"location":"library/#VSFlow.profilerun","page":"Library","title":"VSFlow.profilerun","text":"profilerun(p::Profile, accfunc, isshow=false; is4thorder=true)\n\nSimulate the testcase.\n\nArguments\n\np: Profile simulation to run.\naccfunc: Function giving the force applied on the pivot point of the body.\nisshow = false: Bool indicating whether an animation has to be generated.\n\nKeyword Arguments\n\nis4thorder = true: Bool indicating if the time marching is RK4 or ForwardEuler.\n\n\n\n\n\n","category":"function"},{"location":"library/#Post-processing","page":"Library","title":"Post processing","text":"","category":"section"},{"location":"library/","page":"Library","title":"Library","text":"History\ngetimpulse(p::Profile)\ngetboundpotential","category":"page"},{"location":"library/#VSFlow.History","page":"Library","title":"VSFlow.History","text":"History\n\nA structure allowing to record useful data.\n\nFields\n\nt: timelapse.\nX: position [X, Y, α].\nẊ: velocity [U, V, ω].\nacp: aerodynamic coefficients obtained with pressure integration.\nacn: aerodynamic coefficients obtained with control volume approach.\naci: aerodynamic coefficients obtained with momentum conservation.\nϕs: velocity potential on body surface.\ncps: pressure coefficients on body surface.\nP: impulse [px, py, pm].\nΓ: total body circulation.\nθ: shedding angle (0 is the bissector).\n\n\n\n\n\n","category":"type"},{"location":"library/#VSFlow.getimpulse-Tuple{Profile}","page":"Library","title":"VSFlow.getimpulse","text":"getimpulse(p::Profile)\n\nReturn the linear/angular impulse exerted by the fluid on the body.\n\n\n\n\n\n","category":"method"},{"location":"library/#VSFlow.getboundpotential","page":"Library","title":"VSFlow.getboundpotential","text":"getboundpotential(p::Profile)\n\nIntegrates the velocity on the profile to obtain the velocity potential ϕ.\n\n\n\n\n\n","category":"function"},{"location":"library/#Private-API","page":"Library","title":"Private API","text":"","category":"section"},{"location":"library/#Vortex-Elements","page":"Library","title":"Vortex Elements","text":"","category":"section"},{"location":"library/","page":"Library","title":"Library","text":"VSFlow.VortexPoint\nVSFlow.ConstantPanel\nVSFlow.LinearPanel","category":"page"},{"location":"library/#VSFlow.VortexPoint","page":"Library","title":"VSFlow.VortexPoint","text":"VortexPoint\n\nStructure representing a point vortex.\n\nFields\n\nX: X-position of the point vortex in the inertial frame.\nY: Y-position of the point vortex in the inertial frame.\nΓ: Circulation of the point vortex.\ndv: Velocity correction in case of lumping.\n\n\n\n\n\n","category":"type"},{"location":"library/#VSFlow.ConstantPanel","page":"Library","title":"VSFlow.ConstantPanel","text":"ConstantPanel\n\nStructure representing a vortex panel of uniform strength.\n\nFields\n\nX1: X-position of the start of the panel in the inertial frame.\nY1: Y-position of the start of the panel in the inertial frame.\nX2: X-position of the end of the panel in the inertial frame.\nY2: Y-position of the end of the panel in the inertial frame.\nXc: X-position of the center of the panel in the inertial frame.\nYc: Y-position of the center of the panel in the inertial frame.\nb: half-length of the panel.\nθ: panel orientation with respect to the X axis of the inertial frame.\nprofID: ID of the profile the panel belongs to.\n\n\n\n\n\n","category":"type"},{"location":"library/#VSFlow.LinearPanel","page":"Library","title":"VSFlow.LinearPanel","text":"LinearPanel\n\nStructure representing a vortex panel of linearly varying strength.\n\nFields\n\nX1: X-position of the start of the panel in the inertial frame.\nY1: Y-position of the start of the panel in the inertial frame.\nX2: X-position of the end of the panel in the inertial frame.\nY2: Y-position of the end of the panel in the inertial frame.\nXc: X-position of the center of the panel in the inertial frame.\nYc: Y-position of the center of the panel in the inertial frame.\nb: half-length of the panel.\nθ: panel orientation with respect to the X axis of the inertial frame.\nprofID: ID of the profile the panel belongs to.\n\n\n\n\n\n","category":"type"},{"location":"library/#Airfoil-geometries-2","page":"Library","title":"Airfoil geometries","text":"","category":"section"},{"location":"library/","page":"Library","title":"Library","text":"VSFlow.initiator\nVSFlow.carlingfish\nVSFlow.splinefish\nVSFlow.swimmerthickness\nVSFlow.getbodymass\nVSFlow.getinertia\nVSFlow.getcenterofmass\nVSFlow.getcenterofmassvelocity\nVSFlow.correctvelocity","category":"page"},{"location":"library/#VSFlow.initiator","page":"Library","title":"VSFlow.initiator","text":"initiator(ξ, T)\n\nSmooth increase from 0 to 1 in half a period T and its time derivative.\n\n\n\n\n\n","category":"function"},{"location":"library/#VSFlow.carlingfish","page":"Library","title":"VSFlow.carlingfish","text":"carlingfish(s, t, T)\n\nReturn the curvature of a curve of unit length associated with the motion of a carling fish of period T. Its time derivative is computed too.\n\n\n\n\n\n","category":"function"},{"location":"library/#VSFlow.splinefish","page":"Library","title":"VSFlow.splinefish","text":"splinefish(s, t, T, X, K)\n\nReturn the curvature of a curve of unit length associated with the motion of a fish following a cubic spline.\n\n\n\n\n\n","category":"function"},{"location":"library/#VSFlow.swimmerthickness","page":"Library","title":"VSFlow.swimmerthickness","text":"swimmerthickness(s)\n\nReturn the thickness of the swimming body along its coordinate s.\n\n\n\n\n\n","category":"function"},{"location":"library/#VSFlow.getbodymass","page":"Library","title":"VSFlow.getbodymass","text":"getbodymass(XYm, XYt, XYb)\n\nIntegrate the constant and uniform unit density on the body. Integral computed with 25-point Gauss quadrature.\n\n\n\n\n\n","category":"function"},{"location":"library/#VSFlow.getinertia","page":"Library","title":"VSFlow.getinertia","text":"getinertia(XYm, XYt, XYb)\n\nReturn the moment of inertia of the body through 25-point Gauss quatradure on triangles.\n\n\n\n\n\n","category":"function"},{"location":"library/#VSFlow.getcenterofmass","page":"Library","title":"VSFlow.getcenterofmass","text":"getcenterofmass(xm, ym, xb, yb, xt, yt)\n\nCompute the center of mass of the body. Integrals computed with 25-point Gauss quadrature.\n\n\n\n\n\n","category":"function"},{"location":"library/#VSFlow.getcenterofmassvelocity","page":"Library","title":"VSFlow.getcenterofmassvelocity","text":"getcenterofmassvelocity(xm, ym, xb, yb, xt, yt, ub, vb, ut, vt)\n\nCompute the velocity of the center of mass of the body. Integrals computed with bi-linear assumption on quadrilaterals.\n\n\n\n\n\n","category":"function"},{"location":"library/#VSFlow.correctvelocity","page":"Library","title":"VSFlow.correctvelocity","text":"correctvelocity(xm, ym, xb, yb, xt, yt, ub, vb, ut, vt)\n\nCompute the angular velocity (ratio between angular impulse and scalar inertia moment). Integrals computed with bi-linear assumption on quadrilaterals, and 25-point Gauss quadrature.\n\n\n\n\n\n","category":"function"},{"location":"library/#Maths","page":"Library","title":"Maths","text":"","category":"section"},{"location":"library/","page":"Library","title":"Library","text":"VSFlow.rotate\nVSFlow.simps\nVSFlow.trapz\nVSFlow.backwarddifference\nVSFlow.centereddifference\nVSFlow.checkangle\nVSFlow.gausslegendretriangle\nVSFlow.get∫fdV\nVSFlow.cubicspline","category":"page"},{"location":"library/#VSFlow.rotate","page":"Library","title":"VSFlow.rotate","text":"rotate(X, Y, θ)\n\nRotate frame by angle θ in the clockwise direction.\n\n\n\n\n\n","category":"function"},{"location":"library/#VSFlow.simps","page":"Library","title":"VSFlow.simps","text":"simps(f, x)\n\nIntegrate function f using Simpson's rule.\n\n\n\n\n\n","category":"function"},{"location":"library/#VSFlow.trapz","page":"Library","title":"VSFlow.trapz","text":"trapz(f, x)\n\nIntegrate function f using trapeze rule.\n\n\n\n\n\n","category":"function"},{"location":"library/#VSFlow.backwarddifference","page":"Library","title":"VSFlow.backwarddifference","text":"backwarddifference(x, old_x, dt)\n\nCompute the backward difference of nth order according to the length of old_x.\n\n\n\n\n\n","category":"function"},{"location":"library/#VSFlow.centereddifference","page":"Library","title":"VSFlow.centereddifference","text":"centereddifference(neighbours, dx)\n\nCentered differences. Order depends on length of neighbours.\n\n\n\n\n\n","category":"function"},{"location":"library/#VSFlow.checkangle","page":"Library","title":"VSFlow.checkangle","text":"checkangle(θ)\n\nAssert the continuity of panel angles by subtracting 2π if necessary.\n\n\n\n\n\n","category":"function"},{"location":"library/#VSFlow.gausslegendretriangle","page":"Library","title":"VSFlow.gausslegendretriangle","text":"gausslegendretriangle(f, fx, fy)\n\nComputes the integral of function f(x,y) on the standard triangle abc: a (0; 0) b (1; 0) c (0; 1) using 25 points per dimension. fx, fy are the mapping from the standard triangle to the real one. Coefficients are extracted from Rathod et al (2004) \"Gauss Legendre quadrature over a triangle\"\n\n\n\n\n\n","category":"function"},{"location":"library/#VSFlow.get∫fdV","page":"Library","title":"VSFlow.get∫fdV","text":"get∫fdV(XY, f)\n\nComputes the integral of f on the triangle created by the set of points XY, Y using a Gauss-Legendre quadrature. f can return a vector.\n\n\n\n\n\nget∫fdV(pa::LinearPanel, XY, f)\n\nComputes the integral of f on the triangle created by the panel and point XY using a Gauss-Legendre quadrature. f can return a vector.\n\n\n\n\n\n","category":"function"},{"location":"library/#VSFlow.cubicspline","page":"Library","title":"VSFlow.cubicspline","text":"cubicspline(z, X, U)\n\nComputes the cubic spline at point z ∈ [0, 1], build from interpolation of points (X, U) as well as its derivative.\n\n\n\n\n\n","category":"function"},{"location":"manual/#User-Guide","page":"User Guide","title":"User Guide","text":"","category":"section"},{"location":"manual/#Example","page":"User Guide","title":"Example","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = VSFlow","category":"page"},{"location":"#VSFlow","page":"Home","title":"VSFlow","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"An inviscid planar flow model with vortex shedding for external aerodynamics","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"PotentialFlow can be installed using the Julia package manager. From the Julia REPL, type ] to enter the Pkg REPL mode and run","category":"page"},{"location":"","page":"Home","title":"Home","text":"pkg> add VSFlow","category":"page"},{"location":"","page":"Home","title":"Home","text":"The plots in this documentation are generated using Plots.jl.","category":"page"},{"location":"#Basic-Usage","page":"Home","title":"Basic Usage","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Let's simulate an impulsively started NACA0012 with a 10º pitch angle.","category":"page"},{"location":"","page":"Home","title":"Home","text":"using VSFlow\n\n# Simulation parameters\nN = 100          # number of vortex panels on the body surface\nT = 5            # horizon time\ndt = 5e-2        # time step\nanimate = true   # generate animation\n\n# Body geometry and motion\nairfoil_ID = \"My-naca0012\"\nshape = naca00(12.)\nmotion = uniform(1., 0., 0.)\ninitial_position = [0., 0., deg2rad(10.)]\ninitial_velocity = motion(0.)[2]\n\n# Building the `Profile` object\nairfoil = Profile(id = airfoil_ID,\n                    profileshape = shape,\n                    x0 = initial_position,\n                    ẋ0 = initial_velocity,\n                    N = N,\n                    dt = dt,\n                    T = T)\n\nprofilerun(airfoil, motion, animate)\n\nprintln(airfoil.history.acn)","category":"page"},{"location":"","page":"Home","title":"Home","text":"A history of most useful data in the flow is saved at each time step in a data structure called history.","category":"page"},{"location":"","page":"Home","title":"Home","text":"(Image: Impulsively started NACA0012)","category":"page"},{"location":"","page":"Home","title":"Home","text":"Documentation for VSFlow.","category":"page"}]
}
