var documenterSearchIndex = {"docs":
[{"location":"library/#Public-API","page":"Library","title":"Public API","text":"","category":"section"},{"location":"library/#Airfoil-geometries","page":"Library","title":"Airfoil geometries","text":"","category":"section"},{"location":"library/","page":"Library","title":"Library","text":"naca00\ngaw1\ncircle\nellipse","category":"page"},{"location":"library/#VSFlow.naca00","page":"Library","title":"VSFlow.naca00","text":"naca00(x; kwargs...)\n\nReturn the coordinates of a NACA00ZZ airfoil.\n\nKeyword arguments:\n\nZZ = 12: max thickness (in % of the chord) of the airfoil.\n\n\n\n\n\n","category":"function"},{"location":"library/#VSFlow.gaw1","page":"Library","title":"VSFlow.gaw1","text":"gaw1(x)\n\nReturn the coordinates of a GAW1 airfoil.\n\n\n\n\n\n","category":"function"},{"location":"library/#VSFlow.circle","page":"Library","title":"VSFlow.circle","text":"circle(x; kwargs...)\n\nReturn the coordinates of a circle.\n\nKeyword arguments:\n\nR = .5: radius of the circle.\n\n\n\n\n\n","category":"function"},{"location":"library/#VSFlow.ellipse","page":"Library","title":"VSFlow.ellipse","text":"ellipse(x; kwargs...)\n\nReturn the coordinates of an ellipse.\n\nKeyword arguments:\n\nB = .5: half height of the ellipse.\n\n\n\n\n\n","category":"function"},{"location":"library/#Airfoil-motion","page":"Library","title":"Airfoil motion","text":"","category":"section"},{"location":"library/","page":"Library","title":"Library","text":"heavepitch\ncircularmotion","category":"page"},{"location":"library/#VSFlow.heavepitch","page":"Library","title":"VSFlow.heavepitch","text":"heavepitch(t, h0, αmax, ψ1, ψ2, strouhal)\n\nReturn the vector [x, ẋ, ẍ] corresponding to a heaving & pitching motion where x = [X, Y, α].\n\n\n\n\n\n","category":"function"},{"location":"library/#VSFlow.circularmotion","page":"Library","title":"VSFlow.circularmotion","text":"circularmotion(t, R, strouhal)\n\nReturn the vector [x, ẋ, ẍ] corresponding to a circular motion of radius R where x = [X, Y, α].\n\n\n\n\n\n","category":"function"},{"location":"library/#Running-the-simulation","page":"Library","title":"Running the simulation","text":"","category":"section"},{"location":"library/","page":"Library","title":"Library","text":"Profile(id, profile::Function, N, position, dt, T, δ, ϵ, rng, is_comb, lump, args...)\nsetinitvelocity\nprofilerun","category":"page"},{"location":"library/#VSFlow.Profile-Tuple{Any, Function, Any, Any, Any, Any, Any, Any, Any, Any, Any, Vararg{Any, N} where N}","page":"Library","title":"VSFlow.Profile","text":"Profile(id, profile::Function, N, position, dt, T, δ, ϵ, rng, is_comb, lump, args...)\n\nConstructs a Profile object.\n\nArguments\n\nid: String used to identify the profile.\nprofile: function determining the profile shape.\nN: number of panels on the surface of the profile.\nposition: position of the pivot point of the profile.\ndt: timestep of the simulation.\nT: horizon time.\nδ: kernel cut-off width.\nϵ: percentage of the average panel length to chop the trailing edge.\nrng: random number generator\nis_comb: Bool indicating whether a random comb of dipoles hit the body.\nlump: Bool indicating whether the lumging of vortices should occur.\n\nAdditional Arguments\n\nη: maximum error due to lumping.\nTmin: minimum time between two successive active vortices.\nsheet_size: minimum vortex sheet length in the wake.\nZZ: body maximum thickness (in % of the chord).\n\n\n\n\n\n","category":"method"},{"location":"library/#VSFlow.setinitvelocity","page":"Library","title":"VSFlow.setinitvelocity","text":"setinitvelocity(p::Profile, u, v, α̇)\n\nSet the initial velocities of the body. Note: α̇ represents the rate of change in AoA and is thus opposite to the angular velocity.\n\n\n\n\n\n","category":"function"},{"location":"library/#VSFlow.profilerun","page":"Library","title":"VSFlow.profilerun","text":"profilerun(p::Profile, is_write, is_lumped, is4thorder, forcefun, motion_args, show=false)\n\nSimulate the testcase.\n\nArguments\n\np: Profile simulation to run.\nis_write: Bool indicating if the results are saved into a file.\nis_lumped: Bool indicating if vortices have to be lumped together.\nis4thorder: Bool indicating if the time marching is RK4 or ForwardEuler.\nforcefun: function giving the force applied on the pivot point of the body.\nmotion_args: array containing parameters to feed forcefun with.\nshow: Bool indicating whether an animation has to be generated.\n\n\n\n\n\n","category":"function"},{"location":"library/#Quantities-of-interest","page":"Library","title":"Quantities of interest","text":"","category":"section"},{"location":"library/","page":"Library","title":"Library","text":"getimpulse(p::Profile)\ngetboundpotential\ngetcoeffimpulse!\ngetcoeffpotential!\ngetnoca!","category":"page"},{"location":"library/#VSFlow.getimpulse-Tuple{Profile}","page":"Library","title":"VSFlow.getimpulse","text":"getimpulse(p::Profile)\n\nReturn the linear/angular impulse exerted by the fluid on the body.\n\n\n\n\n\n","category":"method"},{"location":"library/#VSFlow.getboundpotential","page":"Library","title":"VSFlow.getboundpotential","text":"getboundpotential(p::Profile)\n\nIntegrates the velocity on the profile to obtain the velocity potential ϕ.\n\n\n\n\n\n","category":"function"},{"location":"library/#VSFlow.getcoeffimpulse!","page":"Library","title":"VSFlow.getcoeffimpulse!","text":"getcoeffimpulse!()\n\nReturn aerodynamic coefficients computed with impulse conservation in the flow.\n\nKeyword Arguments\n\nupdate = true: indicates whether the impulse history has to be updated.\nβ = 1.0: portion of a timestep.\n\n\n\n\n\n","category":"function"},{"location":"library/#VSFlow.getcoeffpotential!","page":"Library","title":"VSFlow.getcoeffpotential!","text":"getcoeffpotential!(p::Profile)\n\nReturn the aerodynamic coefficients through pressure integration on the body.\n\n\n\n\n\n","category":"function"},{"location":"library/#VSFlow.getnoca!","page":"Library","title":"VSFlow.getnoca!","text":"getnoca!(p::Profile, update=true)\n\nReturn aerodynamic coefficients computed with a control volume approach.\n\nKeyword Arguments\n\nupdate = true: indicates whether the impulse history has to be updated.\n\n\n\n\n\n","category":"function"},{"location":"library/#Private-API","page":"Library","title":"Private API","text":"","category":"section"},{"location":"library/#Vortex-Elements","page":"Library","title":"Vortex Elements","text":"","category":"section"},{"location":"library/","page":"Library","title":"Library","text":"VortexPoint\nConstantPanel\nLinearPanel","category":"page"},{"location":"library/#Airfoil-geometries-2","page":"Library","title":"Airfoil geometries","text":"","category":"section"},{"location":"library/","page":"Library","title":"Library","text":"VSFlow.initiator\nVSFlow.carlingfish\nVSFlow.splinefish\nVSFlow.swimmerthickness\nVSFlow.getbodymass\nVSFlow.getinertia\nVSFlow.getcenterofmass\nVSFlow.getcenterofmassvelocity\nVSFlow.correctvelocity","category":"page"},{"location":"library/#VSFlow.initiator","page":"Library","title":"VSFlow.initiator","text":"initiator(ξ, T)\n\nSmooth increase from 0 to 1 in half a period T and its time derivative.\n\n\n\n\n\n","category":"function"},{"location":"library/#VSFlow.carlingfish","page":"Library","title":"VSFlow.carlingfish","text":"carlingfish(s, t, T)\n\nReturn the curvature of a curve of unit length associated with the motion of a carling fish of period T. Its time derivative is computed too.\n\n\n\n\n\n","category":"function"},{"location":"library/#VSFlow.splinefish","page":"Library","title":"VSFlow.splinefish","text":"splinefish(s, t, T, X, K)\n\nReturn the curvature of a curve of unit length associated with the motion of a fish following a cubic spline.\n\n\n\n\n\n","category":"function"},{"location":"library/#VSFlow.swimmerthickness","page":"Library","title":"VSFlow.swimmerthickness","text":"swimmerthickness(s)\n\nReturn the thickness of the swimming body along its coordinate s.\n\n\n\n\n\n","category":"function"},{"location":"library/#VSFlow.getbodymass","page":"Library","title":"VSFlow.getbodymass","text":"getbodymass(XYm, XYt, XYb)\n\nIntegrate the constant and uniform unit density on the body. Integral computed with 25-point Gauss quadrature.\n\n\n\n\n\n","category":"function"},{"location":"library/#VSFlow.getinertia","page":"Library","title":"VSFlow.getinertia","text":"getinertia(XYm, XYt, XYb)\n\nReturn the moment of inertia of the body through 25-point Gauss quatradure on triangles.\n\n\n\n\n\n","category":"function"},{"location":"library/#VSFlow.getcenterofmass","page":"Library","title":"VSFlow.getcenterofmass","text":"getcenterofmass(xm, ym, xb, yb, xt, yt)\n\nCompute the center of mass of the body. Integrals computed with 25-point Gauss quadrature.\n\n\n\n\n\n","category":"function"},{"location":"library/#VSFlow.getcenterofmassvelocity","page":"Library","title":"VSFlow.getcenterofmassvelocity","text":"getcenterofmassvelocity(xm, ym, xb, yb, xt, yt, ub, vb, ut, vt)\n\nCompute the velocity of the center of mass of the body. Integrals computed with bi-linear assumption on quadrilaterals.\n\n\n\n\n\n","category":"function"},{"location":"library/#VSFlow.correctvelocity","page":"Library","title":"VSFlow.correctvelocity","text":"correctvelocity(xm, ym, xb, yb, xt, yt, ub, vb, ut, vt)\n\nCompute the angular velocity (ratio between angular impulse and scalar inertia moment). Integrals computed with bi-linear assumption on quadrilaterals, and 25-point Gauss quadrature.\n\n\n\n\n\n","category":"function"},{"location":"library/#Maths","page":"Library","title":"Maths","text":"","category":"section"},{"location":"library/","page":"Library","title":"Library","text":"VSFlow.rotate\nVSFlow.simps\nVSFlow.trapz\nVSFlow.backwarddifference\nVSFlow.centereddifference\nVSFlow.checkangle\nVSFlow.gausslegendretriangle\nVSFlow.get∫fdV\nVSFlow.cubicspline","category":"page"},{"location":"library/#VSFlow.rotate","page":"Library","title":"VSFlow.rotate","text":"rotate(X, Y, θ)\n\nRotate frame by angle θ in the clockwise direction.\n\n\n\n\n\n","category":"function"},{"location":"library/#VSFlow.simps","page":"Library","title":"VSFlow.simps","text":"simps(f, x)\n\nIntegrate function f using Simpson's rule.\n\n\n\n\n\n","category":"function"},{"location":"library/#VSFlow.trapz","page":"Library","title":"VSFlow.trapz","text":"trapz(f, x)\n\nIntegrate function f using trapeze rule.\n\n\n\n\n\n","category":"function"},{"location":"library/#VSFlow.backwarddifference","page":"Library","title":"VSFlow.backwarddifference","text":"backwarddifference(x, old_x, dt)\n\nCompute the backward difference of nth order according to the length of old_x.\n\n\n\n\n\n","category":"function"},{"location":"library/#VSFlow.centereddifference","page":"Library","title":"VSFlow.centereddifference","text":"centereddifference(neighbours, dx)\n\nCentered differences. Order depends on length of neighbours.\n\n\n\n\n\n","category":"function"},{"location":"library/#VSFlow.checkangle","page":"Library","title":"VSFlow.checkangle","text":"checkangle(θ)\n\nAssert the continuity of panel angles by subtracting 2π if necessary.\n\n\n\n\n\n","category":"function"},{"location":"library/#VSFlow.gausslegendretriangle","page":"Library","title":"VSFlow.gausslegendretriangle","text":"gausslegendretriangle(f, fx, fy)\n\nComputes the integral of function f(x,y) on the standard triangle abc: a (0; 0) b (1; 0) c (0; 1) using 25 points per dimension. fx, fy are the mapping from the standard triangle to the real one. Coefficients are extracted from Rathod et al (2004) \"Gauss Legendre quadrature over a triangle\"\n\n\n\n\n\n","category":"function"},{"location":"library/#VSFlow.get∫fdV","page":"Library","title":"VSFlow.get∫fdV","text":"get∫fdV(XY, f)\n\nComputes the integral of f on the triangle created by the set of points XY, Y using a Gauss-Legendre quadrature. f can return a vector.\n\n\n\n\n\nget∫fdV(pa::LinearPanel, XY, f)\n\nComputes the integral of f on the triangle created by the panel and point XY using a Gauss-Legendre quadrature. f can return a vector.\n\n\n\n\n\n","category":"function"},{"location":"library/#VSFlow.cubicspline","page":"Library","title":"VSFlow.cubicspline","text":"cubicspline(z, X, U)\n\nComputes the cubic spline at point z ∈ [0, 1], build from interpolation of points (X, U) as well as its derivative.\n\n\n\n\n\n","category":"function"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = VSFlow","category":"page"},{"location":"#VSFlow","page":"Home","title":"VSFlow","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for VSFlow.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"}]
}
