# bigonlight1.0
# Author: Michele Grasso 2019-2021 (Center for Theoretical Physics, Polish Academy of Sciences)
# Copyright: GPL

# Summary:
The BiGONLight is a Mathematica package which encode the the Bi-local Geodesic Operators formalism (BGO) to study light propagation in the geometric optics regime in General Relativit. The package contains a collection of function, including those to compute geodesics, parallel transported frames and solve the BGO’s equation. BiGONLight.m works as an external library that, once is called by a Mathematica notebook, the user can specify the metric g and S and O four-velocity and four-acceleration: in the case that the metric is obtained from a numerical simulation, one can use as input the spatial metric γ, the extrinsic curvature K and the gauge functions lapse α and shift β obtained from the simulation. However, the code also can accept 4D metrics given in analytic form, i.e. with the components of the metric as arbitrarily complicated functions of the spacetime coordinates, and use the powerful Mathematica symbolic algebra manipulation to perform the 3 + 1 splitting of the metric g, to compute the lapse α, the shift β, the spatial metric γ, the extrinsic curvature K and the normal vector n. This is done by the function ADM[] encoded into the package. On top of that, the BiGONLight contains other functions to calculate useful quantities in differential geometry like Christoffel[] and Riemann[] to compute the Christoffel symbols and the Riemann tensor. Together with the powerfull symbolic algebra manipulation, Mathematica also provides a large variety of numerical methods that can be used and customized to adapt them to the particular problem. For instance, this was implemented with SolveGeodesic[], SolveEnergy[], TransportedFrame[] and SolveBGO[] functions in which the user can chose the numerical methods used to solve the system of ODE. 
Another useful feature is the Mathematica’s precision control options, which allows the user to set the precision (and accuracy) of numerical result through the commands WorkingPrecision, SetPrecision and SetAccuracy.

The main achievement of our package is to simulate light propagation in numerical relativity to extract observables. The recipe to obtain the observables using the functions defined in the BiGONLight package can be summarized as follows:
    1. first, one need to specify the S and O kinematics and the metric, which can be given already through its 3 + 1 components (α, β, γ and K) or as a 4D metric tensor and use ADM[] to do the splitting.
    2. set up the null geodesic giving the initial tangent vector or using InitialConditions[]
    3. find the geodesic equations in 3+1 using GeodesicEquations[] and EnergyEquations[] and then solve them with SolveGeodesic[] and SolveEnergy[]
    4. set up the initial conditions for the frame that will be parallel transported along the line of sight (SNF). The parallel transport of the frame is obtained using the function PTransportedFrame[]
    5. use the OpticalTidalMatrix[] to compute the optical tidal matrix projected into the parallel transported frame. OpticalTidalMatrix[] take the 3+1 components of the metric to compute the Riemann tensor in terms of its 3+1 components (i.e. the Gauss relation, the Codazzi relation and the Ricci relation) and then does the contraction with the 3+1 tangent vector and the projection into the parallel transported frame (in terms of its 3+1 components).
    6. compose the geodesic deviation equation for the BGO with BGOequations[] and solve them using SolveBGO[]
    7. use the expression for the observables in terms of the BGO to compute optical observables


# Additional information:
More information about the BGO framework can be found in Phys.Rev.D 99,064038 (2019) (arXiv:1811.10284) or Phys.Rev.D 101,063506 (2020) (arXiv:1912.04988). 
More information about the package can be found in Grasso, Villa (2021) and in Grasso, Villa, Korzynski, Matarrese (2021).

# License: 
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the ollowing conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

The Software is provided “as is”, without warranty of any kind, expressed or implied, including but not limited to the warranties of merchantability, fitness for a particular purpose and noninfringement. In no event shall the authors or copyright holders be liable for any claim, damages or other liability, whether in an action of contract, tort or otherwise, arising from, out of or in connection with the Software or the use or other dealings in the Software.

