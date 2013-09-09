A library for the approximation of the elastodynamics equation 
==============================================================

TSPEED is a library for the high order approximation of the elastodyanmics equation on triangular unstructured meshes.

Installation on Linux
---------------------
The compilation of the library requires CMake. To install, run

    git clone https://github.com/carlomr/tspeed.git  
    cd tspeed  
    mkdir build  
    cmake ..  
    make
    make install  
    make doc  

This will by default install the library in `/usr/local/lib` and the header files in `/usr/local/include/tspeed`. To choose a different installation directory, use

    cmake .. -DCMAKE_INSTALL_PREFIX=/path/to/dir

Two test are generated as `Examples/Lamb` and `Examples/Wedge`. For the source generating them, see `Examples/src/Wedge.cpp` and `Examples/src/Lamb.cpp`.

Problem description
-------------------
The analysis of elastic and acoustic wave propagation phenomena has been
widely studied by mathematicians and scientists since the XIX century.
Elastic waves in solids have also been historically studied, though
analytic solutions are available only for simple domains and settings.
The approximation of the solution to the elastodynamics equation is
therefore of critical importance for the analysis of the propagation of
seismic waves in complex scenarios.

In recent years, seismological, geophysical and technological advances
have allowed for a greater insight into physical seismological events
and this has contributed to the growth of the demand for accurate and
flexible numerical methods. Those, indeed, permit the comparison between
the empirical observations and accurate numerical wave fields in complex
domains.

TSPEED is a C++ library for the
application of the discontinuous spectral element method on meshes made
of simplicial elements to the approximation of the elastodynamic
equation. Spectral element methods were introduced in the computational
fluid dynamics field (A. T. Patera 1984; Maday and Patera 1989) and they
combine the high order accuracy of spectral methods with the flexibility
and computational feasibility of finite elements methods. They are
strictly related with the h-p version of the finite element method and
they have been extensively used for computational geodynamics in the
last two decades (Komatitsch and Vilotte 1998; Komatitsch and Tromp
2002). Spectral element methods have been introduced and are currently
built on mashes made of tensor product elements (i.e. deformed squares
and cubes), since this is the context in which the extension from one
spatial dimension to d dimensions, d=2,3, is more straightforward.
Spectral elements on triangles and tetrahedra have been historically
less widely studied, though different formulations (modal with different
bases (Karniadakis and Sherwin 2005), nodal with different interpolation
nodes (Hesthaven and Warburton 2008)) have been proposed and analyzed in
the last years. Several of these formulations have been employed in
geodynamical applications (Mercerat, Vilotte, and Sánchez-Sesma 2006;
Pelties et al. 2012), but a thorough analysis has not been carried out
for the coupling of modal bases and discontinuous methods. In general,
they provide the flexibility and geometrical adaptability of simplicial
mesh combined with the high order of spectral methods.

Discontinuous methods were introduced for hyperbolic equations, have
been extended to to the elliptic case and developed independently in
both environments. The advantages of discontinuous methods lies in the
fact that they allow for the accurate approximation of sharp gradients
in the solution, that they can be used to develop h-p adaptive strategy
and that the computational cost can be distributed without much
overhead.

Of the spectral bases on triangular elements presented, one (the
Legendre-Dubiner basis (Dubiner 1991; Koornwinder 1975)) can be used
only in the framework of a fully discontinuous approximation, since
there is no way to enforce the continuity of the space between
neighboring elements. The boundary adapted basis functions (Karniadakis
and Sherwin 2005; Dubiner 1991), on the other hand, are modified in
order to be used in a continuous scheme. Therefore, a non-conforming
scheme as in (Antonietti et al. 2012) is possible, and the following
analysis helps to understand the properties of such a scheme.

References
----------
Antonietti, Paola F., Ilario Mazzieri, Alfio Quarteroni, and Francesca
Rapetti. 2012. “Non-conforming high order approximations of the
elastodynamics equation.” *Comput. Methods Appl. Mech. Engrg.* 209/212:
212–238. doi:10.1016/j.cma.2011.11.004.
<http://dx.doi.org/10.1016/j.cma.2011.11.004>.

Dubiner, Moshe. 1991. “Spectral methods on triangles and other domains.”
*Journal of Scientific Computing* 6 (4): 345–390.

Hesthaven, Jan S., and Tim Warburton. 2008. *Nodal discontinuous
Galerkin methods*. *Texts in Applied Mathematics*. Vol. 54. New York:
Springer. <http://dx.doi.org/10.1007/978-0-387-72067-8>.

Karniadakis, George Em, and Spencer J. Sherwin. 2005. *Spectral/\$hp\$
element methods for computational fluid dynamics*. *Numerical
Mathematics and Scientific Computation*. Second.. New York: Oxford
University Press.
<http://dx.doi.org/10.1093/acprof:oso/9780198528692.001.0001>.

Komatitsch, Dimitri, and Jeroen Tromp. 2002. “Spectral-element
simulations of global seismic wave propagation—I. Validation.”
*Geophysical Journal International* 149 (2): 390–412.
doi:10.1046/j.1365-246X.2002.01653.x.
<http://dx.doi.org/10.1046/j.1365-246X.2002.01653.x>.

Komatitsch, Dimitri, and Jean-Pierre Vilotte. 1998. “The spectral
element method: An efficient tool to simulate the seismic response of 2D
and 3D geological structures.” *Bulletin of the Seismological Society of
America* 88 (2): 368–392.
<http://www.bssaonline.org/content/88/2/368.abstract>.

Koornwinder, Tom. 1975. “Two-variable analogues of the classical
orthogonal polynomials.” In *Theory and application of special functions
(Proc. Advanced Sem., Math. Res. Center, Univ. Wisconsin, Madison, Wis.,
1975)*, 435–495. New York: Academic Press.

Maday, Y., and A. T. Patera. 1989. “Spectral element methods for the
incompressible Navier-Stokes equations.” In *State-of-the-art surveys on
computational mechanics (A90-47176 21-64)*, 71–143. American Society of
Mechanical Engineers.

Mercerat, Enrique D., Jean–Pierre Vilotte, and Francisco J.
Sánchez-Sesma. 2006. “Triangular Spectral Element simulation of
two-dimensional elastic wave propagation using unstructured triangular
grids.” *Geophysical Journal International* 166 (2): 679–698.
doi:10.1111/j.1365-246X.2006.03006.x. [\<Go to
ISI\>://WOS:000239004900015](<Go to ISI>://WOS:000239004900015 "<Go to ISI>://WOS:000239004900015").

Patera, Anthony T. 1984. “A spectral element method for fluid dynamics:
Laminar flow in a channel expansion .” *Journal of Computational Physics
*54 (3): 468–488. doi:http://dx.doi.org/10.1016/0021-9991(84)90128-1.
<http://www.sciencedirect.com/science/article/pii/0021999184901281>.

Pelties, Christian, Josep de la Puente, Jean-Paul Ampuero, Gilbert B.
Brietzke, and Martin Käser. 2012. “Three-dimensional dynamic rupture
simulation with a high-order discontinuous Galerkin method on
unstructured tetrahedral meshes.” *Journal of Geophysical Research:
Solid Earth* 117 (B2). doi:10.1029/2011JB008857.
<http://dx.doi.org/10.1029/2011JB008857>.



