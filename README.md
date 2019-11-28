# GH_Solar_V2
Solar irradiation model plug-in for Rhinoceros Grasshopper. Simulates hourly annual time series of solar irradiance on mesh vertices.

Please refer to this publication for citation: [Waibel et al. (2017)](http://www.sciencedirect.com/science/article/pii/S0038092X17309349)

Features:
- [Perez (1987)](https://www.sciencedirect.com/science/article/pii/S0038092X87800312) all-weather solar model
- considers snow-coverage
- multi-layer tree approach according to [Krayenhoff (2014](https://link.springer.com/article/10.1007/s10546-013-9883-1)
- specular inter-reflections max 2 bounces
- diffuse inter-reflections max 1 bounce

Version 2 with better physics than [version 1](https://github.com/christophwaibel/GH_Solar_V1).

To Do:
- [ ] specular inter-reflections refraction coefficients
- [ ] Input points and normals, instead of analysis surface mesh... more control on where to place sensor points
- [ ] precise calculation of equinox & solstice dates (SunVector.cs, int [] GetEquinoxSolstice(...))
- [ ] replace rhino libraries for geometry operations with open source libraries (https://doc.cgal.org/ ?). solar.dll should have no rhino dependency. ghsolar.gha should be the rhino wrapper


Bugs:
- [ ] Memory build up over time... get's slower over time... 
- [ ] mesh surfaces far away from the origin become crumpled when used as analysis surface 
- [ ] values underestimated for geometries close to origin (0,0,0)
- [ ] Out of Memory in Rhino... write results into .txt file, reather than into Rhino directly
- [ ] rhino geometry intersection class not thread safe
- [ ] only works properly in meter as workspace unit in Rhino
