# GH_Solar_V2
Solar irradiation model plug-in for Rhinoceros Grasshopper. Version 2.

Please refer to this publication for citation: [Waibel et al. (2017)](http://www.sciencedirect.com/science/article/pii/S0038092X17309349)

To Do:
- [x] Perez diffuse sky model
- [x] skydome with equal patches
- [x] snow cover
- [x] trees
- [x] annual interpolation trees
- [x] specular inter-reflections max 2 bounces
- [x] annual interpolation specular inter-reflections
- [x] diffuse inter-reflections max 1 bounce
- [x] annual interpolation diffuse inter-reflections
- [ ] specular inter-reflections refraction coefficients
- [x] multi-threading annual
- [ ] multi-threading hourly
- [ ] Input points and normals, instead of analysis surface mesh... more control on where to place sensor points
- [ ] precise calculation of equinox & solstice dates (SunVector.cs, int [] GetEquinoxSolstice(...))


Bugs:
- [ ] mesh surfaces far away from the origin become crumpled when used as analysis surface 
- [ ] values underestimated for geometries close to origin (0,0,0)
- [ ] Out of Memory in Rhino... write results into .txt file, reather than into Rhino directly
- [ ] rhino geometry intersection class not thread safe
- [ ] only works in meter as workspace unit in Rhino
