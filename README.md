# GH_Solar_V2
Solar irradiation model plug-in for Rhinoceros Grasshopper. Simulates hourly annual time series of solar irradiance on mesh vertices.

Please refer to this publication for citation: [Waibel et al. (2017)](http://www.sciencedirect.com/science/article/pii/S0038092X17309349)

Features:
- [Perez (1987)](https://www.sciencedirect.com/science/article/pii/S0038092X87800312) all-weather solar model
- considers snow-coverage
- multi-layer tree approach according to [Krayenhoff (2014)](https://link.springer.com/article/10.1007/s10546-013-9883-1)
- specular inter-reflections max 2 bounces
- diffuse inter-reflections max 1 bounce

Version 2 with better physics than [version 1](https://github.com/christophwaibel/GH_Solar_V1).

How to install:
- Check the releases [here](https://github.com/christophwaibel/GH_Solar_V2/releases) and put the SolarModel.dll and GHSolar.gha into your Rhino Grasshopper components folder.
- For the current built navigate [here](GHSolar/bin) and get the SolarModel.dll and GHSolar.gha
