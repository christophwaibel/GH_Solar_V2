
using System;
using System.Collections.Generic;
using Grasshopper.Kernel;
using Rhino.Geometry;
using System.Linq;

using SolarModel;

/*
 * GHSolarMeshYear.cs
 * Copyright 2017 Christoph Waibel <chwaibel@student.ethz.ch>
 * 
 * This work is licensed under the GNU GPL license version 3.
*/


namespace GHSolar
{
    public class GHSolarMeshYear : GH_Component
    {

        public GHSolarMeshYear()
            : base("IrradiationMeshYear", "IrrMshY",
                "Calculates hourly solar irradiation on a mesh for an entire year.",
                "EnergyHubs", "Solar Simulation")
        {
        }


        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("AnalysisObjects", "AnalysisObjs", "Input obstacle objects for solar irradiation calculation.", GH_ParamAccess.item);
            pManager.AddGenericParameter("ObstacleObjects", "ObstclObjs", "Input obstacle objects.", GH_ParamAccess.list); //each obstacle consists of: (i) a mesh obstacle, (ii) a 8760-timeseries of 0-1 fraction indicating albedo, (iii) a boolean 0-1 for either specular (0) or diffuse (1) reflection.
            pManager[1].Optional = true;
            pManager.AddGenericParameter("TreeObject", "TreeObj", "Input tree objects.", GH_ParamAccess.list);    //each tree object is : (i) a mesh obstacle, (ii) a 8760-timeseries of 0-1 fractions, indicating leave-coverage 1 is full of leaves=full obstruction.
            pManager[2].Optional = true;
            
            
            //________________________________________________________________________________________________________________________________________________________________
            pManager.AddIntegerParameter("////////////////////", "////////////////////", "////////////////////", GH_ParamAccess.item);
            pManager[3].Optional = true;
            //________________________________________________________________________________________________________________________________________________________________
            
            
            pManager.AddNumberParameter("Latitude", "Latitude", "Latitude of the location in [°].", GH_ParamAccess.item);
            pManager.AddNumberParameter("Longitude", "Longitude", "Longitude of the location in [°].", GH_ParamAccess.item);
            pManager.AddIntegerParameter("Year", "Year", "Year", GH_ParamAccess.item);
            pManager.AddNumberParameter("SnwThrsh", "SnwThrsh", "Snow thickness threshold, before which snow does not remain on surface. Default is 1.", GH_ParamAccess.item);
            pManager[7].Optional = true;
            pManager.AddNumberParameter("TiltThrsh", "TiltThrsh", "Tilst threshold, after which snow does not remain on surface. Default is 60.", GH_ParamAccess.item);
            pManager[8].Optional = true;
            
            
            //________________________________________________________________________________________________________________________________________________________________
            pManager.AddIntegerParameter("////////////////////", "////////////////////", "////////////////////", GH_ParamAccess.item);
            pManager[9].Optional = true;
            //________________________________________________________________________________________________________________________________________________________________
            
            
            pManager.AddNumberParameter("DNI", "DNI", "Direct normal irradiation 8760-time series.", GH_ParamAccess.list);
            pManager.AddNumberParameter("DHI", "DHI", "Diffuse horizontal irradiation 8760-time series.", GH_ParamAccess.list);
            pManager.AddNumberParameter("Snow", "Snow", "Snow coverage 8760-time series. Default is 0 throughout the year. Used for coverage of analysis surfaces. Surfaces with an inclination angle higher than the specified threshold are not affected.", GH_ParamAccess.list);
            pManager[12].Optional = true;
            pManager.AddNumberParameter("GrndAlb", "GrndAlb", "Ground albedo, 8760 time series. Default is 0.2 throughout the year.", GH_ParamAccess.list);
            pManager[13].Optional = true;
            pManager.AddNumberParameter("Azimuth", "Azimuth", "Solar Azimuth 8760-time series. Optional. If not data is provided, the solar vectors are calculated according to latitude, longitude, timezone and year.", GH_ParamAccess.list);
            pManager[14].Optional = true;
            pManager.AddNumberParameter("Altitude", "Altitude", "Solar Altitude 8760-time series. Optional. If not data is provided, the solar vectors are calculated according to latitude, longitude, timezone and year.", GH_ParamAccess.list);
            pManager[15].Optional = true;
            
            
            //________________________________________________________________________________________________________________________________________________________________
            pManager.AddIntegerParameter("////////////////////", "////////////////////", "////////////////////", GH_ParamAccess.item);
            pManager[16].Optional = true;
            //________________________________________________________________________________________________________________________________________________________________
            
            
            pManager.AddIntegerParameter("MainSkyRes", "MainSkyRes", "Sykdome resolution for diffuse shading mask. I.e. recursion level of the icosahedron hemisphere. 0: 10 rays; 1: 29 rays; 2: 97 rays; 3: 353 rays. Default is (2).", GH_ParamAccess.item);
            pManager[17].Optional = true;
            pManager.AddIntegerParameter("MainInterpMode", "MainInterpMode", "Interpolation mode for beam obstruction calculations. (0) = 3-days; (1) = 12-days. Default is (1).", GH_ParamAccess.item);
            pManager[18].Optional = true;
            
            
            //________________________________________________________________________________________________________________________________________________________________
            pManager.AddIntegerParameter("////////////////////", "////////////////////", "////////////////////", GH_ParamAccess.item);
            pManager[19].Optional = true;
            //________________________________________________________________________________________________________________________________________________________________
            
            
            pManager.AddIntegerParameter("SpecBounces", "SpecBounces", "Number of specular bounces for inter-reflections. 0 (min) - 2 (max). Default is (1).", GH_ParamAccess.item);
            pManager[20].Optional = true;
            pManager.AddIntegerParameter("SpecInterpMode", "SpecInterpMode", "Interpolation mode for specular inter-reflections obstruction calculations. (0) = 3-days; (1) = 12-days. Default is (0).", GH_ParamAccess.item);
            pManager[21].Optional = true;


            //________________________________________________________________________________________________________________________________________________________________
            pManager.AddIntegerParameter("////////////////////", "////////////////////", "////////////////////", GH_ParamAccess.item);
            pManager[22].Optional = true;
            //________________________________________________________________________________________________________________________________________________________________


            pManager.AddIntegerParameter("DiffIReflSkyRes", "DiffIReflSkyRes", "SkyDome resolution for diffuse shading mask, i.e. rays leaving maing sensor point. (0): 10 rays; (1): 29 rays; (2): 97 rays. Default is (0).", GH_ParamAccess.item);
            pManager[23].Optional = true;
            pManager.AddIntegerParameter("DiffIReflSkyRes2nd", "DiffIReflSkyRes2nd", "SkyDome resolution for diffuse shading mask of secondary sensor points, i.e. rays leaving secondary sensor points. (0): 10 rays; (1): 29 rays; (2): 97 rays. Default is (0).", GH_ParamAccess.item);
            pManager[24].Optional = true;
            pManager.AddIntegerParameter("DiffIReflMode", "DiffIReflMode", "Diffuse inter-reflection mode. (0) = no secondary sensor points, only unobstructed simple ground reflection.; (1) = secondary sensor points without obstruction calculation; (2) = secondary sensor points with 3-days beam interpolation and skydomes of resolution diffIReflSkyRes2nd; Default is (1).", GH_ParamAccess.item);
            pManager[25].Optional = true;


            //________________________________________________________________________________________________________________________________________________________________
            pManager.AddIntegerParameter("////////////////////", "////////////////////", "////////////////////", GH_ParamAccess.item);
            pManager[26].Optional = true;
            //________________________________________________________________________________________________________________________________________________________________

    
            pManager.AddGenericParameter("PrecalcIrefl", "PrecalcIrefl", "Input precalculated interreflected specular and diffuse irradiation.", GH_ParamAccess.item);
            pManager[27].Optional = true;


            //________________________________________________________________________________________________________________________________________________________________
            pManager.AddIntegerParameter("////////////////////", "////////////////////", "////////////////////", GH_ParamAccess.item);
            pManager[28].Optional = true;
            //________________________________________________________________________________________________________________________________________________________________


            pManager.AddBooleanParameter("MT", "MT", "Multi threading. Default is (false).", GH_ParamAccess.item);
            pManager[29].Optional = true;
        }


        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
              pManager.AddGenericParameter("Results", "Results", "Results data of solar irradiation calculation", GH_ParamAccess.item);
        }


        protected override void SolveInstance(IGH_DataAccess DA)
        {
            //___________________________________________________________________
            /////////////////////////////////////////////////////////////////////
            //INPUTS
            Mesh msh = new Mesh();
            CObstacleObject mshobj = null;
            if (!DA.GetData(0, ref mshobj)) { return; }
            msh = mshobj.mesh;

            //should containt analysis surface itself
            List<CObstacleObject> objObst = new List<CObstacleObject>();
            if (!DA.GetDataList(1, objObst)) { return; }

            List<CPermObject> treeObst = new List<CPermObject>();
            DA.GetDataList(2, treeObst);

            double latitude = 0.0;
            if (!DA.GetData(4, ref latitude)) { return; }
            double longitude = 0.0;
            if (!DA.GetData(5, ref longitude)) { return; }

            int year = 0;
            if (!DA.GetData(6, ref year)) { return; }

            double snow_threshold = 1.0;
            if (!DA.GetData(7, ref snow_threshold)) { snow_threshold = 1.0; }
            double tilt_threshold = 60.0;
            if (!DA.GetData(8, ref tilt_threshold)) { tilt_threshold = 60.0; }


            List<double> DNI = new List<double>();
            if (!DA.GetDataList(10, DNI)) { return; }
            List<double> DHI = new List<double>();
            if (!DA.GetDataList(11, DHI)) { return; }
            List<double> SNOW = new List<double>();
            if (!DA.GetDataList(12, SNOW)) { for (int t = 0; t < 8760; t++) SNOW.Add(0.0); }
            List<double> groundalbedo = new List<double>();
            if (!DA.GetDataList(13, groundalbedo)) { for (int t = 0; t < 8760; t++) groundalbedo.Add(0.2); }
            List<double> solarAzimuth = new List<double>();
            DA.GetDataList(14, solarAzimuth);
            List<double> solarAltitude = new List<double>();
            DA.GetDataList(15, solarAltitude);


            int MainSkyRes = 0;
            if (!DA.GetData(17, ref MainSkyRes)) { MainSkyRes = 1; }
            int MainInterpMode = 0;
            if (!DA.GetData(18, ref MainInterpMode)) { MainInterpMode = 0; }


            int SpecBounces = 1;
            if (!DA.GetData(20, ref SpecBounces)) { SpecBounces = 1; }
            SpecBounces = (SpecBounces < 0) ? 0 : SpecBounces;
            int SpecInterpMode = 0;
            if (!DA.GetData(21, ref SpecInterpMode)) { SpecInterpMode = 0; }


            int DiffIReflSkyRes = 0;
            if (!DA.GetData(23, ref DiffIReflSkyRes)) { DiffIReflSkyRes = 0; }
            int DiffIReflSkyRes2nd = 0;
            if (!DA.GetData(24, ref DiffIReflSkyRes2nd)) { DiffIReflSkyRes2nd = 0; }
            int DiffIReflMode = 1;
            if (!DA.GetData(25, ref DiffIReflMode)) { DiffIReflMode = 1; }
            
            CResultsInterreflections ResultsIreflIn = null;
            DA.GetData(27, ref ResultsIreflIn);

            bool mt = false;
            if (!DA.GetData(29, ref mt)) { mt = false; }

            /////////////////////////////////////////////////////////////////////
            //___________________________________________________________________




            #region SIMULATE
            //______________________________________________________________________________________________
            ////////////////////////////////        S I M U L A T E      ///////////////////////////////////
            CCalculateSolarMesh calc = new CCalculateSolarMesh(
                mshobj, objObst, treeObst, latitude, longitude, DNI, DHI, SNOW, groundalbedo, snow_threshold, tilt_threshold,
                year, null, mt, solarAzimuth, solarAltitude);
            if (mt)
            {
                calc.RunAnnualSimulation_MT(mshobj.tolerance,
                    MainSkyRes, MainInterpMode, SpecBounces, SpecInterpMode, DiffIReflSkyRes, DiffIReflSkyRes2nd, DiffIReflMode);
            }
            else
            {
                calc.RunAnnualSimulation(mshobj.tolerance,
                    MainSkyRes, MainInterpMode, SpecBounces, SpecInterpMode, DiffIReflSkyRes, DiffIReflSkyRes2nd, DiffIReflMode);
            }
            CResults results = calc.getResults();
            //cResultsInterreflections resultsIreflOut = calc.getResultsInterreflections();
            ////////////////////////////////////////////////////////////////////////////////////////////////
            //______________________________________________________________________________________________
            #endregion






            #region OUTPUT
            //______________________________________________________________________________________________
            ////////////////////////////////          O U T P U T        ///////////////////////////////////
            DA.SetData(0, results);

            ////////////////////////////////////////////////////////////////////////////////////////////////
            //______________________________________________________________________________________________
            #endregion



            Rhino.RhinoApp.WriteLine("SOLAR MODEL... Done");
        }



        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                //You can add image files to your project resources and access them like this:
                // return Resources.IconForThisComponent;
                return GHSolar.Properties.Resources.pic_irrad_y;
            }
        }



        public override Guid ComponentGuid
        {
            get { return new Guid("{6d5b9658-b84d-47ee-a5a1-1dd080c0caf1}"); }
        }
    }
}
