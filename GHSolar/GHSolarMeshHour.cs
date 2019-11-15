using System;
using System.Collections.Generic;
using Grasshopper.Kernel;
using Rhino.Geometry;

using SolarModel;
using System.Threading.Tasks;


/*
 * GHSOlarMeshHour.cs
 * Copyright 2017 Christoph Waibel <chwaibel@student.ethz.ch>
 * 
 * This work is licensed under the GNU GPL license version 3.
*/

namespace GHSolar
{
    /// <summary>
    /// Calculates solar irradation on a mesh for one hour of the year.
    /// </summary>
    public class GHSolarMeshHour : GH_Component
    {
        public GHSolarMeshHour()
            : base("IrradiationMeshhour", "IrrMshH",
                "Calculates solar irradiation on a mesh for a certain hour of the year.",
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
            pManager.AddIntegerParameter("Month", "Month", "Month ∈ [1, 12]", GH_ParamAccess.item);
            pManager.AddIntegerParameter("Day", "Day", "Day ∈ [1, 31]", GH_ParamAccess.item);
            pManager.AddIntegerParameter("Hour", "Hour", "Hour ∈ [0, 23]", GH_ParamAccess.item);
            pManager.AddNumberParameter("SnwThrsh", "SnwThrsh", "Snow thickness threshold, before which snow does not remain on surface. Default is 1.", GH_ParamAccess.item);
            pManager[10].Optional = true;
            pManager.AddNumberParameter("TiltThrsh", "TiltThrsh", "Tilst threshold, after which snow does not remain on surface. Default is 60.", GH_ParamAccess.item);
            pManager[11].Optional = true;


            //________________________________________________________________________________________________________________________________________________________________
            pManager.AddIntegerParameter("////////////////////", "////////////////////", "////////////////////", GH_ParamAccess.item);
            pManager[12].Optional = true;
            //________________________________________________________________________________________________________________________________________________________________


            pManager.AddNumberParameter("DNI", "DNI", "Direct normal irradiation 8760-time series.", GH_ParamAccess.list);
            pManager.AddNumberParameter("DHI", "DHI", "Diffuse horizontal irradiation 8760-time series.", GH_ParamAccess.list);
            pManager.AddNumberParameter("Snow", "Snow", "Snow coverage 8760-time series. Default is 0 throughout the year. Used for coverage of analysis surfaces. Surfaces with an inclination angle higher than the specified threshold are not affected.", GH_ParamAccess.list);
            pManager[15].Optional = true;
            pManager.AddNumberParameter("GrndAlb", "GrndAlb", "Ground albedo, 8760 time series. Default is 0.2 throughout the year.", GH_ParamAccess.list);
            pManager[16].Optional = true;
            pManager.AddNumberParameter("Azimuth", "Azimuth", "Solar Azimuth 8760-time series. Optional. If not data is provided, the solar vectors are calculated according to latitude, longitude, timezone and year.", GH_ParamAccess.list);
            pManager[17].Optional = true;
            pManager.AddNumberParameter("Altitude", "Altitude", "Solar Altitude 8760-time series. Optional. If not data is provided, the solar vectors are calculated according to latitude, longitude, timezone and year.", GH_ParamAccess.list);
            pManager[18].Optional = true;


            //________________________________________________________________________________________________________________________________________________________________
            pManager.AddIntegerParameter("////////////////////", "////////////////////", "////////////////////", GH_ParamAccess.item);
            pManager[19].Optional = true;
            //________________________________________________________________________________________________________________________________________________________________


            pManager.AddIntegerParameter("MainSkyRes", "MainSkyRes", "Sykdome resolution for diffuse shading mask. I.e. recursion level of the icosahedron hemisphere. 0: 10 rays; 1: 29 rays; 2: 97 rays; 3: 353 rays. Default is (2).", GH_ParamAccess.item);
            pManager[20].Optional = true;
            pManager.AddIntegerParameter("SpecBounces", "SpecBounces", "Number of specular bounces for inter-reflections. 0 (min) - 2 (max). Default is (1).", GH_ParamAccess.item);
            pManager[21].Optional = true;
            pManager.AddIntegerParameter("DiffIReflSkyRes", "DiffIReflSkyRes", "SkyDome resolution for diffuse shading mask, i.e. rays leaving maing sensor point. (0): 10 rays; (1): 29 rays; (2): 97 rays. Default is (0). (-1): simple unobstructed ground reflection only. ", GH_ParamAccess.item);
            pManager[22].Optional = true;
            pManager.AddIntegerParameter("DiffIReflSkyRes2nd", "DiffIReflSkyRes2nd", "SkyDome resolution for diffuse shading mask of secondary sensor points, i.e. rays leaving secondary sensor points. (0): 10 rays; (1): 29 rays; (2): 97 rays. Default is (0).", GH_ParamAccess.item);
            pManager[23].Optional = true;


            //________________________________________________________________________________________________________________________________________________________________
            pManager.AddIntegerParameter("////////////////////", "////////////////////", "////////////////////", GH_ParamAccess.item);
            pManager[24].Optional = true;
            //________________________________________________________________________________________________________________________________________________________________


            pManager.AddGenericParameter("PrecalcIrefl", "PrecalcIrefl", "Input precalculated interreflected specular and diffuse irradiation.", GH_ParamAccess.item);
            pManager[25].Optional = true;


            //________________________________________________________________________________________________________________________________________________________________
            pManager.AddIntegerParameter("////////////////////", "////////////////////", "////////////////////", GH_ParamAccess.item);
            pManager[26].Optional = true;
            //________________________________________________________________________________________________________________________________________________________________


            pManager.AddBooleanParameter("MT", "MT", "Multi threading. Default is (false).", GH_ParamAccess.item);
            pManager[27].Optional = true;
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddLineParameter("Vec", "Vec", "Solar vector of current hour.", GH_ParamAccess.item);
            pManager.AddGenericParameter("Results", "Results", "Results data of solar irradiation calculation", GH_ParamAccess.item);
            pManager.AddGenericParameter("Irefl", "Irefl", "Irradiation from inter-reflections only.", GH_ParamAccess.item);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {       
            #region INPUTS
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
            int month = 0;
            if (!DA.GetData(7, ref month)) { return; }
            int day = 0;
            if (!DA.GetData(8, ref day)) { return; }
            int hour = 0;
            if (!DA.GetData(9, ref hour)) { return; }

            double snow_threshold = 1.0;
            if (!DA.GetData(10, ref snow_threshold)) { snow_threshold = 1.0; }
            double tilt_threshold = 60.0;
            if (!DA.GetData(11, ref tilt_threshold)) { tilt_threshold = 60.0; }


            List<double> DNI = new List<double>();
            if (!DA.GetDataList(13, DNI)) { return; }
            List<double> DHI = new List<double>();
            if (!DA.GetDataList(14, DHI)) { return; }
            List<double> SNOW = new List<double>();
            if (!DA.GetDataList(15, SNOW)) { for (int t = 0; t < 8760; t++) SNOW.Add(0.0); }
            List<double> groundalbedo = new List<double>();
            if (!DA.GetDataList(16, groundalbedo)) { for (int t = 0; t < 8760; t++) groundalbedo.Add(0.2); }
            List<double> solarAzimuth = new List<double>();
            DA.GetDataList(17, solarAzimuth);
            List<double> solarAltitude = new List<double>();
            DA.GetDataList(18, solarAltitude);


            int MainSkyRes = 0;
            if (!DA.GetData(20, ref MainSkyRes)) { MainSkyRes = 1; }
            int SpecBounces = 1;
            if (!DA.GetData(21, ref SpecBounces)) { SpecBounces = 1; }
            SpecBounces = (SpecBounces < 0) ? 0 : SpecBounces;
            int DiffIReflSkyRes = 0;
            if (!DA.GetData(22, ref DiffIReflSkyRes)) { DiffIReflSkyRes = 0; }
            int DiffIReflSkyRes2nd = 0;
            if (!DA.GetData(23, ref DiffIReflSkyRes2nd)) { DiffIReflSkyRes2nd = 0; }
            
            CResultsInterreflections ResultsIreflIn = null;
            DA.GetData(25, ref ResultsIreflIn);

            bool mt = false;
            if (!DA.GetData(27, ref mt)) { mt = false; }

            /////////////////////////////////////////////////////////////////////
            //___________________________________________________________________
            #endregion




            #region SIMULATE
            //______________________________________________________________________________________________
            ////////////////////////////////        S I M U L A T E      ///////////////////////////////////
            CCalculateSolarMesh calc = new CCalculateSolarMesh(
                mshobj, objObst, treeObst, latitude, longitude, DNI, DHI, SNOW, groundalbedo, snow_threshold, tilt_threshold,
                year, null, mt, solarAzimuth, solarAltitude);
            calc.RunHourSimulation(month, day, hour, MainSkyRes, SpecBounces, DiffIReflSkyRes, DiffIReflSkyRes2nd, mt);
            Line ln = calc.getSolarVec();
            CResults results = calc.getResults();
            CResultsInterreflections resultsIreflOut = calc.getResultsInterreflections();
            ////////////////////////////////////////////////////////////////////////////////////////////////
            //______________________________________________________________________________________________
            #endregion



            #region OUTPUT
            //______________________________________________________________________________________________
            ////////////////////////////////          O U T P U T        ///////////////////////////////////
            DA.SetData(0, ln);
            DA.SetData(1, results);
            DA.SetData(2, resultsIreflOut);
            ////////////////////////////////////////////////////////////////////////////////////////////////
            //______________________________________________________________________________________________
            #endregion
        }

        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                // You can add image files to your project resources and access them like this:
                //return Resources.IconForThisComponent;
                return GHSolar.Properties.Resources.pic_irrad_h;
            }
        }

        public override Guid ComponentGuid
        {
            get { return new Guid("{eca57c0f-954f-4db5-801a-68461e81b21c}"); }
        }
    }






}
