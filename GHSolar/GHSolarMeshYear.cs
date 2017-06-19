
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

            pManager.AddNumberParameter("Latitude", "Latitude", "Latitude of the location in [°].", GH_ParamAccess.item);
            pManager.AddNumberParameter("Longitude", "Longitude", "Longitude of the location in [°].", GH_ParamAccess.item);

            pManager.AddNumberParameter("DNI", "DNI", "Direct normal irradiation 8760-time series.", GH_ParamAccess.list);
            pManager.AddNumberParameter("DHI", "DHI", "Diffuse horizontal irradiation 8760-time series.", GH_ParamAccess.list);
            pManager.AddNumberParameter("Snow", "Snow", "Snow coverage 8760-time series. Used for coverage of horizontal analysis surfaces ~ -/+ 45°. More inclined surfaces are not affected.", GH_ParamAccess.list);
            pManager[7].Optional = true;

            pManager.AddIntegerParameter("Year", "Year", "Year", GH_ParamAccess.item);

            pManager.AddIntegerParameter("Spec. bounces", "Spec.bounces", "Number of specular bounces for inter-reflections. 0 (min) - 2 (max).", GH_ParamAccess.item);
            pManager[9].Optional = true;
            pManager.AddIntegerParameter("Skydome Resolution", "SkyRes", "Sykdome resolution for diffuse shading mask. I.e. recursion level of the icosahedron hemisphere. 0: 10 rays; 1: 29 rays; 2: 97 rays; 3: 353 rays.", GH_ParamAccess.item);
            pManager[10].Optional = true;
            pManager.AddIntegerParameter("Interrefl. res.", "InterreflRes", "Precision of interreflections calculations. For values >=0, the value sets the recursion level of the diffuse interreflection hemisphere (0: 10 rays; 1: 29 rays; 2: 97 rays; 3: 353 rays), and 12 days are calculated for specular interreflections and interpolated for all other hours of the year. For -1, 3 days only are used for specular interreflections, and a recursion level of 0 for the skydome of the diffuse interrefelction. For assuming unobstructed diffuse ground reflection only, use: -2.", GH_ParamAccess.item);
            pManager[11].Optional = true;
            pManager.AddIntegerParameter("Interpolation mode", "Interp", "Interpolation mode. Options: 0 = 3 days interpolation, 1 = 12 days interpolation (slower, but more accurate).", GH_ParamAccess.item);
            pManager[12].Optional = true;

            pManager.AddGenericParameter("Precalc Irefl", "Precalc.Irefl", "Input precalculated interreflected irradiation", GH_ParamAccess.item);
            pManager[13].Optional = true;

            pManager.AddBooleanParameter("MT", "MT", "Multi threading?", GH_ParamAccess.item);
            pManager[14].Optional = true;

            pManager.AddNumberParameter("GrndAlb", "GrndAlb", "Ground albedo, 8760 time series", GH_ParamAccess.list);
            pManager[15].Optional = true;
            pManager.AddNumberParameter("SnwThrsh", "SnwThrsh", "Snow thickness threshold, before which snow does not remain on surface", GH_ParamAccess.item);
            pManager[16].Optional = true;
            pManager.AddNumberParameter("TiltThrsh", "TiltThrsh", "Tilst threshold, after which snow does not remain on surface", GH_ParamAccess.item);
            pManager[17].Optional = true;
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
            cObstacleObject mshobj = null;
            if (!DA.GetData(0, ref mshobj)) { return; }
            msh = mshobj.mesh;

            //should containt analysis surface itself
            List<cObstacleObject> objObst = new List<cObstacleObject>();
            if (!DA.GetDataList(1, objObst)) { return; }

            List<cPermObject> treeObst = new List<cPermObject>();
            DA.GetDataList(2, treeObst);

            double latitude = 0.0;
            if (!DA.GetData(3, ref latitude)) { return; }
            double longitude = 0.0;
            if (!DA.GetData(4, ref longitude)) { return; }
            
            List<double> DNI = new List<double>();
            if (!DA.GetDataList(5, DNI)) { return; }
            List<double> DHI = new List<double>();
            if (!DA.GetDataList(6, DHI)) { return; }
            List<double> SNOW = new List<double>();
            if (!DA.GetDataList(7, SNOW)) { return; }

            int year = 0;
            if (!DA.GetData(8, ref year)) { return; }
           
            int bounces = 1;
            if (!DA.GetData(9, ref bounces)) { bounces = 1; }
            bounces = (bounces < 0) ? 0 : bounces;

            int rec = 0;
            if (!DA.GetData(10, ref rec)) { rec = 1; }
            int diffRes = -2;
            if (!DA.GetData(11, ref diffRes)) { diffRes = -2; }

            int interpmode = 0;
            if (!DA.GetData(12, ref interpmode)) { interpmode = 0; }

            cResultsInterreflections ResultsIreflIn = null;
            DA.GetData(13, ref ResultsIreflIn);

            bool mt = false;
            if (!DA.GetData(14, ref mt)) { mt = false; }

            List<double> groundalbedo = new List<double>();
            if (!DA.GetDataList(15, groundalbedo)) { for (int t = 0; t < 8760; t++) groundalbedo.Add(0.3); }

            double snow_threshold = 1.0;
            if (!DA.GetData(16, ref snow_threshold)) { snow_threshold = 1.0; }
            double tilt_threshold = 30.0;
            if (!DA.GetData(17, ref tilt_threshold)) { tilt_threshold = 30.0; }
            /////////////////////////////////////////////////////////////////////
            //___________________________________________________________________




            //___________________________________________________________________
            /////////////////////////////////////////////////////////////////////
            //Status?
            //double status = (100 / Convert.ToDouble(mshvrt.Length)) * Convert.ToDouble(i + 1);
            //Rhino.RhinoApp.WriteLine("SOLAR_V2: (1/6) Shadow calculation... " + Convert.ToString(Math.Round(status,2) + "%"));
            //Rhino.RhinoApp.WriteLine("SOLAR_V2: (2/6) Setting interpolated shadows...");
            //Rhino.RhinoApp.WriteLine("SOLAR_V2: (3/6) Setting snow cover...");
            //Rhino.RhinoApp.WriteLine("SOLAR_V2: (4/6) Calculating inter-reflections...");
            //Rhino.RhinoApp.WriteLine("SOLAR_V2: (5/6) Calculating irradiation...");
            //Rhino.RhinoApp.WriteLine("SOLAR_V2: (6/6) Writing data to Rhino...");
            /////////////////////////////////////////////////////////////////////
            //___________________________________________________________________



            #region SIMULATE
            //______________________________________________________________________________________________
            ////////////////////////////////        S I M U L A T E      ///////////////////////////////////
            cCalculateSolarMesh calc = new cCalculateSolarMesh(
                mshobj, objObst, treeObst, latitude, longitude, DNI, DHI, SNOW, groundalbedo, snow_threshold, tilt_threshold,
                year, bounces, rec, diffRes, null, mt);
            calc.RunAnnualSimulation(0.01);
            cResults results = calc.getResults();
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



            Rhino.RhinoApp.WriteLine("SOLAR_V2... Done");
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
