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
            pManager.AddGenericParameter("TreeObject", "TreeObj", "Input tree objects (generic).", GH_ParamAccess.list);    //each tree object is : (i) a mesh obstacle, (ii) a 8760-timeseries of 0-1 fractions, indicating leave-coverage 1 is full of leaves=full obstruction.
            pManager[2].Optional = true;

            pManager.AddNumberParameter("Latitude", "Latitude", "Latitude of the location in [°].", GH_ParamAccess.item);
            pManager.AddNumberParameter("Longitude", "Longitude", "Longitude of the location in [°].", GH_ParamAccess.item);

            pManager.AddNumberParameter("DNI", "DNI", "Direct normal irradiation 8760-time series.", GH_ParamAccess.list);
            pManager.AddNumberParameter("DHI", "DHI", "Diffuse horizontal irradiation 8760-time series.", GH_ParamAccess.list);
            pManager.AddNumberParameter("Snow", "Snow", "Snow coverage 8760-time series. Used for coverage of horizontal analysis surfaces ~ -/+ 45°. More inclined surfaces are not affected.", GH_ParamAccess.list);
            pManager[7].Optional = true;

            pManager.AddIntegerParameter("Year", "Year", "Year", GH_ParamAccess.item);
            pManager.AddIntegerParameter("Month", "Month", "Month ∈ [1, 12]", GH_ParamAccess.item);
            pManager.AddIntegerParameter("Day", "Day", "Day", GH_ParamAccess.item);
            pManager.AddIntegerParameter("Hour", "Hour", "Hour ∈ [0, 23]", GH_ParamAccess.item);

            ////don´t know if I need these 2
            //pManager.AddBooleanParameter("Simplified Direct", "SimplDir", "Simplified shading mask for direct radiation? Shading mask is then only caculated for the mesh face center of the analysis surface, instead of for each mesh vertex.", GH_ParamAccess.item);
            //pManager[12].Optional = true;
            //pManager.AddBooleanParameter("Simplified Diffuse", "SimplDiff", "Simplified shading mask for diffuse radiation? Shading mask is then only caculated for the mesh face center of the analysis surface, instead of for each mesh vertex.", GH_ParamAccess.item);
            //pManager[13].Optional = true;


            pManager.AddIntegerParameter("Spec. bounces", "Spec.bounces", "Number of specular bounces for inter-reflections. 0 (min) - 2 (max).", GH_ParamAccess.item);
            pManager[12].Optional = true;
            pManager.AddIntegerParameter("Skydome Resolution", "SkyRes", "Sykdome resolution for diffuse shading mask. I.e. recursion level of the icosahedron hemisphere. 0: 10 rays; 1: 29 rays; 2: 97 rays; 3: 353 rays.", GH_ParamAccess.item);
            pManager[13].Optional = true;
            pManager.AddIntegerParameter("Diff. interrefl. res.", "DiffReflRes", "Hemisphere resolution for diffuse interreflections. I.e. recursion level of the icosahedron hemisphere. 0: 10 rays; 1: 29 rays; 2: 97 rays; 3: 353 rays. For deactivating diffuse interreflection, input -1.", GH_ParamAccess.item);
            pManager[14].Optional = true;

            pManager.AddGenericParameter("Precalc Irefl", "Precalc.Irefl", "Input precalculated interreflected irradiation", GH_ParamAccess.item);
            pManager[15].Optional = true;

            pManager.AddBooleanParameter("MT", "MT", "Multi threading?", GH_ParamAccess.item);
            pManager[16].Optional = true;
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddLineParameter("Vec", "Vec", "Solar vector of current hour.", GH_ParamAccess.item);
            pManager.AddGenericParameter("Results", "Results", "Results data of solar irradiation calculation", GH_ParamAccess.item);
            pManager.AddGenericParameter("Irefl", "Irefl", "Irradiation from inter-reflections only.", GH_ParamAccess.item);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {       
            Mesh msh = new Mesh();
            cObstacleObject mshobj = null;
            if (!DA.GetData(0, ref mshobj)) { return; }
            msh = mshobj.mesh;

            //should containt analysis surface itself
            List<cObstacleObject> objObst = new List<cObstacleObject>();
            if (!DA.GetDataList(1, objObst)) { return; }
            Mesh[] obst = new Mesh[objObst.Count];


            double latitude = 0.0;
            if (!DA.GetData(3, ref latitude)) { return; }
            double longitude = 0.0;
            if (!DA.GetData(4, ref longitude)) { return; }

            List<double> DNI = new List<double>();
            if (!DA.GetDataList(5, DNI)) { return; }
            List<double> DHI = new List<double>();
            if (!DA.GetDataList(6, DHI)) { return; }
            List<double> SNOW = new List<double>();

             

            int year = 0;
            if (!DA.GetData(8, ref year)) { return; }
            int month = 0;
            if (!DA.GetData(9, ref month)) { return; }
            int day = 0;
            if (!DA.GetData(10, ref day)) { return; }
            int hour = 0;
            if (!DA.GetData(11, ref hour)) { return; }
            int bounces = 1;
            if (!DA.GetData(12, ref bounces)) { bounces = 1; }
            bounces = (bounces < 0) ? 0 : bounces;  

            int rec = 0;
            if (!DA.GetData(13, ref rec)) { rec = 1; }
            int diffRes = -1;
            if (!DA.GetData(14, ref diffRes)) { diffRes = -1; }

            cResultsInterreflections ResultsIreflIn = null;
            DA.GetData(15, ref ResultsIreflIn);

            bool mt = false;
            if (!DA.GetData(16, ref mt)) { mt = false; }
            



            cCalculateSolarMesh calc = new cCalculateSolarMesh(
                mshobj, objObst, null, latitude, longitude, DNI, DHI, SNOW, 1.0, 30.0, 
                year, month, day, hour, bounces, rec, diffRes, null, mt);
            calc.RunHourSimulation();
            Line ln = calc.getSolarVec();
            cResults results = calc.getResults();
            cResultsInterreflections resultsIreflOut = calc.getResultsInterreflections();


            DA.SetData(0, ln);
            DA.SetData(1, results);
            DA.SetData(2, resultsIreflOut);
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
