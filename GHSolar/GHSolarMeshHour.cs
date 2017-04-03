using System;
using System.Collections.Generic;
using Grasshopper.Kernel;
using Rhino.Geometry;

using SolarModel;
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
            pManager.AddMeshParameter("AnalysisMesh", "AnalysisMsh", "Input analysis meshes for solar irradiation calculation.", GH_ParamAccess.item);
            pManager.AddGenericParameter("ObstacleObjects", "ObstclObjs", "Input obstacle objects (generic).", GH_ParamAccess.list); //each obstacle consists of: (i) a mesh obstacle, (ii) a 8760-timeseries of 0-1 fraction indicating albedo, (iii) a boolean 0-1 for either specular (0) or diffuse (1) reflection.
            pManager[1].Optional = true;
            pManager.AddGenericParameter("TreeObject", "TreeObj", "Input tree objects (generic).", GH_ParamAccess.list);    //each tree object is : (i) a mesh obstacle, (ii) a 8760-timeseries of 0-1 fractions, indicating leave-coverage 1 is full of leaves=full obstruction.
            pManager[2].Optional = true;

            pManager.AddNumberParameter("φ", "φ", "Latitude of the location in [°].", GH_ParamAccess.item);
            pManager.AddNumberParameter("λ", "λ", "Longitude of the location in [°].", GH_ParamAccess.item);

            pManager.AddNumberParameter("DNI", "DNI", "Direct normal irradiation 8760-time series.", GH_ParamAccess.list);
            pManager.AddNumberParameter("DHI", "DHI", "Diffuse horizontal irradiation 8760-time series.", GH_ParamAccess.list);
            pManager.AddNumberParameter("Snow", "Snow", "Snow coverage 8760-time series. Used for coverage of horizontal analysis surfaces ~ -/+ 45°. More inclined surfaces are not affected.", GH_ParamAccess.list);
            pManager[7].Optional = true;

            pManager.AddIntegerParameter("Year", "Year", "Year", GH_ParamAccess.item);
            pManager.AddIntegerParameter("Month", "Month", "Month ∈ [1, 12]", GH_ParamAccess.item);
            pManager.AddIntegerParameter("Day", "Day", "Day", GH_ParamAccess.item);
            pManager.AddIntegerParameter("Hour", "Hour", "Hour ∈ [0, 23]", GH_ParamAccess.item);

            pManager.AddBooleanParameter("Simplified Direct", "SimplDir", "Simplified shading mask for direct radiation? Shading mask is then only caculated for the mesh face center of the analysis surface, instead of for each mesh vertex.", GH_ParamAccess.item);
            pManager[12].Optional = true;
            pManager.AddBooleanParameter("Simplified Diffuse", "SimplDiff", "Simplified shading mask for diffuse radiation? Shading mask is then only caculated for the mesh face center of the analysis surface, instead of for each mesh vertex.", GH_ParamAccess.item);
            pManager[13].Optional = true;
            pManager.AddIntegerParameter("Bounces", "Bounces", "Number of bounces for inter-reflections. 0 (min) - 2 (max).", GH_ParamAccess.item);
            pManager[14].Optional = true;
            pManager.AddIntegerParameter("Skydome Resolution", "SkyRes", "Sykdome resolution for diffuse shading mask. I.e. recursion level of the icosahedron hemisphere. 0: 10 rays; 1: 29 rays; 2: 97 rays; 3: 353 rays.", GH_ParamAccess.item);
            pManager.AddIntegerParameter("Interreflection Resolution", "RefllllRes", "Hemisphere resolution for interreflections. I.e. recursion level of the icosahedron hemisphere. 0: 10 rays; 1: 29 rays; 2: 97 rays; 3: 353 rays.", GH_ParamAccess.item);
            pManager[16].Optional = true;
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {

            pManager.AddNumberParameter("Total I", "I", "I", GH_ParamAccess.list);
            pManager.AddNumberParameter("Ib", "Ib", "Ib", GH_ParamAccess.list);
            pManager.AddNumberParameter("Ih", "Ih", "Ih", GH_ParamAccess.list);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            Mesh msh = new Mesh();
            if (!DA.GetData(0, ref msh)) { return; }


            //List<ObstacleObject> obstacles = new List<ObstacleObject>();
            //if (!DA.GetDataList(1, obstacles)) { return; }
            ////Mesh mesh = obstacles[0].mesh;


            //List<TreeObject> trees = new List<TreeObject>();
            //if (!DA.GetDataList(2, trees)) { return; }
            ////Mesh treemesh = trees[0].mesh;


            double latitude = 0.0;
            if (!DA.GetData(3, ref latitude)) { return; }
            double longitude = 0.0;
            if (!DA.GetData(4, ref longitude)) { return; }

            List<double> DNI = new List<double>();
            if (!DA.GetDataList(5, DNI)) { return; }
            List<double> DHI = new List<double>();
            if (!DA.GetDataList(6, DHI)) { return; }


            int year = 0;
            if (!DA.GetData(8, ref year)) { return; }
            int month = 0;
            if (!DA.GetData(9, ref month)) { return; }
            int day = 0;
            if (!DA.GetData(10, ref day)) { return; }
            int hour = 0;
            if (!DA.GetData(11, ref hour)) { return; }

            int rec = 0;
            if (!DA.GetData(15, ref rec)) { return; }



            double rad = Math.PI / 180;

            List<SunVector> sunvectors = new List<SunVector>();
            Context.Create8760SunVectors(ref sunvectors, longitude, latitude, year);
            Context.cWeatherdata weather;
            weather.DHI = new List<double>(DHI);
            weather.DNI = new List<double>(DNI);

            Context.cLocation location;
            location.dLatitude = latitude;
            location.dLongitude = longitude;
            location.dTgmt = 1;

            int []daysInMonth = new int[12];
            for (int i=0; i< daysInMonth.Length; i++)
            {
                daysInMonth[i] = System.DateTime.DaysInMonth(year, month);
            }
            int DOY = 0;
            for (int i = 0; i < month - 1; i++)
            {
                DOY += (i + 1) * daysInMonth[i];
                //month * 1;            //1-365
            }
            DOY += day;



            List<double> I = new List<double>();
            List<double> Ih = new List<double>();
            List<double> Ib = new List<double>();
            int HOY = (DOY - 1) * 24 + hour;




            List<Sensorpoint> ps = new List<Sensorpoint>();
            Point3d [] mshvrt = msh.Vertices.ToPoint3dArray();
            Vector3f[] mshvrtnorm = new Vector3f[mshvrt.Length];
            msh.FaceNormals.ComputeFaceNormals();
            for (int i=0; i< mshvrt.Length; i++)
            {
                mshvrtnorm[i] = msh.Normals[i];


                //sensor point tilt angle (beta) and azimuth (psi)
                Vector3d betaangle = new Vector3d(0, 0, 1);
                Vector3d psiangle = new Vector3d(0, -1, 0);
                Plane psiplane = new Plane(new Point3d(0, 0, 0), new Vector3d(0, 0, 1));
                double beta = Vector3d.VectorAngle(mshvrtnorm[i], betaangle) / rad;
                double psi = Vector3d.VectorAngle(mshvrtnorm[i], psiangle, psiplane) / rad;


                Sensorpoint p = new Sensorpoint(year, weather, location, sunvectors,beta, psi, rec);
                p.CalcIrradiation(DOY, hour);
                ps.Add(p);


                I.Add(p.I[HOY]);
                Ih.Add(p.Idiff[HOY][0]);
                Ib.Add(p.Ibeam[HOY]);
            }


            DA.SetDataList(0, I);
            DA.SetDataList(1, Ib);
            DA.SetDataList(2, Ih);


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
