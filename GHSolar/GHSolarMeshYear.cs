using System;
using System.Collections.Generic;
using Grasshopper.Kernel;
using Rhino.Geometry;

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
        /// <summary>
        /// Initializes a new instance of the GHSolarMeshYear class.
        /// </summary>
        public GHSolarMeshYear()
            : base("IrradiationMeshYear", "IrrMshY",
                "Calculates hourly solar irradiation on a mesh for an entire year.",
                "EnergyHubs", "Solar Simulation")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddMeshParameter("AnalysisMesh", "AnalysisMsh", "Input analysis meshes for solar irradiation calculation.", GH_ParamAccess.item);
            pManager.AddGenericParameter("ObstacleObjects", "ObstclObjs", "Input obstacle objects (generic).", GH_ParamAccess.list); //each obstacle consists of: (i) a mesh obstacle, (ii) a 8760-timeseries of 0-1 fraction indicating albedo, (iii) a boolean 0-1 for either specular (0) or diffuse (1) reflection.
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

            pManager.AddIntegerParameter("Bounces", "Bounces", "Number of bounces for inter-reflections. 0 (min) - 2 (max).", GH_ParamAccess.item);
            pManager[9].Optional = true;
            pManager.AddIntegerParameter("Skydome Resolution", "SkyRes", "Sykdome resolution for diffuse shading mask. I.e. recursion level of the icosahedron hemisphere. 0: 10 rays; 1: 29 rays; 2: 97 rays; 3: 353 rays.", GH_ParamAccess.item);
            pManager[10].Optional = true;
            pManager.AddIntegerParameter("Interreflection Resolution", "RefllllRes", "Hemisphere resolution for interreflections. I.e. recursion level of the icosahedron hemisphere. 0: 10 rays; 1: 29 rays; 2: 97 rays; 3: 353 rays.", GH_ParamAccess.item);
            pManager[11].Optional = true;
            pManager.AddIntegerParameter("Interpolation mode", "Interp", "Interpolation mode. Options: 0 = 3 days interpolation, 1 = 12 days interpolation (slower, but more accurate).", GH_ParamAccess.item);
            pManager[12].Optional = true;

            pManager.AddBooleanParameter("MT", "MT", "Multi threading?", GH_ParamAccess.item);
            pManager[13].Optional = true;
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            //pManager.AddNumberParameter("I annual", "I annual", "Annual total irradiation in [kWh/a]", GH_ParamAccess.list);
            //pManager.AddNumberParameter("Ib annual", "Ib annual", "Annual beam irradiation in [kWh/a]", GH_ParamAccess.list);
            //pManager.AddNumberParameter("Ih annual", "Ih annual", "Annual diffuse irradiation in [kWh/h]", GH_ParamAccess.list);

            //pManager.AddGenericParameter("I hourly", "I hourly", "Hourly total irradiation in [Wh]", GH_ParamAccess.item);
            //pManager.AddGenericParameter("Ib hourly", "Ib hourly", "Hourly beam irradiation in [Wh]", GH_ParamAccess.item);
            //pManager.AddGenericParameter("Ih hourly", "Ih hourly", "Hourly diffuse irradiation in [Wh]", GH_ParamAccess.item);

            pManager.AddGenericParameter("Results", "Results", "Results data of solar irradiation calculation", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            Mesh msh = new Mesh();
            if (!DA.GetData(0, ref msh)) { return; }

            List<ObstacleObject> objObst = new List<ObstacleObject>();
            if (!DA.GetDataList(1, objObst)) { return; }
            Mesh[] obst = new Mesh[objObst.Count];
            for (int i = 0; i < objObst.Count; i++)
            {
                obst[i] = objObst[i].mesh;
            }

            double latitude = 0.0;
            if (!DA.GetData(3, ref latitude)) { return; }
            double longitude = 0.0;
            if (!DA.GetData(4, ref longitude)) { return; }

            List<double> DNI = new List<double>();
            if (!DA.GetDataList(5, DNI)) { return; }
            List<double> DHI = new List<double>();
            if (!DA.GetDataList(6, DHI)) { return; }
            List<double> SNOW = new List<double>();
            if (!DA.GetDataList(7, SNOW))
                for (int i = 0; i < 8760; i++)
                    SNOW.Add(0.0);

            int year = 0;
            if (!DA.GetData(8, ref year)) { return; }


            int rec = 0;
            if (!DA.GetData(10, ref rec)) { rec = 1; }


            int interpmode = 0;
            if (!DA.GetData(12, ref interpmode)) { interpmode = 0; }

            bool mt = false;
            if (!DA.GetData(13, ref mt)) { mt = false; }


            //these two should be inputs
            double snow_threshold = 1;
            double tilt_treshold = 30;


            double rad = Math.PI / 180;

            List<SunVector> sunvectors = new List<SunVector>();
            SunVector.Create8760SunVectors(ref sunvectors, longitude, latitude, year);
            Context.cWeatherdata weather;
            weather.DHI = new List<double>(DHI);
            weather.DNI = new List<double>(DNI);
            weather.Snow = new List<double>(SNOW);

            Context.cLocation location;
            location.dLatitude = latitude;
            location.dLongitude = longitude;
            location.dTgmt = 1;

            List<double> I = new List<double>();
            List<double> Ih = new List<double>();
            List<double> Ib = new List<double>();

            List<Sensorpoint> ps = new List<Sensorpoint>();
            Point3d[] mshvrt = msh.Vertices.ToPoint3dArray();
            Vector3f[] mshvrtnorm = new Vector3f[mshvrt.Length];
            msh.FaceNormals.ComputeFaceNormals();

            double[] arrbeta = new double[mshvrt.Length];
            double[] arrpsi = new double[mshvrt.Length];

            Vector3d betaangle = new Vector3d(0, 0, 1);
            Vector3d psiangle = new Vector3d(0, 1, 0);
            Plane psiplane = new Plane(new Point3d(0, 0, 0), new Vector3d(0, 0, 1));
            for (int i = 0; i < mshvrt.Length; i++)
            {
                mshvrtnorm[i] = msh.Normals[i];
                //sensor point tilt angle (beta) and azimuth (psi)
                double beta = Vector3d.VectorAngle(mshvrtnorm[i], betaangle) / rad;
                double psi = Vector3d.VectorAngle(mshvrtnorm[i], psiangle, psiplane) / rad;
                if (Double.IsNaN(psi) || Double.IsInfinity(psi))
                {
                    psi = 0;
                }

                arrbeta[i] = beta;
                arrpsi[i] = psi;

            }
            Sensorpoints p = new Sensorpoints(year, weather, location, sunvectors, arrbeta, arrpsi, rec);


            List<bool[]> ShdwBeam_equinox = new List<bool[]>();
            List<bool[]> ShdwBeam_summer = new List<bool[]>();
            List<bool[]> ShdwBeam_winter = new List<bool[]>();

            List<bool[][]> ShdwBeam = new List<bool[][]>();
            int[] startDays = new int[12];
            int[] endDays = new int[12];

            List<bool[]> ShdwSky = new List<bool[]>();


            int[] equsol = SunVector.GetEquinoxSolstice(year);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            int HOYequ = (equsol[0] - 1) * 24;
            int HOYsum = (equsol[1] - 1) * 24;
            int HOYwin = (equsol[3] - 1) * 24;

            List<Point3d> coords = new List<Point3d>();
            for (int i = 0; i < mshvrt.Length; i++)
            {
                double status = (100 / Convert.ToDouble(mshvrt.Length)) * Convert.ToDouble(i + 1);
                Rhino.RhinoApp.WriteLine("SOLAR_V2: (1/6) Shadow calculation... " + Convert.ToString(Math.Round(status,2) + "%"));

                Point3d orig = new Point3d(mshvrt[i].X, mshvrt[i].Y, mshvrt[i].Z);
                coords.Add(orig);
                //Rhino.RhinoDoc.ActiveDoc.Objects.AddPoint(orig);

                //sky dome diffuse
                Vector3d[] vec_sky = new Vector3d[p.sky[i].VerticesHemisphere.Count];
                for (int u = 0; u < vec_sky.Length; u++)
                {
                    vec_sky[u] = new Vector3d(
                        p.sky[i].VertexCoordinatesSphere[p.sky[i].VerticesHemisphere[u]][0],
                        p.sky[i].VertexCoordinatesSphere[p.sky[i].VerticesHemisphere[u]][1],
                        p.sky[i].VertexCoordinatesSphere[p.sky[i].VerticesHemisphere[u]][2]);
                }
                bool[] shdw_sky = new bool[p.sky[i].VerticesHemisphere.Count];
                if (mt)
                    cShadow.CalcShadowMT(orig, mshvrtnorm[i], 0.1, vec_sky, obst, ref shdw_sky);
                else 
                    cShadow.CalcShadow(orig, mshvrtnorm[i], 0.1, vec_sky, obst, ref shdw_sky);

                ShdwSky.Add(shdw_sky);

                if (interpmode == 0)
                {
                    //beam for 24 hours. no! don´t use those vectors, where its night time anyway!
                    // equinox:             march 20
                    // summer solstice:     june 21
                    // winter solstice:     december 21
                    Vector3d[] vec_beam_equ = new Vector3d[24];
                    Vector3d[] vec_beam_sum = new Vector3d[24];
                    Vector3d[] vec_beam_win = new Vector3d[24];
                    bool[] sunshine_equ = new bool[24];
                    bool[] sunshine_sum = new bool[24];
                    bool[] sunshine_win = new bool[24];
                    for (int t = 0; t < 24; t++)
                    {
                        if (sunvectors[HOYequ + t].Sunshine)
                            sunshine_equ[t] = true;
                        if (sunvectors[HOYsum + t].Sunshine)
                            sunshine_sum[t] = true;
                        if (sunvectors[HOYwin + t].Sunshine)
                            sunshine_win[t] = true;
                        vec_beam_equ[t] = new Vector3d(sunvectors[HOYequ + t].udtCoordXYZ.x, sunvectors[HOYequ + t].udtCoordXYZ.y, sunvectors[HOYequ + t].udtCoordXYZ.z);
                        vec_beam_sum[t] = new Vector3d(sunvectors[HOYsum + t].udtCoordXYZ.x, sunvectors[HOYsum + t].udtCoordXYZ.y, sunvectors[HOYsum + t].udtCoordXYZ.z);
                        vec_beam_win[t] = new Vector3d(sunvectors[HOYwin + t].udtCoordXYZ.x, sunvectors[HOYwin + t].udtCoordXYZ.y, sunvectors[HOYwin + t].udtCoordXYZ.z);
                    }

                    bool[] shdw_beam_equ = new bool[24];
                    bool[] shdw_beam_sum = new bool[24];
                    bool[] shdw_beam_win = new bool[24];
                    if (mt)
                    {
                        cShadow.CalcShadowMT(orig, mshvrtnorm[i], 0.1, vec_beam_equ, sunshine_equ, obst, ref shdw_beam_equ);
                        cShadow.CalcShadowMT(orig, mshvrtnorm[i], 0.1, vec_beam_sum, sunshine_sum, obst, ref shdw_beam_sum);
                        cShadow.CalcShadowMT(orig, mshvrtnorm[i], 0.1, vec_beam_win, sunshine_win, obst, ref shdw_beam_win);
                    }
                    else
                    {
                        cShadow.CalcShadow(orig, mshvrtnorm[i], 0.1, vec_beam_equ, sunshine_equ, obst, ref shdw_beam_equ);
                        cShadow.CalcShadow(orig, mshvrtnorm[i], 0.1, vec_beam_sum, sunshine_sum, obst, ref shdw_beam_sum);
                        cShadow.CalcShadow(orig, mshvrtnorm[i], 0.1, vec_beam_win, sunshine_win, obst, ref shdw_beam_win);
                    }
                    ShdwBeam_equinox.Add(shdw_beam_equ);
                    ShdwBeam_summer.Add(shdw_beam_sum);
                    ShdwBeam_winter.Add(shdw_beam_win);
                }
                else
                {
                    bool[][] shdw_beam = new bool[12][];
                    int dmcount = 1;    //days in month counter
                    for (int d = 0; d < 12; d++)
                    {
                        Vector3d[] vec_beam = new Vector3d[24];
                        int dm = System.DateTime.DaysInMonth(year, d + 1);
                        startDays[d] = dmcount;
                        endDays[d] = dm + dmcount;
                        dmcount += dm;

                        bool[] sunshine = new bool[24];
                        int HOY = (startDays[d] - 1) * 24;
                        for (int t = 0; t < 24; t++) 
                        {
                            if (sunvectors[HOY + t].Sunshine)
                                sunshine[t] = true;
                            vec_beam[t] = new Vector3d(sunvectors[HOY + t].udtCoordXYZ.x, sunvectors[HOY + t].udtCoordXYZ.y, sunvectors[HOY + t].udtCoordXYZ.z);
                        }
                        shdw_beam[d] = new bool[24];
                        if (mt)
                            cShadow.CalcShadowMT(orig, mshvrtnorm[i], 0.1, vec_beam, sunshine, obst, ref shdw_beam[d]);
                        else
                            cShadow.CalcShadow(orig, mshvrtnorm[i], 0.1, vec_beam, sunshine, obst, ref shdw_beam[d]);
                    }
                    ShdwBeam.Add(shdw_beam);
                }
            }


            Rhino.RhinoApp.WriteLine("SOLAR_V2: (2/6) Setting interpolated shadows...");
            if (mt)
            {
                if (interpmode == 0)
                    p.SetShadowsInterpolatedMT(ShdwBeam_equinox, ShdwBeam_summer, ShdwBeam_winter, ShdwSky);
                else
                    p.SetShadowsInterpolatedMT(startDays, endDays, ShdwBeam, ShdwSky);
            }
            else
            {
                if (interpmode == 0)
                    p.SetShadowsInterpolated(ShdwBeam_equinox, ShdwBeam_summer, ShdwBeam_winter, ShdwSky);
                else
                    p.SetShadowsInterpolated(startDays, endDays, ShdwBeam, ShdwSky);
            }

            Rhino.RhinoApp.WriteLine("SOLAR_V2: (3/6) Setting snow cover...");
            p.SetSnowcover(snow_threshold, tilt_treshold);

            Rhino.RhinoApp.WriteLine("SOLAR_V2: (4/6) Calculating inter-reflections...");
            //p.SetInterreflection();

            Rhino.RhinoApp.WriteLine("SOLAR_V2: (5/6) Calculating irradiation...");
            if (mt)
                p.CalcIrradiationMT();
            else
                p.CalcIrradiation();


            Rhino.RhinoApp.WriteLine("SOLAR_V2: (6/6) Writing data to Rhino...");

            Matrix I_hourly = new Matrix(mshvrt.Length, 8760);
            Matrix Ib_hourly = new Matrix(mshvrt.Length, 8760);
            Matrix Ih_hourly = new Matrix(mshvrt.Length, 8760);
            for (int i = 0; i < mshvrt.Length; i++)
            {
                double Itot = 0;
                double Ibtot = 0;
                double Idtot = 0;
                for (int t = 0; t < 8760; t++)
                {
                    I_hourly[i, t] = p.I[i][t];
                    Ib_hourly[i, t] = p.Ibeam[i][t];
                    Ih_hourly[i, t] = p.Idiff[i][t];
                    Itot += p.I[i][t];
                    Ibtot += p.Ibeam[i][t];
                    Idtot += p.Idiff[i][t];
                }
                I.Add(Itot / 1000);     //in kWh/a
                Ib.Add(Ibtot / 1000);   //in kWh/a
                Ih.Add(Idtot / 1000);   //in kWh/a
            }



            //DA.SetDataList(0, I);
            //DA.SetDataList(1, Ib);
            //DA.SetDataList(2, Ih);

            //DA.SetData(3, I_hourly);
            //DA.SetData(4, Ib_hourly);
            //DA.SetData(5, Ih_hourly);

            cResults results = new cResults(I, Ib, Ih, I_hourly,Ib_hourly,Ih_hourly,coords);
            DA.SetData(0, results);

            Rhino.RhinoApp.WriteLine("SOLAR_V2... Done");
        }

        /// <summary>
        /// Provides an Icon for the component.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                //You can add image files to your project resources and access them like this:
                // return Resources.IconForThisComponent;
                return GHSolar.Properties.Resources.pic_irrad_y;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("{6d5b9658-b84d-47ee-a5a1-1dd080c0caf1}"); }
        }
    }
}