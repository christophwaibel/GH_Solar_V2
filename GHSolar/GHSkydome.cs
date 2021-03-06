﻿using System;
using System.Collections.Generic;
using Grasshopper.Kernel;
using Rhino.Geometry;
using SolarModel;
using System.Drawing;
using System.Threading.Tasks;
using System.Xml.Serialization;
using Grasshopper.GUI.IconEditor;
using Rhino.DocObjects;

/*
 * GHSkydome.cs
 * Copyright 2019 Christoph Waibel <waibel@arch.ethz.ch>
 * 
 * This work is licensed under the GNU GPL license version 3.
*/

namespace GHSolar
{
    public class GHSkydome : GH_Component
    {
        public GHSkydome()
            : base("GHSkydome", "Skydome",
                "Skydome component that shows the viewfactors of a sensor point in various resolutions, can show a sunpath diagram, and a cumulative sky matrix (given a weather file).",
                "EnergyHubs", "Solar Simulation")
        {
        }


        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            // 0
            pManager.AddIntegerParameter("recursion level", "Rec.", "Recursion level. 0: 10 vertices, 1: 29 vert., 2: 97 vert., 3: 353 vert., ...", GH_ParamAccess.item);
            // 1
            pManager.AddBooleanParameter("Viewfactors", "ViewFac?", "Showing the viewfactors on the sky dome. I.e. how much is the sensor point obstructed?", GH_ParamAccess.item);
            // 2
            pManager.AddBooleanParameter("CumSkyMatrix", "SkyMat?", "Showing the cumulative sky matrix. Requires a weather file.", GH_ParamAccess.item);
            // 3
            pManager.AddBooleanParameter("SunPath", "SunPath?", "Showing the annual sun path diagram.", GH_ParamAccess.item);
            // 4
            pManager.AddBooleanParameter("SolarVec", "SolarVec?", "Showing hourly solar vector?", GH_ParamAccess.item);
            // 5
            pManager.AddIntegerParameter("HOY", "hoy", "Hour of the year. if empty, no hourly vector will be drawn. Needs a list, but could contain only one hour.", GH_ParamAccess.list);
            // 6
            pManager.AddNumberParameter("location", "location", "Location, 2 numbers. 0: Longitude, 1: Latitude.", GH_ParamAccess.list);
            // 7
            pManager.AddNumberParameter("DNI", "DNI", "DNI. 8760 time series.", GH_ParamAccess.list);
            // 8
            pManager.AddNumberParameter("DHI", "DHI", "DHI. 8760 time series.", GH_ParamAccess.list);
            // 9
            pManager.AddNumberParameter("size", "size", "Size of the skydome. By default it is the bounding box of all the context.", GH_ParamAccess.item);
            // 10
            pManager.AddMeshParameter("context", "context", "Context, i.e. adjacent obstacles", GH_ParamAccess.list);
            // 11
            pManager.AddPointParameter("sp", "sp", "Sensor Point, around which a skydome will be constructed.", GH_ParamAccess.item);
            //12
            pManager.AddIntegerParameter("timezone", "timezone", "timezone", GH_ParamAccess.item, 0);

            int[] ilist = new int[8] { 1, 2, 3, 4, 7, 8, 9, 10 };
            foreach (int i in ilist)
            {
                pManager[i].Optional = true;
            }
        }


        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddMeshParameter("mesh", "mesh", "mesh", GH_ParamAccess.list);
            pManager.AddLineParameter("vectors", "vectors", "Hourly solar vectors.", GH_ParamAccess.list);
            pManager.AddCurveParameter("sunpaths", "sunpaths", "sunpaths", GH_ParamAccess.list);
            pManager.AddBrepParameter("sun geo", "sun geo", "sun geo", GH_ParamAccess.item);
        }


        private readonly List<Line> _solar_vectors = new List<Line>();
        private List<PolylineCurve> _sun_paths = new List<PolylineCurve>();
        private List<bool> _night_time = new List<bool>();
        private List<List<Curve>> _txt = new List<List<Curve>>();
        private Rhino.Display.DisplayMaterial _mat;
        //private Mesh _skydomeViewFactors = new Mesh();
        //private List<Mesh> _colouredMesh = new List<Mesh>();
        private Mesh _viewFactors;
        private List<int> _vectorsObstructed = new List<int>();


        protected override void BeforeSolveInstance()
        {
            _solar_vectors.Clear();
            _night_time.Clear();
            _sun_paths.Clear();
            _txt.Clear();
            //_colouredMesh.Clear();
            _vectorsObstructed.Clear();
            _viewFactors = new Mesh();
        }
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            int year = 2017;
            //int tasks = 1;
            //if (this.mt) tasks = Environment.ProcessorCount;
            int tasks = Environment.ProcessorCount;
            ParallelOptions paropts = new ParallelOptions { MaxDegreeOfParallelism = tasks };
            //ParallelOptions paropts_1cpu = new ParallelOptions { MaxDegreeOfParallelism = 1 };


            //////////////////////////////////////////////////////////////////////////////////////////
            /// INPUTS
            int reclvl = 0;
            if (!DA.GetData(0, ref reclvl)) return;

            bool drawviewfactors = false;
            if (!DA.GetData(1, ref drawviewfactors)) drawviewfactors = false;

            bool drawcumskymatrix = false;
            if (!DA.GetData(2, ref drawcumskymatrix)) drawcumskymatrix = false;

            bool draw_sunpath = true;
            if (!DA.GetData(3, ref draw_sunpath)) draw_sunpath = true;

            bool draw_solarvec = true;
            if (!DA.GetData(4, ref draw_solarvec)) draw_solarvec = true;

            List<int> hoy = new List<int>();
            if (!DA.GetDataList(5, hoy)) return;

            List<double> loc = new List<double>();
            if (!DA.GetDataList(6, loc)) return;
            double longitude = loc[0];
            double latitude = loc[1];

            List<double> dni = new List<double>();
            List<double> dhi = new List<double>();
            DA.GetDataList(7, dni);
            DA.GetDataList(8, dhi);

            double domesize = 1.2;
            if (!DA.GetData(9, ref domesize)) domesize = 1.2;

            List<Mesh> context = new List<Mesh>();
            DA.GetDataList(10, context);

            Point3d sp = new Point3d();
            if (!DA.GetData(11, ref sp)) return;

            int timezone = 0;
            DA.GetData(12, ref timezone);


            //////////////////////////////////////////////////////////////////////////////////////////
            /// size of skydome 
            Point3d anchor = sp;
            Point3d bb_furthest_point;    // max distance to furthest corner of context
            double bb_max_distance = double.MinValue;
            if (context.Count > 0)
            {

                Mesh context_joined = new Mesh();
                foreach (Mesh msh in context)
                    context_joined.Append(msh);
                BoundingBox bb_context = context_joined.GetBoundingBox(false);
                if (bb_context.IsDegenerate(-1.0) == 4)
                    bb_furthest_point = new Point3d(sp.X + 1.0, sp.Y + 1.0, sp.Z + 1.0);
                else
                {
                    Point3d[] _pts = bb_context.GetCorners();
                    int _p_index = 0;
                    for (int i = 0; i < _pts.Length; i++)
                    {
                        double _d = sp.DistanceTo(_pts[i]);
                        if (_d > bb_max_distance)
                        {
                            bb_max_distance = _d;
                            _p_index = i;
                        }
                    }
                    bb_furthest_point = _pts[_p_index];
                }
            }
            else
            {
                bb_furthest_point = new Point3d(sp.X + 1.0, sp.Y + 1.0, sp.Z + 1.0);
            }
            Vector3d vec_sp = bb_furthest_point - sp;
            vec_sp = Vector3d.Multiply(vec_sp, domesize);
            double vec_sp_len = vec_sp.Length;


            //////////////////////////////////////////////////////////////////////////////////////////
            /// SKYDOME
            /// View factors and/or Cumulative SkyMatrix
            SkyDome dome = new SkyDome(reclvl);

            // create 2 meshes, one with obstructed views (make it black transparent), and one unobstructed (regular GH_Mesh color)
            Mesh meshObstructed = new Mesh();

            foreach (double[] p in dome.VertexVectorsSphere)
            {
                Vector3d vec = new Vector3d(p[0], p[1], p[2]);
                vec = Vector3d.Multiply(vec_sp_len, vec);
                meshObstructed.Vertices.Add(vec + sp);
            }
            foreach (int[] f in dome.Faces)
                meshObstructed.Faces.AddFace(f[0], f[1], f[2]);
            meshObstructed.UnifyNormals();

            if (drawviewfactors)
            {

                List<Vector3d> vec_sky_list = new List<Vector3d>();
                List<int> vec_int = new List<int>();
                for (int i = 0; i < meshObstructed.Vertices.Count; i++)
                {
                    Vector3d testvec = meshObstructed.Vertices[i] - sp;
                    if (testvec.Z >= 0.0)
                    {
                        vec_sky_list.Add(testvec);
                        vec_int.Add(i);
                    }
                }
                Color[] colors = new Color[meshObstructed.Vertices.Count];
                for (int i = 0; i < meshObstructed.Vertices.Count; i++)
                {
                    colors[i] = Color.FromArgb(100, 255, 255, 255);  //alpha not working
                }
                meshObstructed.VertexColors.SetColors(colors);
                Vector3d[] vec_sky = vec_sky_list.ToArray();
                bool[] shadow = new bool[vec_sky_list.Count];
                if (context.Count > 0) CShadow.CalcShadowMT(sp, new Vector3d(0, 0, 1), 0.001, vec_sky, context.ToArray(), ref shadow, paropts);

                int j = 0;
                foreach (int i in vec_int)
                {
                    Color c = new Color();
                    if (shadow[j])
                    {
                        // Custom material, DisplayMaterial (rhinostyle) rendering material. and make in override DrawViewportMesh
                        c = Color.FromArgb(100, 0, 0, 0);   //alpha not working
                        meshObstructed.VertexColors.SetColor(i, c);
                        _vectorsObstructed.Add(i);
                    }
                    j++;
                }

            }
            else if (drawcumskymatrix)
            {
                // https://www.sciencedirect.com/science/article/pii/S0038092X04001161
                // http://alexandria.tue.nl/openaccess/635611/p1153final.pdf
                // Solarmodel.dll needs new function to compute cumulative sky view matrix (requires obstruction check from drawviewfactors
                // 1. calc perez diffuse for each hour. use that value (hor, circum, dome) and assign it to each mesh face
                // 2. DNI is computed directly onto SP
                // 3. visualize colored dome for diff only.
                // 4. add text to sensorpoint, stating annual irradiation (DNI plus diff)
                //
                // cumskymatrix seperate component!! coz it can be re-used for several sensor points
                // matrix inversion as in robinson stone to compute irradiation on all sensor points with refl.?
                //
                // needs a separate component that uses cumskymatrix on a number of SPs and visualizes that analysis surface. 
                //... or use this component, output the sensorpoints, give it to a surface and make surface evaluate with the points, and recolor that surface

                if (dni.Count == 8760 && dhi.Count == 8760)     // only continue, if solar irradiance time series are provided
                {
                    //Rhino.RhinoApp.WriteLine("Leggo!");
                    //Context.cWeatherdata weather;
                    //weather.DHI = dhi;
                    //weather.DNI = dni;
                    //weather.Snow = new List<double>();
                    //double[] beta = new double[1] { beta_in };
                    //double[] psi = new double[1] { psi_in };
                    //Sensorpoints.p3d[] coord = new Sensorpoints.p3d[1];   //dummy variables. will not be used in this simplified simulation
                    //coord[0].X = 0;
                    //coord[0].Y = 0;
                    //coord[0].Z = 0;
                    //Sensorpoints.v3d[] normal = new Sensorpoints.v3d[1];   //dummy variables. will not be used in this simplified simulation
                    //normal[0].X = 0;
                    //normal[0].Y = 1;
                    //normal[0].Z = 0;


                    //double[] albedo = new double[8760];
                    //for (int t = 0; t < 8760; t++)
                    //{
                    //    albedo[t] = albedo1;
                    //}

                    //Console.WriteLine("Calculating irradiation...");
                    //Sensorpoints p = new Sensorpoints(beta, psi, coord, normal, reclvl);
                    //p.SetSimpleSkyMT(beta, paropts);
                    //p.SetSimpleGroundReflectionMT(beta, albedo, weather, sunvectors.ToArray(), paropts);
                    //p.CalcIrradiationMT(weather, sunvectors.ToArray(), paropts);

                    // hold on... i need a new function in SolarModel for CumSkyMatrix

                }
                else
                {
                    Rhino.RhinoApp.WriteLine("No data for Direct Normal Irradiance and Diffuse Horizontal Irradiance provided... Please provide 8760 time series for each.");
                }
            }

            if (drawviewfactors || drawcumskymatrix)
            {
                //_colouredMesh.Add(meshObstructed);
                _viewFactors = meshObstructed;
            }

            //////////////////////////////////////////////////////////////////////////////////////////
            /// Solar Vectors
            List<Sphere> spheres = new List<Sphere>();
            double fontsize = vec_sp_len / 50.0;
            List<SunVector> sunvectors_list;
            SunVector.Create8760SunVectors(out sunvectors_list, longitude, latitude, year);
            
            //shifting list of sunvectors according to timezone, so it matches weather file data
            SunVector.ShiftSunVectorsByTimezone(ref sunvectors_list, timezone);

            int count = 0;
            if (draw_solarvec)
            {
                foreach (int h in hoy)
                {
                    Vector3d vec = new Vector3d(sunvectors_list[h].udtCoordXYZ.x, sunvectors_list[h].udtCoordXYZ.y, sunvectors_list[h].udtCoordXYZ.z);
                    vec = Vector3d.Multiply(vec_sp_len, vec);
                    Point3d solarpoint = new Point3d(Point3d.Add(sp, vec));
                    Line ln = new Line(sp, solarpoint);
                    ln.Flip();
                    _solar_vectors.Add(ln);
                    if (sunvectors_list[h].udtCoordXYZ.z < 0) _night_time.Add(true);
                    else _night_time.Add(false);

                    int year_now = sunvectors_list[h].udtTime.iYear;
                    int month_now = sunvectors_list[h].udtTime.iMonth;
                    int day_now = sunvectors_list[h].udtTime.iDay;
                    double hour_now = sunvectors_list[h].udtTime.dHours;
                    string hour_now2 = Convert.ToString(hour_now);
                    if (hour_now < 10) hour_now2 = "0" + Convert.ToString(hour_now);
                    string strval = Convert.ToString(year_now) + "/ " + Convert.ToString(month_now) + "/ " + Convert.ToString(day_now) + "/ " + hour_now2 + ":00";
                    Plane pl = new Plane(ln.From, new Vector3d(-1, 0, 0));
                    //Plane pl = new Plane(ln.From, vec);
                    var te = Rhino.RhinoDoc.ActiveDoc.Objects.AddText(strval, pl, fontsize, "Baskerville", false, false);
                    Rhino.DocObjects.TextObject txt = Rhino.RhinoDoc.ActiveDoc.Objects.Find(te) as Rhino.DocObjects.TextObject;
                    _txt.Add(new List<Curve>());
                    if (txt != null)
                    {
                        var tt = txt.Geometry as Rhino.Geometry.TextEntity;
                        Curve[] A = tt.Explode();

                        foreach (Curve crv in A)
                            _txt[count].Add(crv);
                    }
                    count++;
                    Rhino.RhinoDoc.ActiveDoc.Objects.Delete(te, true);

                    Sphere sph = new Sphere(ln.From, vec_sp_len / 30.0);
                    spheres.Add(sph);
                }
            }


            //////////////////////////////////////////////////////////////////////////////////////////
            /// SUN PATH
            /// !!! wierd sun paths at extreme longitudes -> time shift... +/- UCT
            // draw solar paths: curves that connect each month, but for the same hour
            if (draw_sunpath)
            {
                for (int hod = 0; hod < 24; hod++)
                {
                    List<Point3d> pts = new List<Point3d>();
                    for (int d = 0; d < 365; d++)
                    {
                        int h = hod + 24 * d;
                        Vector3d vec = new Vector3d(sunvectors_list[h].udtCoordXYZ.x, sunvectors_list[h].udtCoordXYZ.y, sunvectors_list[h].udtCoordXYZ.z);
                        vec = Vector3d.Multiply(vec_sp_len, vec);
                        if (vec.Z > 0)
                        {
                            Point3d solarpoint = new Point3d(Point3d.Add(sp, vec));
                            pts.Add(solarpoint);
                        }
                    }
                    if (pts.Count > 0)
                    {
                        PolylineCurve crv = new PolylineCurve(pts);
                        _sun_paths.Add(crv);
                    }
                }

                // draw solar paths; curves that connects each hour, but for the same month
                int interv = 365 / 12;
                for (int m = 0; m < 12; m++)
                {
                    List<Point3d> pts = new List<Point3d>();
                    for (int hod = 0; hod < 24; hod++)
                    {
                        int h = hod + ((m * interv + interv / 2) * 24);
                        Vector3d vec = new Vector3d(sunvectors_list[h].udtCoordXYZ.x, sunvectors_list[h].udtCoordXYZ.y, sunvectors_list[h].udtCoordXYZ.z);
                        vec = Vector3d.Multiply(vec_sp_len, vec);
                        if (vec.Z > 0)
                        {
                            Point3d solarpoint = new Point3d(Point3d.Add(sp, vec));
                            pts.Add(solarpoint);
                        }
                    }
                    if (pts.Count > 0)
                    {
                        PolylineCurve crv = new PolylineCurve(pts);
                        _sun_paths.Add(crv);
                    }
                }
            }


            //////////////////////////////////////////////////////////////////////////////////////////
            /// OUTPUT
            DA.SetData(0, _viewFactors);        // this mesh needs to be colored according to view factor or cumulative sky matrix
            DA.SetDataList(1, _solar_vectors);
            DA.SetDataList(2, _sun_paths);
            DA.SetDataList(3, spheres);
        }

        public override void DrawViewportMeshes(IGH_PreviewArgs args)
        {
            //base.DrawViewportMeshes(args);

            //could be a global variable, if needed to be changed by some smiulation
            Material material = new Material();
            material.Transparency = 0.3;
            _mat = new Rhino.Display.DisplayMaterial(material);

            if (_viewFactors != null)
            {
                args.Display.DrawMeshShaded(_viewFactors, _mat);
                //foreach (Mesh msh in _colouredMesh)
                //{
                //    args.Display.DrawMeshShaded(msh, _mat);

                //    //args.Display.DrawMeshFalseColors(msh);
                //}
            }
        }


        public override void DrawViewportWires(IGH_PreviewArgs args)
        {
            //base.DrawViewportWires(args);

            Color col1 = Color.FromArgb(255, 223, 0);
            Color col2 = Color.FromArgb(120, 120, 120);

            for (int i = 0; i < _solar_vectors.Count; i++)
            {
                Color c = col1;
                if (_night_time[i]) c = col2;
                args.Display.DrawLine(_solar_vectors[i], c, 2);
                args.Display.DrawArrow(_solar_vectors[i], c);
                foreach (Curve crv in _txt[i])
                {
                    args.Display.DrawCurve(crv, c, 1);
                }
            }

            foreach (PolylineCurve crv in _sun_paths)
            {
                args.Display.DrawCurve(crv, col2, 1);
            }


            if (_viewFactors != null)
            {
                foreach (int i in _vectorsObstructed)
                {
                    Point3d point = _viewFactors.Vertices[i];
                    if (args.Display.Viewport.CameraLocation.DistanceTo(point) < 20)
                        args.Display.DrawPoint(point, Rhino.Display.PointStyle.ActivePoint, (float)5, Color.DarkRed);
                }
            }
        }


        public override GH_Exposure Exposure => GH_Exposure.tertiary;


        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                //You can add image files to your project resources and access them like this:
                // return Resources.IconForThisComponent;
                return GHSolar.Properties.Resources.shadingMask;
            }
        }


        public override Guid ComponentGuid
        {
            get { return new Guid("{8d6704b8-515b-4231-905f-c33b7af7316b}"); }
        }
    }
}