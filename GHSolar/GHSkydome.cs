﻿using System;
using System.Collections.Generic;
using Grasshopper.Kernel;
using Rhino.Geometry;
using SolarModel;
using System.Drawing;

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
            pManager.AddIntegerParameter("recursion level", "Rec.", "Recursion level. 0: 10 vertices, 1: 29 vert., 2: 97 vert., 3: 353 vert., ...", GH_ParamAccess.item);
            pManager.AddBooleanParameter("Viewfactors", "ViewFac?", "Showing the viewfactors on the sky dome. I.e. how much is the sensor point obstructed?", GH_ParamAccess.item);
            pManager.AddBooleanParameter("CumSkyMatrix", "SkyMat?", "Showing the cumulative sky matrix. Requires a weather file.", GH_ParamAccess.item);
            pManager.AddBooleanParameter("SunPath", "SunPath?", "Showing the annual sun path diagram.", GH_ParamAccess.item);
            pManager.AddIntegerParameter("HOY", "hoy", "Hour of the year. if empty, no hourly vector will be drawn. Needs a list, but could contain only one hour.", GH_ParamAccess.list);
            pManager.AddNumberParameter("location", "location", "Location, 2 numbers. 0: Latitude, 1: Longitude.", GH_ParamAccess.list);
            pManager.AddNumberParameter("DNI", "DNI", "DNI. 8760 time series.", GH_ParamAccess.list);
            pManager.AddNumberParameter("DHI", "DHI", "DHI. 8760 time series.", GH_ParamAccess.list);
            pManager.AddNumberParameter("size", "size", "Size of the skydome. By default it is the bounding box of all the context.", GH_ParamAccess.item);
            pManager.AddMeshParameter("context", "context", "Context, i.e. adjacent obstacles", GH_ParamAccess.list);
            pManager.AddPointParameter("sp", "sp", "Sensor Point, around which a skydome will be constructed.", GH_ParamAccess.item);

            int[] ilist = new int[7] { 1, 2, 3, 6, 7, 8, 9 };
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
        }


        private readonly List<Line> _solar_vectors = new List<Line>();
        List<PolylineCurve> _sun_paths = new List<PolylineCurve>();
        List<bool> _night_time = new List<bool>();
        List<List<Curve>> _txt = new List<List<Curve>>();
        protected override void BeforeSolveInstance()
        {
            _solar_vectors.Clear();
            _night_time.Clear();
            _sun_paths.Clear();
            _txt.Clear();
        }
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            int year = 2017;

            //////////////////////////////////////////////////////////////////////////////////////////
            /// INPUTS
            int reclvl = 0;
            if (!DA.GetData(0, ref reclvl)) return;

            bool drawviewfactors = false;
            if (!DA.GetData(1, ref drawviewfactors)) drawviewfactors = false;

            bool drawcumskymatrix = false;
            if (!DA.GetData(2, ref drawcumskymatrix)) drawcumskymatrix = false;

            List<int> hoy = new List<int>();
            if (!DA.GetDataList(4, hoy)) return;

            List<double> loc = new List<double>();
            if (!DA.GetDataList(5, loc)) return;
            double longitude = loc[0];
            double latitude = loc[1];

            Point3d sp = new Point3d();
            if (!DA.GetData(10, ref sp)) return;


            //////////////////////////////////////////////////////////////////////////////////////////
            /// SKYDOME
            SkyDome dome = new SkyDome(reclvl);
            Mesh mesh = new Mesh();
            List<Mesh> meshlist = new List<Mesh>();
            foreach (double[] p in dome.VertexVectorsSphere)
            {
                mesh.Vertices.Add(p[0] + sp.X, p[1] + sp.Y, p[2] + sp.Z);
            }
            foreach (int[] f in dome.Faces)
            {
                mesh.Faces.AddFace(f[0], f[1], f[2]);
            }
            mesh.UnifyNormals();

            // !!!!!!!!!!
            if (drawviewfactors)
            {
                // using GHSolar, only compute obstructions for the skydome. W/o doing all the perez stuff.
                // how?
            }else if (drawcumskymatrix)
            {
                // Solarmodel.dll needs new function to compute cumulative sky view matrix (requires obstruction check from drawviewfactors
            }
            meshlist.Add(mesh);


            //////////////////////////////////////////////////////////////////////////////////////////
            /// Solar Vectors
            BoundingBox bb = mesh.GetBoundingBox(false);
            double fontsize = (bb.Max.X - bb.Min.X) / 50.0;
            List<SunVector> sunvectors_list;
            SunVector.Create8760SunVectors(out sunvectors_list, longitude, latitude, year);
            int count = 0;
            foreach (int h in hoy)
            {
                Vector3d vec = new Vector3d(sunvectors_list[h].udtCoordXYZ.x, sunvectors_list[h].udtCoordXYZ.y, sunvectors_list[h].udtCoordXYZ.z);
                Point3d solarpoint = new Point3d(Point3d.Add(sp, vec));
                Line ln = new Line(sp, solarpoint);
                ln.Flip();
                _solar_vectors.Add(ln);
                if(sunvectors_list[h].udtCoordXYZ.z < 0) _night_time.Add(true);
                else _night_time.Add(false);

                int year_now = sunvectors_list[h].udtTime.iYear;
                int month_now = sunvectors_list[h].udtTime.iMonth;
                int day_now = sunvectors_list[h].udtTime.iDay;
                double hour_now = sunvectors_list[h].udtTime.dHours;
                string strval = Convert.ToString(year_now) + "; " + Convert.ToString(month_now) + "; " + Convert.ToString(day_now) + "; " + Convert.ToString(hour_now) + ":00";
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
                    {
                        _txt[count].Add(crv);
                    }
                }
                count++;
                Rhino.RhinoDoc.ActiveDoc.Objects.Delete(te, true);
            }


            //////////////////////////////////////////////////////////////////////////////////////////
            /// SUN PATH
            /// !!! wierd sun paths at extreme longitudes -> time shift... +/- UCT
            // draw solar paths: curves that connect each month, but for the same hour
            for (int hod = 0; hod < 24; hod++)
            {   
                List<Point3d> pts = new List<Point3d>();
                for (int d = 0; d < 365; d++)
                {
                    int h = hod + 24 * d;
                    Vector3d vec = new Vector3d(sunvectors_list[h].udtCoordXYZ.x, sunvectors_list[h].udtCoordXYZ.y, sunvectors_list[h].udtCoordXYZ.z);
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
            for (int m = 0; m<12; m++)
            {
                List<Point3d> pts = new List<Point3d>();
                for (int hod=0; hod<24; hod++)
                {
                    int h = hod + ((m * interv + interv / 2) * 24);
                    Vector3d vec = new Vector3d(sunvectors_list[h].udtCoordXYZ.x, sunvectors_list[h].udtCoordXYZ.y, sunvectors_list[h].udtCoordXYZ.z);
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


            //////////////////////////////////////////////////////////////////////////////////////////
            /// OUTPUT
            DA.SetDataList(0, meshlist);        // this mesh needs to be colored according to view factor or cumulative sky matrix
            DA.SetDataList(1, _solar_vectors);
            DA.SetDataList(2, _sun_paths);
        }


        public override void DrawViewportWires(IGH_PreviewArgs args)
        {
            Color col1 = Color.FromArgb(255, 223, 0);
            Color col2 = Color.FromArgb(120, 120, 120);
            

            for (int i = 0; i < _solar_vectors.Count; i++)
            {
                Color c = col1;
                if (_night_time[i]) c = col2;
                args.Display.DrawLine(_solar_vectors[i], c, 2);
                args.Display.DrawArrow(_solar_vectors[i], c);
                foreach(Curve crv in _txt[i])
                {
                    args.Display.DrawCurve(crv, c, 1);
                }
            }

            foreach (PolylineCurve crv in _sun_paths)
            {
                args.Display.DrawCurve(crv, col2, 1);
            }
        }


        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                //You can add image files to your project resources and access them like this:
                // return Resources.IconForThisComponent;
                return null;
            }
        }


        public override Guid ComponentGuid
        {
            get { return new Guid("{8d6704b8-515b-4231-905f-c33b7af7316b}"); }
        }
    }
}