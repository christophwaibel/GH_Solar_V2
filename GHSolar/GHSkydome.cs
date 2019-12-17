using System;
using System.Collections.Generic;
using Grasshopper.Kernel;
using Rhino.Geometry;
using SolarModel;


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


        protected override void SolveInstance(IGH_DataAccess DA)
        {
            int year = 2017;
            int reclvl = 0;
            if (!DA.GetData(0, ref reclvl)) return;

            List<int> hoy = new List<int>();
            if (!DA.GetDataList(4, hoy)) return;

            List<double> loc = new List<double>();
            if (!DA.GetDataList(5, loc)) return;
            double longitude = loc[0];
            double latitude = loc[1];

            //10 for SP
            Point3d sp = new Point3d();
            if (!DA.GetData(10, ref sp)) return;


            //Full Sphere
            /*
            IcoSphere ico = new IcoSphere(reclvl);
            List<int[]> facesico = ico.getFaces();
            List<double[]> vertexcoordsico = ico.getVertexCoordinates();
            Mesh meshico = new Mesh();
            foreach (double[] p in vertexcoordsico)
            {
                meshico.Vertices.Add(p[0], p[1], p[2]);
                points.Add(new Point3d(p[0], p[1], p[2]));
            }
            foreach (int[] f in facesico)
            {
                meshico.Faces.AddFace(f[0], f[1], f[2]);
            }
            meshlist.Add(meshico);
            */

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
            meshlist.Add(mesh);

            DA.SetDataList(0, meshlist);        // this mesh needs to be colored according to view factor or cumulative sky matrix

            List<SunVector> sunvectors_list;
            SunVector.Create8760SunVectors(out sunvectors_list, longitude, latitude, year);

            List<Line> ln = new List<Line>();
            foreach (int h in hoy)
            {
                Vector3d vec = new Vector3d(sunvectors_list[h].udtCoordXYZ.x, sunvectors_list[h].udtCoordXYZ.y, sunvectors_list[h].udtCoordXYZ.z);
                Point3d solarpoint = new Point3d(Point3d.Add(sp, vec));
                ln.Add(new Line(sp, solarpoint));
            }

            //CCalculateSolarMesh calc = new CCalculateSolarMesh(
            //    mshobj, objObst, treeObst, latitude, longitude, DNI, DHI, SNOW, groundalbedo, snow_threshold, tilt_threshold,
            //    year, null, mt, solarAzimuth, solarAltitude);
            //calc.RunHourSimulationMT(month, day, hour, MainSkyRes, SpecBounces, DiffIReflSkyRes, DiffIReflSkyRes2nd);
            //Line[] ln = calc.getSolarVec();
            DA.SetDataList(1, ln);


            List<PolylineCurve> crvs = new List<PolylineCurve>();

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
                    crvs.Add(crv);
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
                    crvs.Add(crv);
                }
            }

            DA.SetDataList(2, crvs);
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