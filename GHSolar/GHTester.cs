using System;
using System.Collections.Generic;
using Grasshopper.Kernel;
using Rhino.Geometry;
using SolarModel;


/*
 * GHTester.cs
 * Copyright 2017 Christoph Waibel <chwaibel@student.ethz.ch>
 * 
 * This work is licensed under the GNU GPL license version 3.
*/

namespace GHSolar
{
    public class GHTester : GH_Component
    {

        public GHTester()
            : base("GHTester", "Nickname",
                "Description",
                "EnergyHubs", "Solar Simulation")
        {
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddIntegerParameter("recursion level", "rec", "recusrion level", GH_ParamAccess.item);
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddPointParameter("points", "points", "points", GH_ParamAccess.list);
            pManager.AddMeshParameter("mesh", "mesh", "mesh", GH_ParamAccess.list);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            int reclvl= 0;
            if (!DA.GetData(0, ref reclvl)) { return; }



            List<Point3d> points = new List<Point3d>();
            List<Mesh> meshlist = new List<Mesh>();

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

            Mesh msh = new Mesh();
            int i = msh.Vertices.Count;
            TextDot dot = new TextDot(i.ToString(), new Point3d(0, 0, 0));



            //Hemisphere
            SkyDome dome = new SkyDome(reclvl,1,1,1);           //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            //List<Point3d> points = new List<Point3d>();
            Mesh mesh = new Mesh();
            foreach (double[] p in dome.VertexCoordinatesSphere)
            {
                mesh.Vertices.Add(p[0], p[1], p[2]);
                points.Add(new Point3d(p[0], p[1], p[2]));
            }
            foreach (int[] f in dome.Faces)
            {
                mesh.Faces.AddFace(f[0], f[1], f[2]);
            }
            mesh.UnifyNormals();


            meshlist.Add(mesh);


            DA.SetDataList(0, points);
            DA.SetDataList(1, meshlist);


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