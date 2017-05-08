using System;
using System.Collections.Generic;

using System.Threading.Tasks;

using Grasshopper.Kernel;
using Rhino.Geometry;

/*
 * GHSolarObstacle.cs
 * Copyright 2017 Christoph Waibel <chwaibel@student.ethz.ch>
 * 
 * This work is licensed under the GNU GPL license version 3.
*/

namespace GHSolar
{
    public class GHSolarObstacle : GH_Component
    {
        public GHSolarObstacle()
            : base("Obstacle object", "ObstObj",
                "Obstacle object used for irradiation calculation.",
                "EnergyHubs", "Solar Simulation")
        {
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddMeshParameter("Mesh", "Mesh", "Obstacle as mesh", GH_ParamAccess.item);
            pManager.AddNumberParameter("Albedo", "Albedo", "8760 time series for albedo values. Value between 0 - 1.", GH_ParamAccess.list);
            pManager.AddIntegerParameter("Reflection", "Reflection", "Reflection type: 0: diffuse; 1: specular.", GH_ParamAccess.item);
            pManager.AddNumberParameter("Tolerance", "Tolerance", "Tolerance, used to offset mesh face centers in normal direction", GH_ParamAccess.item);
            pManager[3].Optional = true;
            pManager.AddTextParameter("Name", "Name", "Name", GH_ParamAccess.item);
            pManager.AddBooleanParameter("MT", "MT", "Multi threading", GH_ParamAccess.item);
            pManager[5].Optional = true;
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Obstacle Object", "ObstObj", "Obstacle object for solar irradiation calculation.", GH_ParamAccess.item);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            //each obstacle consists of: 
            //  (i) a mesh obstacle, 
            //  (ii) a 8760-timeseries of 0-1 fraction indicating albedo, 
            //  (iii) a boolean 0-1 for either diffuse (0) or specular (1) reflection.

            Mesh mesh = new Mesh();
            if (!DA.GetData(0, ref mesh)) { return; }

            List<double> alb = new List<double>();
            if (!DA.GetDataList(1, alb)) { return; }

            int refl = 0;
            if (!DA.GetData(2, ref refl)) { return; }

            double tolerance = 0.01;
            if (!DA.GetData(3, ref tolerance)) { tolerance = 0.01; }

            string name = null;
            if (!DA.GetData(4, ref name)) { return; }

            bool mt = false;
            if (!DA.GetData(5, ref mt)) { mt = false; }

            ObstacleObject obst = new ObstacleObject(mesh, alb, refl, tolerance, name, mt);

            DA.SetData(0, obst);
        }

        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                //You can add image files to your project resources and access them like this:
                // return Resources.IconForThisComponent;
                return GHSolar.Properties.Resources.pic_obstacle3;
            }
        }

        public override Guid ComponentGuid
        {
            get { return new Guid("{423b0270-11c7-42a0-9284-19dabae3541b}"); }
        }
    }

    internal class ObstacleObject
    {
        internal Mesh mesh;
        internal List<double> albedos;  //8760 values
        internal int reflType;
        internal Vector3d[] normals;
        internal Vector3d[] normalsRev;
        internal Point3d[] faceCen;

        internal double tolerance;
        internal string name;

        internal ObstacleObject(Mesh _mesh, List<double> _albedos, int _reflType, double _tolerance, string _name, bool mt)
        {
            mesh = _mesh;
            albedos = new List<double>(_albedos);
            reflType = _reflType;
            name = _name;
            tolerance = _tolerance; 

            mesh.FaceNormals.ComputeFaceNormals();

            normals = new Vector3d[mesh.Faces.Count];
            normalsRev = new Vector3d[mesh.Faces.Count];
            faceCen = new Point3d[mesh.Faces.Count];
            if (!mt)
            {
                for (int k = 0; k < mesh.Faces.Count; k++)
                {
                    normals[k] = mesh.FaceNormals[k];
                    normalsRev[k] = Vector3d.Negate(normals[k]);
                    Point3d cen0 = mesh.Faces.GetFaceCenter(k);
                    faceCen[k] = new Point3d(Point3d.Add(cen0, Vector3d.Multiply(Vector3d.Divide(normals[k], normals[k].Length), tolerance)));
                }
            }
            else
            {
                Parallel.For(0, mesh.Faces.Count, k =>
                {
                    normals[k] = mesh.FaceNormals[k];
                    normalsRev[k] = Vector3d.Negate(normals[k]);
                    Point3d cen0 = mesh.Faces.GetFaceCenter(k);
                    faceCen[k] = new Point3d(Point3d.Add(cen0, Vector3d.Multiply(Vector3d.Divide(normals[k], normals[k].Length), tolerance)));
                });
            }


        }
    }
}