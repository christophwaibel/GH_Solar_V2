using System;
using System.Collections.Generic;

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
            //  (iii) a boolean 0-1 for either specular (0) or diffuse (1) reflection.

            Mesh mesh = new Mesh();
            if (!DA.GetData(0, ref mesh)) { return; }

            List<double> alb = new List<double>();
            if (!DA.GetDataList(1, alb)) { return; }

            int refl = 0;
            if (!DA.GetData(2, ref refl)) { return; }

            ObstacleObject obst = new ObstacleObject(mesh, alb, refl);

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
        internal List<double> albedos;
        internal int reflType;

        internal ObstacleObject(Mesh _mesh, List<double> _albedos, int _reflType)
        {
            mesh = _mesh;
            albedos = new List<double>(_albedos);
            reflType = _reflType;
        }
    }
}