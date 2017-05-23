using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;

/*
 * GHSolarTree.cs
 * Copyright 2017 Christoph Waibel <chwaibel@student.ethz.ch>
 * 
 * This work is licensed under the GNU GPL license version 3.
*/

namespace GHSolar
{
    public class GHSolarTree : GH_Component
    {
        public GHSolarTree()
            : base("Tree Obstacle", "TreeObst",
                "Tree obstacle object used for irradiation calculation.",
                "EnergyHubs", "Solar Simulation")
        {
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddMeshParameter("Mesh", "Mesh", "Obstacle as mesh", GH_ParamAccess.item);
            pManager.AddNumberParameter("Leaves", "Leaves", "8760 time series for leaves coverage of the tree. Value between 0 - 1. 1 means it is full of trees and hence a full obstructor (no light penetrates).", GH_ParamAccess.list);
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Tree Object", "TreeObj", "Tree object for solar irradiation calculation.", GH_ParamAccess.item);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            //each tree object is : 
            //  (i) a mesh obstacle, 
            //  (ii) a 8760-timeseries of 0-1 fractions, indicating leave-coverage 1 is full of leaves=full obstruction.
            Mesh mesh = new Mesh();
            if (!DA.GetData(0, ref mesh)) { return; }

            List<double> leaves = new List<double>();
            if (!DA.GetDataList(1, leaves)) { return; }

            cSemiPermObject tree = new cSemiPermObject(mesh, leaves.ToArray());

            DA.SetData(0, tree);
        }

        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                //You can add image files to your project resources and access them like this:
                // return Resources.IconForThisComponent;
                return GHSolar.Properties.Resources.pic_tree2;
            }
        }

        public override Guid ComponentGuid
        {
            get { return new Guid("{e3ea7e61-2151-4f46-b20f-aeb562b27e81}"); }
        }
    }


  
}