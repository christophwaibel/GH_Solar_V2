using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;

using System.IO;


/*
 * GHReadTxt.cs
 * Copyright 2017 Christoph Waibel <chwaibel@student.ethz.ch>
 * 
 * This work is licensed under the GNU GPL license version 3.
*/

namespace GHSolar
{
    public class GHReadTxt : GH_Component
    {
        public GHReadTxt()
            : base("Read Text", "ReadTxt",
                "Reads a *.txt file, given a path.",
                "EnergyHubs", "Misc")
        {
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddTextParameter("Path", "Path", "Path of the text file to read", GH_ParamAccess.item); 
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddNumberParameter("Data", "Data", "Data as a list of doubles.", GH_ParamAccess.list);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            string path = null;
            if (!DA.GetData(0, ref path)) { return; }

            string[] lines = File.ReadAllLines(@path);
            List<double> data = new List<double>();
            foreach (string line in lines)
            {
                data.Add(Convert.ToDouble(line));
            }

            DA.SetDataList(0, data);
        }

        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                //You can add image files to your project resources and access them like this:
                // return Resources.IconForThisComponent;
                return GHSolar.Properties.Resources.pic_data;
            }
        }

        public override Guid ComponentGuid
        {
            get { return new Guid("{c30d858a-6e15-42ae-b673-525f89eec905}"); }
        }
    }
}