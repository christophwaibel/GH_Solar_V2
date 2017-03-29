using System;
using System.Collections.Generic;
using Grasshopper.Kernel;
using Rhino.Geometry;

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
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
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