using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;

using sm = SolarModel;

namespace GHSolar
{
    public class GHTestSunriseSunset : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the GHTestSunriseSunset class.
        /// </summary>
        public GHTestSunriseSunset()
          : base("GHTestSunriseSunset", "Nickname",
              "Description",
              "Category", "Subcategory")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddNumberParameter("theta", "theta", "theta, latitude", GH_ParamAccess.item);
            pManager.AddIntegerParameter("day", "day", "day, 1 to 365", GH_ParamAccess.item);
            pManager.AddBooleanParameter("sunrise?", "sunrise?", "sunrise?", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddNumberParameter("azimuth", "azmiuth", "azimuth", GH_ParamAccess.item);
            pManager.AddNumberParameter("declination angle", "declination angle", "declination angle", GH_ParamAccess.item);
            pManager.AddNumberParameter("hour angle", "hour angle", "hour angle", GH_ParamAccess.item);
            pManager.AddNumberParameter("altitude", "altitude", "altitude", GH_ParamAccess.item);
            pManager.AddNumberParameter("daylength", "daylength", "daylength", GH_ParamAccess.item);
            pManager.AddNumberParameter("sunrisehour", "sunrisehour", "sunrisehour", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            double theta = 0.0;
            if (!DA.GetData(0, ref theta)) { return; }

            int day = 1;
            if (!DA.GetData(1, ref day)) { return; }

            bool sunrise = true;
            if (!DA.GetData(2, ref sunrise)) { return; }

            double[] results = sm.SunVector.GetSunriseSunsetAzimuth(theta, day, sunrise);
            DA.SetData(0, results[0]);
            DA.SetData(1, results[1]);
            DA.SetData(2, results[2]);
            DA.SetData(3, results[3]);
            DA.SetData(4, results[4]);
            DA.SetData(5, results[5]);
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
                return null;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("a9b770d3-2af2-4da2-9db2-af0da8c19bb0"); }
        }
    }
}