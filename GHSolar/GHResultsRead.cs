using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;

namespace GHSolar
{
    public class GHResultsRead : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the GHResultsTimeseries class.
        /// </summary>
        public GHResultsRead()
            : base("ReadResults", "ReadResults",
                "Read irradiation results. Create timeseries for hourly values.",
                "EnergyHubs", "Solar Simulation")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("I", "I", "Results data from solar irradiation calculation.", GH_ParamAccess.item);
            pManager.AddIntegerParameter("Value", "Value",
                "Select the value to output: [0] = Total annual irradiation [kWh/a], [1] = Beam annual [kWh/a], [2] = Diffuse annual [kWh/a], [3] = Total hourly [W], [4] = Beam hourly [W], [5] = Diffuse hourly [W].",
                GH_ParamAccess.item);
            pManager[1].Optional = true;
            pManager.AddIntegerParameter("SP", "SP", "Select sensor point to read values from.", GH_ParamAccess.item);
            pManager[2].Optional = true;
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddNumberParameter("results", "results", "Read results, output as list of double", GH_ParamAccess.list);
            pManager.AddPointParameter("SP", "SP", "Selected sensor point.", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            List<double> valin = new List<double>();

            cResults results = null;
            if (!DA.GetData(0, ref results)) { return; }

            int outputType = 0;
            if (!DA.GetData(1, ref outputType)) { outputType = 0; }

            int sp = 0; //sensor point, or mesh vertex
            if (!DA.GetData(2, ref sp)) { sp = 0; }

            switch (outputType)
            {
                case 0:
                    valin.Add(results.I_total[sp]);
                    break;
                case 1:
                    valin.Add(results.Ib_total[sp]);
                    break;
                case 2:
                    valin.Add(results.Id_total[sp]);
                    break;
                case 3:
                    for (int t = 0; t < results.I_hourly.ColumnCount; t++)
                        valin.Add(results.I_hourly[sp, t]);
                    break;
                case 4:
                    for (int t = 0; t < results.Ib_hourly.ColumnCount; t++)
                        valin.Add(results.Ib_hourly[sp, t]);
                    break;
                case 5:
                    for (int t = 0; t < results.Id_hourly.ColumnCount; t++)
                        valin.Add(results.Id_hourly[sp, t]);
                    break;
            }

            
            DA.SetDataList(0, valin);
            DA.SetData(1, results.coords[sp]);

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
                return GHSolar.Properties.Resources.pic_graph;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("{82f0b34e-96b1-4c2e-8c8f-de3f9fdd47f7}"); }
        }
    }
}