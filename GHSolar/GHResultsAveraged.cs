using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;
using Rhino.Geometry.Collections;

namespace GHSolar
{
    public class GHResultsAveraged : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the GHResultsAveraged class.
        /// </summary>
        public GHResultsAveraged()
            : base("ReadResultsAvg", "ReadResultsAvg",
                "Read irradiation results, averaged per mesh face. Create timeseries for hourly values.",
                "EnergyHubs", "Solar Simulation")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddMeshParameter("Mesh", "Mesh", "Analysis mesh", GH_ParamAccess.item);
            pManager.AddGenericParameter("I", "I", "Results data from solar irradiation calculation.", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddNumberParameter("avgI", "avgI", "averaged hourly solar potentials in W/m^2", GH_ParamAccess.tree);
            pManager.AddNumberParameter("area", "area", "area per mesh face in m^2", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            Mesh mshin = new Mesh();
            if (!DA.GetData(0, ref mshin)) { return; }

            cResults results = null;
            if (!DA.GetData(1, ref results)) { return; }


            //List<double[]> avgI = new List<double[]>();
            Grasshopper.DataTree<double> avgI = new Grasshopper.DataTree<double>();
                
            List<double> areas = new List<double>();

            MeshFaceList mshfaces = mshin.Faces;


            for(int i=0; i<mshfaces.Count; i++)
            {
                MeshFace fc = mshfaces[i];

                int spA = fc.A;
                int spB = fc.B;
                int spC = fc.C;
                int spD = fc.D;

                Grasshopper.Kernel.Data.GH_Path ghpath = new Grasshopper.Kernel.Data.GH_Path(i);

                List<double> sp_I = new List<double>();
                for (int t = 0; t < results.I_hourly.ColumnCount; t++)
                {
                    avgI.Add((results.I_hourly[spA, t] + results.I_hourly[spB, t] + results.I_hourly[spC, t] + results.I_hourly[spD, t]) / 4, ghpath);

                    //sp_I.Add((results.I_hourly[spA, t] + results.I_hourly[spB, t] + results.I_hourly[spC, t] + results.I_hourly[spD, t]) / 4);
                }

               

                areas.Add(cMisc.getMeshFaceArea(i, mshin));

            }

            

            DA.SetDataTree(0, avgI);
            DA.SetDataList(1, areas);


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
                return GHSolar.Properties.Resources.pic_graph_avg;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("08e54258-9975-43a8-8e53-c443407b8dc1"); }
        }
    }
}