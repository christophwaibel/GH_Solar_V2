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
            pManager.AddIntegerParameter("Value", "Value",
                "Select the value to output: [0] = Total specific annual irradiation [kWh/m^2a], [1] = Specific beam annual [kWh/m^2a], [2] = Specific diffuse annual [kWh/m^2a], [3] = Total specific hourly [W/m^2], [4] = Specific beam hourly [W/m^2], [5] = Specific diffuse hourly [W/m^2], [6] = Solar Acces (defined as if beam radiation for at least 3 of the 4 vertices is not 0; outputs 1 if condition fulfilled, 0 otherwise).",
                GH_ParamAccess.item);

        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddNumberParameter("avgI", "avgI", "averaged hourly solar irradiation or irradiance per mesh face as defined by 'Value' in the input.", GH_ParamAccess.tree);
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

            CResults results = null;
            if (!DA.GetData(1, ref results)) { return; }

            int outputType = 0;
            if (!DA.GetData(2, ref outputType)) { outputType = 0; }


            Grasshopper.DataTree<double> avgI = new Grasshopper.DataTree<double>();

            List<double> areas = new List<double>();

            MeshFaceList mshfaces = mshin.Faces;


            for (int m = 0; m < mshfaces.Count; m++)
            {
                MeshFace fc = mshfaces[m];

                int spA = fc.A;
                int spB = fc.B;
                int spC = fc.C;
                int spD = fc.D;
                int intvert = 3;
                double quadmulti = 0;
                if (fc.IsQuad)
                {
                    intvert = 4;
                    quadmulti = 1;
                }

                Grasshopper.Kernel.Data.GH_Path ghpath = new Grasshopper.Kernel.Data.GH_Path(m);

                List<double> sp_I = new List<double>();

                switch (outputType)
                {
                    case 0:
                        avgI.Add((results.I_total[spA] + results.I_total[spB] + results.I_total[spC] + (results.I_total[spD] * quadmulti)) / intvert, ghpath);
                        break;
                    case 1:
                        avgI.Add((results.Ib_total[spA] + results.Ib_total[spB] + results.Ib_total[spC] + (results.Ib_total[spD] * quadmulti)) / intvert, ghpath);
                        break;
                    case 2:
                        avgI.Add((results.Id_total[spA] + results.Id_total[spB] + results.Id_total[spC] + (results.Id_total[spD] * quadmulti)) / intvert, ghpath);
                        break;
                    case 3:
                        for (int t = 0; t < results.I_hourly.ColumnCount; t++)
                            avgI.Add((results.I_hourly[spA, t] + results.I_hourly[spB, t] + results.I_hourly[spC, t] + (results.I_hourly[spD, t] * quadmulti)) / intvert, ghpath);
                        break;
                    case 4:
                        for (int t = 0; t < results.Ib_hourly.ColumnCount; t++)
                            avgI.Add((results.Ib_hourly[spA, t] + results.Ib_hourly[spB, t] + results.Ib_hourly[spC, t] + (results.Ib_hourly[spD, t] * quadmulti)) / intvert, ghpath);
                        break;
                    case 5:
                        for (int t = 0; t < results.Id_hourly.ColumnCount; t++)
                            avgI.Add((results.Id_hourly[spA, t] + results.Id_hourly[spB, t] + results.Id_hourly[spC, t] + (results.Id_hourly[spD, t] * quadmulti)) / intvert, ghpath);
                        break;
                    case 6:
                        for (int t = 0; t < results.Ib_hourly.ColumnCount; t++)
                        {
                            int shdwcount = 0;
                            if (results.Ib_hourly[spA, t] == 0) shdwcount++;
                            if (results.Ib_hourly[spB, t] == 0) shdwcount++;
                            if (results.Ib_hourly[spC, t] == 0) shdwcount++;
                            if (results.Ib_hourly[spD, t] == 0) shdwcount++;
                            if (shdwcount > 1) avgI.Add(0, ghpath);
                            else avgI.Add(1, ghpath);
                        }
                        break;
                }
                areas.Add(CMisc.getMeshFaceArea(m, mshin));
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