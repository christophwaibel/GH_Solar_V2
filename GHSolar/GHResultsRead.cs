using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;

/*
 * GHResultsRead.cs
 * Copyright 2017 Christoph Waibel <chwaibel@student.ethz.ch>
 * 
 * This work is licensed under the GNU GPL license version 3.
*/

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
            pManager.AddMeshParameter("Mesh", "Mesh", "Analysis mesh", GH_ParamAccess.item);
            pManager.AddGenericParameter("I", "I", "Results data from solar irradiation calculation.", GH_ParamAccess.item);
            pManager.AddIntegerParameter("Value", "Value",
                "Select the value to output: [0] = Total specific annual irradiation [kWh/m^2a], [1] = Specific beam annual [kWh/m^2a], [2] = Specific diffuse annual [kWh/m^2a], " +
                "[3] = Total annual irradiation per mesh [kWh/a], " +
                "[4] = Total specific hourly [W/m^2], [5] = Specific beam hourly [W/m^2], [6] = Specific diffuse hourly [W/m^2], " +
                "[7] = Total hourly per mesh [W], [8] = Beam hourly per mesh [W], [9] = Diffuse hourly per mesh [W].",
                GH_ParamAccess.item);
            pManager[2].Optional = true;
            pManager.AddIntegerParameter("SP", "SP", "Select sensor point to read values from (not for Value type 3 or 7).", GH_ParamAccess.item);
            pManager[3].Optional = true;
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
            List<double[]> valin2 = new List<double[]>();

            Mesh mshin = new Mesh();
            if (!DA.GetData(0, ref mshin)) { return; }

            CResults results = null;
            if (!DA.GetData(1, ref results)) { return; }

            int outputType = 0;
            if (!DA.GetData(2, ref outputType)) { outputType = 0; }

            int sp = 0; //sensor point, or mesh vertex
            if (!DA.GetData(3, ref sp)) { sp = 0; }

            double copyval;
            switch (outputType)
            {
                case 0: // [kWh/m^2a] total
                    copyval = results.I_total[sp];
                    valin.Add(copyval);
                    break;
                case 1: // [kWh/m^2a] beam
                    copyval = results.Ib_total[sp];
                    valin.Add(copyval);
                    break;
                case 2: // [kWh/m^2a] diffuse
                    copyval = results.Id_total[sp];
                    valin.Add(copyval);
                    break;
                case 3: // [kWh/a] total
                    valin = new List<double>(results.I_total);
                    break;
                case 4: // [W/m^2] total
                    for (int t = 0; t < results.I_hourly.ColumnCount; t++)
                        valin.Add(results.I_hourly[sp, t]);
                    break;
                case 5: // [W/m^2] beam
                    for (int t = 0; t < results.Ib_hourly.ColumnCount; t++)
                        valin.Add(results.Ib_hourly[sp, t]);
                    break;
                case 6: // [W/m^2] diffuse
                    for (int t = 0; t < results.Id_hourly.ColumnCount; t++)
                        valin.Add(results.Id_hourly[sp, t]);
                    break;
                case 7: // [W] total
                    for (int i = 0; i < results.I_hourly.RowCount; i++)
                    {
                        double[] val_t = new double[results.I_hourly.ColumnCount];
                        for (int t = 0; t < results.I_hourly.ColumnCount; t++)
                            val_t[t] = results.I_hourly[i, t];
                        valin2.Add(val_t);
                    }
                    break;
                case 8: // but only beam
                    for (int i = 0; i < results.Ib_hourly.RowCount; i++)
                    {
                        double[] val_t = new double[results.Ib_hourly.ColumnCount];
                        for (int t = 0; t < results.Ib_hourly.ColumnCount; t++)
                            val_t[t] = results.Ib_hourly[i, t];
                        valin2.Add(val_t);
                    }
                    break;
                case 9: // but only diffuse
                    for (int i = 0; i < results.Id_hourly.RowCount; i++)
                    {
                        double[] val_t = new double[results.Id_hourly.ColumnCount];
                        for (int t = 0; t < results.Id_hourly.ColumnCount; t++)
                            val_t[t] = results.Id_hourly[i, t];
                        valin2.Add(val_t);
                    }
                    break;
            }

            if (outputType == 3)
            {
                double[] mshFaceAreas = new double[mshin.Faces.Count];

                double totVal = 0;
                for (int i = 0; i < mshin.Faces.Count; i++)
                {
                    mshFaceAreas[i] = CMisc.getMeshFaceArea(i, mshin);

                    double FaceVal;
                    double valVertex1 = valin[mshin.Faces[i].A];
                    double valVertex2 = valin[mshin.Faces[i].B];
                    double valVertex3 = valin[mshin.Faces[i].C];
                    if (mshin.Faces[i].IsQuad)
                    {
                        double valVertex4 = valin[mshin.Faces[i].D];
                        FaceVal = ((valVertex1 + valVertex2 + valVertex3 + valVertex4) / 4) * CMisc.getMeshFaceArea(i, mshin);
                    }
                    else
                    {
                        FaceVal = ((valVertex1 + valVertex2 + valVertex3) / 3) * CMisc.getMeshFaceArea(i, mshin);
                    }
                    totVal += FaceVal;
                }
                valin = new List<double>();
                valin.Add(totVal);
            }
            else if (outputType == 7 || outputType == 8 || outputType == 9)
            {
                List<double> valout = new List<double>();

                double[] mshFaceAreas = new double[mshin.Faces.Count];

                for (int t = 0; t < results.I_hourly.ColumnCount; t++)
                {
                    double totVal = 0;
                    for (int i = 0; i < mshin.Faces.Count; i++)
                    {
                        mshFaceAreas[i] = CMisc.getMeshFaceArea(i, mshin);

                        double FaceVal;
                        double valVertex1 = valin2[mshin.Faces[i].A][t];
                        double valVertex2 = valin2[mshin.Faces[i].B][t];
                        double valVertex3 = valin2[mshin.Faces[i].C][t];
                        if (mshin.Faces[i].IsQuad)
                        {
                            double valVertex4 = valin2[mshin.Faces[i].D][t];
                            FaceVal = ((valVertex1 + valVertex2 + valVertex3 + valVertex4) / 4) * CMisc.getMeshFaceArea(i, mshin);
                        }
                        else
                        {
                            FaceVal = ((valVertex1 + valVertex2 + valVertex3) / 3) * CMisc.getMeshFaceArea(i, mshin);
                        }
                        totVal += FaceVal;
                    }
                    valout.Add(totVal);
                }

                valin = new List<double>();
                valin = valout;
            }

            if (outputType == 0 || outputType == 1 || outputType == 2 || outputType == 3)
            {
                List<double> valin_copy = new List<double>(valin);
                for (int i = 0; i < valin.Count; i++)
                {
                    valin_copy[i] *= 0.001;
                }
                valin = new List<double>(valin_copy);
            }

            DA.SetDataList(0, valin);
            DA.SetData(1, results.coords[sp]);

        }


        public override GH_Exposure Exposure => GH_Exposure.tertiary;


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