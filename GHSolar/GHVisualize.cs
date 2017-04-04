using System;
using System.Collections.Generic;
using Grasshopper.Kernel;
using Rhino.Geometry;

using System.Drawing;
using System.Linq;

/*
 * GHVisualize.cs
 * Copyright 2017 Christoph Waibel <chwaibel@student.ethz.ch>
 * 
 * This work is licensed under the GNU GPL license version 3.
*/

namespace GHSolar
{
    public class GHVisualize : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the GHVisualize class.
        /// </summary>
        public GHVisualize()
            : base("VisualizeMeshSolar", "VisuMshSol",
                "Visualize Irradiation on a Mesh.",
                "EnergyHubs", "Solar Simulation")
        {
        }
        
        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddMeshParameter("Mesh", "Mesh", "Analysis mesh", GH_ParamAccess.item);
            pManager.AddNumberParameter("I", "I", "Irradiation values per mesh vertex", GH_ParamAccess.list);
            pManager.AddNumberParameter("min", "min", "Minimum value for color gradient", GH_ParamAccess.item);
            pManager[2].Optional = true;
            pManager.AddNumberParameter("max", "max", "Maximum value for color gradient", GH_ParamAccess.item);
            pManager[3].Optional = true;
            pManager.AddIntegerParameter("clr", "clr", "colour sheme. 0: Blue (min) - Red - Yellow (max); 1: Blue (min) - Red - Yellow (max); 1: Blue (min) - Green - Red (max); 2: Black (min) - White (max).", GH_ParamAccess.item);
            pManager[4].Optional = true;        
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddMeshParameter("Mesh", "Mesh", "Coloured mesh", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            Mesh mshin = new Mesh();
            if (!DA.GetData(0, ref mshin)) { return; }

            List<double> valin = new List<double>();
            if (!DA.GetDataList(1, valin)) { return; }

            double min = 0;
            if (!DA.GetData(2, ref min)) { min = valin.Min(); }

            double max = 1;
            if (!DA.GetData(3, ref max)) { max = valin.Max(); }

            int clr = 0;
            if(!DA.GetData(4, ref clr)) clr = 0;


            Color c = new Color();
            Mesh mshcol = new Mesh();
            int count = 0;
            for (int i = 0; i < mshin.Faces.Count; i++)
            {
                c = Utilities.GetRGB(clr, valin[mshin.Faces[i].A], max, min);
                mshcol.Vertices.Add(mshin.Vertices[mshin.Faces[i].A]);
                mshcol.VertexColors.SetColor(count, c);

                c = Utilities.GetRGB(clr, valin[mshin.Faces[i].B], max, min);
                mshcol.Vertices.Add(mshin.Vertices[mshin.Faces[i].B]);
                mshcol.VertexColors.SetColor(count + 1, c);

                c = Utilities.GetRGB(clr, valin[mshin.Faces[i].C], max, min);
                mshcol.Vertices.Add(mshin.Vertices[mshin.Faces[i].C]);
                mshcol.VertexColors.SetColor(count + 2, c);

                if (mshin.Faces[i].IsQuad)
                {
                    c = Utilities.GetRGB(clr, valin[mshin.Faces[i].D], max, min);
                    mshcol.Vertices.Add(mshin.Vertices[mshin.Faces[i].D]);
                    mshcol.VertexColors.SetColor(count + 3, c);
                    mshcol.Faces.AddFace(count, count + 1, count + 2, count +3);
                    count += 4;
                }else{
                    mshcol.Faces.AddFace(count, count + 1, count + 2);
                    count += 3;
                }
   

            }


            DA.SetData(0,mshcol);

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
                return GHSolar.Properties.Resources.pic_visu;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("{b36e5547-2a9f-4706-8c83-3f776d9f3e97}"); }
        }
    }




    public static class Utilities
    {
        public static Color GetRGB(int colourSheme, double quantity, double top, double low)
        {
            double RR = 0.0;
            double GG = 0.0;
            double BB = 0.0;


            quantity = (quantity - low) / (top - low);
            double third = 1.0 / 5.0;

            switch (colourSheme)
            {
                case 0:
                    if (quantity > third && quantity <= 2.0 * third)
                    {
                        RR = (quantity - third) * (255.0 / third);
                        GG = 0.0;
                        BB = 255 - ((quantity - third) * (255.0 / third));
                    }
                    else if (quantity > 2.0 * third)
                    {
                        RR = 255.0;
                        GG = (quantity - 2.0 * third) * (255.0 / third);
                        BB = 0.0;
                    }
                    else
                    {
                        RR = 0.0;
                        GG = 0.0;
                        BB = 255.0;
                    }
                    break;
                case 1:
                    third = 1.0 / 3.0;
                    if (quantity > third && quantity <= 2.0 * third)
                    {
                        RR = (quantity - third) * (255.0 / third);
                        GG = 255.0;
                        BB = 255.0 - ((quantity - third) * (255.0 / third));
                    }
                    else if (quantity > 2.0 * third)
                    {
                        RR = 255.0;
                        GG = 255.0 - ((quantity - 2.0 * third) * (255.0 / third));
                        BB = 0.0;
                    }
                    else
                    {
                        RR = 0.0;
                        GG = quantity * (255.0 / third);
                        BB = 255.0;
                    }
                    break;
                case 2:
                    RR = quantity * (255.0 / 1);
                    GG = quantity * (255.0 / 1);
                    BB = quantity * (255.0 / 1);
                    break;
            }

            if (RR > 255) RR = 255;
            else if (RR < 0) RR = 0;
            if (GG > 255) GG = 255;
            else if (GG < 0) GG = 0;
            if (BB > 255) BB = 255;
            else if (BB < 0) BB = 0;
            return Color.FromArgb((int)RR, (int)GG, (int)BB);

        }
    }
}