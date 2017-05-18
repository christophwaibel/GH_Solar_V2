using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using Rhino.Geometry;

using SolarModel;

namespace GHSolar
{
    internal class cResults
    {
        internal Matrix I_hourly;
        internal Matrix Ib_hourly;
        internal Matrix Id_hourly;
        internal List<double> I_total = new List<double>();
        internal List<double> Ib_total = new List<double>();
        internal List<double> Id_total = new List<double>();
        internal List<Point3d> coords = new List<Point3d>();

        internal cResults(List<double> I_total, List<double> Ib_total, List<double> Id_total,
            Matrix I_hourly, Matrix Ib_hourly, Matrix Id_hourly, 
            List<Point3d> coords) 
        {
            this.I_hourly = new Matrix(I_hourly.RowCount, I_hourly.ColumnCount);
            this.I_hourly = I_hourly;

            this.Ib_hourly = new Matrix(Ib_hourly.RowCount, Ib_hourly.ColumnCount);
            this.Ib_hourly = Ib_hourly;

            this.Id_hourly = new Matrix(Id_hourly.RowCount, Id_hourly.ColumnCount);
            this.Id_hourly = Id_hourly;

            this.I_total = I_total;
            this.Ib_total = Ib_total;
            this.Id_total = Id_total;
            this.coords = coords;         
        }
    }


    internal class cResultsInterreflections
    {
        //variables for diffuse inter-reflection
        internal Sensorpoints[] Idiffuse_SPs;
        internal int[][] Idiff_obstacles;
        internal int[][] Idiff_domevertices;
        internal SkyDome[] Idiff_domes;

        /// <param name="Ispecular">Normal irradiation values [i][t][m] for each sensor point i, each solar vector t and each reflected ray m.</param>
        /// <param name="Inormals">Normal vectors [i][t][m] for each sensor point i, each solar vector t and each reflected ray m.</param>
        //variables for specular inter-reflection
        internal double[][][] Ispecular2;// = new double[mshvrt.Length][][];
        internal Vector3d[][][] Inormals2;// = new Vector3d[mshvrt.Length][][];

        internal cResultsInterreflections(Sensorpoints[] Idiffuse_SPs, int[][] Idiff_obstacles, int[][] Idiff_domevertices, SkyDome[] Idiff_domes,
            double[][][] Ispecular2, Vector3d[][][] Inormals2)
        {
            int SPcount = Idiffuse_SPs.Length;

            this.Idiffuse_SPs = new Sensorpoints[SPcount];
            this.Idiff_obstacles = new int[SPcount][];
            this.Idiff_domevertices = new int[SPcount][];
            this.Idiff_domes = new SkyDome[SPcount];

            this.Ispecular2 = new double[SPcount][][];
            this.Inormals2 = new Vector3d[SPcount][][];

            for (int i = 0; i < SPcount; i++)
            {
                this.Idiffuse_SPs[i] = Idiffuse_SPs[i];
                this.Idiff_obstacles[i] = Idiff_obstacles[i];
                this.Idiff_domevertices[i] = Idiff_domevertices[i];
                this.Idiff_domes[i] = Idiff_domes[i];

                this.Ispecular2[i] = new double[Ispecular2[i].Length][];
                this.Inormals2[i] = new Vector3d[Inormals2[i].Length][];
                for (int t = 0; t < Ispecular2[i].Length; t++)
                {
                    this.Ispecular2[i][t] = Ispecular2[i][t];
                    this.Inormals2[i][t] = Inormals2[i][t];
                }
            }
        }
    }
}
