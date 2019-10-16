using System.Collections.Generic;
using Rhino.Geometry;
using SolarModel;

namespace GHSolar
{
    public class CResults
    {
        public Matrix I_hourly;
        public Matrix Ib_hourly;
        public Matrix Id_hourly;
        public List<double> I_total = new List<double>();
        public List<double> Ib_total = new List<double>();
        public List<double> Id_total = new List<double>();
        public List<Point3d> coords = new List<Point3d>();

        public CResults(List<double> I_total, List<double> Ib_total, List<double> Id_total,
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


    public class CResultsInterreflections
    {
        //variables for diffuse inter-reflection
        //internal Sensorpoints[] Idiffuse_SPs;
        public int[][] Idiff_obstacles;
        public int[][] Idiff_domevertices;
        public SkyDome[] Idiff_domes;

        public List<List<double>> diffSP_beta_list = new List<List<double>>();
        public List<List<double>> diffSP_psi_list = new List<List<double>>();
        public List<List<Sensorpoints.v3d>> diffSP_normal_list = new List<List<Sensorpoints.v3d>>();
        public List<List<Sensorpoints.p3d>> diffSP_coord_list = new List<List<Sensorpoints.p3d>>();

        /// <param name="Ispecular">Normal irradiation values [i][t][m] for each sensor point i, each solar vector t and each reflected ray m.</param>
        /// <param name="Inormals">Normal vectors [i][t][m] for each sensor point i, each solar vector t and each reflected ray m.</param>
        //variables for specular inter-reflection
        public double[][][] Ispecular2;// = new double[mshvrt.Length][][];
        public Vector3d[][][] Inormals2;// = new Vector3d[mshvrt.Length][][];

        public CResultsInterreflections(List<List<double>> diffSP_beta_list, List<List<double>> diffSP_psi_list, 
            List<List<Sensorpoints.v3d>> diffSP_normal_list,List<List<Sensorpoints.p3d>> diffSP_coord_list,
            int[][] Idiff_obstacles, int[][] Idiff_domevertices, SkyDome[] Idiff_domes,
            double[][][] Ispecular2, Vector3d[][][] Inormals2)
        {
            int SPcount = diffSP_beta_list.Count;

            this.diffSP_beta_list = diffSP_beta_list;
            this.diffSP_psi_list = diffSP_psi_list;
            this.diffSP_normal_list = diffSP_normal_list;
            this.diffSP_coord_list = diffSP_coord_list;

            //this.Idiffuse_SPs = new Sensorpoints[SPcount];
            this.Idiff_obstacles = new int[SPcount][];
            this.Idiff_domevertices = new int[SPcount][];
            this.Idiff_domes = new SkyDome[SPcount];

            this.Ispecular2 = new double[SPcount][][];
            this.Inormals2 = new Vector3d[SPcount][][];

            for (int i = 0; i < SPcount; i++)
            {
                //this.Idiffuse_SPs[i] = Idiffuse_SPs[i];
                this.Idiff_obstacles[i] = Idiff_obstacles[i];
                this.Idiff_domevertices[i] = Idiff_domevertices[i];
                this.Idiff_domes[i] = Idiff_domes[i];

                if (Misc.IsNullOrEmpty(Ispecular2[i])) continue;
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
