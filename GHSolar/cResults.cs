using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using Rhino.Geometry;
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
}
