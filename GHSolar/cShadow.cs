using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using Rhino.Geometry;

using SolarModel;

/*
 * cShadow.cs
 * Copyright 2017 Christoph Waibel <chwaibel@student.ethz.ch>
 * 
 * This work is licensed under the GNU GPL license version 3.
*/

namespace GHSolar
{

    internal static class cShadow
    {
        //Direct shadow list

        //diffuse shadow list


        //inter-reflections? n


        //MT versions !! is it possible with rhino?

        /// <summary>
        /// 
        /// </summary>
        /// <param name="origin"></param>
        /// <param name="origNormal"></param>
        /// <param name="tolerance">Tolerance to offset ray from sensor point.</param>
        /// <param name="vec"></param>
        /// <param name="obstacles"></param>
        /// <param name="shdw"></param>
        internal static void CalcShadow(Point3d origin, Vector3d origNormal, double tolerance, Vector3d[] vec, Mesh[] obstacles, ref bool[] shdw)
        {
            shdw = new bool[vec.Length];        //by default all elements false. true means it's obstructed
            Point3d origOffset = new Point3d(Point3d.Add(origin, Vector3d.Multiply(Vector3d.Divide(origNormal, origNormal.Length), tolerance)));
            for (int t = 0; t < vec.Length; t++)
            {
                for (int u = 0; u < obstacles.Length; u++)
                {
                    Ray3d ray = new Ray3d(origOffset, vec[t]);
                    double inters = Rhino.Geometry.Intersect.Intersection.MeshRay(obstacles[u], ray);
                    if (inters >= 0)
                    {
                        shdw[t] = true;
                        break;
                    }
                }
            }
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="origin"></param>
        /// <param name="origNormal"></param>
        /// <param name="tolerance">Tolerance to offset ray from sensor point.</param>
        /// <param name="vec"></param>
        /// <param name="obstacles"></param>
        /// <param name="shdw"></param>
        internal static void CalcShadowMT(Point3d origin, Vector3d origNormal, double tolerance, Vector3d[] vec, Mesh[] obstacles, ref bool[] shdw)
        {
            shdw = new bool[vec.Length];        //by default all elements false. true means it's obstructed
            bool[] shdw_mt = new bool[vec.Length];
            Point3d origOffset = new Point3d(Point3d.Add(origin, Vector3d.Multiply(Vector3d.Divide(origNormal, origNormal.Length), tolerance)));
            Parallel.For(0, vec.Length, t =>
            {
                for (int u = 0; u < obstacles.Length; u++)
                {
                    Ray3d ray = new Ray3d(origOffset, vec[t]);
                    double inters = Rhino.Geometry.Intersect.Intersection.MeshRay(obstacles[u], ray);
                    if (inters >= 0)
                    {
                        shdw_mt[t] = true;
                        break;
                    }
                }
            });
            shdw_mt.CopyTo(shdw, 0);
        }


        /// <summary>
        /// 
        /// </summary>
        /// <param name="origin"></param>
        /// <param name="origNormal"></param>
        /// <param name="tolerance">Tolerance to offset ray from sensor point.</param>
        /// <param name="vec"></param>
        /// <param name="obstacles"></param>
        /// <param name="shdw"></param>
        internal static void CalcShadow(Point3d origin, Vector3d origNormal, double tolerance, Vector3d[] vec, bool[] sunshine, Mesh[] obstacles, ref bool[] shdw)
        {
            shdw = new bool[vec.Length];        //by default all elements false. true means it's obstructed
            Point3d origOffset = new Point3d(Point3d.Add(origin, Vector3d.Multiply(Vector3d.Divide(origNormal, origNormal.Length), tolerance)));
            for (int t = 0; t < vec.Length; t++)
            {
                if (sunshine[t])
                {
                    for (int u = 0; u < obstacles.Length; u++)
                    {
                        Ray3d ray = new Ray3d(origOffset, vec[t]);
                        double inters = Rhino.Geometry.Intersect.Intersection.MeshRay(obstacles[u], ray);
                        if (inters >= 0)
                        {
                            shdw[t] = true;
                            break;
                        }
                    }
                }
                else
                {
                    shdw[t] = true;
                }
            }
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="origin"></param>
        /// <param name="origNormal"></param>
        /// <param name="tolerance">Tolerance to offset ray from sensor point.</param>
        /// <param name="vec"></param>
        /// <param name="obstacles"></param>
        /// <param name="shdw"></param>
        internal static void CalcShadowMT(Point3d origin, Vector3d origNormal, double tolerance, Vector3d[] vec, bool[] sunshine, Mesh[] obstacles, ref bool[] shdw)
        {
            shdw = new bool[vec.Length];        //by default all elements false. true means it's obstructed
            bool[] shdw_mt = new bool[vec.Length];
            Point3d origOffset = new Point3d(Point3d.Add(origin, Vector3d.Multiply(Vector3d.Divide(origNormal, origNormal.Length), tolerance)));
            Parallel.For(0, vec.Length, t =>
            {
                if (sunshine[t])
                {
                    for (int u = 0; u < obstacles.Length; u++)
                    {
                        Ray3d ray = new Ray3d(origOffset, vec[t]);
                        double inters = Rhino.Geometry.Intersect.Intersection.MeshRay(obstacles[u], ray);
                        if (inters >= 0)
                        {
                            shdw_mt[t] = true;
                            break;
                        }
                    }
                }
                else
                {
                    shdw_mt[t] = true;
                }
            });

            shdw_mt.CopyTo(shdw, 0);
        }



        internal static void CalcShadowMeshMT(int year, List<SunVector> sunvectors, Mesh msh, Mesh[] obst, int interpmode,
            Context.cWeatherdata weather, Context.cLocation location, int rec,
            ref List<double> I, ref List<double> Ih, ref List <double> Ib,
            ref Matrix I_hourly, ref Matrix Ih_hourly, ref Matrix Ib_hourly)
        {
            //            List<double> I = new List<double>();
            //List<double> Ih = new List<double>();
            //List<double> Ib = new List<double>();

            double snow_threshold = 10;
            double tilt_treshold = 30;


            double rad = Math.PI / 180;


            Point3d[] mshvrt = msh.Vertices.ToPoint3dArray();
            Vector3f[] mshvrtnorm = new Vector3f[mshvrt.Length];
            msh.FaceNormals.ComputeFaceNormals();

            double[] arrbeta = new double[mshvrt.Length];
            double[] arrpsi = new double[mshvrt.Length];

            Vector3d betaangle = new Vector3d(0, 0, 1);
            Vector3d psiangle = new Vector3d(0, -1, 0);
            Plane psiplane = new Plane(new Point3d(0, 0, 0), new Vector3d(0, 0, 1));
            for (int i = 0; i < mshvrt.Length; i++)
            {
                mshvrtnorm[i] = msh.Normals[i];
                //sensor point tilt angle (beta) and azimuth (psi)
                double beta = Vector3d.VectorAngle(mshvrtnorm[i], betaangle) / rad;
                double psi = Vector3d.VectorAngle(mshvrtnorm[i], psiangle, psiplane) / rad;
                if (Double.IsNaN(psi) || Double.IsInfinity(psi))
                {
                    psi = 0;
                }

                arrbeta[i] = beta;
                arrpsi[i] = psi;

            }
            Sensorpoints p = new Sensorpoints(year, weather, location, sunvectors, arrbeta, arrpsi, rec);



            List<bool[]> ShdwBeam_equinox = new List<bool[]>();
            List<bool[]> ShdwBeam_summer = new List<bool[]>();
            List<bool[]> ShdwBeam_winter = new List<bool[]>();

            List<bool[][]> ShdwBeam = new List<bool[][]>();
            int[] startDays = new int[12];
            int[] endDays = new int[12];

            List<bool[]> ShdwSky = new List<bool[]>();


            int[] equsol = SunVector.GetEquinoxSolstice(year);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            int HOYequ = (equsol[0] - 1) * 24;
            int HOYsum = (equsol[1] - 1) * 24;
            int HOYwin = (equsol[3] - 1) * 24;

            if (interpmode == 0)
            {
                Parallel.For(0, mshvrt.Length, i =>
                {
                    mshvrtnorm[i] = msh.Normals[i];
                    Point3d orig = new Point3d(mshvrt[i].X, mshvrt[i].Y, mshvrt[i].Z);

                    //sky dome diffuse
                    Vector3d[] vec_sky = new Vector3d[p.sky[i].VerticesHemisphere.Count];
                    for (int u = 0; u < vec_sky.Length; u++)
                    {
                        vec_sky[u] = new Vector3d(
                            p.sky[i].VertexCoordinatesSphere[p.sky[i].VerticesHemisphere[u]][0],
                            p.sky[i].VertexCoordinatesSphere[p.sky[i].VerticesHemisphere[u]][1],
                            p.sky[i].VertexCoordinatesSphere[p.sky[i].VerticesHemisphere[u]][2]);
                    }
                    bool[] shdw_sky = new bool[p.sky[i].VerticesHemisphere.Count];
                    cShadow.CalcShadow(orig, mshvrtnorm[i], 0.1, vec_sky, obst, ref shdw_sky);

                    ShdwSky.Add(shdw_sky);

                });
            }
            else
            {


                Parallel.For(0, mshvrt.Length, i =>
                {
                    mshvrtnorm[i] = msh.Normals[i];
                    Point3d orig = new Point3d(mshvrt[i].X, mshvrt[i].Y, mshvrt[i].Z);

                    //sky dome diffuse
                    Vector3d[] vec_sky = new Vector3d[p.sky[i].VerticesHemisphere.Count];
                    for (int u = 0; u < vec_sky.Length; u++)
                    {
                        vec_sky[u] = new Vector3d(
                            p.sky[i].VertexCoordinatesSphere[p.sky[i].VerticesHemisphere[u]][0],
                            p.sky[i].VertexCoordinatesSphere[p.sky[i].VerticesHemisphere[u]][1],
                            p.sky[i].VertexCoordinatesSphere[p.sky[i].VerticesHemisphere[u]][2]);
                    }
                    bool[] shdw_sky = new bool[p.sky[i].VerticesHemisphere.Count];
                    cShadow.CalcShadow(orig, mshvrtnorm[i], 0.1, vec_sky, obst, ref shdw_sky);

                    ShdwSky.Add(shdw_sky);

                    bool[][] shdw_beam = new bool[12][];
                    int dmcount = 1;    //days in month counter
                    for (int d = 0; d < 12; d++)
                    {
                        Vector3d[] vec_beam = new Vector3d[24];
                        int dm = System.DateTime.DaysInMonth(year, d + 1);
                        startDays[d] = dmcount;
                        endDays[d] = dm + dmcount;
                        dmcount += dm;

                        bool[] sunshine = new bool[24];
                        int HOY = (startDays[d] - 1) * 24;
                        for (int t = 0; t < 24; t++)
                        {
                            if (sunvectors[HOY + t].Sunshine)
                                sunshine[t] = true;
                            vec_beam[t] = new Vector3d(sunvectors[HOY + t].udtCoordXYZ.x, sunvectors[HOY + t].udtCoordXYZ.y, sunvectors[HOY + t].udtCoordXYZ.z);
                        }
                        shdw_beam[d] = new bool[24];
                        cShadow.CalcShadow(orig, mshvrtnorm[i], 0.1, vec_beam, sunshine, obst, ref shdw_beam[d]);

                    }
                    ShdwBeam.Add(shdw_beam);
                });
            }



           if (interpmode == 0)
              p.SetShadowsInterpolatedMT(ShdwBeam_equinox, ShdwBeam_summer, ShdwBeam_winter, ShdwSky);
           else
              p.SetShadowsInterpolatedMT(startDays, endDays, ShdwBeam, ShdwSky);


           p.SetSnowcover(snow_threshold, tilt_treshold);
            //p.SetInterreflection();
           
            
           p.CalcIrradiationMT();



            I_hourly = new Matrix(mshvrt.Length, 8760);
            Ib_hourly = new Matrix(mshvrt.Length, 8760);
            Ih_hourly = new Matrix(mshvrt.Length, 8760);

            for (int i = 0; i < mshvrt.Length; i++)
            {
                double Itot = 0;
                double Ibtot = 0;
                double Idtot = 0;
                for (int t = 0; t < 8760; t++)
                {
                    I_hourly[i, t] = p.I[i][t];
                    Ib_hourly[i, t] = p.Ibeam[i][t];
                    Ih_hourly[i, t] = p.Idiff[i][t];
                    Itot += p.I[i][t];
                    Ibtot += p.Ibeam[i][t];
                    Idtot += p.Idiff[i][t];
                }
                I.Add(Itot / 1000);     //in kWh/a
                Ib.Add(Ibtot / 1000);   //in kWh/a
                Ih.Add(Idtot / 1000);   //in kWh/a
            }
        }





    }
}
