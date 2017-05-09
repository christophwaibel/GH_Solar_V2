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
                if (Vector3d.VectorAngle(origNormal, vec[t]) > 90)  //assumes a surface. a globe in space could of course get a ray from "behind" 
                {
                    shdw[t] = true;
                }
                else
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
                //if (shdw[t] == false)
                //{
                //    Line ln = new Line(origOffset, Vector3d.Multiply(1000, vec[t]));
                //    var attribs = Rhino.RhinoDoc.ActiveDoc.CreateDefaultAttributes();
                //    attribs.ObjectDecoration = Rhino.DocObjects.ObjectDecoration.BothArrowhead;
                //    Rhino.RhinoDoc.ActiveDoc.Objects.AddLine(ln, attribs);
                //}
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
                if (Vector3d.VectorAngle(origNormal, vec[t]) > 90)  //assumes a surface. a globe in space could of course get a ray from "behind" 
                {
                    shdw_mt[t] = true;
                }
                else
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
                    if (Vector3d.VectorAngle(origNormal, vec[t]) > 90)  //assumes a surface. a globe in space could of course get a ray from "behind" 
                    {
                        shdw[t] = true;
                    }
                    else
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
                    //if (shdw[t] == false)
                    //{
                    //    Line ln = new Line(origOffset, Vector3d.Multiply(1000, vec[t]));
                    //    var attribs = Rhino.RhinoDoc.ActiveDoc.CreateDefaultAttributes();
                    //    attribs.ObjectDecoration = Rhino.DocObjects.ObjectDecoration.BothArrowhead;
                    //    Rhino.RhinoDoc.ActiveDoc.Objects.AddLine(ln, attribs);
                    //}
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
                    if (Vector3d.VectorAngle(origNormal, vec[t]) > 90)  //assumes a surface. a globe in space could of course get a ray from "behind" 
                    {
                        shdw_mt[t] = true;
                    }
                    else
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
                }
                else
                {
                    shdw_mt[t] = true;
                }
            });

            shdw_mt.CopyTo(shdw, 0);
        }


        //mesh operations are not multi-threading safe?
        /*
        internal static void CalcShadowMeshMT(int year, List<SunVector> sunvectors, Mesh msh, Mesh[] obst, int interpmode,
            Context.cWeatherdata weather, Context.cLocation location, int rec,
            ref List<double> I, ref List<double> Ih, ref List<double> Ib,
            ref Matrix I_hourly, ref Matrix Ih_hourly, ref Matrix Ib_hourly)
        {
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

        */



        /// <summary>
        /// Calculates specular interreflections on one sensor point for an array of solar vectors. Irradiation values given normal to reflected ray.
        /// </summary>
        /// <remarks>
        /// Approach 1: Backtrace rays from sensor point, directed to each specular object.
        /// </remarks>
        /// <param name="origin">3D coordinate of sensor point.</param>
        /// <param name="origNormal">Normal of sensor point.</param>
        /// <param name="tolerance">Offset from origin in normal direction, to avoid self-intersection.</param>
        /// <param name="solarvec">Array of solar vectors to check interreflections for.</param>
        /// <param name="sunshine">Array of booleans, indicating if sun shines (vector above ground plane).</param>
        /// <param name="obstacles">Array of mesh obstacles.</param>
        /// <param name="albedo">Double array [u][t] of reflection coefficients (albedos) for each obstacle u and for each vector t.</param>
        /// <param name="reflType">Reflection type for each obstacle. 0: diffuse, 1: specular.</param>
        /// <param name="bounces">Number of bounces. Currently only max 1.</param>
        /// <param name="Ispecular">Normal irradiation values [t][m] for each solar vector t and each reflected ray m.</param>
        /// <param name="Inormals">Normal vectors [t][m] for each solar vector t and each reflected ray m.</param>
        internal static void CalcSpecularNormal1(Point3d origin, Vector3d origNormal, double tolerance, Vector3d[] solarvec, bool[] sunshine,
            List<ObstacleObject> obstacles, double[][] albedo, int[] reflType, int bounces,
            ref double[][] Ispecular, ref Vector3d[][] Inormals)
        {
            //Rhino.RhinoDoc doc = Rhino.RhinoDoc.ActiveDoc;



            //return: [0][n] double Irradiation; [1][n] vector3d normals of interreflections
            if (bounces < 1) return;

            Point3d origOffset = new Point3d(Point3d.Add(origin, Vector3d.Multiply(Vector3d.Divide(origNormal, origNormal.Length), tolerance)));
            List<List<double>> IspecNorm = new List<List<double>>();    //for each solarvec, add a list with elements for each reflected ray hitting the SP
            List<List<Vector3d>> IspecVec = new List<List<Vector3d>>(); //same, just vector of bounce. thats vSrfObst

            List<Line> ln = new List<Line>();
            List<TextDot> dot = new List<TextDot>();


            for (int t = 0; t < solarvec.Length; t++)
            {
                IspecNorm.Add(new List<double>());  //add element for each reflected ray
                IspecVec.Add(new List<Vector3d>());
                if (sunshine[t])
                {
                    for (int u = 0; u < obstacles.Count; u++)
                    {
                        if (reflType[u] == 1 && albedo[u][t] > 0.0)   //specular, and it must be able to reflect ( > 0)
                        {
                            for (int k = 0; k < obstacles[u].mesh.Faces.Count; k++)      //specular obstacles are advised to be clean mesh geometries...
                            {
                                bool obstSun;
                                double vanglesun = Vector3d.VectorAngle(obstacles[u].normals[k], solarvec[t]) * (180.0 / Math.PI);
                                if (vanglesun >= 90)
                                {
                                    obstSun = false;
                                }
                                else
                                {
                                    obstSun = true;
                                }
                                if (bounces == 1)   //sun doesn't reach reflector. 
                                {
                                    if (!obstSun)
                                    {
                                        continue;
                                    }
                                }

                                Vector3d vSrfObst = new Vector3d(Point3d.Subtract(obstacles[u].faceCen[k], origOffset));
                                double vangle = Vector3d.VectorAngle(obstacles[u].normalsRev[k], vSrfObst) * (180.0 / Math.PI);
                                if (vangle >= 90)
                                {
                                    continue;
                                }
                                Ray3d rSrfObst = new Ray3d(origOffset, vSrfObst);
                                Line lSrfObst = new Line(origOffset, obstacles[u].faceCen[k]);
                                //get all obstcles, also non specular ones....
                                bool blnX = false;
                                for (int m = 0; m < obstacles.Count; m++)
                                {
                                    //double inters = Rhino.Geometry.Intersect.Intersection.MeshRay(obstacles[m].mesh, rSrfObst);
                                    int[] f;
                                    Point3d[] inters2 = Rhino.Geometry.Intersect.Intersection.MeshLine(obstacles[m].mesh, lSrfObst, out f);

                                    if (inters2 != null && inters2.Length > 0)
                                    {
                                        blnX = true;
                                        break;
                                    }
                                }
                                if (!blnX)
                                {
                                    if (bounces == 1)    //support only 1 bounce for now.
                                    {
                                        Vector3d refl = Vector3d.Subtract(vSrfObst, Vector3d.Multiply(Vector3d.Multiply(vSrfObst, obstacles[u].normals[k]) * 2.0, obstacles[u].normals[k]));
                                        double tol = 0.1 * (Math.PI / 180.0);


                                        //ln.Add(new Line(obstacles[u].faceCen[k], refl, 100));
                                        //ln.Add(new Line(origOffset, obstacles[u].faceCen[k]));

                                        int par = refl.IsParallelTo(solarvec[t], tol);
                                        if (par == 1)
                                        {
                                            IspecNorm[t].Add(albedo[u][t]);  //add beam irriadiation normal to reflected ray. as a factor (albedos). don't work with DHI yet.
                                            IspecVec[t].Add(Vector3d.Negate(vSrfObst));
                                            //ln.Add(new Line(obstacles[u].faceCen[k], refl, 100));
                                            //ln.Add(new Line(origOffset, obstacles[u].faceCen[k]));
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                Ispecular[t] = IspecNorm[t].ToArray();
                Inormals[t] = IspecVec[t].ToArray();
            }

            //foreach (TextDot d in dot)
            //{
            //    doc.Objects.AddTextDot(d);
            //}
            //foreach (Line l in ln)
            //{
            //    doc.Objects.AddLine(l);
            //}
            //doc.Views.Redraw();

        }


        /// <summary>
        /// Calculates specular interreflections on all sensor points for an array of solar vectors. Irradiation values given normal to reflected ray.
        /// </summary>
        /// <remarks>
        /// Approach 2: Backtrace rays from each specular object. Solar vector hits specular object, is reflected, and this reflected ray is traced.
        /// </remarks>
        /// <param name="SPmesh">Analysis mesh, on which the sensor points are placed.</param>
        /// <param name="SP">Sensor points [i] of the analysis mesh.</param>
        /// <param name="SPnormal">Normal vectors [i] of the i sensor points.</param>
        /// <param name="solarvec">Solar vectors for each time step t.</param>
        /// <param name="sunshine">Indicating sunshine for each respective solar vector.</param>
        /// <param name="obstacles">Obstacle objects.</param>
        /// <param name="albedo">Reflection coefficients [u][t] for each obstacle u and at each time step t.</param>
        /// <param name="refltype">Reflection type for each obstacle. 0: diffuse, 1: specular.</param>
        /// <param name="bounces">Number of bounces. Max 2 recommended.</param>
        /// <param name="Ispecular">Normal irradiation values [i][t][m] for each sensor point i, each solar vector t and each reflected ray m.</param>
        /// <param name="Inormals">Normal vectors [i][t][m] for each sensor point i, each solar vector t and each reflected ray m.</param>
        internal static void CalcSpecularNormal2(ObstacleObject SPmesh, Point3d[] SP, Vector3d[] SPnormal,
            Vector3d[] solarvec, bool[] sunshine,
            List<ObstacleObject> obstacles, double[][] albedo, int[] refltype, int bounces,
            ref double[][][] Ispecular, ref Vector3d[][][] Inormals)
        {
            //Rhino.RhinoDoc doc = Rhino.RhinoDoc.ActiveDoc;

            if (bounces < 1) return;


            List<List<List<double>>> IspecList = new List<List<List<double>>>();
            List<List<List<Vector3d>>> InormList = new List<List<List<Vector3d>>>();
            for (int i = 0; i < SP.Length; i++)
            {
                Ispecular[i] = new double[solarvec.Length][];
                Inormals[i] = new Vector3d[solarvec.Length][];

                IspecList.Add(new List<List<double>>());
                InormList.Add(new List<List<Vector3d>>());
                for (int t = 0; t < solarvec.Length; t++)
                {
                    IspecList[i].Add(new List<double>());
                    InormList[i].Add(new List<Vector3d>());
                }
            }


            //foreach specular object                 
            for (int u = 0; u < obstacles.Count; u++)
            {
                if (refltype[u] == 1)
                {
                    //foreach solar vector t
                    for (int t = 0; t < solarvec.Length; t++)
                    {
                        //foreach face
                        for (int k = 0; k < obstacles[u].faceCen.Length; k++)
                        {
                            //hit with solar vec
                            Ray3d rObstSun = new Ray3d(obstacles[u].faceCen[k], solarvec[t]);

                            //check, solar vec hits front of reflector
                            double vAngle = Vector3d.VectorAngle(obstacles[u].normals[k], solarvec[t]) * (180.0 / Math.PI);

                            //if yes, make reflected ray and trace it
                            if (vAngle < 90)
                            {
                                //first check, if solar vec reaches reflector
                                bool xFirst = false;
                                for (int n = 0; n < obstacles.Count; n++)
                                {
                                    double intersFirst = Rhino.Geometry.Intersect.Intersection.MeshRay(obstacles[n].mesh, rObstSun);
                                    if (intersFirst >= 0)
                                    {
                                        xFirst = true;
                                        break;
                                    }
                                }
                                if (!xFirst)
                                {
                                    Vector3d refl = Vector3d.Subtract(Vector3d.Negate(solarvec[t]), Vector3d.Multiply(Vector3d.Multiply(Vector3d.Negate(solarvec[t]), obstacles[u].normals[k]) * 2.0, obstacles[u].normals[k]));
                                    //doc.Objects.AddLine(new Line(obstacles[u].faceCen[k], refl, 1000));
                                    //doc.Views.ActiveView.Redraw();
                                    Ray3d rObstRefl = new Ray3d(obstacles[u].faceCen[k], refl);

                                    //if it hits the SPmesh, assign the beam to the closest SP
                                    double inters = 0;
                                    bool reflX = false;
                                    int intersP_index = 0;
                                    for (int n = 0; n < obstacles.Count; n++)
                                    {
                                        //!!!!!!!!!!!!!!!!!!!!!
                                        //take out SPmesh!!!
                                        if (!String.Equals(SPmesh.name, obstacles[n].name))
                                        {
                                            inters = Rhino.Geometry.Intersect.Intersection.MeshRay(obstacles[n].mesh, rObstRefl);
                                            if (inters >= 0)
                                            {
                                                intersP_index = n;
                                                reflX = true;
                                                break;
                                            }
                                        }
                                    }
                                    if (!reflX)
                                    {
                                        double inters2 = Rhino.Geometry.Intersect.Intersection.MeshRay(SPmesh.mesh, rObstRefl);
                                        if (inters2 >= 0)
                                        {
                                            Point3d p = rObstRefl.PointAt(inters2);
                                            double min_d = double.MaxValue;
                                            int i_index = 0;
                                            for (int i = 0; i < SP.Length; i++)
                                            {
                                                double d = p.DistanceTo(SP[i]);
                                                if (d < min_d)
                                                {
                                                    i_index = i;
                                                    min_d = d;
                                                }
                                            }
                                            IspecList[i_index][t].Add(albedo[u][t]);
                                            InormList[i_index][t].Add(refl);
                                            //doc.Objects.AddLine(new Line(obstacles[u].faceCen[k], solarvec[t], 1000));
                                            //doc.Objects.AddLine(new Line(obstacles[u].faceCen[k], p));
                                            //doc.Views.ActiveView.Redraw();
                                        }
                                    }
                                    else
                                    {
                                        //if it hits another specular object, and bounce = 2, reflect again and trace that.
                                        if (bounces > 1 && refltype[intersP_index] == 1)
                                        {
                                            Point3d px = rObstRefl.PointAt(inters);
                                            MeshPoint mshp = obstacles[intersP_index].mesh.ClosestMeshPoint(px, 0.0);
                                            Vector3d pnormal = obstacles[intersP_index].mesh.NormalAt(mshp);
                                            Vector3d refl2 = Vector3d.Subtract(refl, Vector3d.Multiply(Vector3d.Multiply(refl, pnormal) * 2.0, pnormal));
                                            Point3d pxoffset = new Point3d(Point3d.Add(px, Vector3d.Multiply(Vector3d.Divide(pnormal, pnormal.Length), obstacles[u].tolerance)));

                                            //check 2nd bounce for obstruction and if not, if it hits the SPmsh
                                            Ray3d rObstRefl2 = new Ray3d(pxoffset, refl2);
                                            double intersBounce2 = 0;
                                            bool reflX2 = false;
                                            //int intersP2_index = 0;
                                            for (int n = 0; n < obstacles.Count; n++)
                                            {
                                                //!!!!!!!!!!!!!!!!!!!!!
                                                //take out SPmesh!!!
                                                if (!String.Equals(SPmesh.name, obstacles[n].name))
                                                {
                                                    intersBounce2 = Rhino.Geometry.Intersect.Intersection.MeshRay(obstacles[n].mesh, rObstRefl2);
                                                    if (intersBounce2 >= 0)
                                                    {
                                                        //intersP2_index = n;
                                                        reflX2 = true;
                                                        break;
                                                    }
                                                }
                                            }
                                            if (!reflX2)
                                            {
                                                double inters2 = Rhino.Geometry.Intersect.Intersection.MeshRay(SPmesh.mesh, rObstRefl2);
                                                if (inters2 >= 0)
                                                {
                                                    Point3d p = rObstRefl2.PointAt(inters2);
                                                    double min_d = double.MaxValue;
                                                    int i_index = 0;
                                                    for (int i = 0; i < SP.Length; i++)
                                                    {
                                                        double d = p.DistanceTo(SP[i]);
                                                        if (d < min_d)
                                                        {
                                                            i_index = i;
                                                            min_d = d;
                                                        }
                                                    }
                                                    IspecList[i_index][t].Add(albedo[u][t]);
                                                    InormList[i_index][t].Add(refl2);

                                                    //doc.Objects.AddLine(new Line(obstacles[u].faceCen[k], solarvec[t], 1000));
                                                    //doc.Objects.AddLine(new Line(obstacles[u].faceCen[k], pxoffset));
                                                    //doc.Objects.AddLine(new Line(pxoffset, p));
                                                    //doc.Views.ActiveView.Redraw();
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }



            for (int i = 0; i < SP.Length; i++)
            {
                for (int t = 0; t < solarvec.Length; t++)
                {
                    Ispecular[i][t] = IspecList[i][t].ToArray();
                    Inormals[i][t] = InormList[i][t].ToArray();
                }
            }
        }


        /// <summary>
        /// Calculates specular interreflections on all sensor points for an array of solar vectors. Irradiation values given normal to reflected ray.
        /// </summary>
        /// <remarks>
        /// Approach 3: Compute rays and reflected rays for each obstacle. Translate them to each sensor point. Check if translated rays reach obstacle and are unobstructed. 1st and 2nd order (1 or 2 bounces).
        /// </remarks>
        /// <param name="SPmesh">Analysis mesh, on which the sensor points are placed.</param>
        /// <param name="SP">Sensor points [i] of the analysis mesh.</param>
        /// <param name="SPnormal">Normal vectors [i] of the i sensor points.</param>
        /// <param name="solarvec">Solar vectors for each time step t.</param>
        /// <param name="sunshine">Indicating sunshine for each respective solar vector.</param>
        /// <param name="obstacles">Obstacle objects.</param>
        /// <param name="albedo">Reflection coefficients [u][t] for each obstacle u and at each time step t.</param>
        /// <param name="refltype">Reflection type for each obstacle. 0: diffuse, 1: specular.</param>
        /// <param name="bounces">Number of bounces. Max 2 recommended.</param>
        /// <param name="Ispecular">Normal irradiation values [i][t][m] for each sensor point i, each solar vector t and each reflected ray m.</param>
        /// <param name="Inormals">Normal vectors [i][t][m] for each sensor point i, each solar vector t and each reflected ray m.</param>
        internal static void CalcSpecularNormal3(ObstacleObject SPmesh, Point3d[] SP, Vector3d[] SPnormal,
            Vector3d[] solarvec, bool[] sunshine,
            List<ObstacleObject> obstacles, double[][] albedo, int[] refltype, int bounces,
            ref double[][][] Ispecular, ref Vector3d[][][] Inormals)
        {
            Rhino.RhinoDoc doc = Rhino.RhinoDoc.ActiveDoc;

            // Approach 3: Compute rays and reflected rays for each obstacle. Translate them to each sensor point. Check if translated rays reach obstacle and are unobstructed. 1st and 2nd order (1 or 2 bounces).
            if (bounces < 1) return;

            List<List<List<double>>> IspecList = new List<List<List<double>>>();
            List<List<List<Vector3d>>> InormList = new List<List<List<Vector3d>>>();
            Point3d[] SPoffset = new Point3d[SP.Length];
            for (int i = 0; i < SP.Length; i++)
            {
                SPoffset[i] = new Point3d(Point3d.Add(SP[i], Vector3d.Multiply(Vector3d.Divide(SPnormal[i], SPnormal[i].Length), obstacles[0].tolerance)));

                Ispecular[i] = new double[solarvec.Length][];
                Inormals[i] = new Vector3d[solarvec.Length][];

                IspecList.Add(new List<List<double>>());
                InormList.Add(new List<List<Vector3d>>());
                for (int t = 0; t < solarvec.Length; t++)
                {
                    IspecList[i].Add(new List<double>());
                    InormList[i].Add(new List<Vector3d>());
                }
            }

            //foreach solar vector t
            for (int t = 0; t < solarvec.Length; t++)
            {
                if (!sunshine[t]) continue;

                //foreach specular object                 
                for (int u = 0; u < obstacles.Count; u++)
                {
                    if (refltype[u] != 1) continue;

                    //foreach face
                    for (int k = 0; k < obstacles[u].faceCen.Length; k++)
                    {
                        //check, if solar vec hits front of reflector
                        //if yes, make reflected ray and trace it. if not, go to next iteration.
                        double vAngle = Vector3d.VectorAngle(obstacles[u].normals[k], solarvec[t]) * (180.0 / Math.PI);
                        if (vAngle >= 90) continue;

                        //compute reflected ray
                        Vector3d solarvec_ref = cMisc.ReflectVec(obstacles[u].normals[k], Vector3d.Negate(solarvec[t]));
                        //doc.Objects.AddLine(new Line(obstacles[u].faceCen[k], solarvec[t], 1000));
                        //doc.Objects.AddLine(new Line(obstacles[u].faceCen[k], solarvec_ref, 1000));

                        if (bounces == 1)
                        {
                            //translate two rays onto each SP and check for obstructions
                            for (int i = 0; i < SP.Length; i++)
                            {
                                vAngle = Vector3d.VectorAngle(SPnormal[i], Vector3d.Negate(solarvec_ref)) * (180.0 / Math.PI);
                                if (vAngle >= 90) continue;

                                //check, if translated ray hits initial obstacle face
                                Ray3d rObstSun_refl = new Ray3d(SP[i], Vector3d.Negate(solarvec_ref));
                                Mesh mshface = new Mesh();
                                mshface.Vertices.Add(obstacles[u].mesh.Vertices[obstacles[u].mesh.Faces.GetFace(k).A]);
                                mshface.Vertices.Add(obstacles[u].mesh.Vertices[obstacles[u].mesh.Faces.GetFace(k).B]);
                                mshface.Vertices.Add(obstacles[u].mesh.Vertices[obstacles[u].mesh.Faces.GetFace(k).C]);
                                if (obstacles[u].mesh.Faces.GetFace(k).IsQuad)
                                {
                                    mshface.Vertices.Add(obstacles[u].mesh.Vertices[obstacles[u].mesh.Faces.GetFace(k).D]);
                                    mshface.Faces.AddFace(0, 1, 2, 3);
                                }
                                else
                                {
                                    mshface.Faces.AddFace(0, 1, 2);
                                }
                                double inters = Rhino.Geometry.Intersect.Intersection.MeshRay(mshface, rObstSun_refl);
                                if (inters < 0) continue;

                                Point3d pX = rObstSun_refl.PointAt(inters);
                                MeshPoint mshpX = obstacles[u].mesh.ClosestMeshPoint(pX, 0.0);
                                Vector3d pnormal = obstacles[u].mesh.NormalAt(mshpX);
                                Point3d pXoff = new Point3d(Point3d.Add(pX, Vector3d.Multiply(Vector3d.Divide(pnormal, pnormal.Length), obstacles[u].tolerance)));
                                Ray3d rObstSun = new Ray3d(pXoff, solarvec[t]);
                                //doc.Objects.AddLine(new Line(pX, SP[i]));

                                //check obstacle to sun and obstacle to SP
                                double inters_a;
                                bool bln_inters = false;
                                for (int n = 0; n < obstacles.Count; n++)
                                {
                                    inters_a = Rhino.Geometry.Intersect.Intersection.MeshRay(obstacles[n].mesh, rObstSun);
                                    if (inters_a >= 0)
                                    {
                                        bln_inters = true;
                                        break;
                                    }
                                    int[] f;
                                    Point3d[] inters_b = Rhino.Geometry.Intersect.Intersection.MeshLine(obstacles[n].mesh, new Line(pXoff, SPoffset[i]), out f);
                                    //doc.Objects.AddLine(new Line(pXoff, SP[i]));

                                    if (inters_b != null && inters_b.Length > 0)
                                    {
                                        bln_inters = true;
                                        break;
                                    }
                                }
                                if (bln_inters) continue;  //obstructed

                                //no obstruction. add vector and albedo
                                IspecList[i][t].Add(albedo[u][t]);
                                InormList[i][t].Add(solarvec_ref);
                                //doc.Objects.AddLine(new Line(pXoff, solarvec[t], 1000));
                                //doc.Objects.AddLine(new Line(pX, SP[i]));
                                //doc.Views.ActiveView.Redraw();
                            }
                        }
                        else if (bounces > 1)
                        {
                            //for 2 bounces, check first, if reflected ray hits another obstacle and if so, bounce off from that
                            Ray3d rObstSun_refl = new Ray3d(obstacles[u].faceCen[k], solarvec_ref);
                            Ray3d rObstSun = new Ray3d(obstacles[u].faceCen[k], solarvec[t]);

                            //check, if it reaches sun from face center
                            double inters = 0.0;
                            bool bln_inters = false;
                            for (int n = 0; n < obstacles.Count; n++)
                            {
                                inters = Rhino.Geometry.Intersect.Intersection.MeshRay(obstacles[n].mesh, rObstSun);
                                if (inters >= 0)
                                {
                                    bln_inters = true;
                                    break;
                                }
                            }
                            if (bln_inters) continue;

                            // check, if it reaches 2nd obstacle from face center.... !!actually should do other way round and loop through all vertices 
                            int refl_index = 0;
                            bln_inters = false;
                            for (int n = 0; n < obstacles.Count; n++)
                            {
                                inters = Rhino.Geometry.Intersect.Intersection.MeshRay(obstacles[n].mesh, rObstSun_refl);
                                if (inters >= 0)
                                {
                                    bln_inters = true;
                                    refl_index = n;
                                    break;
                                }
                            }
                            //its blocked and obstacle is also not reflector and its also not the analysis surface. so stop.
                            if (bln_inters && refltype[refl_index] == 0 && !String.Equals(SPmesh.name, obstacles[refl_index].name)) continue;  

                            //its blocked and obstacle is reflector. 2nd bounce
                            if (bln_inters && refltype[refl_index] == 1 && !String.Equals(SPmesh.name, obstacles[refl_index].name))        
                            {
                                Point3d pX = rObstSun_refl.PointAt(inters);
                                MeshPoint mshp = obstacles[refl_index].mesh.ClosestMeshPoint(pX, 0.0);
                                Vector3d pXnormal = obstacles[refl_index].mesh.NormalAt(mshp);
                                vAngle = Vector3d.VectorAngle(pXnormal, Vector3d.Negate(solarvec_ref)) * (180.0 / Math.PI);
                                if (vAngle >= 90.0) continue;

                                Vector3d solarvec_refl2nd = cMisc.ReflectVec(pXnormal, solarvec_ref);
                                Point3d pXoffset = cMisc.OffsetPt(pX, pXnormal, obstacles[refl_index].tolerance);
                                Ray3d rObstSun_refl2nd = new Ray3d(pXoffset, solarvec_refl2nd);

                                //doc.Objects.AddLine(new Line(obstacles[u].faceCen[k], solarvec[t], 100));
                                //doc.Objects.AddLine(new Line(obstacles[u].faceCen[k], pX));
                                //doc.Objects.AddLine(new Line(pXoffset, solarvec_refl2nd, 100));

                                for (int i = 0; i < SP.Length; i++)
                                {
                                    vAngle = Vector3d.VectorAngle(SPnormal[i], Vector3d.Negate(solarvec_refl2nd)) * (180.0 / Math.PI);
                                    if (vAngle >= 90.0) continue;

                                    //check, if translated ray hits 2nd reflected obstacle
                                    Ray3d r2ndReflSP = new Ray3d(SP[i], Vector3d.Negate(solarvec_refl2nd));
                                    double inters2nd = Rhino.Geometry.Intersect.Intersection.MeshRay(obstacles[refl_index].mesh, r2ndReflSP);
                                    if (inters2nd < 0) continue;
                                    //doc.Objects.AddLine(new Line(SP[i], Vector3d.Negate(solarvec_refl2nd), 1000));

                                    //check, if from there, it hits initial obstacle face
                                    Mesh mshface = new Mesh();
                                    mshface.Vertices.Add(obstacles[u].mesh.Vertices[obstacles[u].mesh.Faces.GetFace(k).A]);
                                    mshface.Vertices.Add(obstacles[u].mesh.Vertices[obstacles[u].mesh.Faces.GetFace(k).B]);
                                    mshface.Vertices.Add(obstacles[u].mesh.Vertices[obstacles[u].mesh.Faces.GetFace(k).C]);
                                    if (obstacles[u].mesh.Faces.GetFace(k).IsQuad)
                                    {
                                        mshface.Vertices.Add(obstacles[u].mesh.Vertices[obstacles[u].mesh.Faces.GetFace(k).D]);
                                        mshface.Faces.AddFace(0, 1, 2, 3);
                                    }
                                    else
                                    {
                                        mshface.Faces.AddFace(0, 1, 2);
                                    }

                                    Point3d pX2nd = r2ndReflSP.PointAt(inters2nd);
                                    MeshPoint mshpX2nd = obstacles[refl_index].mesh.ClosestMeshPoint(pX2nd, 0.0);
                                    Vector3d pX2ndnormal = obstacles[refl_index].mesh.NormalAt(mshpX2nd);
                                    Point3d pX2ndOffset = cMisc.OffsetPt(pX2nd, pX2ndnormal, obstacles[refl_index].tolerance);
                                    Ray3d r1stRefl2ndRefl = new Ray3d(pX2ndOffset, Vector3d.Negate(solarvec_ref));
                                    double inters1st = Rhino.Geometry.Intersect.Intersection.MeshRay(mshface, r1stRefl2ndRefl);
                                    if (inters1st < 0) continue;
                                    //doc.Objects.AddLine(new Line(pX2ndOffset, Vector3d.Negate(solarvec_ref),100));

                                    //yes, it hits the initial face. from there, does it reach the sky?
                                    Point3d pX1st = r1stRefl2ndRefl.PointAt(inters1st);
                                    MeshPoint mshpX1st = mshface.ClosestMeshPoint(pX1st, 0.0);
                                    mshface.Normals.ComputeNormals();
                                    Vector3d pXnormal1st = mshface.NormalAt(mshpX1st);
                                    //MeshPoint mshpX1st = obstacles[u].mesh.ClosestMeshPoint(pX1st, 0.0);
                                    //Vector3d pXnormal1st = obstacles[u].mesh.NormalAt(mshpX1st);
                                    Point3d pXoff1st = cMisc.OffsetPt(pX1st,pXnormal1st, obstacles[u].tolerance);
                                    Ray3d rObstSun1st = new Ray3d(pXoff1st, solarvec[t]);
                                    //doc.Objects.AddLine(new Line(pX2ndOffset, pX1st));
                                    //doc.Objects.AddLine(new Line(pX2ndOffset, solarvec[t]));

                                    //check obstacle to sun and obstacle to SP
                                    double inters_a;
                                    bool bln_inters_1st = false;
                                    for (int n = 0; n < obstacles.Count; n++)
                                    {
                                        inters_a = Rhino.Geometry.Intersect.Intersection.MeshRay(obstacles[n].mesh, rObstSun1st);
                                        if (inters_a >= 0)
                                        {
                                            bln_inters_1st = true;
                                            break;
                                        }
                                        int[] f;
                                        Point3d[] inters_b = Rhino.Geometry.Intersect.Intersection.MeshLine(obstacles[n].mesh, new Line(pXoff1st, pX2ndOffset), out f);
                                        //doc.Objects.AddLine(new Line(pXoff, SP[i]));

                                        if (inters_b != null && inters_b.Length > 0)
                                        {
                                            bln_inters_1st = true;
                                            break;
                                        }
                                    }
                                    if (bln_inters_1st) continue;  //obstructed

                                    //no obstruction. add vector and albedo
                                    IspecList[i][t].Add(albedo[u][t] * albedo[refl_index][t]);
                                    InormList[i][t].Add(solarvec_refl2nd);
                                    //doc.Objects.AddLine(new Line(pXoff1st, solarvec[t], 1000));
                                    //doc.Objects.AddLine(new Line(pX2ndOffset, pXoff1st));
                                    //doc.Objects.AddLine(new Line(pX2ndOffset, SP[i]));
                                    //doc.Views.ActiveView.Redraw();
                                }
                            }
                                //reflected ray either reaches analysis surface, or nothing. either way, translate refletec ray to each SP
                            else
                            {
                                //translate two rays onto each SP and check for obstructions
                                for (int i = 0; i < SP.Length; i++)
                                {
                                    vAngle = Vector3d.VectorAngle(SPnormal[i], Vector3d.Negate(solarvec_ref)) * (180.0 / Math.PI);
                                    if (vAngle >= 90) continue;

                                    //check, if translated ray hits initial obstacle face
                                    Ray3d rObstSun_refl_1b = new Ray3d(SP[i], Vector3d.Negate(solarvec_ref));
                                    Mesh mshface = new Mesh();
                                    mshface.Vertices.Add(obstacles[u].mesh.Vertices[obstacles[u].mesh.Faces.GetFace(k).A]);
                                    mshface.Vertices.Add(obstacles[u].mesh.Vertices[obstacles[u].mesh.Faces.GetFace(k).B]);
                                    mshface.Vertices.Add(obstacles[u].mesh.Vertices[obstacles[u].mesh.Faces.GetFace(k).C]);
                                    if (obstacles[u].mesh.Faces.GetFace(k).IsQuad)
                                    {
                                        mshface.Vertices.Add(obstacles[u].mesh.Vertices[obstacles[u].mesh.Faces.GetFace(k).D]);
                                        mshface.Faces.AddFace(0, 1, 2, 3);
                                    }
                                    else
                                    {
                                        mshface.Faces.AddFace(0, 1, 2);
                                    }
                                    double inters_1b = Rhino.Geometry.Intersect.Intersection.MeshRay(mshface, rObstSun_refl_1b);
                                    if (inters_1b < 0) continue;

                                    Point3d pX = rObstSun_refl_1b.PointAt(inters_1b);
                                    MeshPoint mshpX = obstacles[u].mesh.ClosestMeshPoint(pX, 0.0);
                                    Vector3d pnormal = obstacles[u].mesh.NormalAt(mshpX);
                                    Point3d pXoff = new Point3d(Point3d.Add(pX, Vector3d.Multiply(Vector3d.Divide(pnormal, pnormal.Length), obstacles[u].tolerance)));
                                    Ray3d rObstSun_1b = new Ray3d(pXoff, solarvec[t]);
                                    //doc.Objects.AddLine(new Line(pX, SP[i]));

                                    //check obstacle to sun and obstacle to SP
                                    double inters_a;
                                    bool bln_inters_1b = false;
                                    for (int n = 0; n < obstacles.Count; n++)
                                    {
                                        inters_a = Rhino.Geometry.Intersect.Intersection.MeshRay(obstacles[n].mesh, rObstSun_1b);
                                        if (inters_a >= 0)
                                        {
                                            bln_inters_1b = true;
                                            break;
                                        }
                                        int[] f;
                                        Point3d[] inters_b = Rhino.Geometry.Intersect.Intersection.MeshLine(obstacles[n].mesh, new Line(pXoff, SPoffset[i]), out f);
                                        //doc.Objects.AddLine(new Line(pXoff, SP[i]));

                                        if (inters_b != null && inters_b.Length > 0)
                                        {
                                            bln_inters_1b = true;
                                            break;
                                        }
                                    }
                                    if (bln_inters_1b) continue;  //obstructed

                                    //no obstruction. add vector and albedo
                                    IspecList[i][t].Add(albedo[u][t]);
                                    InormList[i][t].Add(solarvec_ref);
                                    //doc.Objects.AddLine(new Line(pXoff, solarvec[t], 1000));
                                    //doc.Objects.AddLine(new Line(pX, SP[i]));
                                    //doc.Views.ActiveView.Redraw();
                                }
                            }
                        }
                    }
                }
            }



            for (int i = 0; i < SP.Length; i++)
            {
                for (int t = 0; t < solarvec.Length; t++)
                {
                    Ispecular[i][t] = IspecList[i][t].ToArray();
                    Inormals[i][t] = InormList[i][t].ToArray();
                }
            }

        }





        /// <summary>
        /// Calculate beam irradiation incident on a sensor point, considering incidence angles.
        /// </summary>
        /// <param name="origNormal">Normal of sensor point.</param>
        /// <param name="Ispecular">Specular irradiation values [t][i] for each time stept t and for each reflected ray.</param>
        /// <param name="Inormals">Normals for each specular irradiation value [t][i].</param>
        /// <param name="HOY">Hours of the year ∈ [0, 8759] corresponding to each time step t.</param>
        /// <param name="DNI">Direct normal irradiation values for each time step t.</param>
        /// <param name="IspecularIncident">Effective specular reflected irradiation [t] incident on the sensor point, for each time step t.</param>
        internal static void CalcSpecularIncident(Vector3d origNormal, double[][] Ispecular, Vector3d[][] Inormals,
            double[] DNI, ref double[] IspecularIncident)
        {
            //convert Inormals vectors into solar zenith and solar azimuth. coz thats basically my sun.
            //DNI = DNI * Ispecular (here are my albedos)
            if (Misc.IsNullOrEmpty(Ispecular)) return;

            for (int t = 0; t < IspecularIncident.Length; t++)
            {
                IspecularIncident[t] = 0.0;

                if (Misc.IsNullOrEmpty(Ispecular[t]) == false)
                {
                    for (int m = 0; m < Ispecular[t].Length; m++)
                    {
                        double DNI_t = DNI[t] * Ispecular[t][m];
                        IspecularIncident[t] += DNI_t * Math.Sin((90 * Math.PI / 180.0) - Vector3d.VectorAngle(origNormal, Vector3d.Negate(Inormals[t][m])));
                    }
                }
            }
        }



        internal static void CalcDiffuse(Point3d origin, Vector3d origNormal, double tolerance,
            Mesh[] obstacles, double[][] albedo, int difDomeRes, ref double Idiffuse)
        {

            //    // rhino obstructions only once.
            //    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //    /////////////////////////////////////////////////     D I F F U S E     /////////////////////////////////////////////////////////////
            //    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //    //DIFFUSE
            //    // double average_Diffuse_Budget = 0;
            //    // set diffuse hedgehog ray count. using icosahedron vertices, but not those at the "horizon"/perimeter. too flat anyway.
            //    // the more rays, the more precise the average value will be.

            //    // 01:
            //    // for all hedgehog rays:
            //    //      - check, if it hits a diffuse obstacle object. if not, break. else, go to 02:

            //    // 02:
            //    //      - calc I_obstacle = total irradiation on this obstacle. multiply with its diffuse reflection (albedo?) coefficient.
            //    //      - average_Diffuse_Budget += I_obstacle

            //    // 03:
            //    // average_Diffuse_Budget /= hedgehog_ray_Count


            //    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        }


    }
}
