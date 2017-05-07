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
        /// <param name="origin">3D coordinate of sensor point.</param>
        /// <param name="origNormal">Normal of sensor point.</param>
        /// <param name="tolerance">Offset from origin in normal direction, to avoid self-intersection.</param>
        /// <param name="solarvec">Array of solar vectors to check interreflections for.</param>
        /// <param name="sunshine">Array of booleans, indicating if sun shines (vector above ground plane).</param>
        /// <param name="obstacles">Array of mesh obstacles.</param>
        /// <param name="albedo">Double array [u][t] of reflection coefficients (albedos) for each obstacle u and for each vector t.</param>
        /// <param name="reflType">Array indicating reflection type for each obstacle. 0: diffuse, 1: specular.</param>
        /// <param name="bounces">Number of bounces. Max 2 recommended.</param>
        /// <param name="Ispecular">Normal irradiation values [t][m] for each solar vector t and each refleced ray m.</param>
        /// <param name="Inormals">Normal vectors [t][m] for each solar vector t and each reflected ray m.</param>
        internal static void CalcSpecularNormal(Point3d origin, Vector3d origNormal, double tolerance, Vector3d[] solarvec, bool[] sunshine,
            Mesh[] obstacles, double[][] albedo, int[] reflType, int bounces,
            ref double[][] Ispecular, ref Vector3d[][] Inormals)
        {
            //return: [0][n] double Irradiation; [1][n] vector3d normals of interreflections
            if (bounces < 1) return;

            //Rhino.RhinoDoc doc = Rhino.RhinoDoc.ActiveDoc;

            Point3d origOffset = new Point3d(Point3d.Add(origin, Vector3d.Multiply(Vector3d.Divide(origNormal, origNormal.Length), tolerance)));
            List<List<double>> IspecNorm = new List<List<double>>();    //for each solarvec, add a list with elements for each reflected ray hitting the SP
            List<List<Vector3d>> IspecVec = new List<List<Vector3d>>(); //same, just vector of bounce. thats vSrfObst

            //List<Line> ln = new List<Line>();
            //List<TextDot> dot = new List<TextDot>();



            Vector3d[][] obstNorms = new Vector3d[obstacles.Length][];
            Vector3d[][] obstNormsRev = new Vector3d[obstacles.Length][];
            Point3d[][] obstCen = new Point3d[obstacles.Length][];
            for (int u = 0; u < obstacles.Length; u++)
            {
                obstacles[u].FaceNormals.ComputeFaceNormals();
                obstNorms[u] = new Vector3d[obstacles[u].Faces.Count];
                obstNormsRev[u] = new Vector3d[obstacles[u].Faces.Count];
                obstCen[u] = new Point3d[obstacles[u].Faces.Count];
                for (int k = 0; k < obstacles[u].Faces.Count; k++)
                {
                    obstNorms[u][k] = obstacles[u].FaceNormals[k];
                    obstNormsRev[u][k] = Vector3d.Negate(obstNorms[u][k]);
                    Point3d cen0 = obstacles[u].Faces.GetFaceCenter(k);
                    obstCen[u][k] = new Point3d(Point3d.Add(cen0, Vector3d.Multiply(Vector3d.Divide(obstNorms[u][k], obstNorms[u][k].Length), tolerance)));
                }
            }


            for (int t = 0; t < solarvec.Length; t++)
            {
                IspecNorm.Add(new List<double>());  //add element for each reflected ray
                IspecVec.Add(new List<Vector3d>());
                if (sunshine[t])
                {
                    for (int u = 0; u < obstacles.Length; u++)
                    {
                        if (reflType[u] == 1 && albedo[u][t] > 0.0)   //specular, and it must be able to reflect ( > 0)
                        {
                            for (int k = 0; k < obstacles[u].Faces.Count; k++)      //specular obstacles are advised to be clean mesh geometries...
                            {
                                bool obstSun;
                                double vanglesun = Vector3d.VectorAngle(obstNorms[u][k], solarvec[t]) * (180.0 / Math.PI);
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

                                Vector3d vSrfObst = new Vector3d(Point3d.Subtract(obstCen[u][k], origOffset));
                                double vangle = Vector3d.VectorAngle(obstNormsRev[u][k], vSrfObst) * (180.0 / Math.PI);
                                if (vangle >= 90)
                                {
                                    continue;
                                }
                                Ray3d rSrfObst = new Ray3d(origOffset, vSrfObst);

                                //get all obstcles, also non specular ones....
                                bool blnX = false;
                                for (int m = 0; m < obstacles.Length; m++)
                                {
                                    if (m != u)
                                    {
                                        double inters = Rhino.Geometry.Intersect.Intersection.MeshRay(obstacles[m], rSrfObst);
                                        if (inters >= 0)
                                        {
                                            blnX = true;
                                            break;
                                        }
                                    }
                                }
                                if (!blnX)
                                {
                                    if (bounces == 1)    //support only 1 bounce for now.
                                    {
                                        Vector3d refl = Vector3d.Subtract(vSrfObst, Vector3d.Multiply(Vector3d.Multiply(vSrfObst, obstNorms[u][k]) * 2.0, obstNorms[u][k]));
                                        double tol = 360.0 * tolerance * ( Math.PI / 180.0);
                                        int par = refl.IsParallelTo(solarvec[t], tol);
                                        if (par == 1)
                                        {
                                            IspecNorm[t].Add(albedo[u][t]);  //add beam irriadiation normal to reflected ray. as a factor (albedos). don't work with DHI yet.
                                            IspecVec[t].Add(vSrfObst);
                                            //ln.Add(new Line(obstCen[u][k], refl, 100));
                                            //ln.Add(new Line(origOffset, obstCen[u][k]));
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


            //  rhino obstructions forall t. use INTERPOLATION 
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            /////////////////////////////////////////////////     S P E C U L A R   /////////////////////////////////////////////////////////////
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //SPECULAR
            //double total_Specular_Budget = 0.0;
            // forall specular objects:

            // 01:
            //  from sensorpoint (SP), shoot rays ONLY to specular objects. 
            //      - check, if vectorangle(ray, sensorpoint normal) <= 90°. if not, break. else:
            //          - check, if ray unobstructed. if not, break. else, go to 02.

            // 02:
            // if ray <= 90° (otherwise it hits the backside of the specular surface) and unobstructed:
            //      - if bounce > 1 (else go to 03):
            //          - for each bounce > 1, calc reflected ray on the specular obstacle and for new ray, go to 01:
            //              - go to 03, if unobstructed:

            // 03:
            //      - for all t in [0,8759]
            //          - check, if daytime. if not, break. else:
            //              - from specular obstacle object, make ray to sunvector. calc vectorangle(sunvector, obstacle normal). if not <= 90°, break. else:
            //                  - make reflection of sunvector on obstacle. if this reflecting vector *-1 is not coincident with connecting ray to SP (tolerance of +- 1°?), break. else:
            //                      - check if the sunvector (not the reflected) is obstructed. if yes, break. else, go to 03.

            // 03:
            // calc Specular_interreflection = I_direct incident on obstacle object. Multiply with its specular coefficient (albedo?). 
            //      and account for incidence angle on sensorpoint: I_sensorpoint = I_incident sin(90°-vectorangle(sp_normal, incident_ray))
            // total_Specular_Budget += Specular_interreflection;


            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
            double [] DNI, ref double[] IspecularIncident)
        {
            //convert Inormals vectors into solar zenith and solar azimuth. coz thats basically my sun.
            //DNI = DNI * Ispecular (here are my albedos)

            for(int t=0; t<IspecularIncident.Length;t++)
            {
                IspecularIncident[t] = 0.0;
                if (Ispecular[t] != null)
                {
                    for (int m = 0; m < Ispecular[t].Length; m++)
                    {
                        double DNI_t = DNI[t] * Ispecular[t][m];
                        IspecularIncident[t] += DNI_t * Math.Sin((90 * Math.PI / 180.0) - Vector3d.VectorAngle(origNormal, Inormals[t][m]));
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
