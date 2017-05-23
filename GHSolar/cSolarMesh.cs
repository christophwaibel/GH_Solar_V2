using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using Rhino.Geometry;
using SolarModel;

/*
 * cSolarMeshHour.cs
 * Copyright 2017 Christoph Waibel <chwaibel@student.ethz.ch>
 * 
 * This work is licensed under the GNU GPL license version 3.
*/


namespace GHSolar
{
    internal class cCalculateSolarMesh
    {
        private cObstacleObject mshobj;
        private List<cObstacleObject> objObst = new List<cObstacleObject>();
        private List<cSemiPermObject> objTrees = new List<cSemiPermObject>();
        private double latitude;
        private double longitude;
        private List<double> DNI = new List<double>();
        private List<double> DHI = new List<double>();
        private List<double> SNOW = new List<double>();
        private int year;
        private int month;
        private int day;
        private int hour;
        private int bounces;
        private int rec;
        private int diffRes;
        private cResultsInterreflections ResultsIreflIn;
        private bool mt;
        double snow_threshold;
        double tilt_treshold;

        private const double rad = Math.PI / 180;


        private Line ln = new Line();
        private cResults results;
        private cResultsInterreflections resultsIreflOut;




        internal cCalculateSolarMesh(cObstacleObject _mshobj, List<cObstacleObject> _objObst, List<cSemiPermObject> _objTrees,
            double _latitude, double _longitude, List<double> _DNI, List<double> _DHI, List<double> _SNOW, double _snow_threshold, double _tilt_threshold,
            int _year, int _month, int _day, int _hour,
            int _bounces, int _rec, int _diffRes, cResultsInterreflections _ResultsIreflIn,
            bool _mt)
        {
            mshobj = _mshobj;
            objObst = _objObst;
            objTrees = _objTrees;
            latitude = _latitude;
            longitude = _longitude;
            DNI = _DNI;
            DHI = _DHI;
            SNOW = _SNOW;

            year = _year;
            month = _month;
            day = _day;
            hour = _hour;
            bounces = _bounces;
            rec = _rec;
            diffRes = _diffRes;
            ResultsIreflIn = _ResultsIreflIn;
            mt = _mt;

            snow_threshold = _snow_threshold;
            tilt_treshold = _tilt_threshold;
        }


        internal void RunHourSimulation()
        {
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            ///////////////////////   INPUTS   ///////////////////////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //analysis mesh
            Mesh msh = new Mesh();
            msh = mshobj.mesh;

            //should containt analysis surface itself
            Mesh[] obst = new Mesh[objObst.Count];


            if (SNOW.Count == 0)
            {
                for (int i = 0; i < 8760; i++)
                {
                    SNOW.Add(0.0);
                }
            }


            bounces = (bounces < 0) ? 0 : bounces;
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            ///////////////////////////   VARIABLES   ////////////////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            List<SunVector> sunvectors_list = new List<SunVector>();
            SunVector.Create8760SunVectors(ref sunvectors_list, longitude, latitude, year);
            SunVector[] sunvectors = sunvectors_list.ToArray();
            Context.cWeatherdata weather;
            weather.DHI = new List<double>(DHI);
            weather.DNI = new List<double>(DNI);
            weather.Snow = new List<double>(SNOW);

            Context.cLocation location;
            location.dLatitude = latitude;
            location.dLongitude = longitude;
            location.dTgmt = 1;

            int[] daysInMonth = new int[12];
            for (int i = 0; i < daysInMonth.Length; i++)
                daysInMonth[i] = System.DateTime.DaysInMonth(year, i + 1);
            int DOY = 0;
            for (int i = 0; i < month - 1; i++)
                DOY += daysInMonth[i];
            DOY += day;
            int HOY = (DOY - 1) * 24 + hour;



            double[][] albedos = new double[objObst.Count][];
            int[] reflType = new int[objObst.Count];
            for (int i = 0; i < objObst.Count; i++)
            {
                reflType[i] = objObst[i].reflType;
                albedos[i] = new double[1];
                albedos[i][0] = objObst[i].albedos[HOY];
                obst[i] = objObst[i].mesh;
            }
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////   DIRECT AND DIFFUSE OBSTRUCTIONS   ///////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            List<double> I = new List<double>();
            List<double> Ih = new List<double>();
            List<double> Ib = new List<double>();

            Point3d[] mshvrt = msh.Vertices.ToPoint3dArray();
            Vector3d[] mshvrtnorm = new Vector3d[mshvrt.Length];
            Sensorpoints.v3d[] v3dnormals = new Sensorpoints.v3d[mshvrtnorm.Length];
            Sensorpoints.p3d[] p3dcoords = new Sensorpoints.p3d[mshvrt.Length];
            msh.FaceNormals.ComputeFaceNormals();

            double[] arrbeta = new double[mshvrt.Length];
            double[] arrpsi = new double[mshvrt.Length];

            Vector3d betaangle = new Vector3d(0, 0, 1);
            Vector3d psiangle = new Vector3d(0, 1, 0);
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
                v3dnormals[i].X = mshvrtnorm[i].X;
                v3dnormals[i].Y = mshvrtnorm[i].Y;
                v3dnormals[i].Z = mshvrtnorm[i].Z;
                p3dcoords[i].X = mshvrt[i].X;
                p3dcoords[i].Y = mshvrt[i].Y;
                p3dcoords[i].Z = mshvrt[i].Z;
            }
            Sensorpoints p = new Sensorpoints(arrbeta, arrpsi, p3dcoords, v3dnormals, rec);


            List<double> ShdwBeam_hour = new List<double>();
            List<double[]> ShdwSky = new List<double[]>();
            List<Point3d> coords = new List<Point3d>();



            for (int i = 0; i < mshvrt.Length; i++)
            {
                Point3d orig = new Point3d(mshvrt[i].X, mshvrt[i].Y, mshvrt[i].Z);
                coords.Add(orig);
                ShdwSky.Add(new double[] { });
                ShdwBeam_hour.Add(0.0);
            }
            if (!mt)
            {
                for(int i=0;i< mshvrt.Length; i++)
                {
                    /////////////////////////////////////////////////////////////////////
                    //sky dome diffuse
                    Vector3d[] vec_sky = new Vector3d[p.sky[i].VerticesHemisphere.Count];
                    for (int u = 0; u < vec_sky.Length; u++)
                    {
                        vec_sky[u] = new Vector3d(
                            p.sky[i].VertexVectorsSphere[p.sky[i].VerticesHemisphere[u]][0],
                            p.sky[i].VertexVectorsSphere[p.sky[i].VerticesHemisphere[u]][1],
                            p.sky[i].VertexVectorsSphere[p.sky[i].VerticesHemisphere[u]][2]);
                    }
                    bool[] shdw_sky = new bool[p.sky[i].VerticesHemisphere.Count];
                    //if (!mt)
                    cShadow.CalcShadow(coords[i], mshvrtnorm[i], mshobj.tolerance, vec_sky, obst, ref shdw_sky);
                    //else
                    //    cShadow.CalcShadowMT(orig, mshvrtnorm[i], 0.01, vec_sky, obst, ref shdw_sky);
                    if (objTrees.Count > 0)
                    {
                        double[] shdw_sky_dbl = shdw_sky.Select<bool, double>(s => Convert.ToDouble(s)).ToArray<double>();
                        cShadow.CalcSemiPerm(coords[i], mshvrtnorm[i], mshobj.tolerance, vec_sky, objTrees, HOY, ref shdw_sky_dbl);
                        ShdwSky[i] = shdw_sky_dbl;
                    }
                    else
                    {
                        ShdwSky[i] = shdw_sky.Select<bool, double>(s => Convert.ToDouble(s)).ToArray<double>();
                    }
                    /////////////////////////////////////////////////////////////////////




                    /////////////////////////////////////////////////////////////////////
                    //beam for one hour only.
                    Vector3d[] vec_beam = new Vector3d[1];
                    vec_beam[0] = new Vector3d(sunvectors_list[HOY].udtCoordXYZ.x, sunvectors_list[HOY].udtCoordXYZ.y, sunvectors_list[HOY].udtCoordXYZ.z);
                    bool[] shdw_beam = new bool[1];
                    //if (!mt)
                    cShadow.CalcShadow(coords[i], mshvrtnorm[i], mshobj.tolerance, vec_beam, obst, ref shdw_beam);
                    //else
                    //    cShadow.CalcShadowMT(orig, mshvrtnorm[i], 0.01, vec_beam, obst, ref shdw_beam);
                    
                    if (objTrees.Count > 0)
                    {
                        double[] shdw_beam_dbl = new double[1] { Convert.ToDouble(shdw_beam[0]) };
                        cShadow.CalcSemiPerm(coords[i], mshvrtnorm[i], mshobj.tolerance, vec_beam, objTrees, HOY, ref shdw_beam_dbl);
                        ShdwBeam_hour[i] = shdw_beam_dbl[0];
                    }
                    else
                    {
                        ShdwBeam_hour[i] = Convert.ToDouble(shdw_beam[0]);
                    }
                    /////////////////////////////////////////////////////////////////////



                    /////////////////////////////////////////////////////////////////////
                    //one vector for visualization
                    if (i == mshvrt.Length - 1) ln = new Line(coords[i], Vector3d.Multiply(1000, vec_beam[0]));
                    //var attribs = Rhino.RhinoDoc.ActiveDoc.CreateDefaultAttributes();
                    //attribs.ObjectDecoration = Rhino.DocObjects.ObjectDecoration.BothArrowhead;
                    //Rhino.RhinoDoc.ActiveDoc.Objects.AddLine(ln, attribs);
                    /////////////////////////////////////////////////////////////////////
                }
            }
            else
            {
                Parallel.For(0, mshvrt.Length, i =>
                {
                    /////////////////////////////////////////////////////////////////////
                    //sky dome diffuse
                    Vector3d[] vec_sky = new Vector3d[p.sky[i].VerticesHemisphere.Count];
                    for (int u = 0; u < vec_sky.Length; u++)
                    {
                        vec_sky[u] = new Vector3d(
                            p.sky[i].VertexVectorsSphere[p.sky[i].VerticesHemisphere[u]][0],
                            p.sky[i].VertexVectorsSphere[p.sky[i].VerticesHemisphere[u]][1],
                            p.sky[i].VertexVectorsSphere[p.sky[i].VerticesHemisphere[u]][2]);
                    }
                    bool[] shdw_sky = new bool[p.sky[i].VerticesHemisphere.Count];
                    //if (!mt)
                    cShadow.CalcShadow(coords[i], mshvrtnorm[i], mshobj.tolerance, vec_sky, obst, ref shdw_sky);
                    //else
                    //    cShadow.CalcShadowMT(orig, mshvrtnorm[i], 0.01, vec_sky, obst, ref shdw_sky);
                    if (objTrees.Count > 0)
                    {
                        double[] shdw_sky_dbl = shdw_sky.Select<bool, double>(s => Convert.ToDouble(s)).ToArray<double>();
                        cShadow.CalcSemiPerm(coords[i], mshvrtnorm[i], mshobj.tolerance, vec_sky, objTrees, HOY, ref shdw_sky_dbl);
                        ShdwSky[i] = shdw_sky_dbl;
                    }
                    else
                    {
                        ShdwSky[i] = shdw_sky.Select<bool, double>(s => Convert.ToDouble(s)).ToArray<double>();
                    }
                    /////////////////////////////////////////////////////////////////////




                    /////////////////////////////////////////////////////////////////////
                    //beam for one hour only.
                    Vector3d[] vec_beam = new Vector3d[1];
                    vec_beam[0] = new Vector3d(sunvectors_list[HOY].udtCoordXYZ.x, sunvectors_list[HOY].udtCoordXYZ.y, sunvectors_list[HOY].udtCoordXYZ.z);
                    bool[] shdw_beam = new bool[1];
                    //if (!mt)
                    cShadow.CalcShadow(coords[i], mshvrtnorm[i], mshobj.tolerance, vec_beam, obst, ref shdw_beam);
                    //else
                    //    cShadow.CalcShadowMT(orig, mshvrtnorm[i], 0.01, vec_beam, obst, ref shdw_beam);

                    if (objTrees.Count > 0)
                    {
                        double[] shdw_beam_dbl = new double[1] { Convert.ToDouble(shdw_beam[0]) };
                        cShadow.CalcSemiPerm(coords[i], mshvrtnorm[i], mshobj.tolerance, vec_beam, objTrees, HOY, ref shdw_beam_dbl);
                        ShdwBeam_hour[i] = shdw_beam_dbl[0];
                    }
                    else
                    {
                        ShdwBeam_hour[i] = Convert.ToDouble(shdw_beam[0]);
                    }
                    /////////////////////////////////////////////////////////////////////




                    /////////////////////////////////////////////////////////////////////
                    //one vector for visualization
                    if (i == mshvrt.Length - 1) ln = new Line(coords[i], Vector3d.Multiply(1000, vec_beam[0]));
                    //var attribs = Rhino.RhinoDoc.ActiveDoc.CreateDefaultAttributes();
                    //attribs.ObjectDecoration = Rhino.DocObjects.ObjectDecoration.BothArrowhead;
                    //Rhino.RhinoDoc.ActiveDoc.Objects.AddLine(ln, attribs);
                    /////////////////////////////////////////////////////////////////////
                });
            }











            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////   INTER-REFLECTIONS   ///////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            double[] _Idiffuse = new double[mshvrt.Length];
            double[] Ispec_onehour = new double[mshvrt.Length];

            //Sensorpoints[] Idiffuse_SPs = new Sensorpoints[mshvrt.Length];
            List<List<double>> diffSP_beta_list = new List<List<double>>();
            List<List<double>> diffSP_psi_list = new List<List<double>>();
            List<List<Sensorpoints.v3d>> diffSP_normal_list = new List<List<Sensorpoints.v3d>>();
            List<List<Sensorpoints.p3d>> diffSP_coord_list = new List<List<Sensorpoints.p3d>>();
            int[][] Idiff_obstacles = new int[mshvrt.Length][];
            int[][] Idiff_domevertices = new int[mshvrt.Length][];
            SkyDome[] Idiff_domes = new SkyDome[mshvrt.Length];

            Vector3d[] vec_beam2 = new Vector3d[1];
            vec_beam2[0] = new Vector3d(sunvectors_list[HOY].udtCoordXYZ.x, sunvectors_list[HOY].udtCoordXYZ.y, sunvectors_list[HOY].udtCoordXYZ.z);
            double[][][] Ispecular2 = new double[mshvrt.Length][][];
            Vector3d[][][] Inormals2 = new Vector3d[mshvrt.Length][][];


            if (ResultsIreflIn != null) //only calculate angles, not obstruction
            {
                //Idiffuse_SPs = ResultsIreflIn.Idiffuse_SPs;
                diffSP_beta_list = ResultsIreflIn.diffSP_beta_list;
                diffSP_psi_list = ResultsIreflIn.diffSP_psi_list;
                diffSP_normal_list = ResultsIreflIn.diffSP_normal_list;
                diffSP_coord_list = ResultsIreflIn.diffSP_coord_list;
                Idiff_obstacles = ResultsIreflIn.Idiff_obstacles;
                Idiff_domevertices = ResultsIreflIn.Idiff_domevertices;
                Idiff_domes = ResultsIreflIn.Idiff_domes;

                Ispecular2 = ResultsIreflIn.Ispecular2;
                Inormals2 = ResultsIreflIn.Inormals2;

                if (!mt)
                {
                    //cShadow.CalcDiffuse(Idiffuse_SPs, Idiff_obstacles, Idiff_domevertices, Idiff_domes, DOY, hour, weather, sunvectors, obst, objObst, 0.01, snow_threshold, tilt_treshold, out _Idiffuse);
                    cShadow.CalcDiffuse2(diffSP_beta_list, diffSP_psi_list, diffSP_normal_list, diffSP_coord_list, diffRes, Idiff_obstacles, Idiff_domevertices, Idiff_domes, DOY, hour, weather, sunvectors, obst, objObst, 0.01, snow_threshold, tilt_treshold, out _Idiffuse);
                    cShadow.CalcSpecularIncident(mshvrtnorm, Ispecular2, Inormals2, weather.DNI[HOY], ref Ispec_onehour);
                }
                else
                {
                    cShadow.CalcDiffuse2MT(diffSP_beta_list, diffSP_psi_list, diffSP_normal_list, diffSP_coord_list, diffRes, Idiff_obstacles, Idiff_domevertices, Idiff_domes, DOY, hour, weather, sunvectors, obst, objObst, 0.01, snow_threshold, tilt_treshold, out _Idiffuse);
                    cShadow.CalcSpecularIncidentMT(mshvrtnorm, Ispecular2, Inormals2, weather.DNI[HOY], ref Ispec_onehour);
                }
            }
            else //calculate full interreflection (obstruction and angles)
            {
                //interreflections diffuse
                if (diffRes > -1)
                {
                    if (!mt)
                    {
                        //cShadow.CalcDiffuse_GetSPs(mshobj, mshvrt, mshvrtnorm, objObst, diffRes, out Idiffuse_SPs, out Idiff_obstacles, out Idiff_domevertices, out Idiff_domes);
                        //cShadow.CalcDiffuse(Idiffuse_SPs, Idiff_obstacles, Idiff_domevertices, Idiff_domes, DOY, hour, weather, sunvectors, obst, objObst, 0.01, snow_threshold, tilt_treshold, out _Idiffuse);
                        cShadow.CalcDiffuse_GetSPs2(mshobj, mshvrt, mshvrtnorm, objObst, diffRes, out diffSP_beta_list, out diffSP_psi_list, out diffSP_normal_list, out diffSP_coord_list,
                            out Idiff_obstacles, out Idiff_domevertices, out Idiff_domes);
                        cShadow.CalcDiffuse2(diffSP_beta_list, diffSP_psi_list, diffSP_normal_list, diffSP_coord_list, diffRes, Idiff_obstacles, Idiff_domevertices, Idiff_domes, DOY, hour, weather, sunvectors, obst, objObst, 0.01, snow_threshold, tilt_treshold, out _Idiffuse);
                    }
                    else
                    {
                        cShadow.CalcDiffuse_GetSPs2MT(mshobj, mshvrt, mshvrtnorm, objObst, diffRes, out diffSP_beta_list, out diffSP_psi_list, out diffSP_normal_list, out diffSP_coord_list, 
                            out Idiff_obstacles, out Idiff_domevertices, out Idiff_domes);
                        cShadow.CalcDiffuse2MT(diffSP_beta_list, diffSP_psi_list, diffSP_normal_list, diffSP_coord_list, diffRes, Idiff_obstacles, Idiff_domevertices, Idiff_domes, DOY, hour, weather, sunvectors, obst, objObst, 0.01, snow_threshold, tilt_treshold, out _Idiffuse);
                    }
                }

                //interreflections specular
                if (!mt)
                {
                    cShadow.CalcSpecularNormal3(mshvrt, mshvrtnorm, vec_beam2, new bool[1] { true }, HOY, objObst, bounces, ref Ispecular2, ref Inormals2);
                    cShadow.CalcSpecularIncident(mshvrtnorm, Ispecular2, Inormals2, weather.DNI[HOY], ref Ispec_onehour);
                }
                else
                {
                    cShadow.CalcSpecularNormal3MT(mshvrt, mshvrtnorm, vec_beam2, new bool[1] { true }, HOY, objObst, bounces, ref Ispecular2, ref Inormals2);
                    cShadow.CalcSpecularIncidentMT(mshvrtnorm, Ispecular2, Inormals2, weather.DNI[HOY], ref Ispec_onehour);
                }
            }
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////










            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //set shadows, snow, tree and interreflections and calculate irradiation
            if (!mt)
            {
                p.SetShadows(ShdwBeam_hour, ShdwSky, HOY);           //p.SetShadows(ShdwBeam_hour, ShdwSky, ShdwTrees_hour, HOY)
                p.SetSnowcover(snow_threshold, tilt_treshold, weather);
                p.SetInterreflection(HOY, Ispec_onehour, _Idiffuse);
                p.CalcIrradiation(DOY, hour, weather, sunvectors);
            }
            else
            {
                p.SetShadowsMT(ShdwBeam_hour, ShdwSky, HOY);           //p.SetShadows(ShdwBeam_hour, ShdwSky, ShdwTrees_hour, HOY)
                p.SetSnowcoverMT(snow_threshold, tilt_treshold, weather);
                p.SetInterreflectionMT(HOY, Ispec_onehour, _Idiffuse);
                p.CalcIrradiationMT(DOY, hour, weather, sunvectors);
            }
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////













            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////    COPY RESULTS INTO RHINO GH     ///////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            for (int i = 0; i < mshvrt.Length; i++)
            {
                I.Add(p.I[i][HOY]);
                Ib.Add(p.Ibeam[i][HOY]);
                Ih.Add(p.Idiff[i][HOY]);
            }
            Matrix I_hourly = new Matrix(mshvrt.Length, 1);
            Matrix Ib_hourly = new Matrix(mshvrt.Length, 1);
            Matrix Id_hourly = new Matrix(mshvrt.Length, 1);
            for (int i = 0; i < mshvrt.Length; i++)
            {
                I_hourly[i, 0] = I[i];
                Ib_hourly[i, 0] = Ib[i];
                Id_hourly[i, 0] = Ih[i];
            }
            results = new cResults(I, Ib, Ih, I_hourly, Ib_hourly, Id_hourly, coords);
            //resultsIreflOut = new cResultsInterreflections(Idiffuse_SPs, Idiff_obstacles, Idiff_domevertices, Idiff_domes, Ispecular2, Inormals2);
            resultsIreflOut = new cResultsInterreflections(diffSP_beta_list,diffSP_psi_list,diffSP_normal_list,diffSP_coord_list, Idiff_obstacles, Idiff_domevertices, Idiff_domes, Ispecular2, Inormals2);

            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
     
        }

        internal Line getSolarVec()
        {
            return ln;
        }

        internal cResults getResults()
        {
            return results;
        }

        internal cResultsInterreflections getResultsInterreflections()
        {
            return resultsIreflOut;
        }
   

    }

}