using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using Rhino.Geometry;
using SolarModel;
using System.Diagnostics;

/*
 * cSolarMesh.cs
 * Copyright 2017 Christoph Waibel <chwaibel@student.ethz.ch>
 * 
 * This work is licensed under the GNU GPL license version 3.
*/


namespace GHSolar
{
    /// <summary>
    /// Handles solar calculations for a mesh surface, calling different classes from both SolarModel.dll as well as GHSolar.gha (SunVector, Irradiation, cShadow, ...).
    /// </summary>
    public class CCalculateSolarMesh
    {
        private CObstacleObject mshobj;
        private List<CObstacleObject> objObst = new List<CObstacleObject>();
        private List<CPermObject> objTrees = new List<CPermObject>();

        private double latitude;
        private double longitude;
        private int year;
        private double snow_threshold;
        private double tilt_treshold;

        private List<double> DNI = new List<double>();
        private List<double> DHI = new List<double>();
        private List<double> SNOW = new List<double>();
        private List<double> groundalbedo;

        private CResultsInterreflections ResultsIreflIn;

        private bool mt;

        private const double rad = Math.PI / 180;

        private Line ln = new Line();
        private CResults results;
        private CResultsInterreflections resultsIreflOut;

        private Context.cWeatherdata weather;
        private Context.cLocation location;
        private List<SunVector> sunvectors_list;
        private SunVector[] sunvectors;


        public CCalculateSolarMesh(CObstacleObject _mshobj, List<CObstacleObject> _objObst, List<CPermObject> _objTrees,
            double _latitude, double _longitude, List<double> _DNI, List<double> _DHI,
            List<double> _SNOW, List<double> _groundalbedo, double _snow_threshold, double _tilt_threshold,
            int _year, CResultsInterreflections _ResultsIreflIn,
            bool _mt, List<double> solarAzimuth, List<double> solarAltitude,
            int timezone = 0)
        {
            mshobj = _mshobj;
            objObst = _objObst;
            objTrees = _objTrees;
            latitude = _latitude;
            longitude = _longitude;
            DNI = _DNI;
            DHI = _DHI;
            SNOW = _SNOW;

            weather.DHI = new List<double>(DHI);
            weather.DNI = new List<double>(DNI);
            weather.Snow = new List<double>(SNOW);

            location.dLatitude = latitude;
            location.dLongitude = longitude;
            location.dTgmt = 1;


            groundalbedo = _groundalbedo;

            year = _year;
            ResultsIreflIn = _ResultsIreflIn;
            mt = _mt;

            snow_threshold = _snow_threshold;
            tilt_treshold = _tilt_threshold;



            if (!Misc.IsNullOrEmpty(solarAltitude) && !Misc.IsNullOrEmpty(solarAzimuth) && solarAltitude.Count == 8760 && solarAzimuth.Count == 8760)
            {
                SunVector.Create8760SunVectors(out sunvectors_list, longitude, latitude, year, solarAzimuth.ToArray(), solarAltitude.ToArray());
            }
            else
            {
                SunVector.Create8760SunVectors(out sunvectors_list, longitude, latitude, year);
                SunVector.ShiftSunVectorsByTimezone(ref sunvectors_list, timezone);
            }
            sunvectors = sunvectors_list.ToArray();

        }


        /// <summary>
        /// Run one hour simulation. Full program. No simplifications.
        /// </summary>
        /// <param name="month">Month of simulation (1-12).</param>
        /// <param name="day">Day of simulation (1-31).</param>
        /// <param name="hour">Hour of simulation (0-23).</param>
        /// <param name="mainSkyRes">Skydome resolution.</param>
        /// <param name="specBounces">Specular bounces (0-2).</param>
        /// <param name="diffIReflSkyRes">Skydome resolution for diffuse inter-reflection.</param>
        /// <param name="diffIReflSkyRes2nd">Skydome resolution at secondary sensorpoints of diffuse inter-reflection (0 - 2).</param>
        public void RunHourSimulationMT(int month, int day, int hour,
            int mainSkyRes, int specBounces, int diffIReflSkyRes, int diffIReflSkyRes2nd)
        {
            int tasks = 1;
            if (this.mt) tasks = Environment.ProcessorCount;
            ParallelOptions paropts = new ParallelOptions { MaxDegreeOfParallelism = tasks };
            ParallelOptions paropts_1cpu = new ParallelOptions { MaxDegreeOfParallelism = 1 };

            //________________________________________________________________________________________________________________________________________
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            ///////////////////////   INPUTS   ///////////////////////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //analysis mesh
            Mesh msh = new Mesh();
            msh = mshobj.mesh;

            //should containt analysis surface itself
            Mesh[] obst = new Mesh[objObst.Count];
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //________________________________________________________________________________________________________________________________________





            //________________________________________________________________________________________________________________________________________
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            ///////////////////////////   VARIABLES   ////////////////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
            //________________________________________________________________________________________________________________________________________






            //________________________________________________________________________________________________________________________________________
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
            Sensorpoints p = new Sensorpoints(arrbeta, arrpsi, p3dcoords, v3dnormals, mainSkyRes);

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

            for (int i = 0; i < mshvrt.Length; i++)
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
                CShadow.CalcShadowMT(coords[i], mshvrtnorm[i], mshobj.tolerance, vec_sky, obst, ref shdw_sky, paropts);

                if (objTrees.Count > 0)
                {
                    double[] shdw_sky_dbl = shdw_sky.Select<bool, double>(s => Convert.ToDouble(s)).ToArray<double>();
                    CShadow.CalcPermBeamMT(coords[i], mshvrtnorm[i], mshobj.tolerance, vec_sky, objTrees, HOY, paropts, ref shdw_sky_dbl);
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
                CShadow.CalcShadowMT(coords[i], mshvrtnorm[i], mshobj.tolerance, vec_beam, obst, ref shdw_beam, paropts);

                if (objTrees.Count > 0)
                {
                    double[] shdw_beam_dbl = new double[1] { Convert.ToDouble(shdw_beam[0]) };
                    CShadow.CalcPermBeamMT(coords[i], mshvrtnorm[i], mshobj.tolerance, vec_beam, objTrees, HOY, paropts, ref shdw_beam_dbl);
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

            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //________________________________________________________________________________________________________________________________________

            //________________________________________________________________________________________________________________________________________
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

                CShadow.CalcDiffuse2MT(diffSP_beta_list, diffSP_psi_list, diffSP_normal_list, diffSP_coord_list, diffIReflSkyRes,
                    Idiff_obstacles, Idiff_domevertices, Idiff_domes, DOY, hour, weather, sunvectors,
                    obst, objObst, mshobj.tolerance, snow_threshold, tilt_treshold, paropts, out _Idiffuse);
                CShadow.CalcSpecularIncidentMT(mshvrtnorm, Ispecular2, Inormals2, weather.DNI[HOY], paropts, ref Ispec_onehour);

            }
            else //calculate full interreflection (obstruction and angles)
            {
                //interreflections diffuse
                if (diffIReflSkyRes > -1)
                {
                    CShadow.CalcIReflDiff_GetSPs2MT(mshobj, mshvrt, mshvrtnorm, objObst, objTrees, diffIReflSkyRes, paropts, 
                        out diffSP_beta_list, out diffSP_psi_list, out diffSP_normal_list, out diffSP_coord_list,
                        out Idiff_obstacles, out Idiff_domevertices, out Idiff_domes);
                    CShadow.CalcDiffuse2MT(diffSP_beta_list, diffSP_psi_list, diffSP_normal_list, diffSP_coord_list, diffIReflSkyRes,
                        Idiff_obstacles, Idiff_domevertices, Idiff_domes, DOY, hour, weather, sunvectors,
                        obst, objObst, mshobj.tolerance, snow_threshold, tilt_treshold, paropts, out _Idiffuse);
                }
                else
                {
                    p.SetSimpleGroundReflectionMT(arrbeta, groundalbedo.ToArray(), weather, sunvectors, paropts);
                }

                //interreflections specular
                CShadow.CalcSpecularNormal3MT(mshvrt, mshvrtnorm, vec_beam2, new bool[1] { true }, HOY, objObst, objTrees, specBounces, ref Ispecular2, ref Inormals2);
                CShadow.CalcSpecularIncidentMT(mshvrtnorm, Ispecular2, Inormals2, weather.DNI[HOY], paropts, ref Ispec_onehour);
            }
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //________________________________________________________________________________________________________________________________________

            //________________________________________________________________________________________________________________________________________
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////   IRRADIATION CALCULATION   /////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //set shadows, snow, tree and interreflections and calculate irradiation
            p.SetShadowsMT(ShdwBeam_hour, ShdwSky, HOY);           //p.SetShadows(ShdwBeam_hour, ShdwSky, ShdwTrees_hour, HOY)
            p.SetSnowcoverMT(snow_threshold, tilt_treshold, weather, paropts);
            p.SetInterreflectionMT(HOY, Ispec_onehour, _Idiffuse, paropts);
            p.CalcIrradiationMT(DOY, hour, weather, sunvectors, paropts);
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //________________________________________________________________________________________________________________________________________

            //________________________________________________________________________________________________________________________________________
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
            results = new CResults(I, Ib, Ih, I_hourly, Ib_hourly, Id_hourly, coords);
            //resultsIreflOut = new cResultsInterreflections(Idiffuse_SPs, Idiff_obstacles, Idiff_domevertices, Idiff_domes, Ispecular2, Inormals2);
            resultsIreflOut = new CResultsInterreflections(diffSP_beta_list, diffSP_psi_list, diffSP_normal_list, diffSP_coord_list, Idiff_obstacles, Idiff_domevertices, Idiff_domes, Ispecular2, Inormals2);
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
            //________________________________________________________________________________________________________________________________________
        }


        /// <summary>
        /// Run annual irradiation simulation, with some simplifications like interpolation. Multi-threading version.
        /// </summary>
        /// <param name="tolerance"></param>
        public void RunAnnualSimulation_MT(double tolerance,
            int mainSkyRes, int mainInterpMode, int specBounces, int specInterpMode,
            int diffIReflSkyRes, int diffIReflSkyRes2nd, int diffIReflMode)
        {
            int tasks = 1;
            if (this.mt) tasks = Environment.ProcessorCount;
            ParallelOptions paropts = new ParallelOptions { MaxDegreeOfParallelism = tasks };
            ParallelOptions paropts_1cpu = new ParallelOptions { MaxDegreeOfParallelism = 1 };

            Rhino.RhinoApp.WriteLine("SOLAR MODEL. https://github.com/christophwaibel/GH_Solar_V2");
            Stopwatch stopwatch = Stopwatch.StartNew(); //creates and start the instance of Stopwatch
            Rhino.RhinoApp.WriteLine("(1/4): Preparing Data...");
            //________________________________________________________________________________________________________________________________________
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            ///////////////////////   INPUTS   ///////////////////////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //analysis mesh
            Mesh msh = new Mesh();
            msh = mshobj.mesh;

            //should containt analysis surface itself
            Mesh[] obst = new Mesh[objObst.Count];


            //________________________________________________________________________________________________________________________________________
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            ///////////////////////////   VARIABLES   ////////////////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            double[][] albedos = new double[objObst.Count][];
            int[] reflType = new int[objObst.Count];
            Parallel.For(0, objObst.Count, paropts, i =>
            {
                reflType[i] = objObst[i].reflType;
                albedos[i] = new double[8760];
                obst[i] = objObst[i].mesh;
                for (int t = 0; t < 8760; t++)
                {
                    albedos[i][t] = objObst[i].albedos[t];
                }
            });


            //________________________________________________________________________________________________________________________________________
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////   DIRECT AND DIFFUSE OBSTRUCTIONS   ///////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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

            //3days MAIN
            List<bool[]> ShdwBeam_equinox = new List<bool[]>();
            List<bool[]> ShdwBeam_summer = new List<bool[]>();
            List<bool[]> ShdwBeam_winter = new List<bool[]>();
            //12 days MAIN
            List<bool[][]> ShdwBeam_12d = new List<bool[][]>();
            int[] startDays = new int[12];
            int[] endDays = new int[12];

            bool[][] BeamPermIs_equ = new bool[mshvrt.Length][];
            int[][][] BeamPermRefs_equ = new int[mshvrt.Length][][];
            double[][][] BeamPermLength_equ = new double[mshvrt.Length][][];
            bool[][] BeamPermIs_sum = new bool[mshvrt.Length][];
            int[][][] BeamPermRefs_sum = new int[mshvrt.Length][][];
            double[][][] BeamPermLength_sum = new double[mshvrt.Length][][];
            bool[][] BeamPermIs_win = new bool[mshvrt.Length][];
            int[][][] BeamPermRefs_win = new int[mshvrt.Length][][];
            double[][][] BeamPermLength_win = new double[mshvrt.Length][][];

            bool[][][] BeamPermIs_12d = new bool[mshvrt.Length][][];
            int[][][][] BeamPermRefs_12d = new int[mshvrt.Length][][][];
            double[][][][] BeamPermLength_12d = new double[mshvrt.Length][][][];

            List<bool[]> ShdwSky = new List<bool[]>();
            bool[][] SkyPermIs = new bool[mshvrt.Length][];
            int[][][] SkyPermRefs = new int[mshvrt.Length][][];
            double[][][] SkyPermLength = new double[mshvrt.Length][][];


            //________________________________________________________________________________________________________________________________________
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////   INTER-REFLECTIONS   ///////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            double[][] _Idiffuse = new double[mshvrt.Length][];
            double[][] _Ispec_annual = new double[mshvrtnorm.Length][];

            List<List<double>> diffSP_beta_list = new List<List<double>>();
            List<List<double>> diffSP_psi_list = new List<List<double>>();
            List<List<Sensorpoints.v3d>> diffSP_normal_list = new List<List<Sensorpoints.v3d>>();
            List<List<Sensorpoints.p3d>> diffSP_coord_list = new List<List<Sensorpoints.p3d>>();
            int[][] Idiff_obstacles = new int[mshvrt.Length][];
            int[][] Idiff_domevertices = new int[mshvrt.Length][];
            SkyDome[] Idiff_domes = new SkyDome[mshvrt.Length];


            //________________________________________________________________________________________________________________________________________
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////       PRE-LOOP        ///////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            List<Point3d> coords = new List<Point3d>();
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

                Point3d orig = new Point3d(mshvrt[i].X, mshvrt[i].Y, mshvrt[i].Z);
                coords.Add(orig);

                ShdwSky.Add(new bool[] { });
                ShdwBeam_equinox.Add(new bool[] { });
                ShdwBeam_summer.Add(new bool[] { });
                ShdwBeam_winter.Add(new bool[] { });

                ShdwBeam_12d.Add(new bool[][] { });

                BeamPermIs_12d[i] = new bool[12][];
                BeamPermRefs_12d[i] = new int[12][][];
                BeamPermLength_12d[i] = new double[12][][];
            }
            Sensorpoints p = new Sensorpoints(arrbeta, arrpsi, p3dcoords, v3dnormals, mainSkyRes);


            //3 days Main
            // equinox:             march 20
            // summer solstice:     june 21
            // winter solstice:     december 21
            int[] equsol = SunVector.GetEquinoxSolstice(year);
            int HOYequ = (equsol[0] - 1) * 24;
            int HOYsum = (equsol[1] - 1) * 24;
            int HOYwin = (equsol[3] - 1) * 24;

            Vector3d[] vec_beam_equ = new Vector3d[24];
            Vector3d[] vec_beam_sum = new Vector3d[24];
            Vector3d[] vec_beam_win = new Vector3d[24];
            bool[] sunshine_equ = new bool[24];
            bool[] sunshine_sum = new bool[24];
            bool[] sunshine_win = new bool[24];
            Parallel.For(0, 24, paropts, t =>
            {
                if (sunvectors_list[HOYequ + t].Sunshine)
                    sunshine_equ[t] = true;
                if (sunvectors_list[HOYsum + t].Sunshine)
                    sunshine_sum[t] = true;
                if (sunvectors_list[HOYwin + t].Sunshine)
                    sunshine_win[t] = true;
                vec_beam_equ[t] = new Vector3d(sunvectors_list[HOYequ + t].udtCoordXYZ.x, sunvectors_list[HOYequ + t].udtCoordXYZ.y, sunvectors_list[HOYequ + t].udtCoordXYZ.z);
                vec_beam_sum[t] = new Vector3d(sunvectors_list[HOYsum + t].udtCoordXYZ.x, sunvectors_list[HOYsum + t].udtCoordXYZ.y, sunvectors_list[HOYsum + t].udtCoordXYZ.z);
                vec_beam_win[t] = new Vector3d(sunvectors_list[HOYwin + t].udtCoordXYZ.x, sunvectors_list[HOYwin + t].udtCoordXYZ.y, sunvectors_list[HOYwin + t].udtCoordXYZ.z);
            });


            //12 days main
            int dmcount = 1;    //days in month counter
            Vector3d[][] vec_beam_12d = new Vector3d[12][];
            bool[][] sunshine_12d = new bool[12][];
            for (int d = 0; d < 12; d++)
            {
                vec_beam_12d[d] = new Vector3d[24];
                sunshine_12d[d] = new bool[24];
                int dm = System.DateTime.DaysInMonth(year, d + 1);
                startDays[d] = dmcount;
                endDays[d] = dm + dmcount;
                dmcount += dm;
                int HOY = (startDays[d] - 1) * 24;
                Parallel.For(0, 24, paropts, t =>
                {
                    if (sunvectors_list[HOY + t].Sunshine)
                        sunshine_12d[d][t] = true;
                    vec_beam_12d[d][t] = new Vector3d(sunvectors_list[HOY + t].udtCoordXYZ.x, sunvectors_list[HOY + t].udtCoordXYZ.y, sunvectors_list[HOY + t].udtCoordXYZ.z);
                });
            }


            stopwatch.Stop();
            Rhino.RhinoApp.WriteLine("(1/4): " + Convert.ToString(stopwatch.Elapsed));
            stopwatch.Reset();
            stopwatch.Start();
            Rhino.RhinoApp.WriteLine("(2/4): Simulating...");
            //________________________________________________________________________________________________________________________________________
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////       MAIN-LOOP       ///////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            for (int i = 0; i < mshvrt.Length; i++)
            {
                double status = (100 / Convert.ToDouble(mshvrt.Length)) * Convert.ToDouble(i + 1);
                Rhino.RhinoApp.WriteLine("(2/4) Simulating... " + Convert.ToString(Math.Round(status, 2) + "%"));


                /////////////////////////////////////////////////////////////////////
                //sky dome diffuse
                Vector3d[] vec_sky = new Vector3d[p.sky[i].VerticesHemisphere.Count];
                bool[] sunshinesky = new bool[vec_sky.Length];
                Parallel.For(0, vec_sky.Length, paropts, u =>
                {
                    vec_sky[u] = new Vector3d(
                        p.sky[i].VertexVectorsSphere[p.sky[i].VerticesHemisphere[u]][0],
                        p.sky[i].VertexVectorsSphere[p.sky[i].VerticesHemisphere[u]][1],
                        p.sky[i].VertexVectorsSphere[p.sky[i].VerticesHemisphere[u]][2]);
                    sunshinesky[u] = true;
                });
                bool[] shdw_sky = new bool[p.sky[i].VerticesHemisphere.Count];
                CShadow.CalcShadowMT(coords[i], mshvrtnorm[i], tolerance, vec_sky, obst, ref shdw_sky, paropts);
                ShdwSky[i] = shdw_sky;
                if (this.objTrees.Count > 0)
                {
                    CShadow.CalcPermMT(coords[i], mshvrtnorm[i], mshobj.tolerance, vec_sky, sunshinesky, this.objTrees, shdw_sky, paropts,
                        out SkyPermIs[i], out SkyPermRefs[i], out SkyPermLength[i]);
                }
                /////////////////////////////////////////////////////////////////////


                /////////////////////////////////////////////////////////////////////
                //beam MAIN
                if (mainInterpMode == 0)
                {
                    //beam 3 days
                    bool[] shdw_beam_equ = new bool[24];
                    bool[] shdw_beam_sum = new bool[24];
                    bool[] shdw_beam_win = new bool[24];
                    CShadow.CalcShadowMT(coords[i], mshvrtnorm[i], mshobj.tolerance, vec_beam_equ, sunshine_equ, obst, ref shdw_beam_equ, paropts);
                    CShadow.CalcShadowMT(coords[i], mshvrtnorm[i], mshobj.tolerance, vec_beam_sum, sunshine_sum, obst, ref shdw_beam_sum, paropts);
                    CShadow.CalcShadowMT(coords[i], mshvrtnorm[i], mshobj.tolerance, vec_beam_win, sunshine_win, obst, ref shdw_beam_win, paropts);
                    ShdwBeam_equinox[i] = shdw_beam_equ;
                    ShdwBeam_summer[i] = shdw_beam_sum;
                    ShdwBeam_winter[i] = shdw_beam_win;

                    if (this.objTrees.Count > 0)
                    {
                        CShadow.CalcPermMT(coords[i], mshvrtnorm[i], mshobj.tolerance, vec_beam_equ, sunshine_equ, this.objTrees, shdw_beam_equ, paropts, 
                            out BeamPermIs_equ[i], out BeamPermRefs_equ[i], out BeamPermLength_equ[i]);

                        CShadow.CalcPermMT(coords[i], mshvrtnorm[i], mshobj.tolerance, vec_beam_sum, sunshine_sum, this.objTrees, shdw_beam_sum, paropts, 
                            out BeamPermIs_sum[i], out BeamPermRefs_sum[i], out BeamPermLength_sum[i]);

                        CShadow.CalcPermMT(coords[i], mshvrtnorm[i], mshobj.tolerance, vec_beam_win, sunshine_win, this.objTrees, shdw_beam_win, paropts,
                            out BeamPermIs_win[i], out BeamPermRefs_win[i], out BeamPermLength_win[i]);
                    }
                }
                else    //12 days beam Main
                {
                    bool[][] shdw_beam_12d = new bool[12][];
                    for (int d = 0; d < 12; d++)
                    {
                        shdw_beam_12d[d] = new bool[24];
                        CShadow.CalcShadowMT(coords[i], mshvrtnorm[i], mshobj.tolerance, vec_beam_12d[d], sunshine_12d[d], obst, ref shdw_beam_12d[d], paropts);
                    }
                    ShdwBeam_12d[i] = shdw_beam_12d;

                    if (this.objTrees.Count > 0)
                    {
                        for (int d = 0; d < 12; d++)
                        {
                            CShadow.CalcPermMT(coords[i], mshvrtnorm[i], mshobj.tolerance, vec_beam_12d[d], sunshine_12d[d], this.objTrees, shdw_beam_12d[d], paropts,
                                out BeamPermIs_12d[i][d], out BeamPermRefs_12d[i][d], out BeamPermLength_12d[i][d]);
                        }
                    }
                }
                /////////////////////////////////////////////////////////////////////


                /////////////////////////////////////////////////////////////////////
                //Interreflections
                double[][] _Idiffuse_1;
                double[][] _Ispec_annual_1;
                Point3d[] mshvrt_1 = new Point3d[1];
                mshvrt_1[0] = mshvrt[i];
                Vector3d[] mshvrtnorm_1 = new Vector3d[1];
                mshvrtnorm_1[0] = mshvrtnorm[i];


                //calculate full interreflection (obstruction and angles)
                if (diffIReflMode >= 1)   //also used to ignore specular reflections. full diffuse and 3 days specular
                {
                    List<List<double>> diffSP_beta_list_1 = new List<List<double>>();
                    List<List<double>> diffSP_psi_list_1 = new List<List<double>>();
                    List<List<Sensorpoints.v3d>> diffSP_normal_list_1 = new List<List<Sensorpoints.v3d>>();
                    List<List<Sensorpoints.p3d>> diffSP_coord_list_1 = new List<List<Sensorpoints.p3d>>();
                    int[][] Idiff_obstacles_1 = new int[1][];
                    int[][] Idiff_domevertices_1 = new int[1][];

                    SkyDome[] Idiff_domes_1 = new SkyDome[1];
                    if (diffIReflMode == 1)  //simple diffuse interreflections. 3 days specular.
                    {
                        //interreflection sky diffuse
                        CShadow.CalcIReflDiff_GetSPs2MT(mshobj, mshvrt_1, mshvrtnorm_1, objObst, objTrees, diffIReflSkyRes, paropts,
                            out diffSP_beta_list_1, out diffSP_psi_list_1, out diffSP_normal_list_1, out diffSP_coord_list_1,
                            out Idiff_obstacles_1, out Idiff_domevertices_1, out Idiff_domes_1);
                        CShadow.CalcDiffuse_AnnualSimpleMT(diffSP_beta_list_1, diffSP_psi_list_1, diffSP_normal_list_1, diffSP_coord_list_1,
                            Idiff_obstacles_1, Idiff_domevertices_1, Idiff_domes_1, diffIReflSkyRes2nd,
                            weather, sunvectors, objObst, snow_threshold, tilt_treshold, groundalbedo.ToArray(), paropts, 
                            out _Idiffuse_1);
                    }
                    else // 2... full diffuse
                    {
                        //interreflections diffuse
                        CShadow.CalcIReflDiff_GetSPs2MT(mshobj, mshvrt_1, mshvrtnorm_1, objObst, objTrees, diffIReflSkyRes, paropts, 
                            out diffSP_beta_list_1, out diffSP_psi_list_1, out diffSP_normal_list_1, out diffSP_coord_list_1,
                            out Idiff_obstacles_1, out Idiff_domevertices_1, out Idiff_domes_1);
                        CShadow.CalcDiffuse_AnnualMT(diffSP_beta_list_1, diffSP_psi_list_1, diffSP_normal_list_1, diffSP_coord_list_1,
                            Idiff_obstacles_1, Idiff_domevertices_1, Idiff_domes_1, diffIReflSkyRes2nd,
                            year, weather, sunvectors, objObst, objTrees, tolerance, snow_threshold, tilt_treshold, groundalbedo.ToArray(), paropts,
                            out _Idiffuse_1);
                    }
                }
                else
                {
                    _Idiffuse_1 = new double[1][];
                    _Idiffuse_1[0] = new double[8760];
                }

                if (specBounces > 0)
                {
                    if (specInterpMode == 0)
                    {
                        //interreflections specular 3 days
                        int[][][] IObstRef1st_equ, IObstRef1st_win, IObstRef1st_sum;
                        int[][][] IObstRef2nd_equ, IObstRef2nd_win, IObstRef2nd_sum;
                        Vector3d[][][] Inormals_equ, Inormals_win, Inormals_sum;
                        CShadow.CalcSpecularNormal5MT(mshvrt_1, mshvrtnorm_1, vec_beam_equ, sunshine_equ, this.objObst, this.objTrees, specBounces,
                            paropts,
                            out IObstRef1st_equ, out IObstRef2nd_equ, out Inormals_equ);       //equinox with 24 vectors
                        CShadow.CalcSpecularNormal5MT(mshvrt_1, mshvrtnorm_1, vec_beam_win, sunshine_win, this.objObst, this.objTrees, specBounces,
                            paropts,
                            out IObstRef1st_win, out IObstRef2nd_win, out Inormals_win);       //winter with 24 vectors
                        CShadow.CalcSpecularNormal5MT(mshvrt_1, mshvrtnorm_1, vec_beam_sum, sunshine_sum, this.objObst, this.objTrees, specBounces,
                            paropts,
                            out IObstRef1st_sum, out IObstRef2nd_sum, out Inormals_sum);       //summer

                        CShadow.CalcSpecularIncident_AnnualMT(
                            IObstRef1st_equ, IObstRef2nd_equ, Inormals_equ,
                            IObstRef1st_win, IObstRef2nd_win, Inormals_win,
                            IObstRef1st_sum, IObstRef2nd_sum, Inormals_sum,
                            this.objObst, weather.DNI.ToArray(), mshvrtnorm_1, paropts, 
                            out _Ispec_annual_1);
                    }
                    else    //12 days specular interrefl interpolation
                    {
                        int[][][][] IObstRef1st_12d = new int[12][][][];
                        int[][][][] IObstRef2nd_12d = new int[12][][][];
                        Vector3d[][][][] Inormals_12d = new Vector3d[12][][][];

                        Parallel.For(0, 12, paropts, d =>
                        {
                            CShadow.CalcSpecularNormal5MT(mshvrt_1, mshvrtnorm_1, vec_beam_12d[d], sunshine_12d[d], this.objObst, this.objTrees, specBounces,
                                paropts_1cpu,
                                out IObstRef1st_12d[d], out IObstRef2nd_12d[d], out Inormals_12d[d]);       //summer
                        });
                        CShadow.CalcSpecularIncident_AnnualMT(startDays, endDays, IObstRef1st_12d, IObstRef2nd_12d, Inormals_12d,
                            this.objObst, weather.DNI.ToArray(), mshvrtnorm_1, paropts,
                            out _Ispec_annual_1);
                    }
                }
                else
                {
                    _Ispec_annual_1 = new double[1][];
                    _Ispec_annual_1[0] = new double[8760];
                }


                _Idiffuse[i] = _Idiffuse_1[0];
                _Ispec_annual[i] = _Ispec_annual_1[0];
            }
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //________________________________________________________________________________________________________________________________________


            //put all _Ispec_annual_1 and _Idiffuse_1 into _Ispec_annual and _Idiffuse

            if (diffIReflMode == 0)   //no obstruction calculations for interreflections. 
            {
                if (specBounces > 0)    //could still have calculated spec bounces
                {
                    p.SetInterrefl_AnnualMT(_Ispec_annual, _Idiffuse, paropts);
                }
                p.SetSimpleGroundReflectionMT(arrbeta, groundalbedo.ToArray(), weather, sunvectors, paropts);
            }
            else
            {
                p.SetInterrefl_AnnualMT(_Ispec_annual, _Idiffuse, paropts);
            }



            stopwatch.Stop();
            Rhino.RhinoApp.WriteLine("(2/4): " + Convert.ToString(stopwatch.Elapsed));
            stopwatch.Reset();
            stopwatch.Start();
            Rhino.RhinoApp.WriteLine("(3/4): Calculating Perez Solar Model...");
            //________________________________________________________________________________________________________________________________________
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////   IRRADIATION CALCULATION   /////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            List<double[]> extinctCoeff = new List<double[]>();
            foreach (CPermObject perm in this.objTrees)
            {
                extinctCoeff.Add(perm.permeability.ToArray());
            }

            Stopwatch stopwatch2 = Stopwatch.StartNew(); //creates and start the instance of Stopwatch
            stopwatch2.Start();
            //set shadows, snow, tree and interreflections and calculate irradiation
            if (this.objTrees.Count > 0)
            {
                if (mainInterpMode == 0)   //3 days trees
                {
                    p.SetShadows_Annual_PermeablesMT(ShdwBeam_equinox, ShdwBeam_summer, ShdwBeam_winter,
                        BeamPermIs_equ, BeamPermRefs_equ, BeamPermLength_equ,
                        BeamPermIs_sum, BeamPermRefs_sum, BeamPermLength_sum,
                        BeamPermIs_win, BeamPermRefs_win, BeamPermLength_win,
                        ShdwSky, SkyPermIs, SkyPermRefs, SkyPermLength,
                        extinctCoeff, paropts);
                }
                else                       //12 days trees
                {
                    p.SetShadows_Annual_PermeablesMT(startDays, endDays,
                        ShdwBeam_12d, BeamPermIs_12d,
                        BeamPermRefs_12d, BeamPermLength_12d,
                        ShdwSky, SkyPermIs, SkyPermRefs, SkyPermLength,
                        extinctCoeff, paropts);
                }
            }
            else
            {
                if (mainInterpMode == 0)    //3 days no trees
                {
                    p.SetShadowsInterpolatedMT(ShdwBeam_equinox, ShdwBeam_summer, ShdwBeam_winter, ShdwSky, paropts);
                }
                else                        //12 days no trees
                {
                    p.SetShadowsInterpolatedMT(startDays, endDays, ShdwBeam_12d, ShdwSky, paropts);
                }
            }
            stopwatch2.Stop();
            Rhino.RhinoApp.WriteLine("...SetShadows: " + Convert.ToString(stopwatch2.Elapsed));
            stopwatch2.Reset();
            stopwatch2.Start();
            p.SetSnowcoverMT(snow_threshold, tilt_treshold, weather, paropts);
            stopwatch2.Stop();
            Rhino.RhinoApp.WriteLine("...SetSnowcover: " + Convert.ToString(stopwatch2.Elapsed));
            stopwatch2.Reset();
            stopwatch2.Start();
            p.CalcIrradiationMT(weather, sunvectors, paropts);
            stopwatch2.Stop();
            Rhino.RhinoApp.WriteLine("...CalcIrradiation: " + Convert.ToString(stopwatch2.Elapsed));
            stopwatch2.Reset();
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //________________________________________________________________________________________________________________________________________


            stopwatch.Stop();
            Rhino.RhinoApp.WriteLine("(3/4): " + Convert.ToString(stopwatch.Elapsed));
            stopwatch.Reset();
            stopwatch.Start();
            Rhino.RhinoApp.WriteLine("(4/4): Writing data to Rhino...");
            //________________________________________________________________________________________________________________________________________
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////    COPY RESULTS INTO RHINO GH     ///////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            List<double> I = new List<double>();        // kW/sqm a
            List<double> Ih = new List<double>();
            List<double> Ib = new List<double>();

            Matrix I_hourly = new Matrix(mshvrt.Length, 8760);
            Matrix Ib_hourly = new Matrix(mshvrt.Length, 8760);
            Matrix Id_hourly = new Matrix(mshvrt.Length, 8760);

            for (int i = 0; i < mshvrt.Length; i++)
            {
                I.Add(0.0);
                Ib.Add(0.0);
                Ih.Add(0.0);
                for (int HOY = 0; HOY < 8760; HOY++)
                {
                    double temp = p.I[i][HOY];
                    double pINow = temp < 0.0 ? 0.0 : temp;

                    I[i] += pINow;
                    Ib[i] += p.Ibeam[i][HOY];
                    Ih[i] += p.Idiff[i][HOY];

                    I_hourly[i, HOY] = pINow;
                    Ib_hourly[i, HOY] = p.Ibeam[i][HOY];
                    Id_hourly[i, HOY] = p.Idiff[i][HOY];
                }
            }
            results = new CResults(I, Ib, Ih, I_hourly, Ib_hourly, Id_hourly, coords);
            //resultsIreflOut = new cResultsInterreflections(diffSP_beta_list, diffSP_psi_list, diffSP_normal_list, diffSP_coord_list, Idiff_obstacles, Idiff_domevertices, Idiff_domes, Ispecular2, Inormals2);
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
            //________________________________________________________________________________________________________________________________________
            stopwatch.Stop();
            Rhino.RhinoApp.WriteLine(Convert.ToString(stopwatch.Elapsed));
            stopwatch.Reset();
        }


        public Line getSolarVec()
        {
            return ln;
        }


        public CResults getResults()
        {
            return results;
        }


        public CResultsInterreflections getResultsInterreflections()
        {
            return resultsIreflOut;
        }
    }
}