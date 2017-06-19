using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

/*
 * Sensorpoints.cs
 * Copyright 2017 Christoph Waibel <chwaibel@student.ethz.ch>
 * 
 * This work is licensed under the GNU GPL license version 3.
*/

namespace SolarModel
{
    /// <summary>
    /// Sensorpoints class to calculate solar irradiation.
    /// </summary>
    public class Sensorpoints
    {
        /// <summary>
        /// Sky hemisphere. See class "SykDome".
        /// </summary>
        public SkyDome[] sky { get; private set; }
        /// <summary>
        /// Sensor point tilt angle in degree.
        /// </summary>
        public double[] beta { get; private set; }
        /// <summary>
        /// Sensor point azimuth angle in degree.
        /// </summary>
        public double[] psi { get; private set; }
        /// <summary>
        /// Sensor point normal vectors.
        /// </summary>
        public v3d[] normal { get; private set; }
        /// <summary>
        /// Sensor point coordinates.
        /// </summary>
        public p3d[] coord { get; private set; }
        ///// <summary>
        ///// Diffuse irradiation. 
        ///// <para>Jagged array: [x sensor-points] [8760 hours of the year].</para>
        ///// </summary>
        public double[][] Idiff { get; private set; }
        /// <summary>
        /// Beam (direct) irradiation.
        /// <para>Jagged array: [x sensor-points] [8760 hours of the year].</para>
        /// </summary>
        public double[][] Ibeam { get; private set; }
        /// <summary>
        /// Irradiation from specular inter-reflections. 
        /// <para>Jagged array: [x sensor-points] [8760 hours of the year].</para>
        /// </summary>
        public double[][] Irefl_spec { get; private set; }
        /// <summary>
        /// Irradiation from diffuse inter-reflections.
        /// <para>Jagged array: [x sensor-points] [8760 hours of the year].</para>
        /// </summary>
        public double[][] Irefl_diff { get; private set; }
        /// <summary>
        /// Total irradiation (beam + diffuse).
        /// <para>For each sensor point: Array with 8760 doubles (each hour of the year).</para>
        /// </summary>
        public double[][] I { get; private set; }
        /// <summary>
        /// Indicates snow cover for each sensor point and each hour of the year.
        /// <para>this.snowcovered[i][HOY], i ∈ n sensorpoints, HOY ∈ [0, 8759]</para>
        /// </summary>
        public bool[][] snowcovered { get; private set; }
        /// <summary>
        /// Number of sensor points in this object.
        /// </summary>
        public int SPCount { get; private set; }


        /// <summary>
        /// Sensor point object, to calculate solar irradiation.
        /// </summary>
        /// <param name="beta">Tilt angles of sensor points.</param>
        /// <param name="psi">Azimuth angles of sensor points.</param>
        /// <param name="normal">Normal vectors of sensor points.</param>
        /// <param name="reclvlsky">Recursion level for sky hemisphere. 1 or 2 recommended. See class "SkyDome".</param>
        public Sensorpoints(double[] beta, double[] psi, p3d[] coord, v3d[] normal, int reclvlsky) //Context.cWeatherdata weather, Context.cLocation location, List<SunVector> sunvectors, 
        {
            this.SPCount = beta.Length;

            //this.location = location;
            //this.weather = weather;
            //this.sunvectors = new List<SunVector>(sunvectors);

            this.beta = new double[beta.Length];
            this.psi = new double[psi.Length];
            this.coord = new p3d[coord.Length];
            this.normal = new v3d[normal.Length];
            this.sky = new SkyDome[this.SPCount];
            Array.Copy(beta, this.beta, beta.Length);
            Array.Copy(psi, this.psi, psi.Length);
            Array.Copy(coord, this.coord, coord.Length);
            Array.Copy(normal, this.normal, normal.Length);
            for (int i = 0; i < this.SPCount; i++)
            {
                if (i == 0) this.sky[i] = new SkyDome(reclvlsky);
                else this.sky[i] = new SkyDome(this.sky[0]);
            }

            this.Idiff = new double[this.SPCount][];
            this.Ibeam = new double[this.SPCount][];
            this.Irefl_spec = new double[this.SPCount][];
            this.Irefl_diff = new double[this.SPCount][];
            this.I = new double[this.SPCount][];
            this.snowcovered = new bool[this.SPCount][];

            for (int i = 0; i < this.SPCount; i++)
            {
                this.Idiff[i] = new double[8760];
                this.Ibeam[i] = new double[8760];
                this.Irefl_spec[i] = new double[8760];
                this.Irefl_diff[i] = new double[8760];
                this.I[i] = new double[8760];
                this.snowcovered[i] = new bool[8760];
            }
        }

        /// <summary>
        /// 3d vector.
        /// </summary>
        public struct v3d
        {
            public double X;
            public double Y;
            public double Z;
        }
        /// <summary>
        /// 3d point.
        /// </summary>
        public struct p3d
        {
            public double X;
            public double Y;
            public double Z;
        }




        /// <summary>
        /// Calculates diffuse irradiation on the sensor point for one hour of the year.
        /// <para>Access: total diffuse irradiation: this.Idiff[HOY][0]; horizon: this.Idiff[HOY][1]; sky: this.Idiff[HOY][2]; circumsolar: this.Idiff[HOY][3]. HOY ∈ [0, 8759].</para>
        /// </summary>
        /// <param name="DOY">Day of the year, ∈ [1, 365].</param>
        /// <param name="HOY">Hour of the year, ∈ [0, 8759].</param>
        /// <param name="weather">Weather data.</param>
        /// <param name="sunvectors">Sunvectors, [0, 8759].</param>
        private void CalcIdiff(int DOY, int HOY, Context.cWeatherdata weather, SunVector [] sunvectors)
        {
            for (int i = 0; i < this.SPCount; i++)
            {
                if (this.snowcovered[i][HOY])
                    this.Idiff[i][HOY] = 0.0;
                else
                {
                    this.Idiff[i][HOY] = Irradiation.Diffuse(
                        weather.DHI[HOY], weather.DNI[HOY], sunvectors[HOY].udtCoordinates.dZenithAngle,
                        sunvectors[HOY].udtCoordinates.dAzimuth, this.beta[i], this.psi[i], DOY,
                        this.sky[i].ShdwHorizon[HOY], this.sky[i].ShdwDome[HOY], this.sky[i].ShdwBeam[HOY]);
                }
            }
        }

        /// <summary>
        /// Calculates diffuse irradiation on the sensor point for one hour of the year. Multi-Threading version.
        /// <para>Access: total: this.Idiff[HOY][0]; horizon: this.Idiff[HOY][1]; sky: this.Idiff[HOY][2]; circumsolar: this.Idiff[HOY][3]. HOY ∈ [0, 8759].</para>
        /// </summary>
        /// <param name="DOY">Day of the year, ∈ [1, 365].</param>
        /// <param name="HOY">Hour of the year, ∈ [0, 8759].</param>
        /// <param name="weather">Weather data.</param>
        /// <param name="sunvectors">Sunvectors, [0, 8759].</param>
        private void CalcIdiff_MT(int DOY, int HOY, Context.cWeatherdata weather, SunVector[] sunvectors)
        {
            Parallel.For(0, this.SPCount, i =>
            {
                if (snowcovered[i][HOY])
                    this.Idiff[i][HOY] = 0.0;
                else
                {
                    this.Idiff[i][HOY] = Irradiation.Diffuse(
                        weather.DHI[HOY], weather.DNI[HOY], sunvectors[HOY].udtCoordinates.dZenithAngle,
                        sunvectors[HOY].udtCoordinates.dAzimuth, this.beta[i], this.psi[i], DOY,
                        sky[i].ShdwHorizon[HOY], sky[i].ShdwDome[HOY], sky[i].ShdwBeam[HOY]);
                }
            });
        }

        /// <summary>
        /// Calculates hourly diffuse radiation on the sensor point for the entire year.
        /// <para>Access: total: this.Idiff[HOY][0]; horizon: this.Idiff[HOY][1]; sky: this.Idiff[HOY][2]; circumsolar: this.Idiff[HOY][3]. HOY ∈ [0, 8759].</para>
        /// </summary>        
        /// <param name="weather">Weather data.</param>
        /// <param name="sunvectors">Sunvectors, [0, 8759].</param>
        private void CalcIdiff(Context.cWeatherdata weather, SunVector[] sunvectors)
        {
            int HOY = 0;
            for (int i = 1; i < 366; i++)
            {
                for (int u = 0; u < 24; u++)
                {
                    CalcIdiff(i, HOY, weather, sunvectors);
                    HOY++;
                }
            }
        }

        /// <summary>
        /// Calculates hourly diffuse radiation on the sensor point for the entire year. Multi-Threading version.
        /// <para>Access: total: this.Idiff[HOY][0]; horizon: this.Idiff[HOY][1]; sky: this.Idiff[HOY][2]; circumsolar: this.Idiff[HOY][3]. HOY ∈ [0, 8759].</para>
        /// </summary>
        /// <param name="weather">Weather data.</param>
        /// <param name="sunvectors">Sunvectors, [0, 8759].</param>
        private void CalcIdiff_MT(Context.cWeatherdata weather, SunVector[] sunvectors)
        {
            Parallel.For(1, 366, i =>
            {
                for (int u = 0; u < 24; u++)
                {
                    int HOY = (i - 1) * 24 + u;
                    CalcIdiff(i, HOY, weather, sunvectors);
                }
            });

        }


        /// <summary>
        /// Calculates beam (direct) irradiation on the sensor point for one hour of the year.
        /// <para>Access: this.Ibeam[HOY]; HOY ∈ [0, 8759].</para>
        /// </summary>
        /// <param name="DOY">Day of the year, ∈ [1, 365].</param>
        /// <param name="HOY">Hour of the year, ∈ [0, 8759].</param>
        /// <param name="LT">Local time, i.e. hour of the day, ∈ [0, 23].</param>        
        /// <param name="weather">Weather data.</param>
        /// <param name="sunvectors">Sunvectors, [0, 8759].</param>
        private void CalcIbeam(int HOY, Context.cWeatherdata weather, SunVector[] sunvectors)
        {
            if (sunvectors[HOY].Sunshine == true)
            {
                for (int i = 0; i < this.SPCount; i++)
                {
                    if (this.sky[i].ShdwBeam[HOY] + Convert.ToDouble(snowcovered[i][HOY]) > 1.0)
                        this.Ibeam[i][HOY] = 0.0;
                    else
                    {
                        this.Ibeam[i][HOY] = Irradiation.Beam(
                            weather.DNI[HOY], sunvectors[HOY].udtCoordinates.dZenithAngle,
                            sunvectors[HOY].udtCoordinates.dAzimuth, this.beta[i], this.psi[i], 
                            this.sky[i].ShdwBeam[HOY]);
                    }
                }
            }
            else
            {
                for (int i = 0; i < this.SPCount; i++)
                    this.Ibeam[i][HOY] = 0;
            }
        }

        /// <summary>
        /// Calculates beam (direct) irradiation on the sensor point for one hour of the year. Multi-Threading version.
        /// <para>Access: this.Ibeam[HOY]; HOY ∈ [0, 8759].</para>
        /// </summary>
        /// <param name="DOY">Day of the year, ∈ [1, 365].</param>
        /// <param name="HOY">Hour of the year, ∈ [0, 8759].</param>
        /// <param name="LT">Local time, i.e. hour of the day, ∈ [0, 23].</param>      
        /// <param name="weather">Weather data.</param>
        /// <param name="sunvectors">Sunvectors, [0, 8759].</param>
        private void CalcIbeam_MT(int HOY, Context.cWeatherdata weather, SunVector[] sunvectors)
        {
            if (sunvectors[HOY].Sunshine == true)
            {
                Parallel.For(0, this.SPCount, i =>
                {
                    if (this.sky[i].ShdwBeam[HOY] + Convert.ToDouble(snowcovered[i][HOY]) > 1.0)
                        this.Ibeam[i][HOY] = 0.0;
                    else
                    {
                        this.Ibeam[i][HOY] = Irradiation.Beam(
                                weather.DNI[HOY], sunvectors[HOY].udtCoordinates.dZenithAngle,
                                sunvectors[HOY].udtCoordinates.dAzimuth, this.beta[i], this.psi[i],
                                this.sky[i].ShdwBeam[HOY]);
                    }
                });
            }
            else
            {
                Parallel.For(0, this.SPCount, i =>
                {
                    this.Ibeam[i][HOY] = 0;
                });
            }
        }

        /// <summary>
        /// Calculates hourly beam (direct) irradiation on the sensor point for the entire year. 
        /// <para>Access: this.Ibeam[HOY]; HOY ∈ [0, 8759].</para>
        /// </summary>      
        /// <param name="weather">Weather data.</param>
        /// <param name="sunvectors">Sunvectors, [0, 8759].</param>
        private void CalcIbeam(Context.cWeatherdata weather, SunVector[] sunvectors)
        {
            for (int t = 0; t < 8760; t++)
            {
                CalcIbeam(t, weather, sunvectors);
            }
        }

        /// <summary>
        /// Calculates hourly beam (direct) irradiation on the sensor point for the entire year. Multi-Threading version.
        /// <para>Access: this.Ibeam[HOY]; HOY ∈ [0, 8759].</para>
        /// </summary>
        /// <param name="weather">Weather data.</param>
        /// <param name="sunvectors">Sunvectors, [0, 8759].</param>
        private void CalcIbeam_MT(Context.cWeatherdata weather, SunVector[] sunvectors)
        {
            Parallel.For(0, 8760, t =>
            {
                CalcIbeam(t, weather, sunvectors);
            });
        }




        #region PublicFunctions

        /// <summary>
        /// Calculates total solar irradiation on a sensor point for one hour of the year.
        /// <para>Access: this.I[HOY]; HOY ∈ [0, 8759]. HOY = (DOY-1) * 24 + LT.</para>
        /// </summary>
        /// <param name="DOY">Day of the year, ∈ [1, 365].</param>
        /// <param name="LT">Local time, i.e. hour of the day, ∈ [0, 23].</param>
        /// <param name="weather">Weather data.</param>
        /// <param name="sunvectors">Sunvectors, [0, 8759].</param>
        public void CalcIrradiation(int DOY, int LT, Context.cWeatherdata weather, SunVector[] sunvectors)
        {
            int HOY = (DOY - 1) * 24 + LT;
            CalcIbeam(HOY, weather, sunvectors);
            CalcIdiff(DOY, HOY, weather, sunvectors);
            for (int i = 0; i < this.I.Length; i++)
                this.I[i][HOY] = this.Ibeam[i][HOY] + this.Idiff[i][HOY] + this.Irefl_spec[i][HOY] + this.Irefl_diff[i][HOY];
        }

        /// <summary>
        /// Calculates total solar irradiation on a sensor point for one hour of the year. Multi-Threading version.
        /// <para>Access: this.I[HOY]; HOY ∈ [0, 8759]. HOY = (DOY-1) * 24 + LT.</para>
        /// </summary>
        /// <param name="DOY">Day of the year, ∈ [1, 365].</param>
        /// <param name="LT">Local time, i.e. hour of the day, ∈ [0, 23].</param>
        /// <param name="weather">Weather data.</param>
        /// <param name="sunvectors">Sunvectors, [0, 8759].</param>
        public void CalcIrradiationMT(int DOY, int LT, Context.cWeatherdata weather, SunVector[] sunvectors)
        {
            int HOY = (DOY - 1) * 24 + LT;

            CalcIbeam_MT(HOY, weather, sunvectors);
            CalcIdiff_MT(DOY, HOY, weather, sunvectors);

            Parallel.For(0, this.I.Length, i =>
            {
                this.I[i][HOY] = this.Ibeam[i][HOY] + this.Idiff[i][HOY] + this.Irefl_spec[i][HOY] + this.Irefl_diff[i][HOY];
            });
        }

        /// <summary>
        /// Calculates hourly total solar irradiation on the sensor point for the entire year.
        /// <para>Access: this.I[HOY]; HOY ∈ [0, 8759].</para>
        /// </summary>
        /// <param name="weather">Weather data.</param>
        /// <param name="sunvectors">Sunvectors, [0, 8759].</param>
        public void CalcIrradiation(Context.cWeatherdata weather, SunVector[] sunvectors)
        {
            CalcIbeam(weather, sunvectors);
            CalcIdiff(weather, sunvectors);
            for (int i = 0; i < this.I.Length; i++)
                for (int t = 0; t < 8760; t++)
                    this.I[i][t] = this.Ibeam[i][t] + this.Idiff[i][t] + this.Irefl_spec[i][t] + this.Irefl_diff[i][t];
        }

        /// <summary>
        /// Calculates hourly total solar irradiation on the sensor point for the entire year. Multi-Threading version.
        /// <para>Access: this.I[HOY]; HOY ∈ [0, 8759].</para>
        /// </summary>
        /// <param name="weather">Weather data.</param>
        /// <param name="sunvectors">Sunvectors, [0, 8759].</param>
        public void CalcIrradiationMT(Context.cWeatherdata weather, SunVector[] sunvectors)
        {
            CalcIbeam_MT(weather, sunvectors);
            CalcIdiff_MT(weather, sunvectors);
            Parallel.For(0, this.I.Length, i =>
            {
                for (int t = 0; t < 8760; t++)
                    this.I[i][t] = this.Ibeam[i][t] + this.Idiff[i][t] + this.Irefl_spec[i][t] + this.Irefl_diff[i][t]; 
            });
        }


        /// <summary>
        /// Set simple sky dome obstruction by taking into account the sensor points tilt angles.
        /// </summary>
        /// <param name="tiltangle">Sensor points tilt angles, in degree.</param>
        public void SetSimpleSky(double[] tiltangle)
        {
            for (int i = 0; i < this.SPCount; i++)
            {
                double visibleHemisphere = (1 + Math.Cos(tiltangle[i] * (Math.PI / 180))) / 2;
                for (int HOY = 0; HOY < 8760; HOY++)
                {
                    this.sky[i].ShdwDome[HOY] = 1 - visibleHemisphere;
                }   
            }
        }


        /// <summary>
        /// Set simplified ground reflection values, based on albedo value of the ground.
        /// </summary>
        /// <remarks>Reflected irradiation is added to the current diffuse inter-reflection value.</remarks>
        /// <param name="tiltangle">Sensor points tilt angles, in degree. [i] = for each sensor point.</param>
        /// <param name="albedo">Albedo of the ground. 8760 time series, [t] = for each hour of the year.</param>
        /// <param name="weather">Weather data.</param>
        /// <param name="sunvectors">8760 sunvectors, [t] = for each hour of the year.</param>
        public void SetSimpleGroundReflection(double [] tiltangle, double [] albedo, Context.cWeatherdata weather, SunVector[] sunvectors)
        {
            for (int i = 0; i < this.SPCount; i++)
            {
                for (int HOY = 0; HOY < 8760; HOY++)
                {
                    this.Irefl_diff[i][HOY] += (albedo[HOY] *
                        (weather.DNI[HOY] * Math.Cos((sunvectors[HOY].udtCoordinates.dZenithAngle) * (Math.PI / 180))
                        + weather.DHI[HOY])) * ((1 - Math.Cos(tiltangle[i] * (Math.PI / 180))) / 2);
                }
            }
        }

        #endregion






        #region ExternalComputingNecessary

        /// <summary>
        /// Applys shadow / obstruction factors from externally calculated view factor calculations to the sensor points.
        /// </summary>
        /// <param name="ShdwBeam_hour">Indicate for one hour of the year, if a sensor point is obstructed from beam radiation (true), or not (false). The list must have booleans for each sensor point.</param>
        /// <param name="ShdwSky">Indicate for each vertex of a sensor point's skydome, if the view between sensor point and vertex is obstructed (true), or not (false). The list has boolean arrays for each sensor point; an array is of length of the sky dome vertex count (this.sky[i].VerticesHemisphere.Count).</param>
        /// <param name="HOY">Hour of the year ∈ [0, 8759].</param>
        public void SetShadows(List<bool> ShdwBeam_hour, List<bool[]> ShdwSky, int HOY)
        {
            for (int i = 0; i < sky.Length; i++)
            {
                for (int u = 0; u < ShdwSky[i].Length; u++)
                {
                    this.sky[i].VertexShadowSphere[this.sky[i].VerticesHemisphere[u]] = Convert.ToDouble(ShdwSky[i][u]);
                }

                this.sky[i].SetShadow_Dome();
                this.sky[i].SetShadow_Horizon();
                this.sky[i].SetShadow_Beam(HOY, Convert.ToDouble(ShdwBeam_hour[i]));
            }
        }

        /// <summary>
        /// Applys shadow / obstruction factors from externally calculated view factor calculations to the sensor points. Multi-threading version.
        /// </summary>
        /// <param name="ShdwBeam_hour">Indicate for one hour of the year, if a sensor point is obstructed from beam radiation (true), or not (false). The list must have booleans for each sensor point.</param>
        /// <param name="ShdwSky">Indicate for each vertex of a sensor point's skydome, if the view between sensor point and vertex is obstructed (true), or not (false). The list has boolean arrays for each sensor point; an array is of length of the sky dome vertex count (this.sky[i].VerticesHemisphere.Count).</param>
        /// <param name="HOY">Hour of the year ∈ [0, 8759].</param>
        public void SetShadowsMT(List<bool> ShdwBeam_hour, List<bool[]> ShdwSky, int HOY)
        {
            Parallel.For(0, sky.Length, i =>
            {
                for (int u = 0; u < ShdwSky[i].Length; u++)
                {
                    this.sky[i].VertexShadowSphere[this.sky[i].VerticesHemisphere[u]] = Convert.ToDouble(ShdwSky[i][u]);
                }

                this.sky[i].SetShadow_Dome();
                this.sky[i].SetShadow_Horizon();
                this.sky[i].SetShadow_Beam(HOY, Convert.ToDouble(ShdwBeam_hour[i]));
            });
        }

        /// <summary>
        /// Applys shadow / obstruction factors from externally calculated view factor calculations to the sensor points.
        /// </summary>
        /// <param name="ShdwBeam_hour">Indicate for one hour of the year, if a sensor point is obstructed from beam radiation (1.0), or not (0.0). The list must have coefficients (taking a value from 0.0 to 1.0) for each sensor point. A coefficient makes sense, e.g. if a semi-transparent object is obstructing the sensor point, like a tree.</param>
        /// <param name="ShdwSky">Indicate for each vertex of a sensor point's skydome, if the view between sensor point and vertex is obstructed (1.0), or not (0.0). The list has arrays of doubles for each sensor point; an array is of length of the sky dome vertex count (this.sky[i].VerticesHemisphere.Count).</param>
        /// <param name="HOY">Hour of the year ∈ [0, 8759].</param>
        public void SetShadows(List<double> ShdwBeam_hour, List<double[]> ShdwSky, int HOY)
        {
            for (int i = 0; i < sky.Length; i++)
            {
                for (int u = 0; u < ShdwSky[i].Length; u++)
                {
                    this.sky[i].VertexShadowSphere[this.sky[i].VerticesHemisphere[u]] = ShdwSky[i][u];
                }

                this.sky[i].SetShadow_Dome();
                this.sky[i].SetShadow_Horizon();
                this.sky[i].SetShadow_Beam(HOY, ShdwBeam_hour[i]);
            }
        }

        /// <summary>
        /// Applys shadow / obstruction factors from externally calculated view factor calculations to the sensor points. Multi-threading version.
        /// </summary>
        /// <param name="ShdwBeam_hour">Indicate for one hour of the year, if a sensor point is obstructed from beam radiation (1.0), or not (0.0). The list must have coefficients (taking a value from 0.0 to 1.0) for each sensor point. A coefficient makes sense, e.g. if a semi-transparent object is obstructing the sensor point, like a tree.</param>
        /// <param name="ShdwSky">Indicate for each vertex of a sensor point's skydome, if the view between sensor point and vertex is obstructed (1.0), or not (0.0). The list has arrays of doubles for each sensor point; an array is of length of the sky dome vertex count (this.sky[i].VerticesHemisphere.Count).</param>
        /// <param name="HOY">Hour of the year ∈ [0, 8759].</param>
        public void SetShadowsMT(List<double> ShdwBeam_hour, List<double[]> ShdwSky, int HOY)
        {
            Parallel.For(0, sky.Length, i =>
            {
                for (int u = 0; u < ShdwSky[i].Length; u++)
                {
                    this.sky[i].VertexShadowSphere[this.sky[i].VerticesHemisphere[u]] = ShdwSky[i][u];
                }

                this.sky[i].SetShadow_Dome();
                this.sky[i].SetShadow_Horizon();
                this.sky[i].SetShadow_Beam(HOY, ShdwBeam_hour[i]);
            });
        }

        /// <summary>
        /// Applys shadow / obstruction factors from externally calculated view factor calculations to the sensor points.
        /// </summary>
        /// <param name="ShdwBeam_hour">Indicate for all 8760 hours of the year, if a sensor point is obstructed from beam radiation (true), or not (false). The list must have booleans for each sensor point.</param>
        /// <param name="ShdwSky">Indicate for each vertex of a sensor point's skydome, if the view between sensor point and vertex is obstructed (true), or not (false). The list has boolean arrays for each sensor point; an array is of length of the sky dome vertex count (this.sky[i].VerticesHemisphere.Count).</param>
        /// <param name="HOY"></param>
        public void SetShadows(List<bool[]> ShdwBeam_hour, List<bool[]> ShdwSky)
        {
            for (int i = 0; i < sky.Length; i++)
            {
                for (int u = 0; u < ShdwSky[i].Length; u++)
                {
                    this.sky[i].VertexShadowSphere[this.sky[i].VerticesHemisphere[u]] = Convert.ToDouble(ShdwSky[i][u]);
                }
                this.sky[i].SetShadow_Dome();
                this.sky[i].SetShadow_Horizon();
                for (int t = 0; t < 8760; t++)
                {
                    this.sky[i].SetShadow_Beam(t, Convert.ToDouble(ShdwBeam_hour[i][t]));
                }
            }
        }

        /// <summary>
        /// Set annual shadows, using interpolation of several calculated days.
        /// </summary>
        /// <param name="StartDays">Array of start days, as day of year 1-365.</param>
        /// <param name="EndDays">Array of end days, as day of year 1-365.</param>
        /// <param name="ShdwBeam">Indicate, if vector is obstructed (true) or not (false). List [i]: each sensor point. [d][h], d=each day used, h = 24 hours.</param>
        /// <param name="BeamPermIs">Indicate, iff vector is obstructed by permeable object (true) or not (false). List [i]: each sensor point. [d][h], d=each day used, h = 24 hours.</param>
        /// <param name="BeamPermRefs">Iff vector is obstructed by permeable object, reference to permeable object extinction coefficients. List [i]: each sensor point. [d][h][p], d=each day used, h = 24 hours, p = reference to each hit permeable object.</param>
        /// <param name="BeamPermLength">Iff vector is obstructed by permeable object, reference to permeable object extinction coefficients. List [i]: each sensor point. [d][h], d=each day used, h = 24 hours, p = reference to each hit permeable object.</param>
        /// <param name="ShdwSky">List [i]: each sensor point. [u] = each sky dome vertex.</param>
        /// <param name="SkyPermIs">Indicate, if and only if it is obstructed by permeable object. List [i]: each sensor point. [u] = each sky dome vertex.</param>
        /// <param name="SkyPermRefs">Iff obstructed by perm object, indicate the index of the permeable object. List [i]: each SP, [u][p], u = each sky dome vertex, p = index to permobject ext.coeff.</param>
        /// <param name="SkyPermLength">Iff obstructed by perm object, indicate the length of penetration through object. List [i]: each SP, [u][p], u = each sky dome vertex, p = index to permobject ext.coeff.</param>
        ///<param name="extinctCoeff">For each permeable object, 8760 time series of extinction coefficients. List [p] = each permeable object, [hoy] = hour of the year 1-8760.</param>
        public void SetShadows_Annual_Permeables(int[] StartDays, int[] EndDays,
            List<bool[][]> ShdwBeam, List<bool[][]> BeamPermIs, List<int[][][]> BeamPermRefs, List<double[][][]> BeamPermLength,
            List<bool[]> ShdwSky, List<bool[]> SkyPermIs, List<int[][]> SkyPermRefs, List<double[][]> SkyPermLength, List<double[]> extinctCoeff)
        {
            for (int i = 0; i < this.SPCount; i++)
            {
                for (int u = 0; u < ShdwSky[i].Length; u++)
                {
                    this.sky[i].VertexShadowSphere[this.sky[i].VerticesHemisphere[u]] = Convert.ToDouble(ShdwSky[i][u]);
                }
                this.sky[i].SetShadow_DomeAndHorizon(SkyPermIs[i], SkyPermRefs[i], SkyPermLength[i], extinctCoeff);
            }


            int dayStart, dayEnd;
            int dayIntervals = EndDays[0] - StartDays[0];
            double maxPi = 2 * Math.PI / Convert.ToDouble(ShdwBeam[0].Length);

            int daysUsed = StartDays.Length;

            for (int i = 0; i < this.SPCount; i++)    
            {
                int dd;
                for (int d = 0; d < daysUsed; d++)
                {
                    if (d == daysUsed - 1)
                        dd = 0;
                    else
                        dd = d + 1;
                    dayStart = StartDays[d];
                    dayEnd = EndDays[d];
                    for (int n = dayStart; n < dayEnd; n++)
                    {
                        double dist1 = (Convert.ToDouble(dayIntervals) - Math.Abs(dayStart - n)) / Convert.ToDouble(dayIntervals);
                        if (d < (daysUsed / 4))
                        {
                            dist1 = 1 - dist1;
                            dist1 = Math.Cos(dist1 * (0.5 * Math.PI));
                        }
                        else if (d < ((daysUsed / 4) * 2) && d >= ((daysUsed / 4) * 1))
                        {
                            dist1 = Math.Sin(dist1 * (0.5 * Math.PI));
                        }
                        else if (d < ((daysUsed / 4) * 3) && d >= ((daysUsed / 4) * 2))
                        {
                            dist1 = 1 - dist1;
                            dist1 = Math.Cos(dist1 * (0.5 * Math.PI));
                        }
                        else if (d < ((daysUsed / 4) * 4) && d >= ((daysUsed / 4) * 3))
                        {
                            dist1 = Math.Sin(dist1 * (0.5 * Math.PI));
                        }

                        double dist2 = 1 - dist1;
                        for (int h = 0; h < 24; h++)
                        {
                            int HOY = (n - 1) * 24 + h;
                            if (BeamPermIs[i][d][h])
                            {
                                double extCoeffSum1 = 0;
                                double extCoeffSum2 = 0;
                                double factor1 = 0;
                                double factor2 = 0;
                                for (int p = 0; p < BeamPermLength[i][d][h].Length; p++)
                                {
                                    extCoeffSum1 += BeamPermLength[i][d][h][p] * extinctCoeff[BeamPermRefs[i][d][h][p]][HOY];
                                }
                                if (extCoeffSum1 > 1) extCoeffSum1 = 1;
                                factor1 = (1 - Convert.ToDouble(ShdwBeam[i][d][h]) * extCoeffSum1) * dist1;
                                if (BeamPermIs[i][dd][h])
                                {
                                    for (int p = 0; p < BeamPermLength[i][dd][h].Length; p++)
                                    {
                                        extCoeffSum2 += BeamPermLength[i][dd][h][p] * extinctCoeff[BeamPermRefs[i][dd][h][p]][HOY];
                                    }
                                    if (extCoeffSum2 > 1) extCoeffSum2 = 1;
                                    factor2 = (1 - Convert.ToDouble(ShdwBeam[i][dd][h]) * extCoeffSum2) * dist2;
                                }
                                else
                                {
                                    factor2 = (1 - Convert.ToDouble(ShdwBeam[i][dd][h])) * dist2;
                                }
                                double shdw = factor1 + factor2;
                                this.sky[i].SetShadow_Beam(HOY, shdw);
                            }
                            else
                            {
                                double factor = ((1 - Convert.ToDouble(ShdwBeam[i][d][h])) * dist1) +
                                    ((1 - Convert.ToDouble(ShdwBeam[i][dd][h])) * dist2);
                                bool shdw = (factor >= 0.5) ? false : true;
                                this.sky[i].SetShadow_Beam(HOY, Convert.ToDouble(shdw));
                            }
                        }
                    }
                }
            }
        }

        /// <summary>
        /// Set annual shadows, using interpolation of three calculated days: Equinox, summer and winter solstice.
        /// </summary>
        /// <param name="StartDays">Array of start days, as day of year 1-365.</param>
        /// <param name="EndDays">Array of end days, as day of year 1-365.</param>
        /// <param name="ShdwBeam_Equinox">Indicate for equinox day, if vector is obstructed (true) or not (false). List [i] = each sensor point. [h] = 24 hours.</param>
        /// <param name="BeamPermIs_Equi">Indicate for equinox day, iff vector is obstructed by permeable object (true) or not (false). List [i] = each sensor point. [h] = 24 hours.</param>
        /// <param name="BeamPermRefs_Equi">For equinox day, iff vector is obstructed by permeable object, reference to permeable object extinction coefficients. List [i]: each sensor point. [h][p], h = 24 hours, p = reference to each hit permeable object.</param>
        /// <param name="BeamPermLength_Equi">For equinox day, iff vector is obstructed by permeable object, reference to permeable object extinction coefficients. List [i]: each sensor point. [h] = 24 hours, p = reference to each hit permeable object.</param>
        ///<param name="ShdwBeam_Summer">Indicate for summer day, if vector is obstructed (true) or not (false). List [i] = each sensor point. [h] = 24 hours.</param>
        ///<param name="BeamPermIs_Sum">Indicate for summer day, iff vector is obstructed by permeable object (true) or not (false). List [i] = each sensor point. [h] = 24 hours.</param>
        ///<param name="BeamPermRefs_Sum">For summer day, iff vector is obstructed by permeable object, reference to permeable object extinction coefficients. List [i]: each sensor point. [h][p], h = 24 hours, p = reference to each hit permeable object.</param>
        ///<param name="BeamPermLength_Sum">For summer day, iff vector is obstructed by permeable object, reference to permeable object extinction coefficients. List [i]: each sensor point. [h] = 24 hours, p = reference to each hit permeable object.</param>
        ///<param name="ShdwBeam_Winter">Indicate for winter day, if vector is obstructed (true) or not (false). List [i] = each sensor point. [h] = 24 hours.</param>
        ///<param name="BeamPermIs_Win">Indicate for winter day, iff vector is obstructed by permeable object (true) or not (false). List [i] = each sensor point. [h] = 24 hours.</param>
        ///<param name="BeamPermRefs_Win">For winter day, iff vector is obstructed by permeable object, reference to permeable object extinction coefficients. List [i]: each sensor point. [h][p], h = 24 hours, p = reference to each hit permeable object.</param>
        ///<param name="BeamPermLength_Win">For winter day, iff vector is obstructed by permeable object, reference to permeable object extinction coefficients. List [i]: each sensor point. [h] = 24 hours, p = reference to each hit permeable object.</param>
        /// <param name="ShdwSky">List [i]: each sensor point. [u] = each sky dome vertex.</param>
        /// <param name="SkyPermIs">Indicate, if and only if it is obstructed by permeable object. List [i]: each sensor point. [u] = each sky dome vertex.</param>
        /// <param name="SkyPermRefs">Iff obstructed by perm object, indicate the index of the permeable object. List [i]: each SP, [u][p], u = each sky dome vertex, p = index to permobject ext.coeff.</param>
        /// <param name="SkyPermLength">Iff obstructed by perm object, indicate the length of penetration through object. List [i]: each SP, [u][p], u = each sky dome vertex, p = index to permobject ext.coeff.</param>
        ///<param name="extinctCoeff">For each permeable object, 8760 time series of extinction coefficients. List [p] = each permeable object, [hoy] = hour of the year 1-8760.</param>
        public void SetShadows_Annual_Permeables(List<bool[]> ShdwBeam_Equinox, List<bool[]> ShdwBeam_Summer, List<bool[]> ShdwBeam_Winter, 
            bool[][] BeamPermIs_Equ, int[][][] BeamPermRefs_Equ, double[][][] BeamPermLength_Equ,
            bool[][] BeamPermIs_Sum, int[][][] BeamPermRefs_Sum, double[][][] BeamPermLength_Sum,
            bool[][] BeamPermIs_Win, int[][][] BeamPermRefs_Win, double[][][] BeamPermLength_Win,
            List<bool[]> ShdwSky, bool[][] SkyPermIs, int[][][] SkyPermRefs, double[][][] SkyPermLength, List<double[]> extinctCoeff)
        {
            for (int i = 0; i < this.SPCount; i++)
            {
                for (int u = 0; u < ShdwSky[i].Length; u++)
                {
                    this.sky[i].VertexShadowSphere[this.sky[i].VerticesHemisphere[u]] = Convert.ToDouble(ShdwSky[i][u]);
                }
                this.sky[i].SetShadow_DomeAndHorizon(SkyPermIs[i], SkyPermRefs[i], SkyPermLength[i], extinctCoeff);
            }


            //    '100% equinox is 20 march and 23 Sept
            //    '100% summer is 21 june
            //    '100% winter is 22 decemebr
            int y1, y2, y3, y4, y5, y6;
            y1 = -9;    //winter solstice
            y2 = 78;    //equinox spring
            y3 = 171;   // summer solstice
            y4 = 265;   //equinox solstice
            y5 = 355;   //winter solstice
            y6 = 443;   //equinox spring

            int fullF1, fullF2;
            int InterpolInterv;
            double dist1, dist2;

            for (int i = 0; i < this.SPCount; i++)
            {
                fullF1 = y1;    //100% shadowwinter on this day
                fullF2 = y2;    //100% shadowequinox on this day
                InterpolInterv = y2 + y1 * -1;
                for (int d = 0; d < y2; d++)        //from 0.Jan to 19.March
                {
                    dist1 = Convert.ToDouble((InterpolInterv - Math.Abs(fullF1 - d))) / Convert.ToDouble(InterpolInterv);
                    dist1 = 1 - dist1;
                    dist1 = Math.Cos(dist1 * (0.5 * Math.PI));
                    dist2 = 1 - dist1;
                    for (int h = 0; h < 24; h++)
                    {
                        bool permPresent = false;
                        int HOY = d * 24 + h;
                        double factor1 = 0;
                        double factor2 = 0;
                        if (BeamPermIs_Win[i][h])     //if and only if obstructed by trees
                        {
                            permPresent = true;
                            double extCoeffSum1 = 0;
                            for (int p = 0; p < BeamPermLength_Win[i][h].Length; p++)
                            {
                                extCoeffSum1 += BeamPermLength_Win[i][h][p] * extinctCoeff[BeamPermRefs_Win[i][h][p]][HOY];
                            }
                            if (extCoeffSum1 > 1) extCoeffSum1 = 1;
                            factor1 = (Convert.ToDouble(BeamPermIs_Win[i][h]) * extCoeffSum1) * dist1;
                            if (BeamPermIs_Equ[i][h])
                            {
                                double extCoeffSum2 = 0;
                                for (int p = 0; p < BeamPermLength_Equ[i][h].Length; p++)
                                {
                                    extCoeffSum2 += BeamPermLength_Equ[i][h][p] * extinctCoeff[BeamPermRefs_Equ[i][h][p]][HOY];
                                }
                                if (extCoeffSum2 > 1) extCoeffSum2 = 1;
                                factor2 = (Convert.ToDouble(BeamPermIs_Equ[i][h]) * extCoeffSum2) * dist2;
                            }
                            else
                            {
                                factor2 = (Convert.ToDouble(ShdwBeam_Equinox[i][h])) * dist2;
                            }
                        }
                        else
                        {
                            factor1 = Convert.ToDouble(ShdwBeam_Winter[i][h]) * dist1;
                            if (BeamPermIs_Equ[i][h])
                            {
                                permPresent = true;
                                double extCoeffSum2 = 0;
                                for (int p = 0; p < BeamPermLength_Equ[i][h].Length; p++)
                                {
                                    extCoeffSum2 += BeamPermLength_Equ[i][h][p] * extinctCoeff[BeamPermRefs_Equ[i][h][p]][HOY];
                                }
                                if (extCoeffSum2 > 1) extCoeffSum2 = 1;
                                factor2 = (Convert.ToDouble(BeamPermIs_Equ[i][h]) * extCoeffSum2) * dist2;
                            }
                            else
                            {
                                factor2 = (Convert.ToDouble(ShdwBeam_Equinox[i][h])) * dist2;
                            }
                        }
                        double shdw = factor1 + factor2;
                        if (permPresent)
                        {
                            if (shdw > 1.0) shdw = 1.0;
                        }
                        else
                        {
                            shdw = (shdw < 0.5) ? 0 : 1;
                        }
                        this.sky[i].SetShadow_Beam(HOY, shdw);
                    }
                }

                fullF1 = y2;        // 100% shadowequinox on this day
                fullF2 = y3;        // 100% shadowsummer on this day
                InterpolInterv = y3 - y2;
                for (int d = y2; d < y3; d++)       //from 20.March to 20.june
                {
                    dist1 = Convert.ToDouble((InterpolInterv - Math.Abs(fullF1 - d))) / Convert.ToDouble(InterpolInterv);
                    dist1 = Math.Sin(dist1 * (0.5 * Math.PI));
                    dist2 = 1 - dist1;
                    for (int h = 0; h < 24; h++)
                    {
                        bool permPresent = false;
                        int HOY = d * 24 + h;
                        double factor1 = 0;
                        double factor2 = 0;
                        if (BeamPermIs_Equ[i][h])     //if and only if obstructed by trees
                        {
                            permPresent = true;
                            double extCoeffSum1 = 0;
                            for (int p = 0; p < BeamPermLength_Equ[i][h].Length; p++)
                            {
                                extCoeffSum1 += BeamPermLength_Equ[i][h][p] * extinctCoeff[BeamPermRefs_Equ[i][h][p]][HOY];
                            }
                            if (extCoeffSum1 > 1) extCoeffSum1 = 1;
                            factor1 = (Convert.ToDouble(BeamPermIs_Equ[i][h]) * extCoeffSum1) * dist1;
                            if (BeamPermIs_Sum[i][h])
                            {
                                double extCoeffSum2 = 0;
                                for (int p = 0; p < BeamPermLength_Sum[i][h].Length; p++)
                                {
                                    extCoeffSum2 += BeamPermLength_Sum[i][h][p] * extinctCoeff[BeamPermRefs_Sum[i][h][p]][HOY];
                                }
                                if (extCoeffSum2 > 1) extCoeffSum2 = 1;
                                factor2 = (Convert.ToDouble(BeamPermIs_Sum[i][h]) * extCoeffSum2) * dist2;
                            }
                            else
                            {
                                factor2 = (Convert.ToDouble(ShdwBeam_Summer[i][h])) * dist2;
                            }
                        }
                        else
                        {
                            factor1 = Convert.ToDouble(ShdwBeam_Equinox[i][h]) * dist1;
                            if (BeamPermIs_Sum[i][h])
                            {
                                permPresent = true;
                                double extCoeffSum2 = 0;
                                for (int p = 0; p < BeamPermLength_Sum[i][h].Length; p++)
                                {
                                    extCoeffSum2 += BeamPermLength_Sum[i][h][p] * extinctCoeff[BeamPermRefs_Sum[i][h][p]][HOY];
                                }
                                if (extCoeffSum2 > 1) extCoeffSum2 = 1;
                                factor2 = (Convert.ToDouble(BeamPermIs_Sum[i][h]) * extCoeffSum2) * dist2;
                            }
                            else
                            {
                                factor2 = (Convert.ToDouble(ShdwBeam_Summer[i][h])) * dist2;
                            }
                        }
                        double shdw = factor1 + factor2;
                        if (permPresent)
                        {
                            if (shdw > 1.0) shdw = 1.0;
                        }
                        else
                        {
                            shdw = (shdw < 0.5) ? 0 : 1;
                        }
                        this.sky[i].SetShadow_Beam(HOY, shdw);
                    }
                }

                fullF1 = y3;        //100% shadowsummer on this day
                fullF2 = y4;        //100% shadowequinox on this day
                InterpolInterv = y4 - y3;
                for (int d = y3; d < y4; d++)       //from 21.June to 22.Sept
                {
                    dist1 = Convert.ToDouble((InterpolInterv - Math.Abs(fullF1 - d))) / Convert.ToDouble(InterpolInterv);
                    dist1 = 1 - dist1;
                    dist1 = Math.Cos(dist1 * (0.5 * Math.PI));
                    dist2 = 1 - dist1;
                    for (int h = 0; h < 24; h++)
                    {
                        bool permPresent = false;
                        int HOY = d * 24 + h;
                        double factor1 = 0;
                        double factor2 = 0;
                        if (BeamPermIs_Sum[i][h])     //if and only if obstructed by trees
                        {
                            permPresent = true;
                            double extCoeffSum1 = 0;
                            for (int p = 0; p < BeamPermLength_Sum[i][h].Length; p++)
                            {
                                extCoeffSum1 += BeamPermLength_Sum[i][h][p] * extinctCoeff[BeamPermRefs_Sum[i][h][p]][HOY];
                            }
                            if (extCoeffSum1 > 1) extCoeffSum1 = 1;
                            factor1 = (Convert.ToDouble(BeamPermIs_Sum[i][h]) * extCoeffSum1) * dist1;
                            if (BeamPermIs_Equ[i][h])
                            {
                                double extCoeffSum2 = 0;
                                for (int p = 0; p < BeamPermLength_Equ[i][h].Length; p++)
                                {
                                    extCoeffSum2 += BeamPermLength_Equ[i][h][p] * extinctCoeff[BeamPermRefs_Equ[i][h][p]][HOY];
                                }
                                if (extCoeffSum2 > 1) extCoeffSum2 = 1;
                                factor2 = (Convert.ToDouble(BeamPermIs_Equ[i][h]) * extCoeffSum2) * dist2;
                            }
                            else
                            {
                                factor2 = (Convert.ToDouble(ShdwBeam_Equinox[i][h])) * dist2;
                            }
                        }
                        else
                        {
                            factor1 = Convert.ToDouble(ShdwBeam_Summer[i][h]) * dist1;
                            if (BeamPermIs_Equ[i][h])
                            {
                                permPresent = true;
                                double extCoeffSum2 = 0;
                                for (int p = 0; p < BeamPermLength_Equ[i][h].Length; p++)
                                {
                                    extCoeffSum2 += BeamPermLength_Equ[i][h][p] * extinctCoeff[BeamPermRefs_Equ[i][h][p]][HOY];
                                }
                                if (extCoeffSum2 > 1) extCoeffSum2 = 1;
                                factor2 = (Convert.ToDouble(BeamPermIs_Equ[i][h]) * extCoeffSum2) * dist2;
                            }
                            else
                            {
                                factor2 = (Convert.ToDouble(ShdwBeam_Equinox[i][h])) * dist2;
                            }
                        }
                        double shdw = factor1 + factor2;
                        if (permPresent)
                        {
                            if (shdw > 1.0) shdw = 1.0;
                        }
                        else
                        {
                            shdw = (shdw < 0.5) ? 0 : 1;
                        }
                        this.sky[i].SetShadow_Beam(HOY, shdw);
                    }
                }

                fullF1 = y4;        //100% shadowequinox on this day
                fullF2 = y5;        //100% shadowwinter on this day
                InterpolInterv = y5 - y4;
                for (int d = y4; d < y5; d++)       //from 23.Sept to 21.Dec
                {
                    dist1 = Convert.ToDouble((InterpolInterv - Math.Abs(fullF1 - d))) / Convert.ToDouble(InterpolInterv);
                    dist1 = Math.Sin(dist1 * (0.5 * Math.PI));
                    dist2 = 1 - dist1;
                    for (int h = 0; h < 24; h++)
                    {
                        bool permPresent = false;
                        int HOY = d * 24 + h;
                        double factor1 = 0;
                        double factor2 = 0;
                        if (BeamPermIs_Equ[i][h])     //if and only if obstructed by trees
                        {
                            permPresent = true;
                            double extCoeffSum1 = 0;
                            for (int p = 0; p < BeamPermLength_Equ[i][h].Length; p++)
                            {
                                extCoeffSum1 += BeamPermLength_Equ[i][h][p] * extinctCoeff[BeamPermRefs_Equ[i][h][p]][HOY];
                            }
                            if (extCoeffSum1 > 1) extCoeffSum1 = 1;
                            factor1 = (Convert.ToDouble(BeamPermIs_Equ[i][h]) * extCoeffSum1) * dist1;
                            if (BeamPermIs_Win[i][h])
                            {
                                double extCoeffSum2 = 0;
                                for (int p = 0; p < BeamPermLength_Win[i][h].Length; p++)
                                {
                                    extCoeffSum2 += BeamPermLength_Win[i][h][p] * extinctCoeff[BeamPermRefs_Win[i][h][p]][HOY];
                                }
                                if (extCoeffSum2 > 1) extCoeffSum2 = 1;
                                factor2 = (Convert.ToDouble(BeamPermIs_Win[i][h]) * extCoeffSum2) * dist2;
                            }
                            else
                            {
                                factor2 = (Convert.ToDouble(ShdwBeam_Winter[i][h])) * dist2;
                            }
                        }
                        else
                        {
                            factor1 = Convert.ToDouble(ShdwBeam_Equinox[i][h]) * dist1;
                            if (BeamPermIs_Win[i][h])
                            {
                                permPresent = true;
                                double extCoeffSum2 = 0;
                                for (int p = 0; p < BeamPermLength_Win[i][h].Length; p++)
                                {
                                    extCoeffSum2 += BeamPermLength_Win[i][h][p] * extinctCoeff[BeamPermRefs_Win[i][h][p]][HOY];
                                }
                                if (extCoeffSum2 > 1) extCoeffSum2 = 1;
                                factor2 = (Convert.ToDouble(BeamPermIs_Win[i][h]) * extCoeffSum2) * dist2;
                            }
                            else
                            {
                                factor2 = (Convert.ToDouble(ShdwBeam_Winter[i][h])) * dist2;
                            }
                        }
                        double shdw = factor1 + factor2;
                        if (permPresent)
                        {
                            if (shdw > 1.0) shdw = 1.0;
                        }
                        else
                        {
                            shdw = (shdw < 0.5) ? 0 : 1;
                        }
                        this.sky[i].SetShadow_Beam(HOY, shdw);
                    }
                }

                fullF1 = y5;        //100% shadowwinter on this day
                fullF2 = y6;        //100% shadowequinox spring on this day
                InterpolInterv = y6 - y5;
                for (int d = y5; d < 365; d++)      //from 22.Dec to 31.Dec
                {
                    dist1 = Convert.ToDouble((InterpolInterv - Math.Abs(fullF1 - d))) / Convert.ToDouble(InterpolInterv);
                    dist1 = 1 - dist1;
                    dist1 = Math.Cos(dist1 * (0.5 * Math.PI));
                    dist2 = 1 - dist1;
                    for (int h = 0; h < 24; h++)
                    {
                        bool permPresent = false;
                        int HOY = d * 24 + h;
                        double factor1 = 0;
                        double factor2 = 0;
                        if (BeamPermIs_Win[i][h])     //if and only if obstructed by trees
                        {
                            permPresent = true;
                            double extCoeffSum1 = 0;
                            for (int p = 0; p < BeamPermLength_Win[i][h].Length; p++)
                            {
                                extCoeffSum1 += BeamPermLength_Win[i][h][p] * extinctCoeff[BeamPermRefs_Win[i][h][p]][HOY];
                            }
                            if (extCoeffSum1 > 1) extCoeffSum1 = 1;
                            factor1 = (Convert.ToDouble(BeamPermIs_Win[i][h]) * extCoeffSum1) * dist1;
                            if (BeamPermIs_Equ[i][h])
                            {
                                double extCoeffSum2 = 0;
                                for (int p = 0; p < BeamPermLength_Equ[i][h].Length; p++)
                                {
                                    extCoeffSum2 += BeamPermLength_Equ[i][h][p] * extinctCoeff[BeamPermRefs_Equ[i][h][p]][HOY];
                                }
                                if (extCoeffSum2 > 1) extCoeffSum2 = 1;
                                factor2 = (Convert.ToDouble(BeamPermIs_Equ[i][h]) * extCoeffSum2) * dist2;
                            }
                            else
                            {
                                factor2 = (Convert.ToDouble(ShdwBeam_Equinox[i][h])) * dist2;
                            }
                        }
                        else
                        {
                            factor1 = Convert.ToDouble(ShdwBeam_Winter[i][h]) * dist1;
                            if (BeamPermIs_Equ[i][h])
                            {
                                permPresent = true;
                                double extCoeffSum2 = 0;
                                for (int p = 0; p < BeamPermLength_Equ[i][h].Length; p++)
                                {
                                    extCoeffSum2 += BeamPermLength_Equ[i][h][p] * extinctCoeff[BeamPermRefs_Equ[i][h][p]][HOY];
                                }
                                if (extCoeffSum2 > 1) extCoeffSum2 = 1;
                                factor2 = (Convert.ToDouble(BeamPermIs_Equ[i][h]) * extCoeffSum2) * dist2;
                            }
                            else
                            {
                                factor2 = (Convert.ToDouble(ShdwBeam_Equinox[i][h])) * dist2;
                            }
                        }
                        double shdw = factor1 + factor2;
                        if (permPresent)
                        {
                            if (shdw > 1.0) shdw = 1.0;
                        }
                        else
                        {
                            shdw = (shdw < 0.5) ? 0 : 1;
                        }
                        this.sky[i].SetShadow_Beam(HOY, shdw);
                    }
                }
            }
        }

        /// <summary>
        /// 3 days interpolation (equinox, summer, winter solstices)
        /// </summary>
        /// <param name="ShdwBeam_Equinox"></param>
        /// <param name="ShdwBeam_Summer"></param>
        /// <param name="ShdwBeam_Winter"></param>
        /// <param name="ShdwSky"></param>
        public void SetShadowsInterpolated(List<bool[]> ShdwBeam_Equinox, List<bool[]> ShdwBeam_Summer, List<bool[]> ShdwBeam_Winter, List<bool[]> ShdwSky)
        {

            for (int i = 0; i < this.SPCount; i++)
            {
                for (int u = 0; u < ShdwSky[i].Length; u++)
                {
                    this.sky[i].VertexShadowSphere[this.sky[i].VerticesHemisphere[u]] = Convert.ToDouble(ShdwSky[i][u]);
                }

                this.sky[i].SetShadow_Dome();
                this.sky[i].SetShadow_Horizon();
            }


            int y1, y2, y3, y4, y5, y6;
            y1 = -9;    //winter solstice
            y2 = 78;    //equinox spring
            y3 = 171;   // summer solstice
            y4 = 265;   //equinox solstice
            y5 = 355;   //winter solstice
            y6 = 443;   //equinox spring

            int fullF1, fullF2;

            double factor = 1.0;
            int InterpolInterv;
            double dist1, dist2;

            for (int i = 0; i < this.SPCount; i++)
            {
                fullF1 = y1;    //100% shadowwinter on this day
                fullF2 = y2;    //100% shadowequinox on this day
                InterpolInterv = y2 + y1 * -1;
                for (int d = 0; d < y2; d++)  //from 0.Jan to 19.March
                {
                    dist1 = Convert.ToDouble((InterpolInterv - Math.Abs(fullF1 - d))) / Convert.ToDouble(InterpolInterv);
                    dist1 = 1 - dist1;
                    dist1 = Math.Cos(dist1 * (0.5 * Math.PI));
                    dist2 = 1 - dist1;
                    for (int u = 0; u < 24; u++)
                    {
                        factor = ((1 - Convert.ToDouble(ShdwBeam_Winter[i][u])) * dist1) +
                            ((1 - Convert.ToDouble(ShdwBeam_Equinox[i][u])) * dist2);
                        bool shdw = (factor >= 0.5) ? false : true;
                        int HOY = d * 24 + u;
                        this.sky[i].SetShadow_Beam(HOY, Convert.ToDouble(shdw));
                    }
                }

                fullF1 = y2;        // 100% shadowequinox on this day
                fullF2 = y3;        // 100% shadowsummer on this day
                InterpolInterv = y3 - y2;
                for (int d = y2; d < y3; d++)  //from 20.March to 20.june
                {
                    dist1 = Convert.ToDouble((InterpolInterv - Math.Abs(fullF1 - d))) / Convert.ToDouble(InterpolInterv);
                    dist1 = Math.Sin(dist1 * (0.5 * Math.PI));
                    dist2 = 1 - dist1;
                    for (int u = 0; u < 24; u++)
                    {
                        factor = ((1 - Convert.ToDouble(ShdwBeam_Equinox[i][u])) * dist1) +
                            ((1 - Convert.ToDouble(ShdwBeam_Summer[i][u])) * dist2);
                        bool shdw = (factor >= 0.5) ? false : true;
                        int HOY = d * 24 + u;
                        this.sky[i].SetShadow_Beam(HOY, Convert.ToDouble(shdw));
                    }
                }

                fullF1 = y3;        //100% shadowsummer on this day
                fullF2 = y4;        //100% shadowequinox on this day
                InterpolInterv = y4 - y3;
                for (int d = y3; d < y4; d++)  //from 21.June to 22.Sept
                {
                    dist1 = Convert.ToDouble((InterpolInterv - Math.Abs(fullF1 - d))) / Convert.ToDouble(InterpolInterv);
                    dist1 = 1 - dist1;
                    dist1 = Math.Cos(dist1 * (0.5 * Math.PI));
                    dist2 = 1 - dist1;
                    for (int u = 0; u < 24; u++)
                    {
                        factor = ((1 - Convert.ToDouble(ShdwBeam_Summer[i][u])) * dist1) +
                            ((1 - Convert.ToDouble(ShdwBeam_Equinox[i][u])) * dist2);
                        bool shdw = (factor >= 0.5) ? false : true;
                        int HOY = d * 24 + u;
                        this.sky[i].SetShadow_Beam(HOY, Convert.ToDouble(shdw));
                    }
                }

                fullF1 = y4;    //100% shadowequinox on this day
                fullF2 = y5;    //100% shadowwinter on this day
                InterpolInterv = y5 - y4;
                for (int d = y4; d < y5; d++)   //from 23.Sept to 21.Dec
                {
                    dist1 = Convert.ToDouble((InterpolInterv - Math.Abs(fullF1 - d))) / Convert.ToDouble(InterpolInterv);
                    dist1 = Math.Sin(dist1 * (0.5 * Math.PI));
                    dist2 = 1 - dist1;
                    for (int u = 0; u < 24; u++)
                    {
                        factor = ((1 - Convert.ToDouble(ShdwBeam_Equinox[i][u])) * dist1) +
                            ((1 - Convert.ToDouble(ShdwBeam_Winter[i][u])) * dist2);
                        bool shdw = (factor >= 0.5) ? false : true;
                        int HOY = d * 24 + u;
                        this.sky[i].SetShadow_Beam(HOY, Convert.ToDouble(shdw));
                    }
                }

                fullF1 = y5;    //100% shadowwinter on this day
                fullF2 = y6;    //100% shadowequinox spring on this day
                InterpolInterv = y6 - y5;
                for (int d = y5; d < 365; d++)  //from 22.Dec to 31.Dec
                {
                    dist1 = Convert.ToDouble((InterpolInterv - Math.Abs(fullF1 - d))) / Convert.ToDouble(InterpolInterv);
                    dist1 = 1 - dist1;
                    dist1 = Math.Cos(dist1 * (0.5 * Math.PI));
                    dist2 = 1 - dist1;
                    for (int u = 0; u < 24; u++)
                    {
                        factor = ((1 - Convert.ToDouble(ShdwBeam_Winter[i][u])) * dist1) +
                            ((1 - Convert.ToDouble(ShdwBeam_Equinox[i][u])) * dist2);
                        bool shdw = (factor >= 0.5) ? false : true;
                        int HOY = d * 24 + u;
                        this.sky[i].SetShadow_Beam(HOY, Convert.ToDouble(shdw));
                    }
                }

            }

        }

        /// <summary>
        /// Multiple days interpolation (e.g. 12 days: 1st day of each month)
        /// </summary>
        /// <param name="StartDays"></param>
        /// <param name="EndDays"></param>
        /// <param name="ShdwBeam"></param>
        /// <param name="ShdwSky"></param>
        public void SetShadowsInterpolated(int[] StartDays, int[] EndDays, List<bool[][]> ShdwBeam, List<bool[]> ShdwSky)
        {
            for (int i = 0; i < sky.Length; i++)
            {
                for (int u = 0; u < ShdwSky[i].Length; u++)
                {
                    this.sky[i].VertexShadowSphere[this.sky[i].VerticesHemisphere[u]] = Convert.ToDouble(ShdwSky[i][u]);
                }
                this.sky[i].SetShadow_Dome();
                this.sky[i].SetShadow_Horizon();
            }


            int dayStart, dayEnd;
            int dayIntervals = EndDays[0] - StartDays[0];
            double maxPi = 2 * Math.PI / Convert.ToDouble(ShdwBeam[0].Length);

            int daysUsed = StartDays.Length;

            for (int i = 0; i < this.sky.Length; i++)    //foreach sensor point
            {
                int dd;
                for (int d = 0; d < daysUsed; d++)
                {
                    if (d == daysUsed - 1)
                        dd = 0;
                    else
                        dd = d + 1;
                    dayStart = StartDays[d];
                    dayEnd = EndDays[d];
                    for (int n = dayStart; n < dayEnd; n++)
                    {
                        double dist1 = (Convert.ToDouble(dayIntervals) - Math.Abs(dayStart - n)) / Convert.ToDouble(dayIntervals);
                        if (d < (daysUsed / 4))
                        {
                            dist1 = 1 - dist1;
                            dist1 = Math.Cos(dist1 * (0.5 * Math.PI));
                        }
                        else if (d < ((daysUsed / 4) * 2) && d >= ((daysUsed / 4) * 1))
                        {
                            dist1 = Math.Sin(dist1 * (0.5 * Math.PI));
                        }
                        else if (d < ((daysUsed / 4) * 3) && d >= ((daysUsed / 4) * 2))
                        {
                            dist1 = 1 - dist1;
                            dist1 = Math.Cos(dist1 * (0.5 * Math.PI));
                        }
                        else if (d < ((daysUsed / 4) * 4) && d >= ((daysUsed / 4) * 3))
                        {
                            dist1 = Math.Sin(dist1 * (0.5 * Math.PI));
                        }

                        double dist2 = 1 - dist1;
                        for (int u = 0; u < 24; u++)
                        {
                            double factor = ((1 - Convert.ToDouble(ShdwBeam[i][d][u])) * dist1) +
                                ((1 - Convert.ToDouble(ShdwBeam[i][dd][u])) * dist2);
                            bool shdw = (factor >= 0.5) ? false : true;
                            int HOY = (n - 1) * 24 + u;
                            this.sky[i].SetShadow_Beam(HOY, Convert.ToDouble(shdw));
                        }
                    }
                }
            }
        }

        /// <summary>
        /// 3 days interpolation (equinox, summer, winter solstices)
        /// </summary>
        /// <param name="ShdwBeam_Equinox"></param>
        /// <param name="ShdwBeam_Summer"></param>
        /// <param name="ShdwBeam_Winter"></param>
        /// <param name="ShdwSky"></param>
        public void SetShadowsInterpolatedMT(List<bool[]> ShdwBeam_Equinox, List<bool[]> ShdwBeam_Summer, List<bool[]> ShdwBeam_Winter, List<bool[]> ShdwSky)
        {

            Parallel.For(0, sky.Length, i =>
            {
                for (int u = 0; u < ShdwSky[i].Length; u++)
                {
                    this.sky[i].VertexShadowSphere[this.sky[i].VerticesHemisphere[u]] = Convert.ToDouble(ShdwSky[i][u]);
                }
                this.sky[i].SetShadow_Dome();
                this.sky[i].SetShadow_Horizon();
            });

            //    '100% equinox is 20 march and 23 Sept
            //    '100% summer is 21 june
            //    '100% winter is 22 decemebr
            int y1, y2, y3, y4, y5, y6;
            y1 = -9;    //winter solstice
            y2 = 78;    //equinox spring
            y3 = 171;   // summer solstice
            y4 = 265;   //equinox solstice
            y5 = 355;   //winter solstice
            y6 = 443;   //equinox spring


            Parallel.For(0, this.sky.Length, i =>
            {
                int fullF1, fullF2;
                double factor = 1.0;
                int InterpolInterv;
                double dist1, dist2;
                //    fullF1 = y1                    '100% shadowwinter on this day
                //    fullF2 = y2                     '100% shadowequinox on this day
                fullF1 = y1;    //100% shadowwinter on this day
                fullF2 = y2;    //100% shadowequinox on this day
                InterpolInterv = y2 + y1 * -1;
                //    For i = 0 To y2 - 1             'from 0.Jan to 19.March
                for (int d = 0; d < y2; d++)
                {
                    dist1 = Convert.ToDouble((InterpolInterv - Math.Abs(fullF1 - d))) / Convert.ToDouble(InterpolInterv);
                    dist1 = 1 - dist1;
                    dist1 = Math.Cos(dist1 * (0.5 * Math.PI));
                    dist2 = 1 - dist1;
                    for (int u = 0; u < 24; u++)
                    {
                        factor = ((1 - Convert.ToDouble(ShdwBeam_Winter[i][u])) * dist1) +
                            ((1 - Convert.ToDouble(ShdwBeam_Equinox[i][u])) * dist2);
                        bool shdw = (factor >= 0.5) ? false : true;
                        int HOY = d * 24 + u;
                        this.sky[i].SetShadow_Beam(HOY, Convert.ToDouble(shdw));
                    }
                }


                //    fullF1 = y2                      '100% shadowequinox on this day
                //    fullF2 = y3                     '100% shadowsummer on this day
                fullF1 = y2;        // 100% shadowequinox on this day
                fullF2 = y3;        // 100% shadowsummer on this day
                InterpolInterv = y3 - y2;
                //    For i = y2 To y3 - 1            'from 20.March to 20.june
                for (int d = y2; d < y3; d++)
                {
                    dist1 = Convert.ToDouble((InterpolInterv - Math.Abs(fullF1 - d))) / Convert.ToDouble(InterpolInterv);
                    dist1 = Math.Cos(dist1 * (0.5 * Math.PI));
                    dist2 = 1 - dist1;
                    for (int u = 0; u < 24; u++)
                    {
                        factor = ((1 - Convert.ToDouble(ShdwBeam_Equinox[i][u])) * dist1) +
                            ((1 - Convert.ToDouble(ShdwBeam_Summer[i][u])) * dist2);
                        bool shdw = (factor >= 0.5) ? false : true;
                        int HOY = d * 24 + u;
                        this.sky[i].SetShadow_Beam(HOY, Convert.ToDouble(shdw));
                    }
                }

                //    fullF1 = y3                      '100% shadowsummer on this day
                //    fullF2 = y4                     '100% shadowequinox on this day
                fullF1 = y3;        //100% shadowsummer on this day
                fullF2 = y4;        //100% shadowequinox on this day
                InterpolInterv = y4 - y3;
                //    For i = y3 To y4 - 1            'from 21.June to 22.Sept
                for (int d = y3; d < y4; d++)
                {
                    dist1 = Convert.ToDouble((InterpolInterv - Math.Abs(fullF1 - d))) / Convert.ToDouble(InterpolInterv);
                    dist1 = 1 - dist1;
                    dist1 = Math.Cos(dist1 * (0.5 * Math.PI));
                    dist2 = 1 - dist1;
                    for (int u = 0; u < 24; u++)
                    {
                        factor = ((1 - Convert.ToDouble(ShdwBeam_Summer[i][u])) * dist1) +
                            ((1 - Convert.ToDouble(ShdwBeam_Equinox[i][u])) * dist2);
                        bool shdw = (factor >= 0.5) ? false : true;
                        int HOY = d * 24 + u;
                        this.sky[i].SetShadow_Beam(HOY, Convert.ToDouble(shdw));
                    }
                }
                //    fullF1 = y4                      '100% shadowequinox on this day
                //    fullF2 = y5                     '100% shadowwinter on this day
                //    InterpolInterv = y5 - y4
                fullF1 = y4;
                fullF2 = y5;
                InterpolInterv = y5 - y4;
                //    For i = y4 To y5 - 1            'from 23.Sept to 21.Dec
                for (int d = y4; d < y5; d++)
                {
                    dist1 = Convert.ToDouble((InterpolInterv - Math.Abs(fullF1 - d))) / Convert.ToDouble(InterpolInterv);
                    dist1 = Math.Cos(dist1 * (0.5 * Math.PI));
                    dist2 = 1 - dist1;
                    for (int u = 0; u < 24; u++)
                    {
                        factor = ((1 - Convert.ToDouble(ShdwBeam_Equinox[i][u])) * dist1) +
                            ((1 - Convert.ToDouble(ShdwBeam_Winter[i][u])) * dist2);
                        bool shdw = (factor >= 0.5) ? false : true;
                        int HOY = d * 24 + u;
                        this.sky[i].SetShadow_Beam(HOY, Convert.ToDouble(shdw));
                    }
                }

                //    fullF1 = y5                      '100% shadowwinter on this day
                //    fullF2 = y6                     '100% shadowequinox spring on this day
                //    InterpolInterv = y6 - y5
                fullF1 = y5;
                fullF2 = y6;
                InterpolInterv = y6 - y5;
                //    For i = y5 To 364            'from 22.Dec to 31.Dec
                for (int d = y5; d < 365; d++)
                {
                    dist1 = Convert.ToDouble((InterpolInterv - Math.Abs(fullF1 - d))) / Convert.ToDouble(InterpolInterv);
                    dist1 = 1 - dist1;
                    dist1 = Math.Cos(dist1 * (0.5 * Math.PI));
                    dist2 = 1 - dist1;
                    for (int u = 0; u < 24; u++)
                    {
                        factor = ((1 - Convert.ToDouble(ShdwBeam_Winter[i][u])) * dist1) +
                            ((1 - Convert.ToDouble(ShdwBeam_Equinox[i][u])) * dist2);
                        bool shdw = (factor >= 0.5) ? false : true;
                        int HOY = d * 24 + u;
                        this.sky[i].SetShadow_Beam(HOY, Convert.ToDouble(shdw));
                    }
                }
            });
        }


        /// <summary>
        /// Multiple days interpolation (e.g. 12 days: 1st day of each month)
        /// </summary>
        /// <param name="StartDays"></param>
        /// <param name="EndDays"></param>
        /// <param name="ShdwBeam"></param>
        /// <param name="ShdwSky"></param>
        public void SetShadowsInterpolatedMT(int[] StartDays, int[] EndDays, List<bool[][]> ShdwBeam, List<bool[]> ShdwSky)
        {
            Parallel.For(0, sky.Length, i =>
            {
                for (int u = 0; u < ShdwSky[i].Length; u++)
                {
                    this.sky[i].VertexShadowSphere[this.sky[i].VerticesHemisphere[u]] = Convert.ToDouble(ShdwSky[i][u]);
                }
                this.sky[i].SetShadow_Dome();
                this.sky[i].SetShadow_Horizon();
            });



            int dayIntervals = EndDays[0] - StartDays[0];
            double maxPi = 2 * Math.PI / Convert.ToDouble(ShdwBeam[0].Length);

            int daysUsed = StartDays.Length;//ShdwBeam[0].Length;

            Parallel.For(0, this.sky.Length, i =>    //foreach sensor point
            {
                int dd;
                int dayStart, dayEnd;

                for (int d = 0; d < daysUsed; d++)
                {
                    if (d == daysUsed - 1)
                        dd = 0;
                    else
                        dd = d + 1;
                    dayStart = StartDays[d];
                    dayEnd = EndDays[d];
                    for (int n = dayStart; n < dayEnd; n++)
                    {
                        double dist1 = (Convert.ToDouble(dayIntervals) - Math.Abs(dayStart - n)) / Convert.ToDouble(dayIntervals);
                        if (d < (daysUsed / 4))
                        {
                            dist1 = 1 - dist1;
                            dist1 = Math.Cos(dist1 * (0.5 * Math.PI));
                        }
                        else if (d < ((daysUsed / 4) * 2) && d >= ((daysUsed / 4) * 1))
                        {
                            dist1 = Math.Sin(dist1 * (0.5 * Math.PI));
                        }
                        else if (d < ((daysUsed / 4) * 3) && d >= ((daysUsed / 4) * 2))
                        {
                            dist1 = 1 - dist1;
                            dist1 = Math.Cos(dist1 * (0.5 * Math.PI));
                        }
                        else if (d < ((daysUsed / 4) * 4) && d >= ((daysUsed / 4) * 3))
                        {
                            dist1 = Math.Sin(dist1 * (0.5 * Math.PI));
                        }
                        double dist2 = 1 - dist1;
                        for (int u = 0; u < 24; u++)
                        {
                            double factor = ((1 - Convert.ToDouble(ShdwBeam[i][d][u])) * dist1) +
                                ((1 - Convert.ToDouble(ShdwBeam[i][dd][u])) * dist2);
                            bool shdw = (factor >= 0.5) ? false : true;
                            int HOY = (n - 1) * 24 + u;
                            this.sky[i].SetShadow_Beam(HOY, Convert.ToDouble(shdw));
                        }
                    }
                }
            });
        }

        /// <summary>
        /// Goes through each hour of the year and applys snow cover, if surface angle is flat enough and if the weather data indicates snow on that hour.
        /// </summary>
        /// <param name="snow_threshold">Snow threshold, after which no radiation is assumed to reach the sensor point.</param>
        /// <param name="tilt_treshold">Sensor point tilt threshold (degree). More flat angles will not allow irradiation. Steeper angles are assumed to let the snow slide down the sensor point.</param>
        /// <param name="weather">Weather data.</param>
        public void SetSnowcover(double snow_threshold, double tilt_treshold, Context.cWeatherdata weather)
        {
            for (int i = 0; i < this.I.Length; i++)
            {
                for (int t = 0; t < 8760; t++)
                {
                    if (weather.Snow[t] > snow_threshold && beta[i] < tilt_treshold)
                    {
                        this.snowcovered[i][t] = true;
                    }
                }
            }

        }

        /// <summary>
        /// Goes through each hour of the year and applys snow cover, if surface angle is flat enough and if the weather data indicates snow on that hour. Multi-threading version.
        /// </summary>
        /// <param name="snow_threshold">Snow threshold, after which no radiation is assumed to reach the sensor point.</param>
        /// <param name="tilt_treshold">Sensor point tilt threshold (degree). More flat angles will not allow irradiation. Steeper angles are assumed to let the snow slide down the sensor point.</param>
        /// <param name="weather">Weather data.</param>
        public void SetSnowcoverMT(double snow_threshold, double tilt_treshold, Context.cWeatherdata weather)
        {
            Parallel.For(0, this.I.Length, i =>
            {
                for (int t = 0; t < 8760; t++)
                {
                    if (weather.Snow[t] > snow_threshold && beta[i] < tilt_treshold)
                    {
                        this.snowcovered[i][t] = true;
                    }
                }
            });
        }



        /// <summary>
        /// Set irradiation by inter-reflections for one hour of the year. 
        /// Actual irradiation calculation needs to be done externally.
        /// </summary>
        /// <param name="HOY">Hour of the year, HOY ∈ [0, 8759]</param>
        /// <param name="_Ispecular">Irradiation values by specular reflection for one hour of the year and for each sensor point.</param>
        /// <param name="_Idiffuse">Irradiation values by diffuse reflection for each sensor point.</param>        
        public void SetInterreflection(int HOY, double[] _Ispecular, double[] _Idiffuse)
        {
            for (int i = 0; i < this.SPCount; i++)
            {
                if (this.snowcovered[i][HOY])
                {
                    this.Irefl_spec[i][HOY] = 0.0;
                    this.Irefl_diff[i][HOY] = 0.0;
                }
                else
                {
                    this.Irefl_spec[i][HOY] = _Ispecular[i];
                    this.Irefl_diff[i][HOY] = _Idiffuse[i];
                }
            }
        }

        /// <summary>
        /// Set irradiation by inter-reflections for one hour of the year. Multi-threading version.
        /// Actual irradiation calculation needs to be done externally.
        /// </summary>
        /// <param name="HOY">Hour of the year, HOY ∈ [0, 8759]</param>
        /// <param name="_Ispecular">Irradiation values by specular reflection for one hour of the year and for each sensor point.</param>
        /// <param name="_Idiffuse">Irradiation values by diffuse reflection for each sensor point.</param>
        public void SetInterreflectionMT(int HOY, double[] _Ispecular, double[] _Idiffuse)
        {
            Parallel.For(0, this.SPCount, i =>
            {
                if (this.snowcovered[i][HOY])
                {
                    this.Irefl_spec[i][HOY] = 0.0;
                    this.Irefl_diff[i][HOY] = 0.0;
                }
                else
                {
                    this.Irefl_spec[i][HOY] = _Ispecular[i];
                    this.Irefl_diff[i][HOY] = _Idiffuse[i];
                }
            });
        }

        /// <summary>
        /// Set irradiation by inter-reflections for all hours of the year, using interpolation for specular reflection. 
        /// Actual irradiation calculation needs to be done externally.
        /// </summary>
        /// <param name="_IreflSpecular">Irradiation values by specular reflection for each sensor point and each hour of the year.</param>
        /// <param name="_IreflDiffuse">Irradiation values by diffuse reflection for each sensor point and each hour of the year.</param>
        public void SetInterrefl_Annual(double[][] _IreflSpecular, double[][] _IreflDiffuse)
        {
            for (int i = 0; i < this.SPCount; i++)
            {
                for(int t=0; t<8760; t++)
                {
                    this.Irefl_spec[i][t] = _IreflSpecular[i][t];
                    this.Irefl_diff[i][t] = _IreflDiffuse[i][t];
                }
            }
        }



        #endregion




        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        //public void SetInterreflection()
        //{

        //    // decomposed between SPECULAR and DIFFUSE interreflections
        //    // limit to bounce = 1



        //    //  rhino obstructions forall t. use INTERPOLATION 
        //    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //    /////////////////////////////////////////////////     S P E C U L A R   /////////////////////////////////////////////////////////////
        //    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //    //SPECULAR
        //    // double total_Specular_Budget = 0;
        //    // forall specular objects:

        //    // 01:
        //    //  from sensorpoint (SP), shoot rays ONLY to specular objects. 
        //    //      - check, if vectorangle(ray, sensorpoint normal) <= 90°. if not, break. else:
        //    //          - check, if ray unobstructed. if not, break. else, go to 02.

        //    // 02:
        //    // if ray <= 90° (otherwise it hits the backside of the specular surface) and unobstructed:
        //    //      - if bounce > 1 (else go to 03):
        //    //          - for each bounce > 1, calc reflected ray on the specular obstacle and for new ray, go to 01:
        //    //              - go to 03, if unobstructed:

        //    // 03:
        //    //      - for all t in [0,8759]
        //    //          - check, if daytime. if not, break. else:
        //    //              - from specular obstacle object, make ray to sunvector. calc vectorangle(sunvector, obstacle normal). if not <= 90°, break. else:
        //    //                  - make reflection of sunvector on obstacle. if this reflecting vector *-1 is not coincident with connecting ray to SP (tolerance of +- 1°?), break. else:
        //    //                      - check if the sunvector (not the reflected) is obstructed. if yes, break. else, go to 03.
            
        //    // 03:
        //    // calc Specular_interreflection = I_direct incident on obstacle object. Multiply with its specular coefficient (albedo?). 
        //    //      and account for incidence angle on sensorpoint: I_sensorpoint = I_incident sin(90°-vectorangle(sp_normal, incident_ray))
        //    // total_Specular_Budget += Specular_interreflection;


        //    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




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
        //}




     


        ///// <summary>
        ///// Sets the obstruction factors (i.e. shadows) of each sensor point.
        ///// <para>Split into shadows of diffuse radiation (further split into horizon and dome) and beam radiation.</para>
        ///// <para>Using three days for beam radiation shadow calculation, and interpolates for the other days of the year.</para>
        ///// </summary>
        ///// <param name="SunVShdw_equinox">List (for each sensorpoint) of boolean arrays with each 24 elements (for each hour of the day) indicating for the equinox day, wether beam radiation hits the sensorpoint.</param>
        ///// <param name="SunVShdw_summer">List (for each sensorpoint) of boolean arrays with each 24 elements (for each hour of the day) indicating for the summer solstice day, wether beam radiation hits the sensorpoint.</param>
        ///// <param name="SunVShdw_winter">List (for each sensorpoint) of boolean arrays with each 24 elements (for each hour of the day) indicating for the winter solstice day, wether beam radiation hits the sensorpoint.</param>
        ///// <param name="HorizonShdw">List (for each sensorpoint) of boolean arrays with each n elements (for each segment of the horizon, accessible with Sensorpoints[i].sky.VerticesHorizon and .HorizonSegments) indicating, wether the sensorpoint has free view to the respective horizon segment.</param>
        ///// <param name="SkyShdw">List (for each sensorpoint) of boolean arrays with each m elements (for each vertex of the skydome, accessible with Sensorpoints[i].sky.VerticesHemisphere) indicating, wether the sensorpoint has free view to the respective skydome vertex.</param>
        //public void SetShadows(List<bool[]> SunVShdw_equinox, List<bool[]> SunVShdw_summer, List<bool[]> SunVShdw_winter, List<bool[]> HorizonShdw, List<bool[]> SkyShdw)
        //{
        //    //List: foreach sensorpoint in sensorpoints; bool[24] foreach hour of day.

        //    //!!!!!!!!!! quite trivial. rhino will provide list of vectors for obstruction. INTERPOLATION!!! take 3 days! and itnerpolate in-between.
        //    // I wonder if a meta-model would work, which takes 3 calculated days as features (and longi+latitude) and tells me for the rest days, if the rays are also blocked. I cuold try that quite quickly, coz I can calculate the training data in rhino. Just run 8760 obstruction checks. Then, meta model in matlab, python, Accord.net or whatever
        //    //sky[i].SetShadow_SunVector(

        //    //!!!!!!!!!!!!!!!!!!!!   these factors need to be calculated. rhino will give me bools of faces of the dome or horizon obstructed.
        //    //some calc needs to be further done to account for the regions, that are blocked anyway. these should be ommited for the factor calcualtion! only consider those faces, which would also be hit by sun rays 
        //    //sky[i].SetShadow_Dome(
        //    //sky[i].SetShadow_Horizon(






        //    for (int i = 0; i < this.I.Length; i++)
        //    {
        //        for (int t = 0; t < 8760; t++)
        //        {
        //            this.Ibeam[i][t] *= (1 - Convert.ToInt32(sky[i].ShdwBeam[t]));
        //            this.Idiff[i][t][3] *= (1 - Convert.ToInt32(sky[i].ShdwBeam[t]));
        //            this.Idiff[i][t][2] *= (1 - sky[i].ShdwDome);
        //            this.Idiff[i][t][1] *= (1 - sky[i].ShdwHorizon);
        //            this.Idiff[i][t][0] = this.Idiff[i][t][1] + this.Idiff[i][t][2] + this.Idiff[i][t][3];
        //            this.I[i][t] = this.Ibeam[i][t] + this.Idiff[i][t][0];
        //        }
        //    }
        //}

        //!!!!!!!!!!! SetShadows(DOY, LT)
    }
}
