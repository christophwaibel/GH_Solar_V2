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
    public class Sensorpoints
    {
        /// <summary>
        /// Sky hemisphere. See class "SykDome".
        /// </summary>
        public SkyDome[] sky { get; private set; }
        /// <summary>
        /// 8760 hourly annual sun vectors. See class "Sunvectors".
        /// </summary>
        public List<SunVector> sunvectors { get; private set; }
        /// <summary>
        /// Location data. See structure "Context.cLocation".
        /// </summary>
        public Context.cLocation location { get; private set; }
        /// <summary>
        /// Weather data. See "Context.cWeatherdata".
        /// </summary>
        public Context.cWeatherdata weather { get; private set; }
        /// <summary>
        /// Sensor point tilt angle in degree.
        /// </summary>
        public double[] beta { get; private set; }
        /// <summary>
        /// Sensor point azimuth angle in degree.
        /// </summary>
        public double[] psi { get; private set; }
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
        /// <para>Array with 8760 doubles, for each hour of the year.</para>
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
        /// <param name="year">Year to calculate.</param>
        /// <param name="weather">Weather data. See structure "Context.cWeatherdata".</param>
        /// <param name="location">Location data. See structure "Context.cLocation".</param>
        /// <param name="sunvectors">8760 sun vectors, for each hour of the year. See class "SunVector".</param>
        /// <param name="beta">Tilt angle of sensor point.</param>
        /// <param name="psi">Azimuth angle of sensor point.</param>
        /// <param name="reclvlsky">Recursion level for sky hemisphere. 1 or 2 recommended. See class "SkyDome".</param>
        public Sensorpoints(Context.cWeatherdata weather, Context.cLocation location, List<SunVector> sunvectors, double[] beta, double[] psi, int reclvlsky)
        {
            this.SPCount = beta.Length;

            this.location = location;
            this.weather = weather;
            this.sunvectors = new List<SunVector>(sunvectors);

            this.beta = new double[beta.Length];
            this.psi = new double[psi.Length];
            this.sky = new SkyDome[this.SPCount];
            Array.Copy(beta, this.beta, beta.Length);
            Array.Copy(psi, this.psi, psi.Length);
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
        /// Calculates diffuse irradiation on the sensor point for one hour of the year.
        /// <para>Access: total diffuse irradiation: this.Idiff[HOY][0]; horizon: this.Idiff[HOY][1]; sky: this.Idiff[HOY][2]; circumsolar: this.Idiff[HOY][3]. HOY ∈ [0, 8759].</para>
        /// </summary>
        /// <param name="DOY">Day of the year, ∈ [1, 365].</param>
        /// <param name="HOY">Hour of the year, ∈ [0, 8759].</param>
        private void CalcIdiff(int DOY, int HOY)
        {
            for (int i = 0; i < this.SPCount; i++)
            {
                if (this.snowcovered[i][HOY])
                    this.Idiff[i][HOY] = 0.0;
                else
                {
                    this.Idiff[i][HOY] = Irradiation.Diffuse(
                        this.weather.DHI[HOY], this.weather.DNI[HOY], this.sunvectors[HOY].udtCoordinates.dZenithAngle,
                        this.sunvectors[HOY].udtCoordinates.dAzimuth, this.beta[i], this.psi[i], DOY,
                        sky[i].ShdwHorizon, sky[i].ShdwDome, sky[i].ShdwBeam[HOY]);
                }
            }
        }

        /// <summary>
        /// Calculates diffuse irradiation on the sensor point for one hour of the year. Multi-Threading version.
        /// <para>Access: total: this.Idiff[HOY][0]; horizon: this.Idiff[HOY][1]; sky: this.Idiff[HOY][2]; circumsolar: this.Idiff[HOY][3]. HOY ∈ [0, 8759].</para>
        /// </summary>
        /// <param name="DOY">Day of the year, ∈ [1, 365].</param>
        /// <param name="HOY">Hour of the year, ∈ [0, 8759].</param>
        private void CalcIdiff_MT(int DOY, int HOY)
        {
            Parallel.For(0, this.SPCount, i =>
            {
                if (snowcovered[i][HOY])
                    this.Idiff[i][HOY] = 0.0;
                else
                {
                    this.Idiff[i][HOY] = Irradiation.Diffuse(
                        this.weather.DHI[HOY], this.weather.DNI[HOY], this.sunvectors[HOY].udtCoordinates.dZenithAngle,
                        this.sunvectors[HOY].udtCoordinates.dAzimuth, this.beta[i], this.psi[i], DOY,
                        sky[i].ShdwHorizon, sky[i].ShdwDome, sky[i].ShdwBeam[HOY]);
                }
            });
        }

        /// <summary>
        /// Calculates hourly diffuse radiation on the sensor point for the entire year.
        /// <para>Access: total: this.Idiff[HOY][0]; horizon: this.Idiff[HOY][1]; sky: this.Idiff[HOY][2]; circumsolar: this.Idiff[HOY][3]. HOY ∈ [0, 8759].</para>
        /// </summary>
        private void CalcIdiff()
        {
            int HOY = 0;
            for (int i = 1; i < 366; i++)
            {
                for (int u = 0; u < 24; u++)
                {
                    CalcIdiff(i, HOY);
                    HOY++;
                }
            }
        }

        /// <summary>
        /// Calculates hourly diffuse radiation on the sensor point for the entire year. Multi-Threading version.
        /// <para>Access: total: this.Idiff[HOY][0]; horizon: this.Idiff[HOY][1]; sky: this.Idiff[HOY][2]; circumsolar: this.Idiff[HOY][3]. HOY ∈ [0, 8759].</para>
        /// </summary>
        private void CalcIdiff_MT()
        {
            Parallel.For(1, 366, i =>
            {
                for (int u = 0; u < 24; u++)
                {
                    int HOY = (i - 1) * 24 + u;
                    CalcIdiff(i, HOY);
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
        private void CalcIbeam(int DOY, int HOY, int LT)
        {
            if (this.sunvectors[HOY].Sunshine == true)
            {
                for (int i = 0; i < this.SPCount; i++)
                {
                    if (Convert.ToInt32(sky[i].ShdwBeam[HOY]) + Convert.ToInt32(snowcovered[i][HOY]) > 0)
                        this.Ibeam[i][HOY] = 0.0;
                    else
                    {
                        this.Ibeam[i][HOY] = Irradiation.Beam(
                            this.weather.DNI[HOY], this.sunvectors[HOY].udtCoordinates.dZenithAngle,
                            this.sunvectors[HOY].udtCoordinates.dAzimuth, this.beta[i], this.psi[i]);
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
        private void CalcIbeam_MT(int DOY, int HOY, int LT)
        {
            if (this.sunvectors[HOY].Sunshine == true)
            {
                Parallel.For(0, this.SPCount, i =>
                {
                    if (Convert.ToInt32(sky[i].ShdwBeam[HOY]) + Convert.ToInt32(snowcovered[i][HOY]) > 0)
                        this.Ibeam[i][HOY] = 0.0;
                    else
                    {
                        this.Ibeam[i][HOY] = Irradiation.Beam(
                                this.weather.DNI[HOY], this.sunvectors[HOY].udtCoordinates.dZenithAngle,
                                this.sunvectors[HOY].udtCoordinates.dAzimuth, this.beta[i], this.psi[i]);
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
        private void CalcIbeam()
        {
            int HOY = 0;
            for (int i = 1; i < 366; i++)
            {
                for (int u = 0; u < 24; u++)
                {
                    CalcIbeam(i, HOY, u);
                    HOY++;
                }
            }
        }

        /// <summary>
        /// Calculates hourly beam (direct) irradiation on the sensor point for the entire year. Multi-Threading version.
        /// <para>Access: this.Ibeam[HOY]; HOY ∈ [0, 8759].</para>
        /// </summary>
        private void CalcIbeam_MT()
        {
            Parallel.For(1, 366, i =>
            {
                for (int u = 0; u < 24; u++)
                {
                    int HOY = (i - 1) * 24 + u;
                    CalcIbeam(i, HOY, u);
                }
            });
        }



        /// <summary>
        /// Calculates total solar irradiation on a sensor point for one hour of the year.
        /// <para>Access: this.I[HOY]; HOY ∈ [0, 8759]. HOY = (DOY-1) * 24 + LT.</para>
        /// </summary>
        /// <param name="DOY">Day of the year, ∈ [1, 365].</param>
        /// <param name="LT">Local time, i.e. hour of the day, ∈ [0, 23].</param>
        public void CalcIrradiation(int DOY, int LT)
        {
            int HOY = (DOY - 1) * 24 + LT;
            CalcIbeam(DOY, HOY, LT);
            CalcIdiff(DOY, HOY);
            for (int i = 0; i < this.I.Length; i++)
                this.I[i][HOY] = this.Ibeam[i][HOY] + this.Idiff[i][HOY] + this.Irefl_spec[i][HOY] + this.Irefl_diff[i][HOY];
        }

        /// <summary>
        /// Calculates total solar irradiation on a sensor point for one hour of the year. Multi-Threading version.
        /// <para>Access: this.I[HOY]; HOY ∈ [0, 8759]. HOY = (DOY-1) * 24 + LT.</para>
        /// </summary>
        /// <param name="DOY">Day of the year, ∈ [1, 365].</param>
        /// <param name="LT">Local time, i.e. hour of the day, ∈ [0, 23].</param>
        public void CalcIrradiationMT(int DOY, int LT)
        {
            int HOY = (DOY - 1) * 24 + LT;

            CalcIbeam_MT(DOY, HOY, LT);
            CalcIdiff_MT(DOY, HOY);

            Parallel.For(0, this.I.Length, i =>
            {
                this.I[i][HOY] = this.Ibeam[i][HOY] + this.Idiff[i][HOY] + this.Irefl_spec[i][HOY] + this.Irefl_diff[i][HOY];
            });
        }

        /// <summary>
        /// Calculates hourly total solar irradiation on the sensor point for the entire year.
        /// <para>Access: this.I[HOY]; HOY ∈ [0, 8759].</para>
        /// </summary>
        public void CalcIrradiation()
        {
            CalcIbeam();
            CalcIdiff();
            for (int i = 0; i < this.I.Length; i++)
                for (int t = 0; t < 8760; t++)
                    this.I[i][t] = this.Ibeam[i][t] + this.Idiff[i][t] + this.Irefl_spec[i][t] + this.Irefl_diff[i][t];
        }

        /// <summary>
        /// Calculates hourly total solar irradiation on the sensor point for the entire year. Multi-Threading version.
        /// <para>Access: this.I[HOY]; HOY ∈ [0, 8759].</para>
        /// </summary>
        public void CalcIrradiationMT()
        {
            CalcIbeam_MT();
            CalcIdiff_MT();
            Parallel.For(0, this.I.Length, i =>
            {
                for (int t = 0; t < 8760; t++)
                    this.I[i][t] = this.Ibeam[i][t] + this.Idiff[i][t] + this.Irefl_spec[i][t] + this.Irefl_diff[i][t]; 
            });
        }





        /// <summary>
        /// Applys shadow / obstruction factors from externally calculated view factor calculations to the sensor points.
        /// </summary>
        /// <param name="ShdwBeam_hour">Indicate for one hour of the year, if a sensor point is obstructed from beam radiation (true), or not (false). The list must have booleans for each sensor point.</param>
        /// <param name="ShdwSky">Indicate for each vertex of a sensor point's skydome, if the view between sensor point and vertex is obstructed (true), or not (false). The list has boolean arrays for each sensor point; an array is of length of the sky dome vertex count (this.sky[i].VerticesHemisphere.Count).</param>
        /// <param name="HOY"></param>
        public void SetShadows(List<bool> ShdwBeam_hour, List<bool[]> ShdwSky, int HOY)
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
                    this.sky[i].VertexShadowSphere[this.sky[i].VerticesHemisphere[u]] = ShdwSky[i][u];
                }
                this.sky[i].SetShadow_Dome();
                this.sky[i].SetShadow_Horizon();
                for (int t = 0; t < 8760; t++)
                {
                    this.sky[i].SetShadow_Beam(t, ShdwBeam_hour[i][t]);
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

            for (int i = 0; i < sky.Length; i++)
            {
                for (int u = 0; u < ShdwSky[i].Length; u++)
                {
                    this.sky[i].VertexShadowSphere[this.sky[i].VerticesHemisphere[u]] = ShdwSky[i][u];
                }

                this.sky[i].SetShadow_Dome();
                this.sky[i].SetShadow_Horizon();
                //this.sky[i].SetShadow_Beam(HOY, ShdwBeam_hour[i]);
            }


            //shdwbeam_summer, winter ewquinox have different lengths!!! coz only when sunshine true. 
            //but its always 
            //march 20
            //june 21
            //december 21

            //interpolate!   check github solar_v1


            //            'interpolating all 8760 between correct dates equinox, summer, winter
            //Private Function InterpolateShadowDays(radList As List(Of Double), _
            //                                       _shadowsummer As ShadowFactor, _
            //                                       _shadowwinter As ShadowFactor, _
            //                                       _shadowequinox As ShadowFactor) As List(Of Double)

            //    '100% equinox is 20 march and 23 Sept
            //    '100% summer is 21 june
            //    '100% winter is 22 decemebr

            //    Dim y1, y2, y3, y4, y5, y6 As Integer
            //    y1 = -9     'winter solstice
            //    y2 = 78     'equinox spring
            //    y3 = 171    'summer solstice
            //    y4 = 265    'equinox autumn
            //    y5 = 355    'winter solstice
            //    y6 = 443    'equinox spring
            int y1, y2, y3, y4, y5, y6;
            y1 = -9;    //winter solstice
            y2 = 78;    //equinox spring
            y3 = 171;   // summer solstice
            y4 = 265;   //equinox solstice
            y5 = 355;   //winter solstice
            y6 = 443;   //equinox spring

            //    Dim fullF1, fullF2 As Integer
            int fullF1, fullF2;

            //    Dim factor As Double = 1
            //    Dim i, u As Integer
            //    Dim InterpolInterv As Integer
            //    Dim dist1, dist2 As Double
            double factor = 1.0;
            int InterpolInterv;
            double dist1, dist2;


            for (int i = 0; i < this.sky.Length; i++)
            {
                //    fullF1 = y1                    '100% shadowwinter on this day
                //    fullF2 = y2                     '100% shadowequinox on this day
                //    InterpolInterv = y2 + (y1 * -1)
                fullF1 = y1;    //100% shadowwinter on this day
                fullF2 = y2;    //100% shadowequinox on this day
                InterpolInterv = y2 + y1 * -1;
                //    For i = 0 To y2 - 1             'from 0.Jan to 19.March
                for (int d = 0; d < y2; d++)
                {
                    //        dist1 = (InterpolInterv - Math.Abs(fullF1 - i)) / InterpolInterv
                    //        dist1 = 1 - dist1
                    //        dist1 = Math.Cos(dist1 * (0.5 * Math.PI))
                    dist1 = (InterpolInterv - Math.Abs(fullF1 - d)) / InterpolInterv;
                    dist1 = 1 - dist1;
                    dist1 = Math.Cos(dist1 * (0.5 * Math.PI));
                    //        'dist2 = (InterpolInterv - Math.Abs(fullF2 - i)) / InterpolInterv
                    //        dist2 = 1 - dist1
                    dist2 = (InterpolInterv - Math.Abs(fullF2 - d)) / InterpolInterv;
                    dist2 = 1 - dist1;
                    //        For u = 0 To 23
                    //            factor = ((1 - Math.Round(_shadowwinter.ShadowFactors(u), 4)) * dist1) * _shadowwinter.sunshine(u) + _
                    //                ((1 - Math.Round(_shadowequinox.ShadowFactors(u), 4)) * dist2) * _shadowequinox.sunshine(u)

                    //            radList(i * 24 + u) = radList(i * 24 + u) * factor
                    //            If factor > 1 Then
                    //                factor = 1
                    //            End If
                    //        Next
                    for (int u = 0; u < 24; u++)
                    {
                        factor = ((1 - Convert.ToDouble(ShdwBeam_Winter[i][u])) * dist1) +
                            ((1 - Convert.ToDouble(ShdwBeam_Equinox[i][u])) * dist2);
                        bool shdw = (factor >= 0.5) ? false : true;
                        int HOY = d * 24 + u;
                        this.sky[i].SetShadow_Beam(HOY, shdw);
                    }
                    //    Next
                }


                //    fullF1 = y2                      '100% shadowequinox on this day
                //    fullF2 = y3                     '100% shadowsummer on this day
                //    InterpolInterv = y3 - y2
                fullF1 = y2;        // 100% shadowequinox on this day
                fullF2 = y3;        // 100% shadowsummer on this day
                InterpolInterv = y3 - y2;
                //    For i = y2 To y3 - 1            'from 20.March to 20.june
                for (int d = y2; d < y3; d++)
                {
                    //        dist1 = (InterpolInterv - Math.Abs(fullF1 - i)) / InterpolInterv
                    //        'dist1 = 1 - dist1
                    //        dist1 = Math.Sin(dist1 * (0.5 * Math.PI))
                    //        'dist2 = (InterpolInterv - Math.Abs(fullF2 - i)) / InterpolInterv
                    //        dist2 = 1 - dist1
                    dist1 = (InterpolInterv - Math.Abs(fullF1 - d)) / InterpolInterv;
                    dist1 = 1 - dist1;
                    dist1 = Math.Sin(dist1 * (0.5 * Math.PI));
                    dist2 = (InterpolInterv - Math.Abs(fullF2 - d)) / InterpolInterv;
                    dist2 = 1 - dist1;
                    //        For u = 0 To 23
                    //            factor = ((1 - Math.Round(_shadowequinox.ShadowFactors(u), 4)) * dist1) * _shadowequinox.sunshine(u) + _
                    //               ((1 - Math.Round(_shadowsummer.ShadowFactors(u), 4)) * dist2) * _shadowsummer.sunshine(u)
                    //            radList(i * 24 + u) = radList(i * 24 + u) * factor
                    //            If factor > 1 Then
                    //                factor = 1
                    //            End If
                    //        Next
                    for (int u = 0; u < 24; u++)
                    {
                        factor = ((1 - Convert.ToDouble(ShdwBeam_Equinox[i][u])) * dist1) +
                            ((1 - Convert.ToDouble(ShdwBeam_Summer[i][u])) * dist2);
                        bool shdw = (factor >= 0.5) ? false : true;
                        int HOY = d * 24 + u;
                        this.sky[i].SetShadow_Beam(HOY, shdw);
                    }
                    //    Next
                }

                //    fullF1 = y3                      '100% shadowsummer on this day
                //    fullF2 = y4                     '100% shadowequinox on this day
                //    InterpolInterv = y4 - y3
                fullF1 = y3;        //100% shadowsummer on this day
                fullF2 = y4;        //100% shadowequinox on this day
                InterpolInterv = y4 - y3;
                //    For i = y3 To y4 - 1            'from 21.June to 22.Sept
                for (int d = y3; d < y4; d++)
                {
                    //        dist1 = (InterpolInterv - Math.Abs(fullF1 - i)) / InterpolInterv
                    //        dist1 = 1 - dist1
                    //        dist1 = Math.Cos(dist1 * (0.5 * Math.PI))
                    //        'dist2 = (InterpolInterv - Math.Abs(fullF2 - i)) / InterpolInterv
                    //        dist2 = 1 - dist1
                    dist1 = (InterpolInterv - Math.Abs(fullF1 - d)) / InterpolInterv;
                    dist1 = 1 - dist1;
                    dist1 = Math.Cos(dist1 * (0.5 * Math.PI));
                    dist2 = (InterpolInterv - Math.Abs(fullF2 - d)) / InterpolInterv;
                    dist2 = 1 - dist1;
                    //        For u = 0 To 23
                    //            factor = ((1 - Math.Round(_shadowsummer.ShadowFactors(u), 4)) * dist1) * _shadowsummer.sunshine(u) + _
                    //               ((1 - Math.Round(_shadowequinox.ShadowFactors(u), 4)) * dist2) * _shadowequinox.sunshine(u)
                    //            radList(i * 24 + u) = radList(i * 24 + u) * factor
                    //            If factor > 1 Then
                    //                factor = 1
                    //            End If
                    //        Next
                    for (int u = 0; u < 24; u++)
                    {
                        factor = ((1 - Convert.ToDouble(ShdwBeam_Summer[i][u])) * dist1) +
                            ((1 - Convert.ToDouble(ShdwBeam_Equinox[i][u])) * dist2);
                        bool shdw = (factor >= 0.5) ? false : true;
                        int HOY = d * 24 + u;
                        this.sky[i].SetShadow_Beam(HOY, shdw);
                    }
                    //    Next
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
                    //        dist1 = (InterpolInterv - Math.Abs(fullF1 - i)) / InterpolInterv
                    //        'dist1 = 1 - dist1
                    //        dist1 = Math.Sin(dist1 * (0.5 * Math.PI))
                    //        'dist2 = (InterpolInterv - Math.Abs(fullF2 - i)) / InterpolInterv
                    //        dist2 = 1 - dist1
                    dist1 = (InterpolInterv - Math.Abs(fullF1 - d)) / InterpolInterv;
                    dist1 = 1 - dist1;
                    dist1 = Math.Sin(dist1 * (0.5 * Math.PI));
                    dist2 = (InterpolInterv - Math.Abs(fullF2 - d)) / InterpolInterv;
                    dist2 = 1 - dist1;
                    //        For u = 0 To 23
                    //            factor = ((1 - Math.Round(_shadowequinox.ShadowFactors(u), 4)) * dist1) * _shadowequinox.sunshine(u) + _
                    //               ((1 - Math.Round(_shadowwinter.ShadowFactors(u), 4)) * dist2) * _shadowwinter.sunshine(u)
                    //            radList(i * 24 + u) = radList(i * 24 + u) * factor
                    //            If factor > 1 Then
                    //                factor = 1
                    //            End If
                    //        Next
                    for (int u = 0; u < 24; u++)
                    {
                        factor = ((1 - Convert.ToDouble(ShdwBeam_Equinox[i][u])) * dist1) +
                            ((1 - Convert.ToDouble(ShdwBeam_Winter[i][u])) * dist2);
                        bool shdw = (factor >= 0.5) ? false : true;
                        int HOY = d * 24 + u;
                        this.sky[i].SetShadow_Beam(HOY, shdw);
                    }
                    //    Next
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
                    //        dist1 = (InterpolInterv - Math.Abs(fullF1 - i)) / InterpolInterv
                    //        dist1 = 1 - dist1
                    //        dist1 = Math.Cos(dist1 * (0.5 * Math.PI))
                    //        'dist2 = (InterpolInterv - Math.Abs(fullF2 - i)) / InterpolInterv
                    //        dist2 = 1 - dist1
                    dist1 = (InterpolInterv - Math.Abs(fullF1 - d)) / InterpolInterv;
                    dist1 = 1 - dist1;
                    dist1 = Math.Cos(dist1 * (0.5 * Math.PI));
                    dist2 = (InterpolInterv - Math.Abs(fullF2 - d)) / InterpolInterv;
                    dist2 = 1 - dist1;
                    //        For u = 0 To 23
                    //            factor = ((1 - Math.Round(_shadowwinter.ShadowFactors(u), 4)) * dist1) * _shadowwinter.sunshine(u) + _
                    //               ((1 - Math.Round(_shadowequinox.ShadowFactors(u), 4)) * dist2) * _shadowequinox.sunshine(u)
                    //            radList(i * 24 + u) = radList(i * 24 + u) * factor
                    //            If factor > 1 Then
                    //                factor = 1
                    //            End If
                    //        Next
                    for (int u = 0; u < 24; u++)
                    {
                        factor = ((1 - Convert.ToDouble(ShdwBeam_Winter[i][u])) * dist1) +
                            ((1 - Convert.ToDouble(ShdwBeam_Equinox[i][u])) * dist2);
                        bool shdw = (factor >= 0.5) ? false : true;
                        int HOY = d * 24 + u;
                        this.sky[i].SetShadow_Beam(HOY, shdw);
                    }
                    //    Next
                }

            }
            //    Return radList
            //End Function

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
                    this.sky[i].VertexShadowSphere[this.sky[i].VerticesHemisphere[u]] = ShdwSky[i][u];
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
                            this.sky[i].SetShadow_Beam(HOY, shdw);
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
                    this.sky[i].VertexShadowSphere[this.sky[i].VerticesHemisphere[u]] = ShdwSky[i][u];
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
                    dist1 = (InterpolInterv - Math.Abs(fullF1 - d)) / InterpolInterv;
                    dist1 = 1 - dist1;
                    dist1 = Math.Cos(dist1 * (0.5 * Math.PI));
                    dist2 = (InterpolInterv - Math.Abs(fullF2 - d)) / InterpolInterv;
                    dist2 = 1 - dist1;
                    for (int u = 0; u < 24; u++)
                    {
                        factor = ((1 - Convert.ToDouble(ShdwBeam_Winter[i][u])) * dist1) +
                            ((1 - Convert.ToDouble(ShdwBeam_Equinox[i][u])) * dist2);
                        bool shdw = (factor >= 0.5) ? false : true;
                        int HOY = d * 24 + u;
                        this.sky[i].SetShadow_Beam(HOY, shdw);
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
                    dist1 = (InterpolInterv - Math.Abs(fullF1 - d)) / InterpolInterv;
                    dist1 = 1 - dist1;
                    dist1 = Math.Sin(dist1 * (0.5 * Math.PI));
                    dist2 = (InterpolInterv - Math.Abs(fullF2 - d)) / InterpolInterv;
                    dist2 = 1 - dist1;
                    for (int u = 0; u < 24; u++)
                    {
                        factor = ((1 - Convert.ToDouble(ShdwBeam_Equinox[i][u])) * dist1) +
                            ((1 - Convert.ToDouble(ShdwBeam_Summer[i][u])) * dist2);
                        bool shdw = (factor >= 0.5) ? false : true;
                        int HOY = d * 24 + u;
                        this.sky[i].SetShadow_Beam(HOY, shdw);
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
                    dist1 = (InterpolInterv - Math.Abs(fullF1 - d)) / InterpolInterv;
                    dist1 = 1 - dist1;
                    dist1 = Math.Cos(dist1 * (0.5 * Math.PI));
                    dist2 = (InterpolInterv - Math.Abs(fullF2 - d)) / InterpolInterv;
                    dist2 = 1 - dist1;
                    for (int u = 0; u < 24; u++)
                    {
                        factor = ((1 - Convert.ToDouble(ShdwBeam_Summer[i][u])) * dist1) +
                            ((1 - Convert.ToDouble(ShdwBeam_Equinox[i][u])) * dist2);
                        bool shdw = (factor >= 0.5) ? false : true;
                        int HOY = d * 24 + u;
                        this.sky[i].SetShadow_Beam(HOY, shdw);
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
                    dist1 = (InterpolInterv - Math.Abs(fullF1 - d)) / InterpolInterv;
                    dist1 = 1 - dist1;
                    dist1 = Math.Sin(dist1 * (0.5 * Math.PI));
                    dist2 = (InterpolInterv - Math.Abs(fullF2 - d)) / InterpolInterv;
                    dist2 = 1 - dist1;
                    for (int u = 0; u < 24; u++)
                    {
                        factor = ((1 - Convert.ToDouble(ShdwBeam_Equinox[i][u])) * dist1) +
                            ((1 - Convert.ToDouble(ShdwBeam_Winter[i][u])) * dist2);
                        bool shdw = (factor >= 0.5) ? false : true;
                        int HOY = d * 24 + u;
                        this.sky[i].SetShadow_Beam(HOY, shdw);
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
                    dist1 = (InterpolInterv - Math.Abs(fullF1 - d)) / InterpolInterv;
                    dist1 = 1 - dist1;
                    dist1 = Math.Cos(dist1 * (0.5 * Math.PI));
                    dist2 = (InterpolInterv - Math.Abs(fullF2 - d)) / InterpolInterv;
                    dist2 = 1 - dist1;
                    for (int u = 0; u < 24; u++)
                    {
                        factor = ((1 - Convert.ToDouble(ShdwBeam_Winter[i][u])) * dist1) +
                            ((1 - Convert.ToDouble(ShdwBeam_Equinox[i][u])) * dist2);
                        bool shdw = (factor >= 0.5) ? false : true;
                        int HOY = d * 24 + u;
                        this.sky[i].SetShadow_Beam(HOY, shdw);
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
                    this.sky[i].VertexShadowSphere[this.sky[i].VerticesHemisphere[u]] = ShdwSky[i][u];
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
                            this.sky[i].SetShadow_Beam(HOY, shdw);
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
        public void SetSnowcover(double snow_threshold, double tilt_treshold)
        {
            for (int i = 0; i < this.I.Length; i++)
            {
                for (int t = 0; t < 8760; t++)
                {
                    if (this.weather.Snow[t] > snow_threshold && beta[i] < tilt_treshold)
                    {
                        this.snowcovered[i][t] = true;
                    }
                }
            }

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
        /// Set irradiation by inter-reflections for all hours of the year, using interpolation for specular reflection. 
        /// Actual irradiation calculation needs to be done externally.
        /// </summary>
        /// <param name="_Ispecular">Irradiation values by specular reflection for several days and for each sensor point. List = each element is one day, double[SPCount][24]</param>
        /// <param name="_Idiffuse">Irradiation values by diffuse reflection for each sensor point.</param>
        public void SetInterreflInterpolated(List<double[][]> _Ispecular, double[] _Idiffuse)
        {

        }

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
