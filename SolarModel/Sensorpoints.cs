using System;
using System.Collections.Generic;
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
        public SkyDome [] sky { get; private set; }
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
        public double [] beta;
        /// <summary>
        /// Sensor point azimuth angle in degree.
        /// </summary>
        public double [] psi;
        /// <summary>
        /// Diffuse irradiation. 
        /// <para>2D array, [8760][4].</para>
        /// <para>[t][0]: total diffuse irradiation; [t][1]: horizon component; [t][2]: anisotropic sky component; [t][3]: circumsolar component.</para>
        /// </summary>
        public double[][][] Idiff;
        /// <summary>
        /// Beam (direct) irradiation.
        /// <para>Array with 8760 doubles, for each hour of the year.</para>
        /// </summary>
        public double[][] Ibeam;
        /// <summary>
        /// Total irradiation (beam + diffuse).
        /// <para>Array with 8760 doubles, for each hour of the year.</para>
        /// </summary>
        public double[][] I;

        public bool[][] snowcovered { get; private set; }



        /// <summary>
        /// Sensor point object, to calculate solar irradiation.
        /// </summary>
        /// <param name="year">Year to calculate.</param>
        /// <param name="_weather">Weather data. See structure "Context.cWeatherdata".</param>
        /// <param name="_location">Location data. See structure "Context.cLocation".</param>
        /// <param name="_sunvectors">8760 sun vectors, for each hour of the year. See class "SunVector".</param>
        /// <param name="beta">Tilt angle of sensor point.</param>
        /// <param name="psi">Azimuth angle of sensor point.</param>
        /// <param name="reclvlsky">Recursion level for sky hemisphere. 1 or 2 recommended. See class "SkyDome".</param>
        public Sensorpoints(int year, Context.cWeatherdata _weather, Context.cLocation _location, List<SunVector> _sunvectors, double []beta, double []psi, int reclvlsky)
        {
            this.location = _location;
            this.weather = _weather;
            this.sunvectors = new List<SunVector>(_sunvectors);

            this.beta = new double[beta.Length];
            this.psi = new double[psi.Length];
            this.sky = new SkyDome[beta.Length];
            Array.Copy(beta, this.beta, beta.Length);
            Array.Copy(psi, this.psi, psi.Length);
            for (int i = 0; i < beta.Length; i++)
            {
                if (i == 0) this.sky[i] = new SkyDome(reclvlsky);
                else this.sky[i] = new SkyDome(this.sky[0]);
            }


            this.Idiff = new double[beta.Length][][];
            this.Ibeam = new double[beta.Length][];
            this.I = new double[beta.Length][];

            for (int i = 0; i < beta.Length; i++)
            {
                this.Idiff[i] = new double[8760][];
                for (int u = 0; u < Idiff.Length; u++)
                {
                    this.Idiff[i][u] = new double[4];
                }
                this.Ibeam[i] = new double[8760];
                this.I[i] = new double[8760];
            }

            this.snowcovered = new bool[I.Length][];
            for (int i = 0; i < snowcovered.Length; i++)
            {
                snowcovered[i] = new bool[8760];
            }
        }




        /// <summary>
        /// Calculates diffuse irradiation on the sensor point for one hour of the year.
        /// <para>Access: total: this.Idiff[HOY][0]; horizon: this.Idiff[HOY][1]; sky: this.Idiff[HOY][2]; circumsolar: this.Idiff[HOY][3]. HOY ∈ [0, 8759].</para>
        /// </summary>
        /// <param name="DOY">Day of the year, ∈ [1, 365].</param>
        /// <param name="HOY">Hour of the year, ∈ [0, 8759].</param>
        private void CalcIdiff(int DOY, int HOY)
        {
            if (this.sunvectors[HOY].Sunshine == true)
                for (int i = 0; i < this.Idiff.Length; i++)
                {
                    if (snowcovered[i][HOY])
                        this.Idiff[i][HOY] = new double[4] { 0, 0, 0, 0 };
                    else
                    {
                        this.Idiff[i][HOY] = Irradiation.Diffuse(
                            this.weather.DHI[HOY], this.weather.DNI[HOY], this.sunvectors[HOY].udtCoordinates.dZenithAngle,
                            this.sunvectors[HOY].udtCoordinates.dAzimuth, this.beta[i], this.psi[i], DOY);
                        this.Idiff[i][HOY][3] *= (1 - Convert.ToDouble(sky[i].ShdwBeam[HOY]));
                        this.Idiff[i][HOY][2] *= (1 - sky[i].ShdwDome);
                        this.Idiff[i][HOY][1] *= (1 - sky[i].ShdwHorizon);
                        this.Idiff[i][HOY][0] = this.Idiff[i][HOY][1] + this.Idiff[i][HOY][2] + this.Idiff[i][HOY][3];
                    }
                    //this.Idiff[i][HOY] = Irradiation.Diffuse(
                    //    this.weather.DHI[HOY], this.weather.DNI[HOY], this.sunvectors[HOY].udtCoordinates.dZenithAngle,
                    //    this.sunvectors[HOY].udtCoordinates.dAzimuth, this.beta[i], this.psi[i], DOY);
                    //this.Idiff[i][HOY][3] *= (1 - Convert.ToInt32(sky[i].ShdwBeam[HOY]));
                    //this.Idiff[i][HOY][2] *= (1 - sky[i].ShdwDome);
                    //this.Idiff[i][HOY][1] *= (1 - sky[i].ShdwHorizon);
                    //this.Idiff[i][HOY][0] = this.Idiff[i][HOY][1] + this.Idiff[i][HOY][2] + this.Idiff[i][HOY][3];
                }
            else
                for (int i = 0; i < this.Idiff.Length; i++)
                    this.Idiff[i][HOY] = new double[4] { 0, 0, 0, 0 };
        }

        /// <summary>
        /// Calculates diffuse irradiation on the sensor point for one hour of the year. Multi-Threading version.
        /// <para>Access: total: this.Idiff[HOY][0]; horizon: this.Idiff[HOY][1]; sky: this.Idiff[HOY][2]; circumsolar: this.Idiff[HOY][3]. HOY ∈ [0, 8759].</para>
        /// </summary>
        /// <param name="DOY">Day of the year, ∈ [1, 365].</param>
        /// <param name="HOY">Hour of the year, ∈ [0, 8759].</param>
        private void CalcIdiff_MT(int DOY, int HOY)
        {
            if (this.sunvectors[HOY].Sunshine == true)
            {
                Parallel.For(0, this.Idiff.Length, i =>
                {
                    if (snowcovered[i][HOY])
                        this.Idiff[i][HOY] = new double[4] { 0, 0, 0, 0 };
                    else
                    {
                        this.Idiff[i][HOY] = Irradiation.Diffuse(
                            this.weather.DHI[HOY], this.weather.DNI[HOY], this.sunvectors[HOY].udtCoordinates.dZenithAngle,
                            this.sunvectors[HOY].udtCoordinates.dAzimuth, this.beta[i], this.psi[i], DOY);
                        this.Idiff[i][HOY][3] *= (1 - Convert.ToDouble(sky[i].ShdwBeam[HOY]));
                        this.Idiff[i][HOY][2] *= (1 - sky[i].ShdwDome);
                        this.Idiff[i][HOY][1] *= (1 - sky[i].ShdwHorizon);
                        this.Idiff[i][HOY][0] = this.Idiff[i][HOY][1] + this.Idiff[i][HOY][2] + this.Idiff[i][HOY][3];
                    }
                    //this.Idiff[i][HOY] = Irradiation.Diffuse(
                    //    this.weather.DHI[HOY], this.weather.DNI[HOY], this.sunvectors[HOY].udtCoordinates.dZenithAngle, 
                    //    this.sunvectors[HOY].udtCoordinates.dAzimuth, this.beta[i], this.psi[i], DOY);
                    //this.Idiff[i][HOY][3] *= (1 - Convert.ToInt32(sky[i].ShdwBeam[HOY]));
                    //this.Idiff[i][HOY][2] *= (1 - sky[i].ShdwDome);
                    //this.Idiff[i][HOY][1] *= (1 - sky[i].ShdwHorizon);
                    //this.Idiff[i][HOY][0] = this.Idiff[i][HOY][1] + this.Idiff[i][HOY][2] + this.Idiff[i][HOY][3];
                });
            }
            else
            {
                Parallel.For(0, this.Idiff.Length, i =>
                {
                    this.Idiff[i][HOY] = new double[4] { 0, 0, 0, 0 };
                });
            }
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
                    HOY++;
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
                for (int i = 0; i < this.Ibeam.Length; i++)
                {
                    if (Convert.ToInt32(sky[i].ShdwBeam[HOY]) + Convert.ToInt32(snowcovered[i][HOY]) > 0)
                        this.Ibeam[i][HOY] = 0.0;
                    else
                        this.Ibeam[i][HOY] = Irradiation.Beam(
                            this.weather.DNI[HOY], this.beta[i], this.location.dLatitude,
                            this.location.dLongitude, this.psi[i], DOY, LT, this.location.dTgmt);
                    //this.Ibeam[i][HOY] = Irradiation.Beam(
                    //    this.weather.DNI[HOY], this.beta[i], this.location.dLatitude,
                    //    this.location.dLongitude, this.psi[i], DOY, LT, this.location.dTgmt) *
                    //    (1 - Convert.ToDouble(sky[i].ShdwBeam[HOY]));
                }
            else
                for (int i = 0; i < this.Ibeam.Length; i++)
                    this.Ibeam[i][HOY] = 0;
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
                Parallel.For(0, this.Ibeam.Length, i =>
                {
                    if (Convert.ToInt32(sky[i].ShdwBeam[HOY]) + Convert.ToInt32(snowcovered[i][HOY]) > 0)
                        this.Ibeam[i][HOY] = 0.0;
                    else
                        this.Ibeam[i][HOY] = Irradiation.Beam(
                            this.weather.DNI[HOY], this.beta[i], this.location.dLatitude,
                            this.location.dLongitude, this.psi[i], DOY, LT, this.location.dTgmt);
                    //this.Ibeam[i][HOY] = Irradiation.Beam(
                    //    this.weather.DNI[HOY], this.beta[i], this.location.dLatitude, 
                    //    this.location.dLongitude, this.psi[i], DOY, LT, this.location.dTgmt)*
                    //    (1 - Convert.ToDouble(sky[i].ShdwBeam[HOY]))*
                    //    (1- Convert.ToDouble(this.snowcovered[i][HOY]));
                });
            }
            else
            {
                Parallel.For(0, this.Ibeam.Length, i =>
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
                this.I[i][HOY] = this.Ibeam[i][HOY] + this.Idiff[i][HOY][0];
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
                this.I[i][HOY] = this.Ibeam[i][HOY] + this.Idiff[i][HOY][0];
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
                for (int u = 0; u < 8760; u++)
                    this.I[i][u] = this.Ibeam[i][u] + this.Idiff[i][u][0];
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
                for (int u = 0; u < 8760; u++)
                    this.I[i][u] = this.Ibeam[i][u] + this.Idiff[i][u][0];
            });
        }






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
        /// Goes through each hour and applys snow blockage, if surface angle is flat enough and if the weather data indicates snow on that hour.
        /// </summary>
        /// <param name="snow_threshold">Snow threshold (mm), after which no radiation is assumed to reach the sensor point.</param>
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




        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        /// <summary>
        /// 
        /// </summary>
        public void SetInterreflection()
        {

        }



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
