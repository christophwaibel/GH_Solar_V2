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
        }




        /// <summary>
        /// Calculates diffuse irradiation on the sensor point for one hour of the year.
        /// <para>Access: total: this.Idiff[HOY][0]; horizon: this.Idiff[HOY][1]; sky: this.Idiff[HOY][2]; circumsolar: this.Idiff[HOY][3]. HOY ∈ [0, 8759].</para>
        /// </summary>
        /// <param name="DOY">Day of the year, ∈ [1, 365].</param>
        /// <param name="HOY">Hour of the year, ∈ [0, 8759].</param>
        private void CalcIdiff_unobstr(int DOY, int HOY)
        {
            if (this.sunvectors[HOY].Sunshine == true)
                for (int i = 0; i < this.Idiff.Length; i++)
                    this.Idiff[i][HOY] = Irradiation.Diffuse(this.weather.DHI[HOY], this.weather.DNI[HOY], this.sunvectors[HOY].udtCoordinates.dZenithAngle, this.sunvectors[HOY].udtCoordinates.dAzimuth, this.beta[i], this.psi[i], DOY);
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
        private void CalcIdiff_unobstrMT(int DOY, int HOY)
        {
            if (this.sunvectors[HOY].Sunshine == true)
            {
                Parallel.For(0, this.Idiff.Length, i =>
                {
                    this.Idiff[i][HOY] = Irradiation.Diffuse(this.weather.DHI[HOY], this.weather.DNI[HOY], this.sunvectors[HOY].udtCoordinates.dZenithAngle, this.sunvectors[HOY].udtCoordinates.dAzimuth, this.beta[i], this.psi[i], DOY);
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
        private void CalcIdiff_unobstr()
        {
            int HOY = 0;
            for (int i = 1; i < 366; i++)
            {
                for (int u = 0; u < 24; u++)
                {
                    CalcIdiff_unobstr(i, HOY);
                    HOY++;
                }
            }
        }

        /// <summary>
        /// Calculates hourly diffuse radiation on the sensor point for the entire year. Multi-Threading version.
        /// <para>Access: total: this.Idiff[HOY][0]; horizon: this.Idiff[HOY][1]; sky: this.Idiff[HOY][2]; circumsolar: this.Idiff[HOY][3]. HOY ∈ [0, 8759].</para>
        /// </summary>
        private void CalcIdiff_unobstrMT()
        {
            Parallel.For(1, 366, i =>
            {
                for (int u = 0; u < 24; u++)
                {
                    int HOY = (i - 1) * 24 + u;
                    CalcIdiff_unobstr(i, HOY);
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
        private void CalcIbeam_unobstr(int DOY, int HOY, int LT)
        {
            if (this.sunvectors[HOY].Sunshine == true)
                for(int i=0; i<this.Ibeam.Length; i++)
                    this.Ibeam[i][HOY] = Irradiation.Beam(this.weather.DNI[HOY], this.beta[i], this.location.dLatitude, this.location.dLongitude, this.psi[i], DOY, LT, this.location.dTgmt);
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
        private void CalcIbeam_unobstrMT(int DOY, int HOY, int LT)
        {
            if (this.sunvectors[HOY].Sunshine == true)
            {
                Parallel.For(0, this.Ibeam.Length, i =>
                {
                    this.Ibeam[i][HOY] = Irradiation.Beam(this.weather.DNI[HOY], this.beta[i], this.location.dLatitude, this.location.dLongitude, this.psi[i], DOY, LT, this.location.dTgmt);
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
        private void CalcIbeam_unobstr()
        {
            int HOY = 0;
            for (int i = 1; i < 366; i++)
            {
                for (int u = 0; u < 24; u++)
                {
                    CalcIbeam_unobstr(i, HOY, u);
                    HOY++;
                }
            }

        }

        /// <summary>
        /// Calculates hourly beam (direct) irradiation on the sensor point for the entire year. Multi-Threading version.
        /// <para>Access: this.Ibeam[HOY]; HOY ∈ [0, 8759].</para>
        /// </summary>
        private void CalcIbeam_unobstrMT()
        {
            Parallel.For(1, 366, i =>
            {
                for (int u = 0; u < 24; u++)
                {
                    int HOY = (i - 1) * 24 + u;
                    CalcIbeam_unobstr(i, HOY, u);
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
            CalcIbeam_unobstr(DOY, LT, HOY);
            CalcIdiff_unobstr(DOY, HOY);
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

            CalcIbeam_unobstrMT(DOY, LT, HOY);
            CalcIdiff_unobstrMT(DOY, HOY);

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
            CalcIbeam_unobstr();
            CalcIdiff_unobstr();
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
            CalcIbeam_unobstrMT();
            CalcIdiff_unobstrMT();
            Parallel.For(0, this.I.Length, i =>
            {
                for (int u = 0; u < 8760; u++)
                    this.I[i][u] = this.Ibeam[i][u] + this.Idiff[i][u][0];
            });
        }




        /// <summary>
        /// Sets the obstruction factors (i.e. shadows) of each sensor point.
        /// <para>Split into shadows of diffuse radiation (further split into horizon and dome) and beam radiation.</para>
        /// </summary>
        public void SetShadows()
        {

            //!!!!!!!!!! quite trivial. rhino will provide list of vectors for obstruction. INTERPOLATION!!! take 3 days! and itnerpolate in-between.
            // I wonder if a meta-model would work, which takes 3 calculated days as features (and longi+latitude) and tells me for the rest days, if the rays are also blocked. I cuold try that quite quickly, coz I can calculate the training data in rhino. Just run 8760 obstruction checks. Then, meta model in matlab, python, Accord.net or whatever
            //sky[i].SetShadow_SunVector(

            //!!!!!!!!!!!!!!!!!!!!   these factors need to be calculated. rhino will give me bools of faces of the dome or horizon obstructed.
            //some calc needs to be further done to account for the regions, that are blocked anyway. these should be ommited for the factor calcualtion! only consider those faces, which would also be hit by sun rays 
            //sky[i].SetShadow_Dome(
            //sky[i].SetShadow_Horizon(

            for (int i = 0; i < this.I.Length; i++)
            {
                for (int t = 0; t < 8760; t++)
                {
                    this.Ibeam[i][t] *= (1 - Convert.ToInt32(sky[i].ShdwSunVector[t]));
                    this.Idiff[i][t][3] *= (1 - Convert.ToInt32(sky[i].ShdwSunVector[t]));
                    this.Idiff[i][t][2] *= (1 - sky[i].ShdwDome);
                    this.Idiff[i][t][1] *= (1 - sky[i].ShdwHorizon);
                    this.Idiff[i][t][0] = this.Idiff[i][t][1] + this.Idiff[i][t][2] + this.Idiff[i][t][3];
                    this.I[i][t] = this.Ibeam[i][t] + this.Idiff[i][t][0];
                }
            }
        }

        //!!!!!!!!!!! SetShadows(DOY, LT)


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
                        this.I[i][t] = 0;
                        this.Ibeam[i][t] = 0;
                        this.Idiff[i][t][0] = 0;
                        this.Idiff[i][t][1] = 0;
                        this.Idiff[i][t][2] = 0;
                        this.Idiff[i][t][3] = 0;
                    }
                }
            }

        }

        //!!!!!!!!!!! SetSnowcover(DOY, LT)


        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        /// <summary>
        /// 
        /// </summary>
        public void SetInterreflection()
        {

        }
    }
}
