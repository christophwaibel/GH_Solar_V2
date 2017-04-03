using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SolarModel
{
    public class Sensorpoint
    {
        /// <summary>
        /// Sky hemisphere. See class "SykDome".
        /// </summary>
        public SkyDome sky { get; private set; }
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
        public double beta;
        /// <summary>
        /// Sensor point azimuth angle in degree.
        /// </summary>
        public double psi;    
        /// <summary>
        /// Diffuse irradiation. 
        /// <para>2D array, [8760][4].</para>
        /// <para>[t][0]: total diffuse irradiation; [t][1]: horizon component; [t][2]: anisotropic sky component; [t][3]: circumsolar component.</para>
        /// </summary>
        public double[][] Idiff;
        /// <summary>
        /// Beam (direct) irradiation.
        /// <para>Array with 8760 doubles, for each hour of the year.</para>
        /// </summary>
        public double[] Ibeam;
        /// <summary>
        /// Total irradiation (beam + diffuse).
        /// <para>Array with 8760 doubles, for each hour of the year.</para>
        /// </summary>
        public double[] I;



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
        public Sensorpoint(int year, Context.cWeatherdata _weather, Context.cLocation _location, List<SunVector> _sunvectors, double beta, double psi, int reclvlsky)
        {
            location = _location;
            weather = _weather;
            sunvectors = new List<SunVector>(_sunvectors);

            sky = new SkyDome(reclvlsky, year, location.dLatitude, location.dLongitude);

            this.beta = beta;
            this.psi = psi;

            this.Idiff = new double[8760][];
            for (int i = 0; i < Idiff.Length; i++)
            {
                this.Idiff[i] = new double[4];
            }
            this.Ibeam = new double[8760];
            this.I = new double[8760];
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
                this.Idiff[HOY] = Irradiation.Diffuse(this.weather.DHI[HOY], this.weather.DNI[HOY], this.sunvectors[HOY].udtCoordinates.dZenithAngle, this.sunvectors[HOY].udtCoordinates.dAzimuth, this.beta, this.psi, DOY);
            else
                this.Idiff[HOY] = new double[4] { 0, 0, 0, 0 };
 
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
                for (int u = 0; u < 24;  u++)
                {
                    CalcIdiff_unobstr(i, HOY);
                    HOY++;
                }
            }
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
                this.Ibeam[HOY] = Irradiation.Beam(this.weather.DNI[HOY], this.beta, this.location.dLatitude, this.location.dLongitude, this.psi, DOY, LT, this.location.dTgmt);
            else
                this.Ibeam[HOY] = 0;
        }

        /// <summary>
        /// Calculates hourly beam (direct) irradiation on the sensor point for the entire year.
        /// <para>Access: this.Ibeam[HOY]; HOY ∈ [0, 8759].</para>
        /// </summary>
        private void CalcIbeam_unobstr()
        {
            //int lt = 0;  // 0-23
            //for (int i = 0; i < 8760; i++)
            //{
            //    CalcIbeam_unobstr(i, lt);
            //    lt++;
            //    if (lt % 24 == 0)
            //        lt = 0;
            //}
            int HOY = 0;
            for (int i = 1; i < 366; i++)
            {
                for (int u = 0; u < 24; u++)
                {
                    CalcIbeam_unobstr(i, u, HOY);
                    HOY++;
                }
            }

        }

        /// <summary>
        /// Calculates total solar irradiation on a sensor point for one hour of the year.
        /// <para>Access: this.I[HOY]; HOY ∈ [0, 8759]. HOY = (DOY-1) * 24 + LT.</para>
        /// </summary>
        /// <param name="DOY">Day of the year, ∈ [1, 365].</param>
        /// <param name="LT">Local time, i.e. hour of the day, ∈ [0, 23].</param>
        public void CalcIrradiation(int DOY, int LT)
        {
            int HOY = (DOY-1) * 24 + LT;
            CalcIbeam_unobstr(DOY, LT, HOY);
            CalcIdiff_unobstr(DOY, HOY);
            this.I[HOY] = this.Ibeam[HOY] + this.Idiff[HOY][0];
        }

        /// <summary>
        /// Calculates hourly total solar irradiation on the sensor point for the entire year.
        /// <para>Access: this.I[HOY]; HOY ∈ [0, 8759].</para>
        /// </summary>
        public void CalcIrradiation()
        {
            CalcIbeam_unobstr();
            CalcIdiff_unobstr();
            for (int i = 0; i < 8760; i++)
            {
                this.I[i] = this.Ibeam[i] + this.Idiff[i][0];
            }

        }

    }
}
