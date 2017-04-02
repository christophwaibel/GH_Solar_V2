using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SolarModel
{
    public class Sensorpoint
    {
        public SkyDome sky;
        public List<SunVector> sunvectors;
        public Context.cLocation location;
        public Context.cWeatherdata weather;
        public double beta;     // sensor point tilt angle
        public double psi;      // sensor point azimuth

        public double[] Idiff;
        public double[] Ibeam;
        public double[] I;

        public Sensorpoint(int year, Context.cWeatherdata _weather, Context.cLocation _location, List<SunVector> _sunvectors, double beta, double psi, int reclvlsky)
        {
            location = _location;
            weather = _weather;
            sunvectors = new List<SunVector>(_sunvectors);

            sky = new SkyDome(reclvlsky, year, location.dLatitude, location.dLongitude);

            this.beta = beta;
            this.psi = psi;

            this.Idiff = new double[8760];
            this.Ibeam = new double[8760];
            this.I = new double[8760];
        }




        /// <summary>
        /// Calculates diffuse irradiation on the sensor point for one hour of the year.
        /// </summary>
        /// <param name="DHI"></param>
        /// <param name="DNI"></param>
        /// <param name="DOY">1-365</param>
        /// <param name="HOY">0-8759</param>
        /// <param name="sunvector">Sun vector for the current location and time.</param>
        private void CalcIdiff_unobstr(int DOY, int HOY)
        {
            if(this.sunvectors[HOY].Sunshine == true)
                this.Idiff[HOY] = Irradiation.Diffuse(this.weather.DHI[HOY], this.weather.DNI[HOY], this.sunvectors[HOY].udtCoordinates.dZenithAngle, this.sunvectors[HOY].udtCoordinates.dAzimuth, this.beta, this.psi, DOY)[0];
            else
                this.Idiff[HOY] = 0;
 
        }

        /// <summary>
        /// Calculates hourly diffuse radiation on the sensor point for the entire year.
        /// Access: this.Idiff[DOY]
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
        /// </summary>
        /// <param name="DOY"></param>
        private void CalcIbeam_unobstr(int DOY, int LT, int HOY)
        {
            if (this.sunvectors[HOY].Sunshine == true)
                this.Ibeam[HOY] = Irradiation.Beam(this.weather.DNI[HOY], this.beta, this.location.dLatitude, this.location.dLongitude, this.psi, DOY, LT, this.location.dTgmt);
            else
                this.Ibeam[HOY] = 0;
        }

        /// <summary>
        /// Calculates hourly beam (direct) irradiation on the sensor point for the entire year.
        /// </summary>
        /// <returns></returns>
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
        /// Calculates solar irradiation on a sensor point for one hour of the year.
        /// </summary>
        /// <param name="DOY">1-365</param>
        /// <param name="LT">0-23</param>
        public void CalcIrradiation(int DOY, int LT)
        {
            int HOY = (DOY-1) * 24 + LT;
            CalcIbeam_unobstr(DOY, LT, HOY);
            CalcIdiff_unobstr(DOY, HOY);
            this.I[HOY] = this.Ibeam[HOY] + this.Idiff[HOY];

        }

        /// <summary>
        /// Calculates hourly irradiation on the sensor point for the entire year.
        /// </summary>
        /// <returns></returns>
        public void CalcIrradiation()
        {
            CalcIbeam_unobstr();
            CalcIdiff_unobstr();
            for (int i = 0; i < 8760; i++)
            {
                this.I[i] = this.Ibeam[i] + this.Idiff[i];
            }

        }

    }
}
