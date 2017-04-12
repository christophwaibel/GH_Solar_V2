using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

/*
 * Context.cs
 * Copyright 2017 Christoph Waibel <chwaibel@student.ethz.ch>
 * Modified after: http://blog.andreaskahler.com/search/label/c%23
 * 
 * This work is licensed under the GNU GPL license version 3.
*/

namespace SolarModel
{
    /// <summary>
    /// Class to define the geographic, climatological and time context.
    /// </summary>
    public static class Context
    {
        /// <summary>
        /// Time context.
        /// </summary>
        public struct cTime
        {
            public int iYear;
            public int iMonth;
            public int iDay;
            public double dHours;
            public double dMinutes;
            public double dSeconds;
        }
       
        /// <summary>
        /// Geographic context.
        /// </summary>
        public struct cLocation
        {
            /// <summary>
            /// Longitude in degree.
            /// </summary>
            public double dLongitude;
            /// <summary>
            /// Latitude in degree.
            /// </summary>
            public double dLatitude;
            /// <summary>
            /// Time difference to GMT (Greenwich Mean Time), in hours.
            /// </summary>
            public int dTgmt;
        }

        /// <summary>
        /// Climatological context, i.e. weather information.
        /// </summary>
        public struct cWeatherdata
        {
            /// <summary>
            /// Diffuse horizontal irradiation, in W/m2, ∀ hours of the year.
            /// </summary>
            public List<double> DHI;
            /// <summary>
            /// Direct normal irradiation, in W/m2, ∀ hours of the year.
            /// </summary>
            public List<double> DNI;
            /// <summary>
            /// Snow cover, in cm, ∀ hours of the year.
            /// </summary>
            public List<double> Snow;
            //public List<double> Rain;           // mm. ∀ hours of the year
            //public List<double> Windspeed;      // m/s. ∀ hours of the year
            //public List<double> Winddirection;  // degree. ∀ hours of the year
        }
    }
}
