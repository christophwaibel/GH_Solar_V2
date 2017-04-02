using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SolarModel
{
    public static class Context
    {
        public struct cTime
        {
            public int iYear;
            public int iMonth;
            public int iDay;
            public double dHours;
            public double dMinutes;
            public double dSeconds;
        }
       
        public struct cLocation
        {
            public double dLongitude;
            public double dLatitude;
            public int dTgmt;
        }

        public struct cWeatherdata
        {
            public List<double> DHI;            // W/m2. ∀ hour of the year
            public List<double> DNI;            // W/m2. ∀ hour of the year
            //public List<double> Snow;           // mm. ∀ hour of the year
            //public List<double> Rain;           // mm. ∀ hour of the year
            //public List<double> Windspeed;      // m/s. ∀ hour of the year
            //public List<double> Winddirection;  // degree. ∀ hour of the year
        }

    }
}
