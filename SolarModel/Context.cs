using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SolarModel
{
    public static class Context
    {

        public static void Create8760SunVectors(ref List<SunVector> sunvectors, double longitude, double latitude, int year)
        {
            for (int m = 1; m <= 12; m++)
            {
                int daysInMonth = System.DateTime.DaysInMonth(year, m);
                for (int d = 1; d <= daysInMonth; d++)
                {
                    for (int i = 1; i <= 24; i++)
                    {
                        SunVector sunvec = new SunVector(year, m, d, i, 0, 0, longitude, latitude);
                        sunvectors.Add(sunvec);
                    }
                }
            }
        }


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
            public List<double> DHI;            // W/m2. ∀ hours of the year
            public List<double> DNI;            // W/m2. ∀ hours of the year
            public List<double> Snow;           // mm. ∀ hours of the year
            //public List<double> Rain;           // mm. ∀ hours of the year
            //public List<double> Windspeed;      // m/s. ∀ hours of the year
            //public List<double> Winddirection;  // degree. ∀ hours of the year
        }

    }
}
