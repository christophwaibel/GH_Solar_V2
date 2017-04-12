using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

/*
 * SunVector.cs
 * Copyright 2017 Christoph Waibel <chwaibel@student.ethz.ch>
 * Source modified after: http://www.psa.es/sdg/sunpos.htm
 * 
 * This work is licensed under the GNU GPL license version 3.
*/

namespace SolarModel
{
    /// <summary>
    /// Calculates the solar vector.
    /// <para></para>
    /// </summary>
    /// <remarks>
    /// Source: Blanco-Muriel, M., Alarcon-Padilla, D.C., Lopez-Moratalla, T., Lara-Coira, M. (2001). Computing the solar vector. Solar Energy. vol. 70, issue 5, pp.431-441.
    /// <para>Source Code: http://www.psa.es/sdg/sunpos.htm </para>
    /// </remarks>
    public class SunVector
    {
        const double pi = Math.PI;
        const double twopi = 2 * pi;
        const double rad = pi / 180;
        const double dEarthMeanRadius = 6371.01;
        const double dAstronomicalUnit = 149597890.0;

        public Context.cTime udtTime;// = new Context.cTime();
        public Context.cLocation udtLocation;
        public cSunCoordinates udtCoordinates;
        public cSunXYZ udtCoordXYZ;

        public bool Sunshine;

        /// <summary>
        /// Create a Solar Vector for a specific time and location.
        /// </summary>
        /// <param name="_year"></param>
        /// <param name="_month"></param>
        /// <param name="_day"></param>
        /// <param name="_hours"></param>
        /// <param name="_minutes"></param>
        /// <param name="_seconds"></param>
        /// <param name="_longitude">Longitude in degree.</param>
        /// <param name="_latitude">Latitude in degree.</param>
        public SunVector(int _year, int _month, int _day, double _hours, double _minutes, double _seconds, double _longitude, double _latitude)
        {

            udtTime.iYear = _year;
            udtTime.iMonth = _month;
            udtTime.iDay = _day;
            udtTime.dHours = _hours;
            udtTime.dMinutes = _minutes;
            udtTime.dSeconds = _seconds;
            udtLocation.dLatitude = _latitude;
            udtLocation.dLongitude = _longitude;

            udtCoordinates = sunpos();
            udtCoordXYZ = SunposXYZ();

            if (udtCoordinates.dZenithAngle <= 90)
                Sunshine = true;
            else
                Sunshine = false;
        }


        /// <summary>
        /// Identify the equinox and solstice days of a year.
        /// <para>TO DO: Algorithm, calculating precise equinox and solstice days.</para>
        /// </summary>
        /// <param name="year"></param>
        /// <returns>[0]: spring equinox, [1]: summer solstice, [2]: autumn equinox, [3]: winter solstice. ∈ [1, 365]</returns>
        public static int [] GetEquinoxSolstice(int year)
        {
            //http://farside.ph.utexas.edu/Books/Syntaxis/Almagest/node36.html
            //view-source:https://stellafane.org/misc/equinox.html
            //https://ch.mathworks.com/matlabcentral/fileexchange/39356-a-matlab-script-for-predicting-equinoxes-and-solstices?requestedDomain=www.mathworks.com
            //

            int[] dm = new int[12];
            int[] equsol = new int[4];
            for (int i = 0; i < 12; i++)
            {
                dm[i] = System.DateTime.DaysInMonth(year, i + 1);
                if (i < 2)
                {
                    equsol[0] += dm[i];
                }
                if (i < 5)
                {
                    equsol[1] += dm[i];
                }
                if (i < 8)
                {
                    equsol[2] += dm[i];
                }
                if (i < 11)
                {
                    equsol[3] += dm[i];
                }
            }

            equsol[0] += 20;
            equsol[1] += 21;
            equsol[2] += 22;
            equsol[3] += 21;
            // spring equinox:      march 20
            // summer solstice:     june 21
            // autumn equinox:      september 22
            // winter solstice:     december 21     

            return equsol;
        }



    //PSA position sun algorithm
    //http://www.psa.es/sdg/sunpos.htm
    //http://www.sciencedirect.com/science/article/pii/S0038092X00001560
    private cSunCoordinates sunpos() 
    { 
        //Main variables
        double dElapsedJulianDays;
        double dDecimalHours;
        double dEclipticLongitude;
        double dEclipticObliquity;
        double dRightAscension;
        double dDeclination;

        //Auxuliary variables
        double dY;
        double dX;

        //Calculate difference in days between the current Julian Day
        //and JD 2451545.0, which is noon 1 January 2000 Universal Time
        double dJulianDate;
        long liAux1;
        long liAux2;
        //Calculate time of the day in UT decimal hours
        dDecimalHours = udtTime.dHours + (udtTime.dMinutes + udtTime.dSeconds / 60.0) / 60.0;
        //Calculate current Julian Day
        liAux1 = (udtTime.iMonth - 14) / 12;
        liAux2 = (1461 * (udtTime.iYear + 4800 + liAux1)) / 4 + (367 * (udtTime.iMonth
                    - 2 - 12 * liAux1)) / 12 - (3 * ((udtTime.iYear + 4900
                + liAux1) / 100)) / 4 + udtTime.iDay - 32075;
        dJulianDate = liAux2 - 0.5 + dDecimalHours / 24.0;
        //Calculate difference between current Julian Day and JD 2451545.0 
        dElapsedJulianDays = dJulianDate - 2451545.0;

        //Calculate ecliptic coordinates (ecliptic longitude and obliquity of the 
        //ecliptic in radians but without limiting the angle to be less than 2*Pi 
        //(i.e., the result may be greater than 2*Pi)
        double dMeanLongitude;
        double dMeanAnomaly;
        double dOmega;
        dOmega = 2.1429 - 0.0010394594 * dElapsedJulianDays;
        dMeanLongitude = 4.895063 + 0.017202791698 * dElapsedJulianDays;   //Radians
        dMeanAnomaly = 6.24006 + 0.0172019699 * dElapsedJulianDays;
        dEclipticLongitude = dMeanLongitude + 0.03341607 * Math.Sin(dMeanAnomaly) 
            + 0.00034894 * Math.Sin(2 * dMeanAnomaly) - 0.0001134 
            - 0.0000203 * Math.Sin(dOmega);
        dEclipticObliquity = 0.4090928 - 0.000000006214 * dElapsedJulianDays 
            + 0.0000396 * Math.Cos(dOmega);


        //Calculate celestial coordinates ( right ascension and declination ) in radians 
        //but without limiting the angle to be less than 2*Pi (i.e., the result may be 
        //greater than 2*Pi)
        double dSin_EclipticLongitude;
        dSin_EclipticLongitude = Math.Sin(dEclipticLongitude);
        dY = Math.Cos(dEclipticObliquity) * dSin_EclipticLongitude;
        dX = Math.Cos(dEclipticLongitude);
        dRightAscension = Math.Atan2(dY, dX);
        if (dRightAscension < 0.0) dRightAscension = dRightAscension + twopi;
        dDeclination = Math.Asin(Math.Sin(dEclipticObliquity) * dSin_EclipticLongitude);


        //Calculate local coordinates ( azimuth and zenith angle ) in degrees
        double dGreenwichMeanSiderealTime;
        double dLocalMeanSiderealTime;
        double dLatitudeInRadians;
        double dHourAngle;
        double dCos_Latitude;
        double dSin_Latitude;
        double dCos_HourAngle;
        double dParallax;
        dGreenwichMeanSiderealTime = 6.6974243242 + 0.0657098283 * dElapsedJulianDays + dDecimalHours;
        dLocalMeanSiderealTime = (dGreenwichMeanSiderealTime * 15 + udtLocation.dLongitude) * rad;
        dHourAngle = dLocalMeanSiderealTime - dRightAscension;
        dLatitudeInRadians = udtLocation.dLatitude * rad;
        dCos_Latitude = Math.Cos(dLatitudeInRadians);
        dSin_Latitude = Math.Sin(dLatitudeInRadians);
        dCos_HourAngle = Math.Cos(dHourAngle);

        cSunCoordinates sunposout;
        sunposout.dZenithAngle = (Math.Acos(dCos_Latitude * dCos_HourAngle 
            * Math.Cos(dDeclination) + Math.Sin(dDeclination) * dSin_Latitude));
        dY = -Math.Sin(dHourAngle);
        dX = Math.Tan(dDeclination) * dCos_Latitude - dSin_Latitude * dCos_HourAngle;
        sunposout.dAzimuth = Math.Atan2(dY, dX);
        if (sunposout.dAzimuth < 0.0) sunposout.dAzimuth = sunposout.dAzimuth + twopi;
        sunposout.dAzimuth = sunposout.dAzimuth / rad;
        //Parallax Correction
        dParallax = (dEarthMeanRadius / dAstronomicalUnit) * Math.Sin(sunposout.dZenithAngle);
        sunposout.dZenithAngle = (sunposout.dZenithAngle + dParallax) / rad;
        return sunposout;
    }


    //translate Spherical Coordinate System (Azimuth and Zenith) to Cartesian (XYZ)
    private cSunXYZ SunposXYZ() 
    {
        //http://ch.mathworks.com/help/matlab/ref/sph2cart.html?requestedDomain=www.mathworks.com
        double r = 1;
        //don't know, why I have to do this 90°- stuff and especially why I have to make *-1 for X...
        cSunXYZ SunposXYZout;
        SunposXYZout.x = (r * Math.Cos(rad * 90.0 - rad * udtCoordinates.dZenithAngle) * Math.Cos((rad * udtCoordinates.dAzimuth) + (rad * 90.0))) * -1.0;
        SunposXYZout.y = r * Math.Cos(rad * 90.0 - rad * udtCoordinates.dZenithAngle) * Math.Sin((rad * udtCoordinates.dAzimuth) + (rad * 90.0));
        SunposXYZout.z = r * Math.Sin(rad * 90.0 - rad * udtCoordinates.dZenithAngle);
        return SunposXYZout;
     }

        
        public struct cSunCoordinates
        {
            public double dZenithAngle;
            public double dAzimuth;
        }
        public struct cSunXYZ
        {
            public double x;
            public double y;
            public double z;
        }





        public static void Create8760SunVectors(ref List<SunVector> sunvectors, double longitude, double latitude, int year)
        {
            for (int m = 1; m <= 12; m++)
            {
                int daysInMonth = System.DateTime.DaysInMonth(year, m);
                for (int d = 1; d <= daysInMonth; d++)
                {
                    for (int i = 0; i <= 23; i++)
                    {
                        SunVector sunvec = new SunVector(year, m, d, i, 0, 0, longitude, latitude);
                        sunvectors.Add(sunvec);
                    }
                }
            }
        }
    }

}
