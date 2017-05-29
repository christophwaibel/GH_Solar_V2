using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

/*
 * Irradiation.cs
 * Copyright 2017 Christoph Waibel <chwaibel@student.ethz.ch>
 * 
 * This work is licensed under the GNU GPL license version 3.
*/

namespace SolarModel
{
    /// <summary>
    /// Solar Irradiation on a sensor point.
    /// <para>Decomposed into beam and diffuse radiation. </para>
    /// <para>Diffuse radiation is further decomposed into components for (i) horizon, (ii) circumsolar and (iii) sky.</para>
    /// </summary>
    /// <remarks>
    /// Source: Perez, R., Ineichen, P., Seals, R., Michalsky, J., Stewart, R. (1990). Modeling Daylight Availabiltiy and Irradiance Components from Direct and Global Irradiance. Solar Energy Vol. 44, No. 5, pp. 271-289
    /// Additionally, coefficients from EnergyPlus V.8.5.0. Engineering Reference, p.184-185, table 4.3 have been used, instead of those from Perez et al. (1990). 
    /// </remarks>
    public static class Irradiation
    {
        private const double rad = Math.PI / 180;
        private const double pi = Math.PI;
        private const double twopi = 2 * pi;
        private const double dEarthMeanRadius = 6371.01;        //in km
        private const double dAstronomicalUnit = 149597890.0;   //in km
        private const double Isc = 1367.0;                      //solar constant in W/m2



        /// <summary>
        /// Diffuse radiation, considering anisotropic sky. See Perez 1990, or EnergyPlus Engineering Documentation.
        /// </summary>
        /// <param name="DHI">Diffuse Horizontal Irradiation (DHI). From a weatherfile, e.g. *.epw.</param>
        /// <param name="DNI">Direct Normal Irradiation (DNI). From a weatherfile, e.g. *.epw.</param>
        /// <param name="θZ">Solar Zenith, in degree.</param>
        /// <param name="θA">Solar Azimuth, in degree.</param>
        /// <param name="θβ">Analysis surface tilt angle.</param>
        /// <param name="θAsrf">Analysis surface azimuth.</param>
        /// <param name="DOY">Day of year (DOY).</param>
        ///<returns>Diffuse radiation of a sensor point for one moment (e.g. hour) of the year.
        ///[0] is total diffuse, [1]: horizon, [2]: dome, [3]: circumsolar.</returns>
        public static double[] Diffuse_old(double DHI, double DNI, double θZ, double θA, double θβ, double θAsrf, int DOY)
        {
            //relative optical air mass
            double m = 1.0 / Math.Cos(rad * θZ);
            if (m > 40) m = 40; //because the simple model is too simple. value could become too high with high zenith angles

            //extraterrestrial irradiance
            //double I0 = 1353.0;           //constant assumed by energyplus
            double b = 2 * pi * (Convert.ToDouble(DOY) / 365.0);// *rad;
            double I0 = Isc * (1.00011 + 0.034221 * Math.Cos(b) + 0.00128 * Math.Sin(b) + 0.000719 * Math.Cos(2 * b) + 0.000077 * Math.Sin(2 * b));

            //Sky Brightness factor
            double Δ = DHI * m / I0;

            //Sky clearness factor (ε)
            double ε = ((DHI + DNI) / DHI + 1.041 * Math.Pow((θZ * rad), 3)) / (1 + 1.041 * Math.Pow((θZ * rad), 3));

            Func<int, double, double> Fij = (ij, εcheck) =>
            {
                //Perez coefficients, from EnergyPlus 6.5 Engineering Documentation p.185. More precise than in Perez 1990.
                double[,] Fijmatrix = new double[,] 
                { 
                { -0.0083117, 0.1299457, 0.3296958, 0.5682053, 0.8730280, 1.1326077, 1.0601591, 0.6777470},             //F11
                { 0.5877285, 0.6825954, 0.4868735, 0.1874525, -0.3920403, -1.2367284, -1.5999137, -0.3272588},          //F12
                { -0.0620636, -0.1513752, -0.2210958, -0.2951290, -0.3616149, -0.4118494, -0.3589221, -0.2504286},      //F13
                { -0.0596012, -0.0189325, 0.0554140, 0.1088631, 0.2255647, 0.2877813, 0.2642124, 0.1561313},            //F21
                { 0.0721249, 0.0659650, -0.0639588, -0.1519229, -0.4620442, -0.8230357, -1.1272340, -1.3765031},        //F22
                { -0.0220216, -0.0288748, -0.0260542, -0.0139754, 0.0012448, 0.0558651, 0.1310694, 0.2506212}           //F23
                };
                double Fijvalue = 0.0;
                if (εcheck < 1.065) Fijvalue = Fijmatrix[ij, 0];
                else if (εcheck >= 1.065 && εcheck < 1.230) Fijvalue = Fijmatrix[ij, 1];
                else if (εcheck >= 1.230 && εcheck < 1.500) Fijvalue = Fijmatrix[ij, 2];
                else if (εcheck >= 1.500 && εcheck < 1.950) Fijvalue = Fijmatrix[ij, 3];
                else if (εcheck >= 1.950 && εcheck < 2.800) Fijvalue = Fijmatrix[ij, 4];
                else if (εcheck >= 2.800 && εcheck < 4.500) Fijvalue = Fijmatrix[ij, 5];
                else if (εcheck >= 4.500 && εcheck < 6.200) Fijvalue = Fijmatrix[ij, 6];
                else Fijvalue = Fijmatrix[ij, 7];//(ε >= 6.200)
                return Fijvalue;
            };

            double F1 = Math.Max(0, Fij(0, ε) + Fij(1, ε) * Δ + Fij(2, ε) * (θZ * rad));
            double F2 = Fij(3, ε) + Fij(4, ε) * Δ + Fij(5, ε) * (θZ * rad);

            //angle of incidence
            double AOI = Math.Acos(Math.Cos(θZ * rad) * Math.Cos(θβ * rad) + Math.Sin(θZ * rad) * Math.Sin(θβ * rad) * Math.Cos((θA - θAsrf) * rad));
            double a = Math.Max(0, Math.Cos(AOI));
            b = Math.Max(Math.Cos(85 * rad), Math.Cos(θZ * rad));

            // Dhorizon + Ddome + Dcircumsolar
            //double D = DHI * ((1 - F1) *((1 + Math.Cos(θβ * rad)) / 2) + F1 * (a / b) + F2 * Math.Sin(θβ * rad));
            double Dhorizon = DHI * F2 * Math.Sin(θβ * rad);
            double Ddome = DHI * (1 - F1) * (1 + Math.Cos(θβ * rad)) / 2;
            double Dcircum = DHI * F1 * (a / b);
            double D = Dhorizon + Ddome + Dcircum;
            return new double[4] { D, Dhorizon, Ddome, Dcircum };


        }

        /// <summary>
        /// Diffuse radiation, considering anisotropic sky. See Perez 1990, or EnergyPlus Engineering Documentation.
        /// </summary>
        /// <param name="DHI">Diffuse Horizontal Irradiation (DHI). From a weatherfile, e.g. *.epw.</param>
        /// <param name="DOY">Direct Normal Irradiation (DNI). From a weatherfile, e.g. *.epw.</param>
        /// <param name="θZ">Solar Zenith, in degree.</param>
        /// <param name="θA">Solar Azimuth, in degree.</param>
        /// <param name="θβ">Analysis surface tilt angle.</param>
        /// <param name="θAsrf">Analysis surface azimuth.</param>
        /// <param name="DNI">Day of year (DOY).</param>
        /// <param name="horizonshdw">0-1 fraction of obstructed horizon. (1 = fully obstructed; 0 = no obstruction)</param>
        /// <param name="domeshdw">0-1 fraction of the obstructed skydome. (1 = fully obstructed; 0 = no obstruction)</param>
        /// <param name="circumsolshdw">Boolean, indicating if solar vector of that moment is obstructed. (true = obstructed; false = no obstruction)</param>
        ///<returns>Diffuse radiation of a sensor point for one moment (e.g. hour) of the year.</returns>
        public static double[] Diffuse_old(double DHI, double DNI, double θZ, double θA, double θβ, double θAsrf, int DOY, double horizonshdw, double domeshdw, bool circumsolshdw)
        {
            double[] D4 = Diffuse_old(DHI, DNI, θZ, θA, θβ, θAsrf, DOY);
            double Dhorizon = D4[1] * (1 - horizonshdw);
            double Ddome = D4[2] * (1 - domeshdw);
            double Dcircumsolar = D4[3] * (1 - Convert.ToInt32(circumsolshdw));

            return new double[4] { Dhorizon + Ddome + Dcircumsolar, Dhorizon, Ddome, Dcircumsolar };
        }



        /// <summary>
        /// Diffuse radiation, considering anisotropic sky. See Perez 1990, or EnergyPlus Engineering Documentation.
        /// </summary>
        /// <param name="DHI">Diffuse Horizontal Irradiation (DHI). From a weatherfile, e.g. *.epw.</param>
        /// <param name="DOY">Direct Normal Irradiation (DNI). From a weatherfile, e.g. *.epw.</param>
        /// <param name="θZ">Solar Zenith, in degree.</param>
        /// <param name="θA">Solar Azimuth, in degree.</param>
        /// <param name="θβ">Analysis surface tilt angle.</param>
        /// <param name="θAsrf">Analysis surface azimuth.</param>
        /// <param name="DNI">Day of year (DOY).</param>
        /// <param name="horizonshdw">0-1 fraction of obstructed horizon. (1 = fully obstructed; 0 = no obstruction)</param>
        /// <param name="domeshdw">0-1 fraction of the obstructed skydome. (1 = fully obstructed; 0 = no obstruction)</param>
        /// <param name="circumsolshdw">Value indicating how much solar vector of that moment is obstructed. (1.0 = fully obstructed; 0.0 = no obstruction)</param>
        ///<returns>Diffuse radiation of a sensor point for one moment (e.g. hour) of the year.</returns>
        public static double Diffuse(double DHI, double DNI, double θZ, double θA, double θβ, double θAsrf, int DOY, double horizonshdw, double domeshdw, double circumsolshdw)
        {
            double D = 0.0;
            if (θZ > 90.0)    //after sunset, could still be an hour of diffuse light. assume isotropic sky
            {
                D = DHI * (1.0 - domeshdw);
            }
            else
            {
                //relative optical air mass
                double m = 1.0 / Math.Cos(rad * θZ);
                if (m > 40.0) m = 40.0; //because the simple model is too simple. value could become too high with high zenith angles

                //extraterrestrial irradiance
                //double I0 = 1353.0;           //constant assumed by energyplus
                double b = 2.0 * pi * (Convert.ToDouble(DOY) / 365.0);// *rad;
                double I0 = Isc * (1.00011 + 0.034221 * Math.Cos(b) + 0.00128 * Math.Sin(b) + 0.000719 * Math.Cos(2.0 * b) + 0.000077 * Math.Sin(2.0 * b));

                //Sky Brightness factor
                double Δ = DHI * m / I0;

                //Sky clearness factor (ε)
                double ε = ((DHI + DNI) / DHI + 1.041 * Math.Pow((θZ * rad), 3)) / (1 + 1.041 * Math.Pow((θZ * rad), 3));

                Func<int, double, double> Fij = (ij, εcheck) =>
                {
                    //Perez coefficients, from EnergyPlus 6.5 Engineering Documentation p.185. More precise than in Perez 1990.
                    double[,] Fijmatrix = new double[,] 
                { 
                { -0.0083117, 0.1299457, 0.3296958, 0.5682053, 0.8730280, 1.1326077, 1.0601591, 0.6777470},             //F11
                { 0.5877285, 0.6825954, 0.4868735, 0.1874525, -0.3920403, -1.2367284, -1.5999137, -0.3272588},          //F12
                { -0.0620636, -0.1513752, -0.2210958, -0.2951290, -0.3616149, -0.4118494, -0.3589221, -0.2504286},      //F13
                { -0.0596012, -0.0189325, 0.0554140, 0.1088631, 0.2255647, 0.2877813, 0.2642124, 0.1561313},            //F21
                { 0.0721249, 0.0659650, -0.0639588, -0.1519229, -0.4620442, -0.8230357, -1.1272340, -1.3765031},        //F22
                { -0.0220216, -0.0288748, -0.0260542, -0.0139754, 0.0012448, 0.0558651, 0.1310694, 0.2506212}           //F23
                };
                    double Fijvalue = 0.0;
                    if (εcheck < 1.065) Fijvalue = Fijmatrix[ij, 0];
                    else if (εcheck >= 1.065 && εcheck < 1.230) Fijvalue = Fijmatrix[ij, 1];
                    else if (εcheck >= 1.230 && εcheck < 1.500) Fijvalue = Fijmatrix[ij, 2];
                    else if (εcheck >= 1.500 && εcheck < 1.950) Fijvalue = Fijmatrix[ij, 3];
                    else if (εcheck >= 1.950 && εcheck < 2.800) Fijvalue = Fijmatrix[ij, 4];
                    else if (εcheck >= 2.800 && εcheck < 4.500) Fijvalue = Fijmatrix[ij, 5];
                    else if (εcheck >= 4.500 && εcheck < 6.200) Fijvalue = Fijmatrix[ij, 6];
                    else Fijvalue = Fijmatrix[ij, 7];//(ε >= 6.200)
                    return Fijvalue;
                };

                double F1 = Math.Max(0, Fij(0, ε) + Fij(1, ε) * Δ + Fij(2, ε) * (θZ * rad));
                double F2 = Fij(3, ε) + Fij(4, ε) * Δ + Fij(5, ε) * (θZ * rad);

                //angle of incidence
                double AOI = Math.Acos(Math.Cos(θZ * rad) * Math.Cos(θβ * rad) + Math.Sin(θZ * rad) * Math.Sin(θβ * rad) * Math.Cos((θA - θAsrf) * rad));
                double a = Math.Max(0, Math.Cos(AOI));
                b = Math.Max(Math.Cos(85 * rad), Math.Cos(θZ * rad));

                // Dhorizon + Ddome + Dcircumsolar
                //double D = DHI * ((1 - F1) * ((1 + Math.Cos(θβ * rad)) / 2) + F1 * (a / b) + F2 * Math.Sin(θβ * rad));
                double Dhorizon = (DHI * F2 * Math.Sin(θβ * rad)) * (1.0 - horizonshdw);
                double Ddome = (DHI * (1 - F1) * (1 + Math.Cos(θβ * rad)) / 2.0) * (1.0 - domeshdw);
                double Dcircum = (DHI * F1 * (a / b)) * (1.0 - circumsolshdw);
                D = Dhorizon + Ddome + Dcircum;
            }
            return Math.Max(0.0, D);
        }









        /// <summary>
        /// Unobstructed Beam (or direct) radiation, depending on day of year, surface tilt angle, local time, and direct normal irradiation from a weather file.
        /// </summary>
        /// <param name="DNI">Direct normal irradiation for a local time (e.g. hour). Data e.g. from a weather file.</param>
        /// <param name="β">Analysis surface tilt angle.</param>
        /// <param name="φ">Latitude of the location.</param>
        /// <param name="λ">Longitude of the location.</param>
        /// <param name="ψ">Analysis surface azimuth (orientation measured from South to West).</param>
        /// <param name="DOY">Day of year</param>
        /// <param name="LT">Local time.</param>
        /// <param name="δTgmt">difference Local Time (LT) from Greenwich Mean Tmie (GMT).</param>
        /// <returns>Beam radiation in [W/m2] on an analysis surface for one moment (e.g. hour) of the year.</returns>
        public static double Beam(double DNI, double β, double φ, double λ, double ψ, int DOY, int LT, int δTgmt)
        {
            //Equation of Time (EoT)
            double Bangle = (360 / 365) * (DOY - 81);      //in degrees
            double EoT = 9.87 * Math.Sin(rad * (2 * Bangle)) - 7.53 * Math.Cos(rad * Bangle) - 1.5 * Math.Sin(rad * Bangle);

            ///Local Standard Time Meridian (LSTM)
            double LSTM = 15.0 * Convert.ToDouble(δTgmt);

            //Time Correction Factor (TC)
            double TC = 4.0 * (λ - LSTM) + EoT;

            //Local Solar Time (LST)
            double LST = Convert.ToDouble(LT) + (TC / 60.0);

            //Hour angle (HRA)
            double HRA = 15 * (LST - 12);

            //declination angle
            double δ = Math.Asin(Math.Sin(rad * 23.45) * Math.Sin(rad * ((360 / 365) * (DOY - 81))));

            double B = DNI * (Math.Sin(rad * δ) * Math.Sin(rad * φ) * Math.Cos(rad * β) -
                                     Math.Sin(rad * δ) * Math.Cos(rad * φ) * Math.Sin(rad * β) * Math.Cos(rad * ψ) +
                                     Math.Cos(rad * δ) * Math.Cos(rad * φ) * Math.Cos(rad * β) * Math.Cos(rad * HRA) +
                                     Math.Cos(rad * δ) * Math.Sin(rad * φ) * Math.Sin(rad * β) * Math.Cos(rad * ψ) * Math.Cos(rad * HRA) +
                                     Math.Cos(rad * δ) * Math.Sin(rad * ψ) * Math.Sin(rad * HRA) * Math.Sin(rad * β));
            return Math.Max(0, B);
        }



        /// <summary>
        /// Unobstructed Beam (or direct) radiation, depending on angle of incidence (solar zenith and azimuth, surface tilt and azimuth), and direct normal irradiation from a weather file.
        /// </summary>
        /// <param name="DNI">Direct Normal Irradiance (DNI). From a weatherfile, e.g. *.epw.</param>
        /// <param name="θZ">Solar Zenith, in degree.</param>
        /// <param name="θA">Solar Azimuth, in degree.</param>
        /// <param name="θβ">Analysis surface tilt angle.</param>
        /// <param name="θAsrf">Analysis surface azimuth.</param>
        /// <param name="shadow">Value indicating how much solar vector of that moment is obstructed. (1.0 = fully obstructed; 0.0 = no obstruction)</param>
        /// <returns></returns>
        public static double Beam(double DNI, double θZ, double θA, double θβ, double θAsrf, double shadow)
        {
            double B = DNI * (Math.Cos(θZ * rad) * Math.Cos(θβ * rad) + Math.Sin(θZ * rad) * Math.Sin(θβ * rad) * Math.Cos((θA - θAsrf) * rad));
            return Math.Max(0, B * (1.0 - shadow)); 
        }





    }




}
