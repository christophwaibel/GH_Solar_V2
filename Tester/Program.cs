using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

//using Rhino.Geometry;
//using Rhino;

using SolarModel;
using System.Diagnostics;




namespace Tester
{
    class Program
    {
        static void Main(string[] args)
        {
            Console.WriteLine("SolarModel V2.0");
            Console.WriteLine();
            Console.WriteLine("You will need: 2 txt files with hourly values for DHI and DNI respectively. They need to be named DNI.txt and DHI.txt");
            Console.WriteLine("Please define path for reading inputs and writing result files:");
            string path = Console.ReadLine(); //"C:\Users\wach\Desktop
            Console.WriteLine();

            //load in weather file DHI
            List<double> DHI = new List<double>();
            List<double> DNI = new List<double>();
            int counter = 0;
            string line;

            // Read the file and display it line by line.
            Console.WriteLine("Reading inputs...");
            System.IO.StreamReader file =
               new System.IO.StreamReader(path + "\\DHI.txt");
            while ((line = file.ReadLine()) != null)
            {
                DHI.Add(Convert.ToDouble(line));
                //Console.WriteLine(line);
                counter++;
            }
            file.Close();

            System.IO.StreamReader file2 = new System.IO.StreamReader(path + "\\DNI.txt");
            while ((line = file2.ReadLine()) != null)
            {
                DNI.Add(Convert.ToDouble(line));
                //Console.WriteLine(line);
                counter++;
            }
            file.Close();
            Console.WriteLine();



            int recursion = 2;      //resolution of skydome

            int year;
            double longitude, latitude;
            Console.WriteLine("Longitude: ");
            if (!double.TryParse(Console.ReadLine(), out longitude))
            {
                Console.WriteLine("You need to input a real number. Hit any key to close.");
                Console.ReadLine();
            }
            Console.WriteLine("Latitude: ");
            if (!double.TryParse(Console.ReadLine(), out latitude))
            {
                Console.WriteLine("You need to input a real number. Hit any key to close.");
                Console.ReadLine();
            }
            Console.WriteLine("Year: ");
            if (!int.TryParse(Console.ReadLine(), out year))
            {
                Console.WriteLine("You need to input an integer number. Hit any key to close.");
                Console.ReadLine();
            }
            Console.WriteLine();
            //double longitude = 8.539;
            //double latitude = 47.370;
            //int year = 2005;



            List<SunVector> sunvectors; 
            SunVector.Create8760SunVectors(out sunvectors, longitude, latitude, year);
            Context.cWeatherdata weather;
            weather.DHI = new List<double>();
            weather.DNI = new List<double>();
            weather.Snow = new List<double>();
            for (int i = 0; i < 8760; i++)
            {
                weather.DHI.Add(DHI[i]);
                weather.DNI.Add(DNI[i]);
            }

            Context.cLocation location;
            location.dLatitude = latitude;
            location.dLongitude = longitude;
            location.dTgmt = 1;


            Dictionary<string, double> albedos = new Dictionary<string, double>();
            albedos.Add("LAWN", 0.205);
            albedos.Add("UNTILTEDFIELD", 0.26);
            albedos.Add("NAKEDGROUND", 0.17);
            albedos.Add("CONCRETE", 0.3);
            albedos.Add("SNOW", 0.85);
            albedos.Add("OLDSNOW", 0.58);

            Console.WriteLine("Ground albedo: allowed inputs 'LAWN', 'UNTILTEDFIELD', 'NAKEDGROUND', 'CONCRETE', 'SNOW', 'OLDSNOW'.");
            string albedo_string = Console.ReadLine();
            double albedo1 = albedos[albedo_string];
            double[] albedo = new double[8760];
            for (int t = 0; t < 8760; t++)
            {
                albedo[t] = albedo1;
            }

            double beta_in, psi_in;
            Console.WriteLine("Tilt angle in degree: ");
            if (!double.TryParse(Console.ReadLine(), out beta_in))
            {
                Console.WriteLine("You need to input a real number. Hit any key to close.");
                Console.ReadLine(); 
                return;
            }
            Console.WriteLine("Azimuth angle in degree (North is 0, South is 180): ");
            if (!double.TryParse(Console.ReadLine(), out psi_in))
            {
                Console.WriteLine("You need to input a real number. Hit any key to close.");
                Console.ReadLine();
                return;
            }
            Console.WriteLine();

            //double []beta = new double[1]{20};
            //double [] psi=new double[1]{180};
            double[] beta = new double[1]{beta_in};
            double [] psi=new double[1]{psi_in};
            Sensorpoints.p3d [] coord  = new Sensorpoints.p3d[1];   //dummy variables. will not be used in this simplified simulation
            coord[0].X=0;
            coord[0].Y=0;
            coord[0].Z=0;
            Sensorpoints.v3d[] normal = new Sensorpoints.v3d[1];   //dummy variables. will not be used in this simplified simulation
            normal[0].X=0;
            normal[0].Y=1;
            normal[0].Z=0;

            Console.WriteLine("Calculating irradiation...");
            Sensorpoints p = new Sensorpoints(beta, psi, coord, normal, recursion);
            p.SetSimpleSky(beta);
            p.SetSimpleGroundReflection(beta, albedo, weather, sunvectors.ToArray());
            p.CalcIrradiation(weather, sunvectors.ToArray());

            Console.WriteLine("Writing to path...");
            System.IO.StreamWriter write = new System.IO.StreamWriter(path + "\\calc.txt");
            System.IO.StreamWriter write2 = new System.IO.StreamWriter(path + "\\calcbeam.txt");
            System.IO.StreamWriter write3 = new System.IO.StreamWriter(path + "\\calcdiff.txt");
            System.IO.StreamWriter write4 = new System.IO.StreamWriter(path + "\\calcgroundrefl.txt");
            for (int i = 0; i < p.I[0].Count(); i++)
            {
                //Console.WriteLine(p.I[0][i]);
                write.WriteLine(p.I[0][i]);
                write2.WriteLine(p.Ibeam[0][i]);
                write3.WriteLine(p.Idiff[0][i]);
                write4.WriteLine(p.Irefl_diff[0][i]);
            }
            write.Close();
            write2.Close();
            write3.Close();
            write4.Close();
            Console.WriteLine();
            Console.WriteLine("Done. Press any key to quit");
            Console.ReadKey();
        }
    }
}
