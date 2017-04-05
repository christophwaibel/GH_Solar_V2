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
            ////Mesh mesh = new Mesh();
            ////mesh.Vertices.Add(0.0, 0.0, 1.0);   //0
            ////mesh.Vertices.Add(1.0, 0.0, 1.0);   //1
            ////mesh.Vertices.Add(1.0, 1.0, 1.0);   //2
            ////mesh.Vertices.Add(0.0, 1.0, 1.0);   //3

            ////mesh.Faces.AddFace(0, 1, 2, 3);

            //Ray3d ray = new Ray3d(new Point3d(0.5,0.5,0), new Vector3d(0,0,1));
            //Point3d p = ray.PointAt(0.1);

            ////double inters = Rhino.Geometry.Intersect.Intersection.MeshRay(mesh, ray);

            ////Console.WriteLine(inters);
            //Console.WriteLine("x: {0}, y: {1}, z: {2}", p.X, p.Y, p.Z);
            //Console.ReadKey();

            Console.WriteLine("hi");

            //IcoSphere ico = new IcoSphere(0);
            //List<int[]> test = ico.getFaces();
            //List<double[]> coords = ico.getVertexCoordinates();

            ////coords.RemoveAt(1);
            


            //Stopwatch watch = new Stopwatch();
            //watch.Start();
            ////SkyDome dome = new SkyDome(2, 2013, 47.3673,8.55 );        //!!!!!!!!!!!!!!!
            //Console.WriteLine(Convert.ToInt32(false));
            //Console.WriteLine("hi");
            //watch.Stop();
            //Console.WriteLine(watch.Elapsed.TotalMilliseconds);
            //Console.ReadKey();




            //load in weather file DHI
            List<double> DHI = new List<double>();
            List<double> DNI = new List<double>();
            int counter = 0;
            string line;

            // Read the file and display it line by line.
            System.IO.StreamReader file =
               new System.IO.StreamReader("C:\\Users\\wach\\Desktop\\SolarV2\\DHI.txt");
            while ((line = file.ReadLine()) != null)
            {
                DHI.Add(Convert.ToDouble(line));
                //Console.WriteLine(line);
                counter++;
            }
            file.Close();

            System.IO.StreamReader file2 = new System.IO.StreamReader("C:\\Users\\wach\\Desktop\\SolarV2\\DNI.txt");
            while ((line = file2.ReadLine()) != null)
            {
                DNI.Add(Convert.ToDouble(line));
                //Console.WriteLine(line);
                counter++;
            }
            file.Close();





            //Sunvectors are always the same for the location and year
            int recursion = 2;
            double longitude = 8.5500025;
            double latitude = 47.367347;
            int year = 2016;



            List<SunVector> sunvectors = new List<SunVector>(); 
            Context.Create8760SunVectors(ref sunvectors, longitude, latitude, year);
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

            Sensorpoint p = new Sensorpoint(year, weather, location, sunvectors,90,30,recursion);
            p.CalcIrradiation();

            //for (int i = 0; i < p.I.Count(); i++)
            //{
            //    Console.WriteLine(p.I[i]);
            //}

     


            Console.ReadKey();
        }
    }
}
