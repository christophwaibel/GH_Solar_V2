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
               new System.IO.StreamReader("C:\\Users\\Chris\\Desktop\\WORK\\DHI_sarah.txt");
            while ((line = file.ReadLine()) != null)
            {
                DHI.Add(Convert.ToDouble(line));
                //Console.WriteLine(line);
                counter++;
            }
            file.Close();

            System.IO.StreamReader file2 = new System.IO.StreamReader("C:\\Users\\Chris\\Desktop\\WORK\\DNI_sarah.txt");
            while ((line = file2.ReadLine()) != null)
            {
                DNI.Add(Convert.ToDouble(line));
                //Console.WriteLine(line);
                counter++;
            }
            file.Close();





            //Sunvectors are always the same for the location and year
            int recursion = 2;
            double longitude = 8.539;
            double latitude = 47.370;
            int year = 2000;



            List<SunVector> sunvectors = new List<SunVector>(); 
            SunVector.Create8760SunVectors(ref sunvectors, longitude, latitude, year);
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

            //Sensorpoint p = new Sensorpoint(year, weather, location, sunvectors, 90, 30, recursion);
            //p.CalcIrradiation();

           







            //int[] a = SunVector.GetEquinoxSolstice(2010);

            double []beta = new double[1]{0};
            double [] psi=new double[1]{0};
            Sensorpoints.p3d [] coord  = new Sensorpoints.p3d[1];
            coord[0].X=0;
            coord[0].Y=0;
            coord[0].Z=0;
            Sensorpoints.v3d[] normal = new Sensorpoints.v3d[1];
            normal[0].X=0;
            normal[0].Y=1;
            normal[0].Z=0;

            Sensorpoints p = new Sensorpoints(beta, psi, coord, normal, recursion);
            p.CalcIrradiation(weather, sunvectors.ToArray());

            System.IO.StreamWriter write = new System.IO.StreamWriter("C:\\Users\\Chris\\Desktop\\sara.txt");
            System.IO.StreamWriter write2 = new System.IO.StreamWriter("C:\\Users\\Chris\\Desktop\\sarabeam.txt");
            System.IO.StreamWriter write3 = new System.IO.StreamWriter("C:\\Users\\Chris\\Desktop\\saradiff.txt");
            for (int i = 0; i < p.I[0].Count(); i++)
            {
                //Console.WriteLine(p.I[0][i]);
                write.WriteLine(p.I[0][i]);
                write2.WriteLine(p.Ibeam[0][i]);
                write3.WriteLine(p.Idiff[0][i]);
            }
            write.Close();
            write2.Close();
            write3.Close();

            Console.ReadKey();
        }
    }
}
