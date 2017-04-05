using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using Rhino.Geometry;

/*
 * cShadow.cs
 * Copyright 2017 Christoph Waibel <chwaibel@student.ethz.ch>
 * 
 * This work is licensed under the GNU GPL license version 3.
*/

namespace GHSolar
{

    internal static class cShadow
    {
        //Direct shadow list

        //diffuse shadow list


        //inter-reflections? n


        /// <summary>
        /// 
        /// </summary>
        /// <param name="origin"></param>
        /// <param name="origNormal"></param>
        /// <param name="tolerance">Tolerance to offset ray from sensor point.</param>
        /// <param name="vec"></param>
        /// <param name="obstacles"></param>
        /// <param name="shdw"></param>
        internal static void CalcShadow(Point3d origin, Vector3d origNormal, double tolerance, Vector3d[] vec, Mesh[] obstacles, ref bool[] shdw)
        {
            //foreach vector

            //foreach obstacle

            shdw = new bool[vec.Length];        //by default all elements false. true means it's obstructed
            Point3d origOffset = new Point3d(Point3d.Add(origin, Vector3d.Multiply(Vector3d.Divide(origNormal, origNormal.Length), tolerance)));
            for (int i = 0; i < vec.Length; i++)
            {

                for (int u = 0; u < obstacles.Length; u++)
                {
                    Ray3d ray = new Ray3d(origOffset, vec[i]);
                    double inters = Rhino.Geometry.Intersect.Intersection.MeshRay(obstacles[u], ray);
                    if (inters >= 0)
                    {
                        shdw[i] = true;
                        break;
                    }
                }

            }


        }


    }
}
