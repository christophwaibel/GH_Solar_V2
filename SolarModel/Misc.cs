using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

/*
 * Misc.cs
 * Copyright 2017 Christoph Waibel <chwaibel@student.ethz.ch>
 * 
 * This work is licensed under the GNU GPL license version 3.
*/

namespace SolarModel
{
    public static class Misc
    {
        public static bool IsNullOrEmpty<T>(this ICollection<T> collection)
        {
            return collection == null || collection.Count == 0;
        }

        /// <summary>
        /// Calculates the area of a 3D triangle.
        /// </summary>
        /// <param name="p1">Corner point 1. double[2] {x,y,z}.</param>
        /// <param name="p2">Corner point 2. double[2] {x,y,z}.</param>
        /// <param name="p3">Corner point 3. double[2] {x,y,z}.</param>
        /// <returns>The area of the 3d triangle.</returns>
        public static double getTriangleArea(double[] p1, double[] p2, double[] p3)
        {
            double area = 0.5 * Math.Sqrt(
                Math.Pow((p2[0] * p1[1]) - (p3[0] * p1[1]) - (p1[0] * p2[1]) + (p3[0] * p2[1]) + (p1[0] * p3[1]) - (p2[0] * p3[1]), 2) +
                Math.Pow((p2[0] * p1[2]) - (p3[0] * p1[2]) - (p1[0] * p2[2]) + (p3[0] * p2[2]) + (p1[0] * p3[2]) - (p2[0] * p3[2]), 2) +
                Math.Pow((p2[1] * p1[2]) - (p3[1] * p1[2]) - (p1[1] * p2[2]) + (p3[1] * p2[2]) + (p1[1] * p3[2]) - (p2[1] * p3[2]), 2));
            return area;
        }
    }
}
