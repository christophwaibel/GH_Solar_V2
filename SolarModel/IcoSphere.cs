﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Media.Media3D;

/*
 * IcoSphere.cs
 * Copyright 2017 Christoph Waibel <chwaibel@student.ethz.ch>
 * Modified after: http://blog.andreaskahler.com/search/label/c%23
 * 
 * This work is licensed under the GNU GPL license version 3.
*/

namespace SolarModel
{
    /// <summary>
    /// Create an icosphere with varying resolution.
    /// Modified after: http://blog.andreaskahler.com/search/label/c%23
    /// </summary>
    public class IcoSphere
    {
        /// <summary>
        /// Mesh Geometry of the IcoSphere as Windows.Media.Media3D object.
        /// </summary>
        public MeshGeometry3D IcoSphere3D;

        /// <summary>
        /// Create instance of IcoSphere.
        /// </summary>
        /// <param name="recursionLevel">Recursion Level, i.e. resolution of the icosphere.</param>
        public IcoSphere(int recursionLevel)
        {
            IcoSphere3D = CreateGeometry(recursionLevel);
        }

        /// <summary>
        /// Get a list of faces of the IcoSphere. 
        /// </summary>
        /// <returns>List of integer arrays of length 3, indicating vertex indices of each face respectively.</returns>
        public List<int[]> getFaces()
        {
            MeshGeometry3D geo = this.IcoSphere3D;
            List<int[]> facesOut = new List<int[]>();
            for (int i = 0; i < geo.TriangleIndices.Count; i+=3)
            {
                facesOut.Add(new int[3] { geo.TriangleIndices[i], geo.TriangleIndices[i+1],geo.TriangleIndices[i+2]});
            }
            return facesOut;
        }

        /// <summary>
        /// Get a list of 3d vertex coordinates.
        /// </summary>
        /// <returns>List of double array of length 3 containing x,y,z -coordinates of each vertex.</returns>
        public List<double[]> getVertexCoordinates()
        {
            MeshGeometry3D geo = this.IcoSphere3D;
            List<double[]> coordsOut = new List<double[]>();
            for (int i = 0; i < geo.Positions.Count; i++)
            {
                coordsOut.Add(new double[3] { geo.Positions[i].X, geo.Positions[i].Y, geo.Positions[i].Z});
            }
            return coordsOut;
        }

        private MeshGeometry3D CreateGeometry(int recursionLevel)
        {
            this.geometry = new MeshGeometry3D();
            this.middlePointIndexCache = new Dictionary<long, int>();
            this.index = 0;

            // create 12 vertices of a icosahedron
            var t = (1.0 + Math.Sqrt(5.0)) / 2.0;

            addVertex(new Point3D(-1, t, 0));
            addVertex(new Point3D(1, t, 0));
            addVertex(new Point3D(-1, -t, 0));
            addVertex(new Point3D(1, -t, 0));

            addVertex(new Point3D(0, -1, t));
            addVertex(new Point3D(0, 1, t));
            addVertex(new Point3D(0, -1, -t));
            addVertex(new Point3D(0, 1, -t));

            addVertex(new Point3D(t, 0, -1));
            addVertex(new Point3D(t, 0, 1));
            addVertex(new Point3D(-t, 0, -1));
            addVertex(new Point3D(-t, 0, 1));


            // create 20 triangles of the icosahedron
            var faces = new List<TriangleIndices>();

            // 5 faces around point 0
            faces.Add(new TriangleIndices(0, 11, 5));
            faces.Add(new TriangleIndices(0, 5, 1));
            faces.Add(new TriangleIndices(0, 1, 7));
            faces.Add(new TriangleIndices(0, 7, 10));
            faces.Add(new TriangleIndices(0, 10, 11));

            // 5 adjacent faces 
            faces.Add(new TriangleIndices(1, 5, 9));
            faces.Add(new TriangleIndices(5, 11, 4));
            faces.Add(new TriangleIndices(11, 10, 2));
            faces.Add(new TriangleIndices(10, 7, 6));
            faces.Add(new TriangleIndices(7, 1, 8));

            // 5 faces around point 3
            faces.Add(new TriangleIndices(3, 9, 4));
            faces.Add(new TriangleIndices(3, 4, 2));
            faces.Add(new TriangleIndices(3, 2, 6));
            faces.Add(new TriangleIndices(3, 6, 8));
            faces.Add(new TriangleIndices(3, 8, 9));

            // 5 adjacent faces 
            faces.Add(new TriangleIndices(4, 9, 5));
            faces.Add(new TriangleIndices(2, 4, 11));
            faces.Add(new TriangleIndices(6, 2, 10));
            faces.Add(new TriangleIndices(8, 6, 7));
            faces.Add(new TriangleIndices(9, 8, 1));


            // refine triangles
            for (int i = 0; i < recursionLevel; i++)
            {
                var faces2 = new List<TriangleIndices>();
                foreach (var tri in faces)
                {
                    // replace triangle by 4 triangles
                    int a = getMiddlePoint(tri.v1, tri.v2);
                    int b = getMiddlePoint(tri.v2, tri.v3);
                    int c = getMiddlePoint(tri.v3, tri.v1);

                    faces2.Add(new TriangleIndices(tri.v1, a, c));
                    faces2.Add(new TriangleIndices(tri.v2, b, a));
                    faces2.Add(new TriangleIndices(tri.v3, c, b));
                    faces2.Add(new TriangleIndices(a, b, c));
                }
                faces = faces2;
            }

            // done, now add triangles to mesh
            foreach (var tri in faces)
            {
                this.geometry.TriangleIndices.Add(tri.v1);
                this.geometry.TriangleIndices.Add(tri.v2);
                this.geometry.TriangleIndices.Add(tri.v3);
            }

            return this.geometry;
        }


        private struct TriangleIndices
        {
            public int v1;
            public int v2;
            public int v3;

            public TriangleIndices(int v1, int v2, int v3)
            {
                this.v1 = v1;
                this.v2 = v2;
                this.v3 = v3;
            }
        }

        private MeshGeometry3D geometry;
        private int index;
        private Dictionary<Int64, int> middlePointIndexCache;

        // add vertex to mesh, fix position to be on unit sphere, return index
        private int addVertex(Point3D p)
        {
            double length = Math.Sqrt(p.X * p.X + p.Y * p.Y + p.Z * p.Z);
            geometry.Positions.Add(new Point3D(p.X / length, p.Y / length, p.Z / length));
            return index++;
        }

        // return index of point in the middle of p1 and p2
        private int getMiddlePoint(int p1, int p2)
        {
            // first check if we have it already
            bool firstIsSmaller = p1 < p2;
            Int64 smallerIndex = firstIsSmaller ? p1 : p2;
            Int64 greaterIndex = firstIsSmaller ? p2 : p1;
            Int64 key = (smallerIndex << 32) + greaterIndex;

            int ret;
            if (this.middlePointIndexCache.TryGetValue(key, out ret))
            {
                return ret;
            }

            // not in cache, calculate it
            Point3D point1 = this.geometry.Positions[p1];
            Point3D point2 = this.geometry.Positions[p2];
            Point3D middle = new Point3D(
                (point1.X + point2.X) / 2.0,
                (point1.Y + point2.Y) / 2.0,
                (point1.Z + point2.Z) / 2.0);

            // add vertex makes sure point is on unit sphere
            int i = addVertex(middle);

            // store it, return index
            this.middlePointIndexCache.Add(key, i);
            return i;
        }

    }
}
