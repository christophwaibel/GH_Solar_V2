using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Media.Media3D;
using System.Windows;

/*
 * SkyDome.cs
 * Copyright 2017 Christoph Waibel <chwaibel@student.ethz.ch>
 * 
 * This work is licensed under the GNU GPL license version 3.
*/

namespace SolarModel
{
    /// <summary>
    /// Constructs a skydome (hemisphere). Used in Sensorpoints.cs.
    /// <para>Each sensor point will have one skydome.</para>
    /// </summary>
    public class SkyDome
    {
        /// <summary>
        /// Icosphere, constructed using the IcoSphere class.
        /// </summary>
        private IcoSphere ico;
        /// <summary>
        /// Faces of the hemisphere. Indices, referencing to VertexCoordinatesSphere.
        /// </summary>
        public List<int[]> Faces { get; private set; } 
        /// <summary>
        /// Surface area of each face. Most of them are identical, but since some had to be split when turning the sphere to a hemisphere, they are not all the same
        /// </summary>
        public List<double> FaceAreas { get; private set; }    
        /// <summary>
        /// All vertex-vectors (points) of the complete sphere.
        /// </summary>
        public List<double[]> VertexVectorsSphere { get; private set; }
        /// <summary>
        /// Indicating for each vertex of the hemisphere, if it is obstructed (true).
        /// </summary>
        public double[] VertexShadowSphere { get; internal set; }
        /// <summary>
        /// Vertices of the hemisphere. Referencing to VertexCoordinatesSphere.
        /// </summary>
        public List<int> VerticesHemisphere { get; private set; }
        /// <summary>
        /// Vertices lying on z=0, i.e. the Horizon. Referencing to VertexCoordinatesSphere
        /// </summary>
        public List<int> VerticesHorizon { get; private set; }
        /// <summary>
        /// Unfortunately, the horizon vertices are not spaced equally. so this is the contribution of each vertex of the entire circle. In degree. Should sum up to 720 (2 circles).
        /// </summary>
        public List<double> HorizonSegments { get; private set; }

        /// <summary>
        /// Weighted fraction of the horizon, which is obstructed. 1 = full obstruction.
        /// </summary>
        public double[] ShdwHorizon { get; private set; }
        /// <summary>
        /// Weighted fraction of the hemisphere, which is obstructed. 1 = full obstruction.
        /// </summary>
        public double[] ShdwDome { get; private set; }
        /// <summary>
        /// 8760 list (for each hour of the year) of booleans, indicating if the sun vector is obstructed or not (1=obstructed).
        /// </summary>
        public double[] ShdwBeam { get; private set; }




        /// <summary>
        /// Creates a sky dome (hemisphere) as a halfed icosahedron. 
        /// </summary>
        /// <param name="resolution">Resolution level of the sky dome. 0: 12 faces, 10 vertices; 1: 44 faces, 29 vertices; 2: 168 faces, 97 vertices; 3: 656 faces, 353 vertices. 1 or 2 recommended.</param>
        public SkyDome(int resolution)
        {
            ico = new IcoSphere(resolution);
            Faces = ico.getFaces();
            VertexVectorsSphere = ico.getVertexCoordinates();

            VerticesHemisphere = new List<int>();
            VerticesHorizon = new List<int>();
            HorizonSegments = new List<double>();
            FaceAreas = new List<double>();

            ShdwHorizon = new double[8760];
            ShdwDome = new double[8760];
            ShdwBeam = new double[8760];

            CalcHemisphere();
            CalcHorizonSegmentWeights();

            VertexShadowSphere = new double[VertexVectorsSphere.Count];
            //create a list of size of the facaes of the dome. use this list for shadow factors...
        }

        /// <summary>
        /// Copy object.
        /// </summary>
        /// <param name="copy"></param>
        public SkyDome(SkyDome copy)
        {
            ico = copy.ico;
            Faces = copy.Faces;
            FaceAreas = copy.FaceAreas;
            ShdwHorizon = new double[8760];
            ShdwDome = new double[8760];
            ShdwBeam = new double[8760];
            VertexVectorsSphere = copy.VertexVectorsSphere;
            VerticesHemisphere = copy.VerticesHemisphere;
            VerticesHorizon = copy.VerticesHorizon;
            HorizonSegments = copy.HorizonSegments;
            VertexShadowSphere = new double[VertexVectorsSphere.Count];
            //ShdwHorizon, ShdwDome, ShdwSunVector must be re-evaluated for new sensor point
        }

        /// <summary>
        /// Replace the vertex vectors of the hemisphere.
        /// <para>E.g. for externally rotating the hemisphere.</para>
        /// </summary>
        /// <param name="NewVectors">New vectors.</param>
        public void ReplaceVertexVectors(double[][] NewVectors)
        {
            for (int i = 0; i < this.VertexVectorsSphere.Count; i++)
            {
                for (int u = 0; u < this.VertexVectorsSphere[i].Length; u++)
                {
                    double copy =  NewVectors[i][u];
                    this.VertexVectorsSphere[i][u] = copy;
                }
            }
        }

        /// <summary>
        /// Rotate vertex vectors of the hemisphere with a rotation matrix.
        /// </summary>
        /// <param name="rotMatrix">Symmetric rotation matrix.</param>
        public void RotateVertexVectors(double[,] R)
        {
            if (R == null) return;
            for (int i = 0; i < this.VertexVectorsSphere.Count; i++)
            {
                double[] rotated = new double[this.VertexVectorsSphere[i].Length];
                for (int j = 0; j < this.VertexVectorsSphere[i].Length; j++)
                {
                    rotated[j] = R[j, 0] * this.VertexVectorsSphere[i][0] + R[j, 1] * this.VertexVectorsSphere[i][1] + R[j, 2] * this.VertexVectorsSphere[i][2];
                }
                for (int j = 0; j < this.VertexVectorsSphere[i].Length; j++)
                {
                    this.VertexVectorsSphere[i][j] = rotated[j];
                }
            }
        }


        /// <summary>
        /// Set the fraction of the dome, which is obstructed.
        /// <para>Needs some obstruction calculation, which you have to run in another program (e.g. Rhinoceros)</para>
        /// </summary>
        internal void SetShadow_Dome()
        {
            double totalArea = 0.0;
            double weights = 0.0;
            for(int i=0; i<this.Faces.Count;i++)        
            {
                double w = 0;
                int wl = this.Faces[i].Length;
                for (int u = 0; u < wl; u++)
                {
                    w += this.VertexShadowSphere[this.Faces[i][u]];
                }
                w = w / Convert.ToDouble(wl);
                weights += w * this.FaceAreas[i];
                totalArea += this.FaceAreas[i];
            }
            double val = weights / totalArea;
            for (int t=0; t<8760;t++)
                this.ShdwDome[t] = val;
        }


        /// <summary>
        /// Set the fraction of the horizon, which is obstructed.
        /// <para>Needs some obstruction calculation, which you have to run in another program (e.g. Rhinoceros)</para>
        /// </summary>
        internal void SetShadow_Horizon()
        {
            double weights = 0.0;
            for (int i = 0; i < this.VerticesHorizon.Count; i++)
            {
                weights += (this.VertexShadowSphere[this.VerticesHorizon[i]]) * this.HorizonSegments[i];
            }
            double val = weights / 720.0;//720 is 2 circles : the sum angle of all horizon segments
            for (int t = 0; t < 8760; t++)
                this.ShdwHorizon[t] = val;
        }


        /// <summary>
        /// Set the fraction of the dome, which is obstructed. With permeable objects which have a time-series as extinction coefficient.
        /// <para>Needs some obstruction calculation, which you have to run in another program (e.g. Rhinoceros)</para>
        /// </summary>
        /// <param name="permObj">[u] indicate for each dome vertex, if and only if it is obstructed by permeable object.</param>
        /// <param name="permObjInd">[u][p] for each dome vertex u, array of indices to permObjlen[p]</param>
        /// <param name="permObjLen">[u][p] for each element in permObjInd, it's corresponding length.</param>
        /// <param name="extCoeff">[p][t] for each perm obcject p and time step t, extinction coefficient.</param>
        internal void SetShadow_DomeAndHorizon(bool[] permObj, int [][] permObjInd, double [][] permObjLen, List<double[]> extCoeff)
        {
            double[][] extCoeffSum = new double[this.VertexVectorsSphere.Count][];  //[u][t] extinction coefficients for each dome vertex [u] and for each time step [t].
            for (int u = 0; u < this.VertexVectorsSphere.Count; u++)
            {
                extCoeffSum[u] = new double[8760];
                for (int t = 0; t < 8760; t++)
                {
                    extCoeffSum[u][t] = 0.0;
                }
            }

            bool [] permObjSphere = new bool[this.VertexVectorsSphere.Count];
            for (int u = 0; u < permObj.Length; u++)
            {
                if (permObj[u])
                {
                    permObjSphere[this.VerticesHemisphere[u]] = true;
                    for (int t = 0; t < 8760; t++)
                    {
                        for (int k = 0; k < permObjInd[u].Length; k++)
                        {
                            extCoeffSum[this.VerticesHemisphere[u]][t] += permObjLen[u][k] * extCoeff[permObjInd[u][k]][t];
                        }
                        if (extCoeffSum[this.VerticesHemisphere[u]][t] > 1.0) extCoeffSum[this.VerticesHemisphere[u]][t] = 1.0;
                    }
                }
            }


            double totalArea = 0.0;
            double[] weights = new double[8760];
            double[] weightsHorizon = new double[8760];
            for (int t = 0; t < 8760; t++)
            {
                weights[t] = 0.0;
                weightsHorizon[t] = 0.0;
            }

            for (int i = 0; i < this.Faces.Count; i++)
            {
                int wl = this.Faces[i].Length;
                for (int t = 0; t < 8760; t++)
                {
                    double w = 0.0;
                    for (int u = 0; u < wl; u++)
                    {
                        if (permObjSphere[this.Faces[i][u]])
                        {
                            w += extCoeffSum[this.Faces[i][u]][t];
                        }
                        else
                        {
                            w += this.VertexShadowSphere[this.Faces[i][u]];
                        }
                    }
                    w = w / Convert.ToDouble(wl);
                    weights[t] += w * this.FaceAreas[i];
                }
                totalArea += this.FaceAreas[i];
            }



            for (int i = 0; i < this.VerticesHorizon.Count; i++)
            {
                for (int t = 0; t < 8760; t++)
                {
                    double wh = 0.0;
                    if (permObjSphere[this.VerticesHorizon[i]])
                    {
                        wh =  extCoeffSum[this.VerticesHorizon[i]][t];
                    }
                    else
                    {
                        wh = this.VertexShadowSphere[this.VerticesHorizon[i]];
                    }
                    weightsHorizon[t] += wh * this.HorizonSegments[i];
                }
            }



            for (int t = 0; t < 8760; t++)
            {
                this.ShdwDome[t] = weights[t] / totalArea;
                this.ShdwHorizon[t] = weightsHorizon[t] / 720.0;    //720 is 2 circles : the sum angle of all horizon segments
            }
        }





        /// <summary>
        /// Indicate for a certain hour of the year, if the center of the skydome lies in shadow. I.e. has no beam radiation.
        /// </summary>
        /// <param name="HOY">Hour of year, HOY ∈ [0,8759].</param>
        /// <param name="shadow">shadow = true.</param>
        internal void SetShadow_Beam(int HOY, double shadow)
        {
            this.ShdwBeam[HOY] = shadow;
        }


        /// <summary>
        /// Takes a sphere and halfs it, so we have a hemisphere.
        /// Sets:
        /// <para>Faces: Faces of the sphere. Each face has three vertex indices.</para>
        /// <para>FaceAreas: Areas of each face of the hemisphere mesh.</para>
        /// <para>VertexCoordinatesSphere: X,Y,Z -coordinates of all vertices.</para>
        /// <para>VerticesHemisphere: Vertices of the hemisphere.</para>
        /// </summary>
        private void CalcHemisphere()
        {
            Dictionary<int, bool> delVertex = new Dictionary<int, bool>();
            for (int i = 0; i < this.VertexVectorsSphere.Count; i++) 
                delVertex.Add(i, true);

            Dictionary<int, bool> splitFace = new Dictionary<int, bool>();

            List<int> removeFace = new List<int>();
            List<int[]> newFaces = new List<int[]>();
            List<double[]> newPoints = new List<double[]>();    //new projected points. check for duplicate
            List<int> newIndices = new List<int>();
            int newpointscreated = 0;
            int StartIndOfNewPts = this.VertexVectorsSphere.Count;
            for (int i = 0; i < this.Faces.Count; i++)
            {
                splitFace.Add(i, false);
                bool deleteface = false;
                int countMinusZ = 0;
                int countZeroZ = 0;
                int vminus = 0;
                for (int u = 0; u < this.Faces[i].Length; u++)
                {
                    //conditions to delete face, and indicate which face to split
                    if (this.VertexVectorsSphere[this.Faces[i][u]][2] < 0)
                    {
                        deleteface = true;
                        vminus = u;
                        countMinusZ++;
                    }
                    else if (this.VertexVectorsSphere[this.Faces[i][u]][2] == 0)
                    {
                        countZeroZ++;
                    }
                }

                // construct a new face and add coordinates to coordinates list
                // later store new vertex into vertices list
                if (countMinusZ == 1 && countZeroZ != 2)
                {
                    splitFace[i] = true;

                    int[] v = new int[2];
                    int count = 0;
                    for (int u = 0; u < this.Faces[i].Length; u++)
                    {
                        if (u != vminus)
                        {
                            v[count] = u;
                            count++;
                        }
                    }

                    //this point will be twice...
                    double[] newP = new double[3];
                    newP[0] = this.VertexVectorsSphere[this.Faces[i][vminus]][0];
                    newP[1] = this.VertexVectorsSphere[this.Faces[i][vminus]][1];
                    newP[2] = 0;

                    int index = 0;
                    if (newpointscreated == 0)
                    {
                        newPoints.Add(newP);
                        index = StartIndOfNewPts;
                        StartIndOfNewPts++;
                        newpointscreated++;
                        newIndices.Add(index);
                    }
                    else
                    {
                        //check for duplicate in newPoints list
                        bool exists = false;
                        for (int j = 0; j < newPoints.Count; j++)
                        {
                            if (newP[0] == newPoints[j][0])
                            {
                                if (newP[1] == newPoints[j][1])
                                {
                                    index = j + this.VertexVectorsSphere.Count;
                                    exists = true;
                                    break;
                                }
                            }
                        }
                        if (!exists)
                        {
                            newPoints.Add(newP);
                            index = StartIndOfNewPts;
                            StartIndOfNewPts++;
                            newpointscreated++;
                            newIndices.Add(index);
                        }
                    }

                    //delVertex.Add(_vertcoord.Count-1, true);        // i still need to delete the vertex point, but the order is wrong.
                    newFaces.Add(new int[] { this.Faces[i][v[0]], this.Faces[i][v[1]], index });

                }

                // just a flag to delete. but actually deletes later
                if (deleteface)
                    removeFace.Add(i);
            }


            this.Faces.AddRange(newFaces);

            //remove Faces
            foreach (int index in removeFace.OrderByDescending(i => i)) 
                this.Faces.RemoveAt(index);

            //indicating unused vertices
            foreach (int[] face in this.Faces)
                foreach (int f in face)
                    delVertex[f] = false;

            //put only used vertices into the Vertex list
            for (int i = 0; i < this.VertexVectorsSphere.Count; i++)
                if (delVertex[i] == false) this.VerticesHemisphere.Add(i);
            this.VertexVectorsSphere.AddRange(newPoints);
            this.VerticesHemisphere.AddRange(newIndices);

            //calculating face areas
            foreach(int[]face in this.Faces)
            {
                this.FaceAreas.Add(Misc.getTriangleArea(this.VertexVectorsSphere[face[0]], this.VertexVectorsSphere[face[1]], this.VertexVectorsSphere[face[2]]));
            }

        }

        /// <summary>
        /// Calculate the segment weights of the horizon vertices. 
        /// <para>We need this, because unfortunately the vertices on the horizon are not equally distributed.</para>
        /// <para>Sets:</para>
        /// <para>HorizonSegments: Weight for each vector. It is the angle between its two surrounding vectors.</para>
        /// <para>VerticesHorizon: Vertices of the horizon. Index reference to the hemisphere vertex coordinates list.</para>
        /// <para>VerticesHemisphere: All vertices of the hemisphere.</para>
        /// <para>VertexCoordinatesSphere: All coordinates of the vertices of the entire sphere.</para>
        /// </summary>
        private void CalcHorizonSegmentWeights()
        {
            List<double> anglescenter = new List<double>();

            for (int i = 0; i < this.VerticesHemisphere.Count; i++)
            {
                if (this.VertexVectorsSphere[this.VerticesHemisphere[i]][2] == 0)
                {
                    this.VerticesHorizon.Add(this.VerticesHemisphere[i]); //needs to be sorted according to anglescenter
                    anglescenter.Add(Vector.AngleBetween(new Vector(0, 1), new Vector(this.VertexVectorsSphere[this.VerticesHemisphere[i]][0], this.VertexVectorsSphere[this.VerticesHemisphere[i]][1])));
                    this.HorizonSegments.Add(0);
                }
            }

            int[] items = this.VerticesHorizon.ToArray();
            double[] order = anglescenter.ToArray();
            Array.Sort(order, items);
            this.VerticesHorizon = new List<int>(items.ToList()); 

            for (int i = 0; i < this.VerticesHorizon.Count; i++)
            {
                int a;
                int b;
                if (i == 0)
                {
                    a = this.VerticesHorizon.Count - 1;
                    b = i + 1;
                }
                else if (i == this.VerticesHorizon.Count - 1)
                {
                    a = i - 1;
                    b = 0;
                }
                else
                {
                    a = i - 1;
                    b = i + 1;
                }

                Vector v1 = new Vector(this.VertexVectorsSphere[this.VerticesHorizon[a]][0], this.VertexVectorsSphere[this.VerticesHorizon[a]][1]);
                Vector v2 = new Vector(this.VertexVectorsSphere[this.VerticesHorizon[b]][0], this.VertexVectorsSphere[this.VerticesHorizon[b]][1]);
                double angle = Vector.AngleBetween(v1, v2);
                this.HorizonSegments[i] = angle;
            }
        }

    }
}
