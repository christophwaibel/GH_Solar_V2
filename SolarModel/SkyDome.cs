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
    //skycumulative values?? do i need?

    //shadow days for direct?

    //shadow fraction for skydome (isotropic compnoetn)
    //shadow fraction for horizon component
    //shadow hours for circumsolar ..... (3 days interpolation`???)




    /// <summary>
    /// Constructs a skydome with patches, each holding value for solar radiation.
    /// Each sensor point will have one SkyDome.
    /// </summary>
    public class SkyDome
    {
        //public sPatches patches; //maybe I don't need this.

        /// <summary>
        /// Icosphere, constructed using the IcoSphere class.
        /// </summary>
        private IcoSphere ico;


        ///// <summary>
        ///// A list of sun vectors for each hour of the year. 8760 in total. Using the SunVector class.
        ///// </summary>
        //public List<SunVector> SunVectors;
        
        ///// <summary>
        ///// Boolean to indicate wether an hour of the year is day-time (true = there is sunshine) or night-time (false = no sunshine).
        ///// </summary>
        //public List<bool> Sunshine;

        /// <summary>
        /// Faces of the hemisphere. Indices, referencing to VertexCoordinatesSphere.
        /// </summary>
        public List<int[]> Faces;
        /// <summary>
        /// Surface area of each face. Most of them are identical, but since some had to be split when turning the sphere to a hemisphere, they are not all the same
        /// </summary>
        public List<double> FaceAreas;    
        /// <summary>
        /// All points (coordinates) of the complete sphere.
        /// </summary>
        public List<double[]> VertexCoordinatesSphere;

        public bool [] VertexShadowSphere;

        /// <summary>
        /// Vertices of the hemisphere. Referencing to VertexCoordinatesSphere.
        /// </summary>
        public List<int> VerticesHemisphere;
        /// <summary>
        /// Vertices lying on z=0, i.e. the Horizon. Referencing to VertexCoordinatesSphere
        /// </summary>
        public List<int> VerticesHorizon;
        /// <summary>
        /// Unfortunately, the horizon vertices are not spaced equally. so this is the contribution of each vertex of the entire circle. In degree. Should sum up to 720 (2 circles).
        /// </summary>
        public List<double> HorizonSegments;

        /// <summary>
        /// Weighted fraction of the horizon, which is obstructed. 1 = full obstruction.
        /// </summary>
        public double ShdwHorizon { get; private set; }
        /// <summary>
        /// Weighted fraction of the hemisphere, which is obstructed. 1 = full obstruction.
        /// </summary>
        public double ShdwDome { get; private set; }
        /// <summary>
        /// 8760 list (for each hour of the year) of booleans, indicating if the sun vector is obstructed or not (1=obstructed).
        /// </summary>
        public bool[] ShdwBeam { get; private set; }




        /// <summary>
        /// Creates a sky dome (hemisphere) as a halfed icosahedron. 
        /// </summary>
        /// <param name="resolution">Resolution level of the sky dome. 0: 12 faces, 10 vertices; 1: 44 faces, 29 vertices; 2: 168 faces, 97 vertices; 3: 656 faces, 353 vertices. 1 or 2 recommended.</param>
        public SkyDome(int resolution)
        {
            ico = new IcoSphere(resolution);
            Faces = ico.getFaces();
            VertexCoordinatesSphere = ico.getVertexCoordinates();

            VerticesHemisphere = new List<int>();
            VerticesHorizon = new List<int>();
            HorizonSegments = new List<double>();
            FaceAreas = new List<double>();

            ShdwBeam = new bool[8760];

            CalcHalfSphere(ref Faces, ref VertexCoordinatesSphere, ref VerticesHemisphere, ref FaceAreas);
            CalcHorizonSegmentWeights(ref HorizonSegments, ref VerticesHorizon, VerticesHemisphere, VertexCoordinatesSphere);

            VertexShadowSphere = new bool[VertexCoordinatesSphere.Count];
            //create a list of size of the facaes of the dome. use this list for shadow factors...

        }

        public SkyDome(SkyDome copy)
        {
            ico = copy.ico;
            Faces = copy.Faces;
            FaceAreas = copy.FaceAreas;
            ShdwBeam = new bool[8760];
            VertexCoordinatesSphere = copy.VertexCoordinatesSphere;
            VerticesHemisphere = copy.VerticesHemisphere;
            VerticesHorizon = copy.VerticesHorizon;
            HorizonSegments = copy.HorizonSegments;
            VertexShadowSphere = new bool[VertexCoordinatesSphere.Count];
            //ShdwHorizon, ShdwDome, ShdwSunVector must be re-evaluated for new sensor point
        }




        /// <summary>
        /// Set the fraction of the dome, which is obstructed.
        /// <para>Needs some obstruction calculation, which you have to run in another program (e.g. Rhinoceros)</para>
        /// </summary>
        public void SetShadow_Dome()
        {
            double totalArea = 0.0;
            double weights = 0.0;
            for(int i=0; i<this.Faces.Count;i++)        
            {
                double w = 0;
                int wl = this.Faces[i].Length;
                for (int u = 0; u < wl; u++)
                {
                    w += Convert.ToDouble(this.VertexShadowSphere[this.Faces[i][u]]);
                }
                w = w / Convert.ToDouble(wl);
                weights += w * this.FaceAreas[i];
                totalArea += this.FaceAreas[i];
            }
            this.ShdwDome = weights / totalArea;
        }

        /// <summary>
        /// Set the fraction of the horizon, which is obstructed.
        /// <para>Needs some obstruction calculation, which you have to run in another program (e.g. Rhinoceros)</para>
        /// </summary>
        /// <param name="obstructedHorizonVectors">Indicates, which segment of the horizon is obstructed (true) and which is not (false).</param>
        public void SetShadow_Horizon()
        {
            double weights = 0.0;
            for (int i = 0; i < this.VerticesHorizon.Count; i++)
            {
                weights += (1 - Convert.ToInt32(this.VertexShadowSphere[this.VerticesHorizon[i]])) * this.HorizonSegments[i];
            }
            this.ShdwHorizon = Convert.ToDouble(weights) / 720.0;     //720 is 2 circles : the sum angle of all horizon segments
        }




        public void SetShadow_Beam(int HOY, bool shadow)
        {
            this.ShdwBeam[HOY] = shadow;
        }


        ///// <summary>
        ///// Input has 8760 values. Just copy.
        ///// </summary>
        ///// <param name="obstructedSunVectors"></param>
        //public void SetShadow_Beam(List<bool> obstructedSunVectors)
        //{
        //    //this.ShdwBeam = new List<bool>(obstructedSunVectors);
        //}




        ////////////////////////////////////////////////////////////////////
        //// interpolation here!!!!!!!!!?????
        ////!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        ///// <summary>
        ///// list of list of bools. each list is one day (24) of values. the provided days will be used to fill the entire year.
        ///// </summary>
        ///// <param name="obstructedsunvectors">keys: days of the year [1,365]. items: list of length 24, for each hour of the day, true=obstructed, false=no obstruction.</param>
        //public void setshadow_beam(dictionary<int, list<bool>> obstructedsunvectors)
        //{

        //    // do some interpolation

        //}











        /// <summary>
        /// Takes a sphere and halfs it, so we have a hemisphere.
        /// </summary>
        /// <param name="_faces">Faces of the sphere. Each face has three vertex indices.</param>
        /// <param name="_vertcoord">X,Y,Z -coordinates of all vertices.</param>
        /// <param name="_faceAreas">Areas of each face of the hemisphere mesh.</param>
        /// <param name="_vertHemi">Vertices of the hemisphere.</param>
        private void CalcHalfSphere(ref List<int[]> _faces, ref List<double[]> _vertcoord, ref List<int> _vertHemi, ref List<double> _faceAreas)
        {
            Dictionary<int, bool> delVertex = new Dictionary<int, bool>();
            for (int i = 0; i < _vertcoord.Count; i++) 
                delVertex.Add(i, true);

            Dictionary<int, bool> splitFace = new Dictionary<int, bool>();

            List<int> removeFace = new List<int>();
            List<int[]> newFaces = new List<int[]>();
            List<double[]> newPoints = new List<double[]>();    //new projected points. check for duplicate
            List<int> newIndices = new List<int>();
            int newpointscreated = 0;
            int StartIndOfNewPts = _vertcoord.Count;
            for(int i=0; i<_faces.Count; i++)
            {
                splitFace.Add(i, false);
                bool deleteface = false;
                int countMinusZ = 0;
                int countZeroZ = 0;
                int vminus = 0;
                for (int u = 0; u < _faces[i].Length; u++)
                {
                    //conditions to delete face, and indicate which face to split
                    if (_vertcoord[_faces[i][u]][2] < 0)
                    {
                        deleteface = true;
                        vminus = u;
                        countMinusZ++;
                    }
                    else if (_vertcoord[_faces[i][u]][2] == 0)
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
                    for (int u = 0; u < _faces[i].Length; u++)
                    {
                        if (u != vminus)
                        {
                            v[count] = u;
                            count++;
                        }
                    }

                    //this point will be twice...
                    double[] newP = new double[3];
                    newP[0] = _vertcoord[_faces[i][vminus]][0];
                    newP[1] = _vertcoord[_faces[i][vminus]][1];
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
                        for (int j=0; j < newPoints.Count; j++)
                        {
                            if (newP[0] == newPoints[j][0])
                            {
                                if (newP[1] == newPoints[j][1])
                                {
                                    index = j + _vertcoord.Count;
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
                    newFaces.Add(new int[] { _faces[i][v[0]], _faces[i][v[1]], index });

                }

                // just a flag to delete. but actually deletes later
                if (deleteface) 
                    removeFace.Add(i);
            }


            _faces.AddRange(newFaces);

            //remove Faces
            foreach (int index in removeFace.OrderByDescending(i => i)) 
                _faces.RemoveAt(index);

            //indicating unused vertices
            foreach (int[] face in _faces)
                foreach (int f in face)
                    delVertex[f] = false;

            //put only used vertices into the Vertex list
            for (int i = 0; i < _vertcoord.Count; i++)
                if (delVertex[i] == false) _vertHemi.Add(i);
            _vertcoord.AddRange(newPoints);
            _vertHemi.AddRange(newIndices);

            //calculating face areas
            foreach(int[]face in _faces)
            {
                _faceAreas.Add(Misc.getTriangleArea(_vertcoord[face[0]], _vertcoord[face[1]], _vertcoord[face[2]]));
            }

        }

        /// <summary>
        /// Calculate the segment weights of the horizon vertices. 
        /// <para>We need this, because unfortunately the vertices on the horizon are not equally distributed.</para>
        /// </summary>
        /// <param name="_weights">Weight for each vector. It is the angle between its two surrounding vectors.</param>
        /// <param name="_verthorizon">Vertices of the horizon. Index reference to the hemisphere vertex coordinates list.</param>
        /// <param name="_vertHemisphere">All vertices of the hemisphere.</param>
        /// <param name="_allcoords">All coordinates of the vertices of the entire sphere.</param>
        private void CalcHorizonSegmentWeights(ref List<double> _weights, ref List<int> _verthorizon, List<int> _vertHemisphere, List<double[]> _allcoords)
        {
            List<double> anglescenter = new List<double>();

            for (int i = 0; i < _vertHemisphere.Count; i++)
            {
                if (_allcoords[_vertHemisphere[i]][2] == 0)
                {
                    _verthorizon.Add(_vertHemisphere[i]); //needs to be sorted according to anglescenter
                    anglescenter.Add(Vector.AngleBetween(new Vector(0, 1), new Vector(_allcoords[_vertHemisphere[i]][0], _allcoords[_vertHemisphere[i]][1])));
                    _weights.Add(0);
                }
            }

            int[] items = _verthorizon.ToArray();
            double[] order = anglescenter.ToArray();
            Array.Sort(order, items);
            _verthorizon = new List<int>(items.ToList()); 

            for (int i = 0; i < _verthorizon.Count; i++)
            {
                int a;
                int b;
                if (i == 0)
                {
                    a = _verthorizon.Count - 1;
                    b = i + 1;
                }
                else if (i == _verthorizon.Count - 1)
                {
                    a = i - 1;
                    b = 0;
                }
                else
                {
                    a = i - 1;
                    b = i + 1;
                }

                Vector v1 = new Vector(_allcoords[_verthorizon[a]][0], _allcoords[_verthorizon[a]][1]);
                Vector v2 = new Vector(_allcoords[_verthorizon[b]][0], _allcoords[_verthorizon[b]][1]);
                double angle = Vector.AngleBetween(v1, v2);
                _weights[i] = angle;
            }
        }





        //public struct sPatches
        //{
        //    public int count { get; private set; }                       //number of patches
        //    public List<sVector> vertices { get; private set; }
        //    public List<double> obstructed { get; private set; }         //overshadowed 0-1 , 1 = no sun

        //}


        //public struct sVector
        //{
        //    public double x { get; private set; }
        //    public double y { get; private set; }
        //    public double z { get; private set; }  
        //}



    }
}
