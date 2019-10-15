using System.Collections.Generic;
using System.Threading.Tasks;
using Rhino.Geometry;

/*
 * cObstacleObject.cs
 * Copyright 2017 Christoph Waibel <chwaibel@student.ethz.ch>
 * 
 * This work is licensed under the GNU GPL license version 3.
*/

namespace GHSolar
{
    public class CObstacleObject
    {
        public Mesh mesh;
        public List<double> albedos;  //8760 values for diffuse
        public List<double> specCoeff; //8760 values for specular
        public int reflType;
        public Vector3d[] normals;
        public Vector3d[] normalsRev;
        public Point3d[] faceCen;

        public double tolerance;
        public string name;

        /// <summary>
        /// Create an obstacle object, used for solar calculations.
        /// </summary>
        /// <param name="_mesh">Mesh object. Face normals are important.</param>
        /// <param name="_albedos">8760 albedo values.</param>
        /// <param name="_specCoeff">8760 time series for specular reflection coefficient values. Value between 0 - 1.</param>
        /// <param name="_reflType">Reflection type. 0 = diffuse, 1 = specular (no refraction currently considered), 2 = specular and diffuse (snow), all other numbers = blind (no inter-reflections)</param>
        /// <param name="_tolerance">Tolerance, used to offset point from actual face center point, to avoid self obstruction.</param>
        /// <param name="_name">Name of the obstacle. E.g. use to indicate an analysis surface.</param>
        /// <param name="mt">Multi-threading.</param>
        public CObstacleObject(Mesh _mesh, List<double> _albedos, List<double> _specCoeff, int _reflType, double _tolerance, string _name, bool mt)
        {
            mesh = _mesh;
            albedos = new List<double>(_albedos);
            specCoeff = new List<double>(_specCoeff);
            reflType = _reflType;
            name = _name;
            tolerance = _tolerance;

            mesh.FaceNormals.ComputeFaceNormals();

            normals = new Vector3d[mesh.Faces.Count];
            normalsRev = new Vector3d[mesh.Faces.Count];
            faceCen = new Point3d[mesh.Faces.Count];
            if (!mt)
            {
                for (int k = 0; k < mesh.Faces.Count; k++)
                {
                    normals[k] = mesh.FaceNormals[k];
                    normalsRev[k] = Vector3d.Negate(normals[k]);
                    Point3d cen0 = mesh.Faces.GetFaceCenter(k);
                    faceCen[k] = new Point3d(Point3d.Add(cen0, Vector3d.Multiply(Vector3d.Divide(normals[k], normals[k].Length), tolerance)));
                }
            }
            else
            {
                Parallel.For(0, mesh.Faces.Count, k =>
                {
                    normals[k] = mesh.FaceNormals[k];
                    normalsRev[k] = Vector3d.Negate(normals[k]);
                    Point3d cen0 = mesh.Faces.GetFaceCenter(k);
                    faceCen[k] = new Point3d(Point3d.Add(cen0, Vector3d.Multiply(Vector3d.Divide(normals[k], normals[k].Length), tolerance)));
                });
            }


        }
    }
}
