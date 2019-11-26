using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using Rhino.Geometry;
using SolarModel;

/*
 * cShadow.cs
 * Copyright 2017 Christoph Waibel <chwaibel@student.ethz.ch>
 * 
 * This work is licensed under the GNU GPL license version 3.
*/

namespace GHSolar
{
    /// <summary>
    /// Shadow and obstruction calculations class, using the Rhino geometry library.
    /// </summary>
    public static class CShadow
    {
        const double rad2deg = 180.0 / Math.PI;
        const double deg2rad = Math.PI / 180.0;


        /// <summary>
        /// Calculates obstruction of vectors to a sensor point. Multi-threading version.
        /// </summary>
        /// <param name="SP">Sensor point.</param>
        /// <param name="SPnormal">Sensor point normal.</param>
        /// <param name="offset">Offset vectors from sensor point to avoid self-obstruction.</param>
        /// <param name="vec">Vectors to be checked for obstruction to sensor point.</param>
        /// <param name="obstacles">Obstacles.</param>
        /// <param name="shadow">Indicates for each input vector, if it is obstructed (true) or if it reaches the sensor point (false).</param>
        public static void CalcShadowMT(Point3d SP, Vector3d SPnormal, double offset, Vector3d[] vec, Mesh[] obstacles, ref bool[] shadow,
            ParallelOptions paropts)
        {
            shadow = new bool[vec.Length];        //by default all elements false. true means it's obstructed
            bool[] shdw_mt = new bool[vec.Length];
            Point3d origOffset = new Point3d(Point3d.Add(SP, Vector3d.Multiply(Vector3d.Divide(SPnormal, SPnormal.Length), offset)));
            Parallel.For(0, vec.Length, paropts, t =>
            {
                if (Math.Round(Vector3d.VectorAngle(SPnormal, vec[t]) * rad2deg, 0) > 90)  //assumes a surface. a globe in space could of course get a ray from "behind" 
                {
                    shdw_mt[t] = true;
                }
                else
                {
                    for (int u = 0; u < obstacles.Length; u++)
                    {
                        Ray3d ray = new Ray3d(origOffset, vec[t]);
                        double inters = Rhino.Geometry.Intersect.Intersection.MeshRay(obstacles[u], ray);
                        if (inters >= 0)
                        {
                            shdw_mt[t] = true;
                            break;
                        }
                    }
                }
            });
            shdw_mt.CopyTo(shadow, 0);
        }


        /// <summary>
        /// Calculates obstruction of vectors to a sensor point. Multi-threading version.
        /// </summary>
        /// <param name="SP">Sensor point.</param>
        /// <param name="SPnormal">Sensor point normal.</param>
        /// <param name="offset">Offset vectors from sensor point to avoid self-obstruction.</param>
        /// <param name="solarVec">Solar vectors to be checked for obstruction to sensor point.</param>
        /// <param name="sunshine">Indicates, if a vector is during daytime (true) or not (false).</param>
        /// <param name="obstacles">Obstacles.</param>
        /// <param name="shadow">Indicates for each input vector, if it is obstructed (true) or if it reaches the sensor point (false).</param>
        public static void CalcShadowMT(Point3d SP, Vector3d SPnormal, double offset, Vector3d[] vec, bool[] sunshine, Mesh[] obstacles, ref bool[] shadow,
            ParallelOptions paropts)
        {
            shadow = new bool[vec.Length];        //by default all elements false. true means it's obstructed
            bool[] shdw_mt = new bool[vec.Length];
            Point3d origOffset = new Point3d(Point3d.Add(SP, Vector3d.Multiply(Vector3d.Divide(SPnormal, SPnormal.Length), offset)));
            Parallel.For(0, vec.Length, paropts, t =>
            {
                if (sunshine[t])
                {
                    if (Math.Round(Vector3d.VectorAngle(SPnormal, vec[t]) * rad2deg, 0) > 90)  //assumes a surface. a globe in space could of course get a ray from "behind" 
                    {
                        shdw_mt[t] = true;
                    }
                    else
                    {
                        for (int u = 0; u < obstacles.Length; u++)
                        {
                            Ray3d ray = new Ray3d(origOffset, vec[t]);
                            double inters = Rhino.Geometry.Intersect.Intersection.MeshRay(obstacles[u], ray);
                            if (inters >= 0)
                            {
                                shdw_mt[t] = true;
                                break;
                            }
                        }
                    }
                }
                else
                {
                    shdw_mt[t] = true;
                }
            });

            shdw_mt.CopyTo(shadow, 0);
        }


        /// <summary>
        /// Calculates obstruction of vectors via semi-permeable objects to a sensor point. 
        /// </summary>
        /// <param name="SP">Sensor point.</param>
        /// <param name="SPnormal">Sensor point normal.</param>
        /// <param name="offset">Offset vectors from sensor point to avoid self-obstruction.</param>
        /// <param name="vec">Vectors to be checked for obstruction to sensor point.</param>
        /// <param name="obstacles">Semi-permeable obstacle objects.</param>
        /// <param name="HOY">Hour of the year ∈ [0, 8759].</param>
        /// <param name="shadow">Indicates for each input vector, how much it permeates through the obstacles. 1.0 : no obstruction, 0.0 : full obstruction.</param>
        public static void CalcPermBeamMT(Point3d SP, Vector3d SPnormal, double offset, Vector3d[] vec, List<CPermObject> obstacles, int HOY, 
            ParallelOptions paropts, ref double[] shadow)
        {
            Point3d origOffset = new Point3d(Point3d.Add(SP, Vector3d.Multiply(Vector3d.Divide(SPnormal, SPnormal.Length), offset)));

            double[] shadowtmp = new double[shadow.Length];
            shadow.CopyTo(shadowtmp, 0);

            Parallel.For(0, shadow.Length, paropts, t =>
            {
                if (shadowtmp[t] != 1.0)    // only look at previously unobstructed rays
                {
                    double permeates = 1.0;
                    for (int u = 0; u < obstacles.Count; u++)
                    {
                        if (obstacles[u].permeability[HOY] <= 0.0) continue;

                        Ray3d ray = new Ray3d(origOffset, vec[t]);
                        double inters = Rhino.Geometry.Intersect.Intersection.MeshRay(obstacles[u].mesh, ray);
                        if (inters >= 0)
                        {
                            Point3d pFarAway = Point3d.Add(origOffset, Vector3d.Multiply(999999999, vec[t]));
                            Ray3d rayOpposite = new Ray3d(pFarAway, Vector3d.Negate(vec[t]));
                            double inters2 = Rhino.Geometry.Intersect.Intersection.MeshRay(obstacles[u].mesh, rayOpposite);
                            if (inters2 >= 0)
                            {
                                Point3d p1 = ray.PointAt(inters);
                                Point3d p2 = rayOpposite.PointAt(inters2);
                                double length = Misc.Distance2Pts(new double[3] { p1.X, p1.Y, p1.Z }, new double[3] { p2.X, p2.Y, p2.Z });   //units must be in meters
                                if (length > 0.0)
                                {
                                    permeates -= length * (1 - obstacles[u].permeability[HOY]);
                                }
                                if (permeates <= 0.0)
                                {
                                    permeates = 0.0;
                                    break;
                                }
                            }
                        }
                    }
                    shadowtmp[t] = 1.0 - permeates;
                }
            });

            shadowtmp.CopyTo(shadow, 0);
        }


        /// <summary>
        /// Calculates obstruction of vectors via semi-permeable objects to a sensor point. Multi-threading version.
        /// </summary>
        /// <param name="SP">Sensor point.</param>
        /// <param name="SPnormal">Sensor point normal.</param>
        /// <param name="offset">Offset vectors from sensor point to avoid self-obstruction.</param>
        /// <param name="vec">Vectors to be checked for obstruction to sensor point.</param>
        /// <param name="sunshine">Indicates, if a vector is during daytime (true) or not (false).</param>
        /// <param name="permObst">Semi-permeable obstacle objects.</param>
        /// <param name="shadow">Indicates for each input vector, how much it permeates through the obstacles. 1.0 : no obstruction, 0.0 : full obstruction.</param>
        /// <param name="solidObstruction">[t] Indicates for each input vector, if it is obstructed by a solid object. Hence, in the annual interpolation, the intersection remains a bool, but only if for all hours simulated, there is solid obstruction and no permeable obstruction.</param>
        ///<param name="permObstruction">[t] indicates for each input vector, if it is obstructed by a permeable object ONLY.</param>
        ///<param name="permInd">[t][u] for each vector, an array of indices to the permeable obstacles.</param>
        ///<param name="permLength">[t][u] for each vector, an array of length values of object penetration.</param>
        public static void CalcPermMT(Point3d SP, Vector3d SPnormal, double offset, Vector3d[] vec, bool[] sunshine, List<CPermObject> permObst,
            bool[] shadow, ParallelOptions paropts,
            out bool[] permObstruction, out int[][] permInd, out double[][] permLength)
        {
            //out bool[] solidObstruction,
            Point3d origOffset = new Point3d(Point3d.Add(SP, Vector3d.Multiply(Vector3d.Divide(SPnormal, SPnormal.Length), offset)));

            //solidObstruction = new bool[vec.Length];
            permObstruction = new bool[vec.Length];
            permInd = new int[vec.Length][];
            permLength = new double[vec.Length][];

            bool[] permObstructionPar = new bool[vec.Length];
            int[][] permIndPar = new int[vec.Length][];
            double[][] permLengthPar = new double[vec.Length][];

            Parallel.For(0, vec.Length, paropts, t =>
            {
                if (sunshine[t] == false || shadow[t] == true) return;

                List<int> permIndList = new List<int>();
                List<double> permLenList = new List<double>();

                //if (shadow[t] == true)
                //{
                //    //solidObstruction[t] = true;
                //    continue; //only look at previously unobstructed rays
                //}

                //double permeates = 1.0;
                for (int u = 0; u < permObst.Count; u++)
                {
                    Ray3d ray = new Ray3d(origOffset, vec[t]);
                    double inters = Rhino.Geometry.Intersect.Intersection.MeshRay(permObst[u].mesh, ray);
                    if (inters >= 0)
                    {
                        Point3d pFarAway = Point3d.Add(origOffset, Vector3d.Multiply(999999999, vec[t]));
                        Ray3d rayOpposite = new Ray3d(pFarAway, Vector3d.Negate(vec[t]));
                        double inters2 = Rhino.Geometry.Intersect.Intersection.MeshRay(permObst[u].mesh, rayOpposite);
                        if (inters2 >= 0)
                        {
                            Point3d p1 = ray.PointAt(inters);
                            Point3d p2 = rayOpposite.PointAt(inters2);
                            double length = Misc.Distance2Pts(new double[3] { p1.X, p1.Y, p1.Z }, new double[3] { p2.X, p2.Y, p2.Z });   //units must be in meters
                            if (length > 0.0)
                            {
                                permObstructionPar[t] = true;
                                permIndList.Add(u);
                                permLenList.Add(length);
                            }
                        }
                    }
                }
                if (permObstructionPar[t] == true)
                {
                    permIndPar[t] = permIndList.ToArray();
                    permLengthPar[t] = permLenList.ToArray();
                }
            });

            for (int t = 0; t < vec.Length; t++)
            {
                permObstruction[t] = permObstructionPar[t];
                permInd[t] = permIndPar[t];
                permLength[t] = permLengthPar[t];
            }
        }


        /// <summary>
        /// Calculates specular interreflections on all sensor points for an array of solar vectors. 
        /// Irradiation values given normal to reflected ray. Multi-threading version.
        /// </summary>
        /// <remarks>
        /// Approach 3: Compute rays and reflected rays for each obstacle. Translate them to each sensor point. 
        /// Check if translated rays reach obstacle and are unobstructed. 1st and 2nd order (1 or 2 bounces).
        /// </remarks>
        /// <param name="SP">Sensor points [i] of the analysis mesh.</param>
        /// <param name="SPnormal">Normal vectors [i] of the i sensor points.</param>
        /// <param name="solarvec">Solar vectors for each time step t.</param>
        /// <param name="sunshine">Indicating sunshine for each respective solar vector.</param>
        /// <param name="HOY">Hour of the year ∈ [0, 8759].</param>
        /// <param name="_obstacles">Obstacle objects.</param>
        /// <param name="permeables">Permeable obstacle objects. Are treated as solids.</param>
        /// <param name="bounces">Number of bounces. Max 2 recommended.</param>
        /// <param name="Ispecular">Normal irradiation values [i][t][m] for each sensor point i, each solar vector t and each reflected ray m.</param>
        /// <param name="Inormals">Normal vectors [i][t][m] for each sensor point i, each solar vector t and each reflected ray m.</param>
        public static void CalcSpecularNormal3MT(Point3d[] SP, Vector3d[] SPnormal,
            Vector3d[] solarvec, bool[] sunshine, int HOY,
            List<CObstacleObject> obstacles, List<CPermObject> permeables, int bounces,
            ref double[][][] Ispecular, ref Vector3d[][][] Inormals)
        {
            //Rhino.RhinoDoc doc = Rhino.RhinoDoc.ActiveDoc;

            // Approach 3: Compute rays and reflected rays for each obstacle. Translate them to each sensor point. Check if translated rays reach obstacle and are unobstructed. 1st and 2nd order (1 or 2 bounces).
            if (bounces < 1) return;

            //add permeables to normal obstacles
            List<CObstacleObject> _obstacles = new List<CObstacleObject>();
            _obstacles = obstacles;

            List<double> _alb = new List<double>();
            List<double> _spec = new List<double>();
            for (int t = 0; t < 8760; t++)
            {
                _alb.Add(0);
                _spec.Add(0);
            }
            double _tol = (obstacles.Count > 0) ? obstacles[0].tolerance : 0.01;
            foreach (CPermObject perm in permeables)
            {
                CObstacleObject obst = new CObstacleObject(perm.mesh, _alb, _spec, 3, _tol, "perm", true);
                _obstacles.Add(obst);
            }

            List<List<List<double>>> IspecList = new List<List<List<double>>>();
            List<List<List<Vector3d>>> InormList = new List<List<List<Vector3d>>>();
            Point3d[] SPoffset = new Point3d[SP.Length];
            for (int i = 0; i < SP.Length; i++)
            {
                SPoffset[i] = new Point3d(Point3d.Add(SP[i], Vector3d.Multiply(Vector3d.Divide(SPnormal[i], SPnormal[i].Length), _obstacles[0].tolerance)));

                Ispecular[i] = new double[solarvec.Length][];
                Inormals[i] = new Vector3d[solarvec.Length][];

                IspecList.Add(new List<List<double>>());
                InormList.Add(new List<List<Vector3d>>());
                for (int t = 0; t < solarvec.Length; t++)
                {
                    IspecList[i].Add(new List<double>());
                    InormList[i].Add(new List<Vector3d>());
                }
            }


            Func<Vector3d, Point3d, Vector3d, bool> ObstructionCheck_Pt2Ray =
                (_vecincident, _pt, _ptNormal) =>
                {
                    double _vAngle = Vector3d.VectorAngle(_ptNormal, _vecincident) * (180.0 / Math.PI);
                    if (_vAngle >= 90) return true; //behind surface

                    Ray3d _ray = new Ray3d(_pt, _vecincident);
                    //check obstacle to sun
                    double _inters;
                    bool bln_inters = false;
                    for (int n = 0; n < _obstacles.Count; n++)
                    {
                        _inters = Rhino.Geometry.Intersect.Intersection.MeshRay(_obstacles[n].mesh, _ray);
                        if (_inters >= 0)
                        {
                            bln_inters = true;
                            break;
                        }
                    }
                    if (bln_inters) return true;  //obstructed

                    return false;
                };


            Func<Vector3d, Point3d, Vector3d, int, int, Vector3d, bool> ObstructionCheck_Pt2Face2Sun =
                (_vecincident, _pt, _ptNormal, _uObst, _kFace, _sun) =>
                {
                    double _vAngle = Vector3d.VectorAngle(_ptNormal, Vector3d.Negate(_vecincident)) * (180.0 / Math.PI);
                    if (_vAngle >= 90) return true; //behind surface (no intersection)

                    //check, if translated ray hits initial obstacle face
                    Ray3d rObstFace = new Ray3d(_pt, Vector3d.Negate(_vecincident));
                    Mesh mshface = new Mesh();
                    mshface.Vertices.Add(_obstacles[_uObst].mesh.Vertices[_obstacles[_uObst].mesh.Faces.GetFace(_kFace).A]);
                    mshface.Vertices.Add(_obstacles[_uObst].mesh.Vertices[_obstacles[_uObst].mesh.Faces.GetFace(_kFace).B]);
                    mshface.Vertices.Add(_obstacles[_uObst].mesh.Vertices[_obstacles[_uObst].mesh.Faces.GetFace(_kFace).C]);
                    if (_obstacles[_uObst].mesh.Faces.GetFace(_kFace).IsQuad)
                    {
                        mshface.Vertices.Add(_obstacles[_uObst].mesh.Vertices[_obstacles[_uObst].mesh.Faces.GetFace(_kFace).D]);
                        mshface.Faces.AddFace(0, 1, 2, 3);
                    }
                    else
                    {
                        mshface.Faces.AddFace(0, 1, 2);
                    }
                    double inters = Rhino.Geometry.Intersect.Intersection.MeshRay(mshface, rObstFace);
                    if (inters < 0) return true;

                    //check, if on the way to the face, its obstructed
                    mshface.Normals.ComputeNormals();
                    Point3d pX = rObstFace.PointAt(inters);
                    MeshPoint mshp = mshface.ClosestMeshPoint(pX, 0.0);
                    Vector3d pNormal = mshface.NormalAt(mshp);
                    Point3d pXoffset = CMisc.OffsetPt(pX, pNormal, _obstacles[_uObst].tolerance);
                    bool bln_inters = false;
                    for (int n = 0; n < _obstacles.Count; n++)
                    {
                        int[] f;
                        Point3d[] _inters = Rhino.Geometry.Intersect.Intersection.MeshLine(_obstacles[n].mesh, new Line(pXoffset, _pt), out f);
                        if (_inters != null && _inters.Length > 0)
                        {
                            bln_inters = true;
                            break;
                        }
                    }
                    if (bln_inters) return true;  //obstructed. no intersection point with face

                    //pXoffset to sun
                    if (ObstructionCheck_Pt2Ray(_sun, pXoffset, pNormal)) return true;
                    else return false;      //it reaches the sun
                };


            Func<Vector3d, Vector3d, Point3d, Vector3d, int, int, int, int, Vector3d, bool> ObstructionCheck_Pt2Face2Face2Sun =
                (_vec1storder, _vec2ndorder, _pt, _ptNormal, _u1st, _k1st, _u2nd, _k2nd, _sun) =>
                {
                    //check sp 2 2ndorder face
                    double _vAngle = Vector3d.VectorAngle(_ptNormal, Vector3d.Negate(_vec2ndorder)) * (180.0 / Math.PI);
                    if (_vAngle >= 90) return true; //behind surface (no intersection)

                    //check, if translated ray hits initial obstacle face
                    Ray3d rPt2Face = new Ray3d(_pt, Vector3d.Negate(_vec2ndorder));
                    Mesh mshface = new Mesh();
                    Point3f A = _obstacles[_u2nd].mesh.Vertices[_obstacles[_u2nd].mesh.Faces.GetFace(_k2nd).A];
                    Point3f B = _obstacles[_u2nd].mesh.Vertices[_obstacles[_u2nd].mesh.Faces.GetFace(_k2nd).B];
                    Point3f C = _obstacles[_u2nd].mesh.Vertices[_obstacles[_u2nd].mesh.Faces.GetFace(_k2nd).C];
                    mshface.Vertices.Add(A);
                    mshface.Vertices.Add(B);
                    mshface.Vertices.Add(C);
                    if (_obstacles[_u2nd].mesh.Faces.GetFace(_k2nd).IsQuad)
                    {
                        Point3f D = _obstacles[_u2nd].mesh.Vertices[_obstacles[_u2nd].mesh.Faces.GetFace(_k2nd).D];
                        mshface.Vertices.Add(D);
                        mshface.Faces.AddFace(0, 1, 2, 3);
                    }
                    else
                    {
                        mshface.Faces.AddFace(0, 1, 2);
                    }
                    double inters = Rhino.Geometry.Intersect.Intersection.MeshRay(mshface, rPt2Face);
                    if (inters < 0.0) return true;

                    //check, if on the way to the face, its obstructed
                    mshface.Normals.ComputeNormals();
                    Point3d pX = rPt2Face.PointAt(inters);
                    MeshPoint mshp = mshface.ClosestMeshPoint(pX, 0.0);
                    Vector3d pNormal = mshface.NormalAt(mshp);
                    Point3d pXoffset = CMisc.OffsetPt(pX, pNormal, _obstacles[_u2nd].tolerance);
                    bool bln_inters = false;
                    for (int n = 0; n < _obstacles.Count; n++)
                    {
                        int[] f;
                        Point3d[] _inters = Rhino.Geometry.Intersect.Intersection.MeshLine(_obstacles[n].mesh, new Line(pXoffset, _pt), out f);
                        if (_inters != null && _inters.Length > 0)
                        {
                            bln_inters = true;
                            break;
                        }
                    }
                    if (bln_inters) return true;  //obstructed. no intersection point with face


                    //check 2ndorder face to 1storder face                 //check 1storder face to sun
                    return ObstructionCheck_Pt2Face2Sun(_vec1storder, pXoffset, pNormal, _u1st, _k1st, _sun);
                };


            //List<Line> ln = new List<Line>();


            //foreach solar vector t
            for (int t = 0; t < solarvec.Length; t++)
            {
                if (!sunshine[t]) continue;

                //foreach specular object   
                for (int u = 0; u < _obstacles.Count; u++)
                {
                    if (_obstacles[u].reflType != 1 && _obstacles[u].reflType != 2) continue;

                    //foreach face
                    for (int k = 0; k < _obstacles[u].faceCen.Length; k++)
                    {
                        double vAngle_1 = Vector3d.VectorAngle(_obstacles[u].normals[k], solarvec[t]) * (180.0 / Math.PI);
                        if (vAngle_1 >= 90) continue;

                        //1st order reflection
                        Vector3d refl_1 = CMisc.ReflectVec(_obstacles[u].normals[k], Vector3d.Negate(solarvec[t]));

                        object sync = new object();
                        Parallel.For(0, SP.Length, i =>
                        {
                            if (ObstructionCheck_Pt2Face2Sun(refl_1, SPoffset[i], SPnormal[i], u, k, solarvec[t])) return;
                            lock (sync)
                            {
                                IspecList[i][t].Add(_obstacles[u].specCoeff[HOY]);
                                InormList[i][t].Add(refl_1);
                            }
                        });

                        if (bounces > 1)
                        {
                            for (int n = 0; n < _obstacles.Count; n++)
                            {
                                if (_obstacles[n].reflType != 1 && _obstacles[n].reflType != 2) continue;

                                for (int q = 0; q < _obstacles[n].faceCen.Length; q++)
                                {
                                    if (n == u && k == q) continue;

                                    double vAngle_2 = Vector3d.VectorAngle(_obstacles[n].normals[q], Vector3d.Negate(refl_1)) * (180.0 / Math.PI);
                                    if (vAngle_2 >= 90) continue;

                                    //2nd order reflection
                                    Vector3d refl_2 = CMisc.ReflectVec(_obstacles[n].normals[q], refl_1);

                                    Parallel.For(0, SP.Length, i =>
                                    {
                                        if (Vector3d.Equals(SPnormal[i], _obstacles[n].normals[q])) return;

                                        if (ObstructionCheck_Pt2Face2Face2Sun(refl_1, refl_2, SPoffset[i], SPnormal[i], u, k, n, q, solarvec[t])) return;

                                        lock (sync)
                                        {
                                            IspecList[i][t].Add(_obstacles[u].specCoeff[HOY] * _obstacles[n].specCoeff[HOY]);
                                            InormList[i][t].Add(refl_2);
                                        }
                                    });
                                }
                            }
                        }
                    }
                }
            }

            //foreach (Line l in ln)
            //{
            //    doc.Objects.AddLine(l);
            //}

            for (int i = 0; i < SP.Length; i++)
            {
                for (int t = 0; t < solarvec.Length; t++)
                {
                    Ispecular[i][t] = IspecList[i][t].ToArray();
                    Inormals[i][t] = InormList[i][t].ToArray();
                }
            }
        }


        /// <summary>
        /// Calculates specular interreflections on all sensor points for an array of solar vectors. 
        /// Reference to specular objects as well as reflected rays given.
        /// </summary>
        /// <remarks>
        /// Approach 3: Compute rays and reflected rays for each obstacle. Translate them to each sensor point. 
        /// Check if translated rays reach obstacle and are unobstructed. 1st and 2nd order (1 or 2 bounces).
        /// </remarks>
        /// <param name="SP">Sensor points [i] of the analysis mesh.</param>
        /// <param name="SPnormal">Normal vectors [i] of the i sensor points.</param>
        /// <param name="solarvec">Solar vectors for each time step t.</param>
        /// <param name="sunshine">Indicating sunshine for each respective solar vector.</param>
        /// <param name="HOY">Hour of the year ∈ [0, 8759].</param>
        /// <param name="_obstacles">Obstacle objects.</param>
        /// <param name="permeables">Permeable obstacle objects. Are treated as solids.</param>
        /// <param name="bounces">Number of bounces. Max 2 recommended.</param>
        /// <param name="IObstRef1stOrder">Normal irradiation reflection coefficient values [i][t][m] for each sensor point i, each solar vector t and each reflected ray m.</param>
        /// <param name="IObstRef2ndOrder"></param>
        /// <param name="Inormals">Normal vectors [i][t][m] for each sensor point i, each solar vector t and each reflected ray m.</param>
        public static void CalcSpecularNormal5MT(Point3d[] SP, Vector3d[] SPnormal,
            Vector3d[] solarvec, bool[] sunshine,
            List<CObstacleObject> obstacles, List<CPermObject> permeables, int bounces,
            ParallelOptions paropts,
            out int[][][] IObstRef1stOrder, out int[][][] IObstRef2ndOrder, out Vector3d[][][] Inormals)
        {
            if (bounces < 1)
            {
                IObstRef1stOrder = null;
                IObstRef2ndOrder = null;
                Inormals = null;
                return;
            }

            //add permeables to normal obstacles
            List<CObstacleObject> _obstacles = new List<CObstacleObject>();
            _obstacles = obstacles;

            List<double> _alb = new List<double>();
            List<double> _spec = new List<double>();
            for (int t = 0; t < 8760; t++)
            {
                _alb.Add(0);
                _spec.Add(0);
            }
            double _tol = (obstacles.Count > 0) ? obstacles[0].tolerance : 0.01;
            foreach (CPermObject perm in permeables)
            {
                CObstacleObject obst = new CObstacleObject(perm.mesh, _alb, _spec, 3, _tol, "perm", false);
                _obstacles.Add(obst);
            }


            //___________________________________________________________________
            /////////////////////////////////////////////////////////////////////
            ////////////////////   Obstruction check functions   ////////////////
            Func<Vector3d, Point3d, Vector3d, bool> ObstructionCheck_Pt2Ray =
                (_vecincident, _pt, _ptNormal) =>
                {
                    double _vAngle = Vector3d.VectorAngle(_ptNormal, _vecincident) * (180.0 / Math.PI);
                    if (_vAngle >= 90) return true; //behind surface

                    Ray3d _ray = new Ray3d(_pt, _vecincident);
                    //check obstacle to sun
                    double _inters;
                    bool bln_inters = false;
                    for (int n = 0; n < _obstacles.Count; n++)
                    {
                        _inters = Rhino.Geometry.Intersect.Intersection.MeshRay(_obstacles[n].mesh, _ray);
                        if (_inters >= 0)
                        {
                            bln_inters = true;
                            break;
                        }
                    }
                    if (bln_inters) return true;  //obstructed

                    return false;
                };


            Func<Vector3d, Point3d, Vector3d, int, int, Vector3d, bool> ObstructionCheck_Pt2Face2Sun =
                (_vecincident, _pt, _ptNormal, _uObst, _kFace, _sun) =>
                {
                    double _vAngle = Vector3d.VectorAngle(_ptNormal, Vector3d.Negate(_vecincident)) * (180.0 / Math.PI);
                    if (_vAngle >= 90) return true; //behind surface (no intersection)

                    //check, if translated ray hits initial obstacle face
                    Ray3d rObstFace = new Ray3d(_pt, Vector3d.Negate(_vecincident));
                    Mesh mshface = new Mesh();
                    mshface.Vertices.Add(_obstacles[_uObst].mesh.Vertices[_obstacles[_uObst].mesh.Faces.GetFace(_kFace).A]);
                    mshface.Vertices.Add(_obstacles[_uObst].mesh.Vertices[_obstacles[_uObst].mesh.Faces.GetFace(_kFace).B]);
                    mshface.Vertices.Add(_obstacles[_uObst].mesh.Vertices[_obstacles[_uObst].mesh.Faces.GetFace(_kFace).C]);
                    if (_obstacles[_uObst].mesh.Faces.GetFace(_kFace).IsQuad)
                    {
                        mshface.Vertices.Add(_obstacles[_uObst].mesh.Vertices[_obstacles[_uObst].mesh.Faces.GetFace(_kFace).D]);
                        mshface.Faces.AddFace(0, 1, 2, 3);
                    }
                    else
                    {
                        mshface.Faces.AddFace(0, 1, 2);
                    }
                    double inters = Rhino.Geometry.Intersect.Intersection.MeshRay(mshface, rObstFace);
                    if (inters < 0) return true;

                    //check, if on the way to the face, its obstructed
                    mshface.Normals.ComputeNormals();
                    Point3d pX = rObstFace.PointAt(inters);
                    MeshPoint mshp = mshface.ClosestMeshPoint(pX, 0.0);
                    Vector3d pNormal = mshface.NormalAt(mshp);
                    Point3d pXoffset = CMisc.OffsetPt(pX, pNormal, _obstacles[_uObst].tolerance);
                    bool bln_inters = false;
                    for (int n = 0; n < _obstacles.Count; n++)
                    {
                        int[] f;
                        Point3d[] _inters = Rhino.Geometry.Intersect.Intersection.MeshLine(_obstacles[n].mesh, new Line(pXoffset, _pt), out f);
                        if (_inters != null && _inters.Length > 0)
                        {
                            bln_inters = true;
                            break;
                        }
                    }
                    if (bln_inters) return true;  //obstructed. no intersection point with face

                    //pXoffset to sun
                    if (ObstructionCheck_Pt2Ray(_sun, pXoffset, pNormal)) return true;
                    else return false;      //it reaches the sun
                };


            Func<Vector3d, Vector3d, Point3d, Vector3d, int, int, int, int, Vector3d, bool> ObstructionCheck_Pt2Face2Face2Sun =
                (_vec1storder, _vec2ndorder, _pt, _ptNormal, _u1st, _k1st, _u2nd, _k2nd, _sun) =>
                {
                    //check sp 2 2ndorder face
                    double _vAngle = Vector3d.VectorAngle(_ptNormal, Vector3d.Negate(_vec2ndorder)) * (180.0 / Math.PI);
                    if (_vAngle >= 90) return true; //behind surface (no intersection)

                    //check, if translated ray hits initial obstacle face
                    Ray3d rPt2Face = new Ray3d(_pt, Vector3d.Negate(_vec2ndorder));
                    Mesh mshface = new Mesh();
                    mshface.Vertices.Add(_obstacles[_u2nd].mesh.Vertices[_obstacles[_u2nd].mesh.Faces.GetFace(_k2nd).A]);
                    mshface.Vertices.Add(_obstacles[_u2nd].mesh.Vertices[_obstacles[_u2nd].mesh.Faces.GetFace(_k2nd).B]);
                    mshface.Vertices.Add(_obstacles[_u2nd].mesh.Vertices[_obstacles[_u2nd].mesh.Faces.GetFace(_k2nd).C]);
                    if (_obstacles[_u2nd].mesh.Faces.GetFace(_k2nd).IsQuad)
                    {
                        mshface.Vertices.Add(_obstacles[_u2nd].mesh.Vertices[_obstacles[_u2nd].mesh.Faces.GetFace(_k2nd).D]);
                        mshface.Faces.AddFace(0, 1, 2, 3);
                    }
                    else
                    {
                        mshface.Faces.AddFace(0, 1, 2);
                    }
                    double inters = Rhino.Geometry.Intersect.Intersection.MeshRay(mshface, rPt2Face);
                    if (inters < 0.0) return true;

                    //check, if on the way to the face, its obstructed
                    mshface.Normals.ComputeNormals();
                    Point3d pX = rPt2Face.PointAt(inters);
                    MeshPoint mshp = mshface.ClosestMeshPoint(pX, 0.0);
                    Vector3d pNormal = mshface.NormalAt(mshp);
                    Point3d pXoffset = CMisc.OffsetPt(pX, pNormal, _obstacles[_u2nd].tolerance);
                    bool bln_inters = false;
                    for (int n = 0; n < _obstacles.Count; n++)
                    {
                        int[] f;
                        Point3d[] _inters = Rhino.Geometry.Intersect.Intersection.MeshLine(_obstacles[n].mesh, new Line(pXoffset, _pt), out f);
                        if (_inters != null && _inters.Length > 0)
                        {
                            bln_inters = true;
                            break;
                        }
                    }
                    if (bln_inters) return true;  //obstructed. no intersection point with face


                    //check 2ndorder face to 1storder face                 //check 1storder face to sun
                    return ObstructionCheck_Pt2Face2Sun(_vec1storder, pXoffset, pNormal, _u1st, _k1st, _sun);
                };
            /////////////////////////////////////////////////////////////////////
            //___________________________________________________________________


            //___________________________________________________________________
            /////////////////////////////////////////////////////////////////////
            //////////////////////////   Approach 3    //////////////////////////
            // Approach 3: Compute rays and reflected rays for each obstacle. 
            // Translate them to each sensor point. Check if translated rays reach 
            //   obstacle and are unobstructed. 1st and 2nd order (1 or 2 bounces).
            // save normals and references to specular obstacles
            IObstRef1stOrder = new int[SP.Length][][];
            IObstRef2ndOrder = new int[SP.Length][][];
            Inormals = new Vector3d[SP.Length][][];
 
            List<List<List<int>>> IObstRefList1storder = new List<List<List<int>>>();
            List<List<List<int>>> IObstRefList2ndorder = new List<List<List<int>>>();
            List<List<List<Vector3d>>> InormList = new List<List<List<Vector3d>>>();
            Point3d[] SPoffset = new Point3d[SP.Length];

            for (int i = 0; i < SP.Length; i++)
            {
                SPoffset[i] = new Point3d(Point3d.Add(SP[i], Vector3d.Multiply(Vector3d.Divide(SPnormal[i], SPnormal[i].Length), _obstacles[0].tolerance)));

                IObstRef1stOrder[i] = new int[solarvec.Length][];
                IObstRef2ndOrder[i] = new int[solarvec.Length][];
                Inormals[i] = new Vector3d[solarvec.Length][];

                IObstRefList1storder.Add(new List<List<int>>());
                IObstRefList2ndorder.Add(new List<List<int>>());
                InormList.Add(new List<List<Vector3d>>());
                for (int t=0; t < solarvec.Length; t++)
                {
                    IObstRefList1storder[i].Add(new List<int>());
                    IObstRefList2ndorder[i].Add(new List<int>());
                    InormList[i].Add(new List<Vector3d>());
                }
            }


            //foreach solar vector t
            for (int t = 0; t < solarvec.Length; t++)
            {
                if (!sunshine[t]) continue;

                //foreach specular object                 
                for (int u = 0; u < _obstacles.Count; u++)
                {
                    if (_obstacles[u].reflType != 1 && _obstacles[u].reflType != 2) continue;

                    //foreach face
                    for (int k = 0; k < _obstacles[u].faceCen.Length; k++)
                    {
                        double vAngle_1 = Vector3d.VectorAngle(_obstacles[u].normals[k], solarvec[t]) * (180.0 / Math.PI);
                        if (vAngle_1 >= 90) continue;

                        //1st order reflection
                        Vector3d refl_1 = CMisc.ReflectVec(_obstacles[u].normals[k], Vector3d.Negate(solarvec[t]));

                        Object sync = new Object();
                        Parallel.For(0, SP.Length, paropts, i =>
                        {
                            if (!ObstructionCheck_Pt2Face2Sun(refl_1, SPoffset[i], SPnormal[i], u, k, solarvec[t]))
                            {
                                lock (sync)
                                {
                                    IObstRefList1storder[i][t].Add(u);
                                    IObstRefList2ndorder[i][t].Add(-1); //dummy, indicating it has no 2nd order refl.
                                    InormList[i][t].Add(refl_1);
                                }
                            }
                        });

                        if (bounces > 1)
                        {
                            for (int n = 0; n < _obstacles.Count; n++)
                            {
                                if (_obstacles[n].reflType == 1 || _obstacles[n].reflType == 2)
                                {
                                    for (int q = 0; q < _obstacles[n].faceCen.Length; q++)
                                    {
                                        if (!(n == u && k == q))
                                        {
                                            double vAngle_2 = Vector3d.VectorAngle(_obstacles[n].normals[q], Vector3d.Negate(refl_1)) * (180.0 / Math.PI);
                                            if (vAngle_2 >= 90) continue;
                                            //2nd order reflection
                                            Vector3d refl_2 = CMisc.ReflectVec(_obstacles[n].normals[q], refl_1);
                                            Parallel.For(0, SP.Length, paropts, i =>
                                            {
                                                if (!Vector3d.Equals(SPnormal[i], _obstacles[n].normals[q]))
                                                {
                                                    if (!ObstructionCheck_Pt2Face2Face2Sun(refl_1, refl_2, SPoffset[i], SPnormal[i], u, k, n, q, solarvec[t]))
                                                    {
                                                        lock (sync)
                                                        {
                                                            IObstRefList1storder[i][t].Add(u);
                                                            IObstRefList2ndorder[i][t].Add(n);
                                                            InormList[i][t].Add(refl_2);
                                                        }
                                                    }
                                                }
                                            });
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }


            for (int i = 0; i < SP.Length; i++)
            {
                for (int t = 0; t < solarvec.Length; t++)
                {
                    IObstRef1stOrder[i][t] = IObstRefList1storder[i][t].ToArray();
                    IObstRef2ndOrder[i][t] = IObstRefList2ndorder[i][t].ToArray();
                    Inormals[i][t] = InormList[i][t].ToArray();
                }
            }
        }


        /// <summary>
        /// Calculates incident beam radiation on each SP based on interpolation of interreflected rays of three days: summer solstice, winter solstice and equinox. Multi-threading version.
        /// </summary>
        /// <param name="IObstRef1st_equ">Equinox day: References to 1st order specular object of the reflected ray, before it reaches the sensor point. [i][h][m], i=each SP, h=24 hours, m=rays hitting that SP for that timestep.</param>
        /// <param name="IObstRef2nd_equ">Equinox day: References to 2nd order specular object of the reflected ray, before it reaches the sensor point. [i][h][m], i=each SP, h=24 hours, m=rays hitting that SP for that timestep.</param>
        /// <param name="Inormals_equ">Equinox day: Reflected rays reaching a sensor point. [i][h][m], i=each SP, h=24 hours, m=rays hitting that SP for that timestep.</param>
        /// <param name="IObstRef1st_win">Winter solstice day: References to 1st order specular object of the reflected ray, before it reaches the sensor point. [i][h][m], i=each SP, h=24 hours, m=rays hitting that SP for that timestep.</param>
        /// <param name="IObstRef2nd_win">Winter solstice day: References to 2nd order specular object of the reflected ray, before it reaches the sensor point. [i][h][m], i=each SP, h=24 hours, m=rays hitting that SP for that timestep.</param>
        /// <param name="Inormals_win">Winter solstice day: Reflected rays reaching a sensor point. [i][h][m], i=each SP, h=24 hours, m=rays hitting that SP for that timestep.</param>
        /// <param name="IObstRef1st_sum">Summer solstice day: References to 1st order specular object of the reflected ray, before it reaches the sensor point. [i][h][m], i=each SP, h=24 hours, m=rays hitting that SP for that timestep.</param>
        /// <param name="IObstRef2nd_sum">Summer solstice day: References to 2nd order specular object of the reflected ray, before it reaches the sensor point. [i][h][m], i=each SP, h=24 hours, m=rays hitting that SP for that timestep.</param>
        /// <param name="Inormals_sum">Summer solstice day: Reflected rays reaching a sensor point. [i][h][m], i=each SP, h=24 hours, m=rays hitting that SP for that timestep.</param>
        /// <param name="obstacles">Obstacle objects, containing beam reflection coefficient values for each hour of the year.</param>
        /// <param name="DNI">[t] 8760 DNI values, for each hour of the year.</param>
        /// <param name="origNormal">[i] Normal vectors for each sensor point.</param>
        /// <param name="Ispecular_annual">Interreflected beam irradiation values on each sensor point, for each hour of the year. [i][t], i=each SP, t=8760 hours</param>
        public static void CalcSpecularIncident_AnnualMT(
            int[][][] IObstRef1st_equ, int[][][] IObstRef2nd_equ, Vector3d[][][] Inormals_equ,
            int[][][] IObstRef1st_win, int[][][] IObstRef2nd_win, Vector3d[][][] Inormals_win,
            int[][][] IObstRef1st_sum, int[][][] IObstRef2nd_sum, Vector3d[][][] Inormals_sum,
            List<CObstacleObject> obstacles, double[] DNI, Vector3d[] origNormal, ParallelOptions paropts,
            out double[][] Ispecular_annual)
        {
            //using interpolation of several days. 3 or x. start with three days.
            int SPCount = IObstRef1st_equ.Length;
            Ispecular_annual = new double[SPCount][];

            //    '100% equinox is 20 march and 23 Sept
            //    '100% summer is 21 june
            //    '100% winter is 22 decemebr

            int y1, y2, y3, y4, y5, y6;
            y1 = -9;    //winter solstice
            y2 = 78;    //equinox spring
            y3 = 171;   // summer solstice
            y4 = 265;   //equinox solstice
            y5 = 355;   //winter solstice
            y6 = 443;   //equinox spring

            int fullF1, fullF2;

            //double factor1 = 1.0;
            //double factor2;
            int InterpolInterv;

            for (int i = 0; i < SPCount; i++)
            {
                Ispecular_annual[i] = new double[8760];
                double [] Ispecular_annual_Par = new double[8760];

                //converting arrays which may contain nulls to 0/1 arrays:
                bool[] dayWin = new bool[24];
                bool[] daySum = new bool[24];
                bool[] dayEqu = new bool[24];
                Parallel.For(0, 24, paropts, t =>
                {
                    if (Misc.IsNullOrEmpty(Inormals_win[i][t]))
                        dayWin[t] = false;
                    else
                        dayWin[t] = true;

                    if (Misc.IsNullOrEmpty(Inormals_equ[i][t]))
                        dayEqu[t] = false;
                    else
                        dayEqu[t] = true;

                    if (Misc.IsNullOrEmpty(Inormals_sum[i][t]))
                        daySum[t] = false;
                    else
                        daySum[t] = true;
                });


                fullF1 = y1;    //100% shadowwinter on this day
                fullF2 = y2;    //100% shadowequinox on this day
                InterpolInterv = y2 + y1 * -1;
                Parallel.For(0, y2, paropts, d =>  //from 0.Jan to 19.March
                {
                    double dist1, dist2, factor1, factor2;
                    dist1 = Convert.ToDouble((InterpolInterv - Math.Abs(fullF1 - d))) / Convert.ToDouble(InterpolInterv);
                    dist1 = 1 - dist1;
                    dist1 = Math.Cos(dist1 * (0.5 * Math.PI));
                    dist2 = 1 - dist1;
                    for (int h = 0; h < 24; h++)
                    {
                        factor1 = Convert.ToDouble(dayWin[h]) * dist1;
                        factor2 = Convert.ToDouble(dayEqu[h]) * dist2;
                        int HOY = d * 24 + h;

                        //calc Iincident_win and Iincident_equ and multiply both with their distances. And thats the interpolated value for HOY
                        //first, calc Iincident_win... loop through all rays, which hit that SP on winter day, time step t. 
                        // check for first and second order... can it be both? yeah sure.
                        double IreflWin = 0.0;
                        double IreflEqu = 0.0;
                        for (int m = 0; m < Inormals_win[i][h].Length; m++)
                        {
                            //first order always. otherwise there would be no normal vector
                            double coeff = obstacles[IObstRef1st_win[i][h][m]].specCoeff[HOY];
                            //but second order, check if it exists. -1 means, there is no 2nd order reflection.
                            if (IObstRef2nd_win[i][h][m] >= 0)
                            {
                                coeff *= obstacles[IObstRef2nd_win[i][h][m]].specCoeff[HOY];
                            }
                            double DNI_t = DNI[HOY] * coeff;
                            IreflWin += DNI_t * Math.Sin((90 * Math.PI / 180.0) - Vector3d.VectorAngle(origNormal[i], Vector3d.Negate(Inormals_win[i][h][m])));
                        }

                        for (int m = 0; m < Inormals_equ[i][h].Length; m++)
                        {
                            //first order always. otherwise there would be no normal vector
                            double coeff = obstacles[IObstRef1st_equ[i][h][m]].specCoeff[HOY];
                            //but second order, check if it exists. -1 means, there is no 2nd order reflection.
                            if (IObstRef2nd_equ[i][h][m] >= 0)
                            {
                                coeff *= obstacles[IObstRef2nd_equ[i][h][m]].specCoeff[HOY];
                            }
                            double DNI_t = DNI[HOY] * coeff;
                            IreflEqu += DNI_t * Math.Sin((90 * Math.PI / 180.0) - Vector3d.VectorAngle(origNormal[i], Vector3d.Negate(Inormals_equ[i][h][m])));
                        }

                        Ispecular_annual_Par[HOY] = IreflWin * factor1 + IreflEqu * factor2;
                    }
                });



                fullF1 = y2;        // 100% shadowequinox on this day
                fullF2 = y3;        // 100% shadowsummer on this day
                InterpolInterv = y3 - y2;
                Parallel.For(y2, y3, paropts, d => //from 20.March to 20.june
                {
                    double dist1, dist2, factor1, factor2;
                    dist1 = Convert.ToDouble((InterpolInterv - Math.Abs(fullF1 - d))) / Convert.ToDouble(InterpolInterv);
                    dist1 = Math.Sin(dist1 * (0.5 * Math.PI));
                    dist2 = 1 - dist1;
                    for (int h = 0; h < 24; h++)
                    {
                        factor1 = Convert.ToDouble(dayEqu[h]) * dist1;
                        factor2 = Convert.ToDouble(daySum[h]) * dist2;
                        int HOY = d * 24 + h;

                        double IreflEqu = 0.0;
                        double IreflSum = 0.0;
                        for (int m = 0; m < Inormals_equ[i][h].Length; m++)
                        {
                            //first order always. otherwise there would be no normal vector
                            double coeff = obstacles[IObstRef1st_equ[i][h][m]].specCoeff[HOY];
                            //but second order, check if it exists. -1 means, there is no 2nd order reflection.
                            if (IObstRef2nd_equ[i][h][m] >= 0)
                            {
                                coeff *= obstacles[IObstRef2nd_equ[i][h][m]].specCoeff[HOY];
                            }
                            double DNI_t = DNI[HOY] * coeff;
                            IreflEqu += DNI_t * Math.Sin((90 * Math.PI / 180.0) - Vector3d.VectorAngle(origNormal[i], Vector3d.Negate(Inormals_equ[i][h][m])));
                        }

                        for (int m = 0; m < Inormals_sum[i][h].Length; m++)
                        {
                            //first order always. otherwise there would be no normal vector
                            double coeff = obstacles[IObstRef1st_sum[i][h][m]].specCoeff[HOY];
                            //but second order, check if it exists. -1 means, there is no 2nd order reflection.
                            if (IObstRef2nd_sum[i][h][m] >= 0)
                            {
                                coeff *= obstacles[IObstRef2nd_sum[i][h][m]].specCoeff[HOY];
                            }
                            double DNI_t = DNI[HOY] * coeff;
                            IreflSum += DNI_t * Math.Sin((90 * Math.PI / 180.0) - Vector3d.VectorAngle(origNormal[i], Vector3d.Negate(Inormals_sum[i][h][m])));
                        }

                        Ispecular_annual_Par[HOY] = IreflEqu * factor1 + IreflSum * factor2;
                    }
                });



                fullF1 = y3;        //100% shadowsummer on this day
                fullF2 = y4;        //100% shadowequinox on this day
                InterpolInterv = y4 - y3;
                //    For i = y3 To y4 - 1            'from 21.June to 22.Sept
                Parallel.For(y3, y4, paropts, d =>
                {
                    double dist1, dist2, factor1, factor2;
                    dist1 = Convert.ToDouble((InterpolInterv - Math.Abs(fullF1 - d))) / Convert.ToDouble(InterpolInterv);
                    dist1 = 1 - dist1;
                    dist1 = Math.Cos(dist1 * (0.5 * Math.PI));
                    dist2 = 1 - dist1;
                    for (int h = 0; h < 24; h++)
                    {
                        factor1 = Convert.ToDouble(daySum[h]) * dist1;
                        factor2 = Convert.ToDouble(dayEqu[h]) * dist2;
                        int HOY = d * 24 + h;

                        double IreflSum = 0.0;
                        double IreflEqu = 0.0;
                        for (int m = 0; m < Inormals_sum[i][h].Length; m++)
                        {
                            //first order always. otherwise there would be no normal vector
                            double coeff = obstacles[IObstRef1st_sum[i][h][m]].specCoeff[HOY];
                            //but second order, check if it exists. -1 means, there is no 2nd order reflection.
                            if (IObstRef2nd_sum[i][h][m] >= 0)
                            {
                                coeff *= obstacles[IObstRef2nd_sum[i][h][m]].specCoeff[HOY];
                            }
                            double DNI_t = DNI[HOY] * coeff;
                            IreflSum += DNI_t * Math.Sin((90 * Math.PI / 180.0) - Vector3d.VectorAngle(origNormal[i], Vector3d.Negate(Inormals_sum[i][h][m])));
                        }

                        for (int m = 0; m < Inormals_equ[i][h].Length; m++)
                        {
                            //first order always. otherwise there would be no normal vector
                            double coeff = obstacles[IObstRef1st_equ[i][h][m]].specCoeff[HOY];
                            //but second order, check if it exists. -1 means, there is no 2nd order reflection.
                            if (IObstRef2nd_equ[i][h][m] >= 0)
                            {
                                coeff *= obstacles[IObstRef2nd_equ[i][h][m]].specCoeff[HOY];
                            }
                            double DNI_t = DNI[HOY] * coeff;
                            IreflEqu += DNI_t * Math.Sin((90 * Math.PI / 180.0) - Vector3d.VectorAngle(origNormal[i], Vector3d.Negate(Inormals_equ[i][h][m])));
                        }

                        Ispecular_annual_Par[HOY] = IreflSum * factor1 + IreflEqu * factor2;
                    }
                });



                fullF1 = y4;    //100% shadowequinox on this day
                fullF2 = y5;    //100% shadowwinter on this day
                InterpolInterv = y5 - y4;
                Parallel.For(y4, y5, paropts, d => // from 23.Sept to 21.Dec
                {
                    double dist1, dist2, factor1, factor2;
                    dist1 = Convert.ToDouble((InterpolInterv - Math.Abs(fullF1 - d))) / Convert.ToDouble(InterpolInterv);
                    dist1 = Math.Sin(dist1 * (0.5 * Math.PI));
                    dist2 = 1 - dist1;
                    for (int h = 0; h < 24; h++)
                    {
                        factor1 = Convert.ToDouble(dayEqu[h]) * dist1;
                        factor2 = Convert.ToDouble(dayWin[h]) * dist2;
                        int HOY = d * 24 + h;

                        double IreflEqu = 0.0;
                        double IreflWin = 0.0;
                        for (int m = 0; m < Inormals_equ[i][h].Length; m++)
                        {
                            //first order always. otherwise there would be no normal vector
                            double coeff = obstacles[IObstRef1st_equ[i][h][m]].specCoeff[HOY];
                            //but second order, check if it exists. -1 means, there is no 2nd order reflection.
                            if (IObstRef2nd_equ[i][h][m] >= 0)
                            {
                                coeff *= obstacles[IObstRef2nd_equ[i][h][m]].specCoeff[HOY];
                            }
                            double DNI_t = DNI[HOY] * coeff;
                            IreflEqu += DNI_t * Math.Sin((90 * Math.PI / 180.0) - Vector3d.VectorAngle(origNormal[i], Vector3d.Negate(Inormals_equ[i][h][m])));
                        }

                        for (int m = 0; m < Inormals_win[i][h].Length; m++)
                        {
                            //first order always. otherwise there would be no normal vector
                            double coeff = obstacles[IObstRef1st_win[i][h][m]].specCoeff[HOY];
                            //but second order, check if it exists. -1 means, there is no 2nd order reflection.
                            if (IObstRef2nd_win[i][h][m] >= 0)
                            {
                                coeff *= obstacles[IObstRef2nd_win[i][h][m]].specCoeff[HOY];
                            }
                            double DNI_t = DNI[HOY] * coeff;
                            IreflWin += DNI_t * Math.Sin((90 * Math.PI / 180.0) - Vector3d.VectorAngle(origNormal[i], Vector3d.Negate(Inormals_win[i][h][m])));
                        }

                        Ispecular_annual_Par[HOY] = IreflEqu * factor1 + IreflWin * factor2;
                    }
                });


                fullF1 = y5;   // 100% shadowwinter on this day
                fullF2 = y6;   // 100% shadowequinox spring on this day
                InterpolInterv = y6 - y5;
                Parallel.For(y5, 365, paropts, d =>    //from 22.Dec to 31.Dec
                {
                    double dist1, dist2, factor1, factor2;
                    dist1 = Convert.ToDouble((InterpolInterv - Math.Abs(fullF1 - d))) / Convert.ToDouble(InterpolInterv);
                    dist1 = 1 - dist1;
                    dist1 = Math.Cos(dist1 * (0.5 * Math.PI));
                    dist2 = 1 - dist1;
                    for (int h = 0; h < 24; h++)
                    {
                        factor1 = Convert.ToDouble(dayWin[h]) * dist1;
                        factor2 = Convert.ToDouble(dayEqu[h]) * dist2;
                        int HOY = d * 24 + h;

                        double IreflWin = 0.0;
                        double IreflEqu = 0.0;
                        for (int m = 0; m < Inormals_win[i][h].Length; m++)
                        {
                            //first order always. otherwise there would be no normal vector
                            double coeff = obstacles[IObstRef1st_win[i][h][m]].specCoeff[HOY];
                            //but second order, check if it exists. -1 means, there is no 2nd order reflection.
                            if (IObstRef2nd_win[i][h][m] >= 0)
                            {
                                coeff *= obstacles[IObstRef2nd_win[i][h][m]].specCoeff[HOY];
                            }
                            double DNI_t = DNI[HOY] * coeff;
                            IreflWin += DNI_t * Math.Sin((90 * Math.PI / 180.0) - Vector3d.VectorAngle(origNormal[i], Vector3d.Negate(Inormals_win[i][h][m])));
                        }

                        for (int m = 0; m < Inormals_equ[i][h].Length; m++)
                        {
                            //first order always. otherwise there would be no normal vector
                            double coeff = obstacles[IObstRef1st_equ[i][h][m]].specCoeff[HOY];
                            //but second order, check if it exists. -1 means, there is no 2nd order reflection.
                            if (IObstRef2nd_equ[i][h][m] >= 0)
                            {
                                coeff *= obstacles[IObstRef2nd_equ[i][h][m]].specCoeff[HOY];
                            }
                            double DNI_t = DNI[HOY] * coeff;
                            IreflEqu += DNI_t * Math.Sin((90 * Math.PI / 180.0) - Vector3d.VectorAngle(origNormal[i], Vector3d.Negate(Inormals_equ[i][h][m])));
                        }

                        Ispecular_annual_Par[HOY] = IreflWin * factor1 + IreflEqu * factor2;
                    }
                });

                Ispecular_annual[i] = Ispecular_annual_Par;
            }
        }


        /// <summary>
        /// Calculates incident beam radiation on each SP based on interpolation of interreflected rays of multiple days. Multi-threading version.
        /// </summary>
        /// <param name="IObstRef1st">Multiple days: References to 1st order specular object of the reflected ray, before it reaches the sensor point. [i][h][m], i=each SP, h=24 hours, m=rays hitting that SP for that timestep.</param>
        /// <param name="IObstRef2nd">Multiple days: References to 2nd order specular object of the reflected ray, before it reaches the sensor point. [i][h][m], i=each SP, h=24 hours, m=rays hitting that SP for that timestep.</param>
        /// <param name="Inormals">Multiple days: Reflected rays reaching a sensor point. [i][h][m], i=each SP, h=24 hours, m=rays hitting that SP for that timestep.</param>
        /// <param name="obstacles">Obstacle objects, containing beam reflection coefficient values for each hour of the year.</param>
        /// <param name="DNI">[t] 8760 DNI values, for each hour of the year.</param>
        /// <param name="origNormal">[i] Normal vectors for each sensor point.</param>
        /// <param name="Ispecular_annual">Interreflected beam irradiation values on each sensor point, for each hour of the year. [i][t], i=each SP, t=8760 hours</param>
        public static void CalcSpecularIncident_AnnualMT(int[] StartDays, int[] EndDays,
            int[][][][] IObstRef1st, int[][][][] IObstRef2nd, Vector3d[][][][] Inormals,
            List<CObstacleObject> obstacles, double[] DNI, Vector3d[] origNormal, ParallelOptions paropts,
            out double[][] Ispecular_annual)
        {
            //using interpolation of several days. 3 or x. start with three days.
            int SPCount = IObstRef1st[0].Length;
            Ispecular_annual = new double[SPCount][];

            double maxPi = 2 * Math.PI / Convert.ToDouble(IObstRef1st.Length);
            int daysUsed = StartDays.Length;

            for (int i = 0; i < SPCount; i++)
            {
                Ispecular_annual[i] = new double[8760];
                double [] Ispecular_annual_Par = new double[8760];

                Parallel.For(0, daysUsed, paropts, d =>
                {
                    int dayStart, dayEnd;
                    int dd;
                    if (d == daysUsed - 1)
                        dd = 0;
                    else
                        dd = d + 1;
                    dayStart = StartDays[d];
                    dayEnd = EndDays[d];
                    for (int n = dayStart; n < dayEnd; n++)
                    {
                        double xprime = 0;
                        if (n < (365 / 4))
                        {
                            int dayStart_qrt = 1;
                            int dayEnd_qrt = 365 / 4;
                            int dayIntervals = dayEnd_qrt - dayStart_qrt;
                            double dist_qrt = (Convert.ToDouble(dayIntervals) - Math.Abs(dayStart_qrt - n)) / Convert.ToDouble(dayIntervals);
                            double xprime_qrt = Math.Sin(dist_qrt * (0.5 * Math.PI));
                            double dist1_max = (Convert.ToDouble(dayIntervals) - Math.Abs(dayStart - dayStart_qrt)) / Convert.ToDouble(dayIntervals);
                            double dist1_min = (Convert.ToDouble(dayIntervals) - Math.Abs(dayEnd - 1 - dayStart_qrt)) / Convert.ToDouble(dayIntervals);
                            double xpr_min = Math.Sin(dist1_min * (0.5 * Math.PI));
                            double xpr_max = Math.Sin(dist1_max * (0.5 * Math.PI));
                            double dist1 = dist_qrt;
                            xprime = (xprime_qrt - xpr_min) / (xpr_max - xpr_min);
                        }
                        else if (n < ((365 / 4) * 2) && n >= ((365 / 4) * 1))
                        {
                            int dayStart_qrt = 365 / 4;
                            int dayEnd_qrt = 365 / 4 * 2;
                            int dayIntervals = dayEnd_qrt - dayStart_qrt;
                            double dist_qrt = (Convert.ToDouble(dayIntervals) - Math.Abs(dayStart_qrt - n)) / Convert.ToDouble(dayIntervals);
                            double xprime_qrt = 1 - Math.Cos(dist_qrt * (0.5 * Math.PI));
                            double dist1_max = (Convert.ToDouble(dayIntervals) - Math.Abs(dayStart - dayStart_qrt)) / Convert.ToDouble(dayIntervals);
                            double dist1_min = (Convert.ToDouble(dayIntervals) - Math.Abs(dayEnd - 1 - dayStart_qrt)) / Convert.ToDouble(dayIntervals);
                            double xpr_min = 1 - Math.Cos(dist1_min * (0.5 * Math.PI));
                            double xpr_max = 1 - Math.Cos(dist1_max * (0.5 * Math.PI));
                            double dist1 = dist_qrt;
                            xprime = (xprime_qrt - xpr_min) / (xpr_max - xpr_min);
                        }
                        else if (n < ((365 / 4) * 3) && n >= ((365 / 4) * 2))
                        {
                            int dayStart_qrt = 365 / 4 * 2;
                            int dayEnd_qrt = 365 / 4 * 3;
                            int dayIntervals = dayEnd_qrt - dayStart_qrt;
                            double dist_qrt = (Convert.ToDouble(dayIntervals) - Math.Abs(dayStart_qrt - n)) / Convert.ToDouble(dayIntervals);
                            double xprime_qrt = Math.Sin(dist_qrt * (0.5 * Math.PI));
                            double dist1_max = (Convert.ToDouble(dayIntervals) - Math.Abs(dayStart - dayStart_qrt)) / Convert.ToDouble(dayIntervals);
                            double dist1_min = (Convert.ToDouble(dayIntervals) - Math.Abs(dayEnd - 1 - dayStart_qrt)) / Convert.ToDouble(dayIntervals);
                            double xpr_min = Math.Sin(dist1_min * (0.5 * Math.PI));
                            double xpr_max = Math.Sin(dist1_max * (0.5 * Math.PI));
                            double dist1 = dist_qrt;
                            xprime = (xprime_qrt - xpr_min) / (xpr_max - xpr_min);
                        }
                        else if (n >= ((365 / 4) * 3))
                        {
                            int dayStart_qrt = 365 / 4 * 3;
                            int dayEnd_qrt = 365;
                            int dayIntervals = dayEnd_qrt - dayStart_qrt;
                            double dist_qrt = (Convert.ToDouble(dayIntervals) - Math.Abs(dayStart_qrt - n)) / Convert.ToDouble(dayIntervals);
                            double xprime_qrt = 1 - Math.Cos(dist_qrt * (0.5 * Math.PI));
                            double dist1_max = (Convert.ToDouble(dayIntervals) - Math.Abs(dayStart - dayStart_qrt)) / Convert.ToDouble(dayIntervals);
                            double dist1_min = (Convert.ToDouble(dayIntervals) - Math.Abs(dayEnd - 1 - dayStart_qrt)) / Convert.ToDouble(dayIntervals);
                            double xpr_min = 1 - Math.Cos(dist1_min * (0.5 * Math.PI));
                            double xpr_max = 1 - Math.Cos(dist1_max * (0.5 * Math.PI));
                            double dist1 = dist_qrt;
                            xprime = (xprime_qrt - xpr_min) / (xpr_max - xpr_min);
                        }

                        double xprime2 = 1 - xprime;

                        //converting arrays which may contain nulls to 0/1 arrays:
                        bool[] day_d = new bool[24];
                        bool[] day_dd = new bool[24];
                        for (int t = 0; t < 24; t++)
                        {
                            if (Misc.IsNullOrEmpty(Inormals[d][i][t]))
                                day_d[t] = false;
                            else
                                day_d[t] = true;

                            if (Misc.IsNullOrEmpty(Inormals[dd][i][t]))
                                day_dd[t] = false;
                            else
                                day_dd[t] = true;
                        }

                        for (int h = 0; h < 24; h++)
                        {
                            double factor1 = Convert.ToDouble(day_d[h]) * xprime;
                            double factor2 = Convert.ToDouble(day_dd[h]) * xprime2;
                            int HOY = (n - 1) * 24 + h;
                            double Irefl_d = 0.0;
                            double Irefl_dd = 0.0;
                            for (int m = 0; m < Inormals[d][i][h].Length; m++)
                            {
                                //first order always. otherwise there would be no normal vector
                                double coeff = obstacles[IObstRef1st[d][i][h][m]].specCoeff[HOY];
                                //but second order, check if it exists. -1 means, there is no 2nd order reflection.
                                if (IObstRef2nd[d][i][h][m] >= 0)
                                {
                                    coeff *= obstacles[IObstRef2nd[d][i][h][m]].specCoeff[HOY];
                                }
                                double DNI_t = DNI[HOY] * coeff;
                                Irefl_d += DNI_t * Math.Sin((90 * Math.PI / 180.0) - Vector3d.VectorAngle(origNormal[i], Vector3d.Negate(Inormals[d][i][h][m])));
                            }

                            for (int m = 0; m < Inormals[dd][i][h].Length; m++)
                            {
                                //first order always. otherwise there would be no normal vector
                                double coeff = obstacles[IObstRef1st[dd][i][h][m]].specCoeff[HOY];
                                //but second order, check if it exists. -1 means, there is no 2nd order reflection.
                                if (IObstRef2nd[dd][i][h][m] >= 0)
                                {
                                    coeff *= obstacles[IObstRef2nd[dd][i][h][m]].specCoeff[HOY];
                                }
                                double DNI_t = DNI[HOY] * coeff;
                                Irefl_dd += DNI_t * Math.Sin((90 * Math.PI / 180.0) - Vector3d.VectorAngle(origNormal[i], Vector3d.Negate(Inormals[dd][i][h][m])));
                            }

                            Ispecular_annual_Par[HOY] = Irefl_d * factor1 + Irefl_dd * factor2;
                        }
                    }
                });

                Ispecular_annual[i] = Ispecular_annual_Par;
            }
        }


        /// <summary>
        /// Calculate beam irradiation incident on all sensor points for one time step, considering incidence angles. Multi-threading version.
        /// </summary>
        /// <param name="SPNormal">Normals of sensor points.</param>
        /// <param name="Ispecular">Specular irradiation values [i][t][m] for each sensor point i, for each time stept t and for each reflected ray m.</param>
        /// <param name="Inormals">Normals for each specular irradiation value [i][t][m].</param>
        /// <param name="DNI">Direct normal irradiation values for one time step t.</param>
        /// <param name="IspecularIncident">Effective specular reflected irradiation [i] incident on the sensor point i, for one time step t.</param>
        public static void CalcSpecularIncidentMT(Vector3d[] SPNormal, double[][][] Ispecular, Vector3d[][][] Inormals,
            double DNI, ParallelOptions paropts, ref double[] IspecularIncident)
        {
            //convert Inormals vectors into solar zenith and solar azimuth. coz thats basically my sun.
            //DNI = DNI * Ispecular (here are my albedos)
            if (Misc.IsNullOrEmpty(Ispecular)) return;

            double[] Ispec_ = new double[SPNormal.Length];
            Parallel.For(0, SPNormal.Length, paropts, i =>
            {
                if (Misc.IsNullOrEmpty(Ispecular[i])) return;
                Ispec_[i] = 0.0;
                if (Misc.IsNullOrEmpty(Ispecular[i][0]) == false)
                {
                    for (int m = 0; m < Ispecular[i][0].Length; m++)
                    {
                        double DNI_t = DNI * Ispecular[i][0][m];
                        Ispec_[i] += DNI_t * Math.Sin((90 * Math.PI / 180.0) - Vector3d.VectorAngle(SPNormal[i], Vector3d.Negate(Inormals[i][0][m])));
                    }
                }
            });

            Ispec_.CopyTo(IspecularIncident, 0);
        }


        /// <summary>
        /// For each sensor point of the analysis mesh, identify new sensor points which need to be calculated for diffuse interreflection. Multi-threading version.
        /// </summary>
        /// <remarks>Actual diffuse interreflection for the new sensor points need to be calculated in a separate routine.</remarks>
        /// <param name="SPmesh">Analysis mesh, on which the sensor points are placed.</param>
        /// <param name="SP">Sensor points [i] of the analysis mesh.</param>
        /// <param name="SPnormal">Normal vectors [i] of the i sensor points.</param>
        /// <param name="_obstacles">Obstacle objects.</param>
        /// <param name="permeables">Permeable obstacle objects. Are treated as solids.</param>
        /// <param name="difDomeRes">Resolution of the hemisphere spanned over each sensor point for diffuse interreflection.</param>
        /// <param name="Idiffuse_SPs">New sensor points [i]: for each sensor point i, diffuse sensor points j which need evaluation.</param>
        /// <param name="Idiff_obstacles">For each sensor point i, indices of obstacles that are hit by interreflected diffuse rays.</param>
        /// <param name="Idiff_domevertexindex">For each sensor point i, indices of dome faces that are will emit diffuse interreflected radiation.</param>
        /// <param name="Idiff_domes">For each sensorpoint i, dome objects which are spanned to calculate itnerreflected diffuse radiation.</param>
        public static void CalcIReflDiff_GetSPs2MT(CObstacleObject SPmesh, Point3d[] SP, Vector3d[] SPnormal,
            List<CObstacleObject> obstacles, List<CPermObject> permeables, int difDomeRes, ParallelOptions paropts,
            out List<List<double>> diffSP_beta_list,
            out List<List<double>> diffSP_psi_list,
            out List<List<Sensorpoints.v3d>> diffSP_normal_list,
            out List<List<Sensorpoints.p3d>> diffSP_coord_list,
            out int[][] Idiff_obstacles,
            out int[][] Idiff_domevertexindex,
            out SkyDome[] Idiff_domes)
        {
            //add permeables to normal obstacles
            List<CObstacleObject> _obstacles = new List<CObstacleObject>();
            _obstacles = obstacles;

            List<double> _alb = new List<double>();
            List<double> _spec = new List<double>();
            for (int t = 0; t < 8760; t++)
            {
                _alb.Add(0);
                _spec.Add(0);
            }
            double _tol = (obstacles.Count > 0) ? obstacles[0].tolerance : 0.01;
            foreach (CPermObject perm in permeables)
            {
                CObstacleObject obst = new CObstacleObject(perm.mesh, _alb, _spec, 3, _tol, "perm", false);
                _obstacles.Add(obst);
            }


            //for each SP: return a list of sensor points, which need to be evaluted in a separate routine. 
            //  that routine calcs Global irradiation for each of the new sensor points. and assigns diffuse irradiation to the SP
            Idiff_obstacles = new int[SP.Length][];
            Idiff_domevertexindex = new int[SP.Length][];
            Idiff_domes = new SkyDome[SP.Length];

            SkyDome[] Idiff_domes_temp = new SkyDome[SP.Length];
            List<List<double>> diffSP_beta_list_temp = new List<List<double>>();
            List<List<double>> diffSP_psi_list_temp = new List<List<double>>();
            List<List<Sensorpoints.v3d>> diffSP_normal_list_temp = new List<List<Sensorpoints.v3d>>();
            List<List<Sensorpoints.p3d>> diffSP_coord_list_temp = new List<List<Sensorpoints.p3d>>();

            //temporary. I don't know yet, how many diffuse sensorpoints per sensorpoint I'll have
            diffSP_beta_list = new List<List<double>>();
            diffSP_psi_list = new List<List<double>>();
            diffSP_normal_list = new List<List<Sensorpoints.v3d>>();
            diffSP_coord_list = new List<List<Sensorpoints.p3d>>();
            List<List<int>> diffobstacleindex_list = new List<List<int>>();
            List<List<int>> diffdomevertindex_list = new List<List<int>>();

            for (int i = 0; i < SP.Length; i++)
            {
                diffSP_beta_list.Add(new List<double>());
                diffSP_psi_list.Add(new List<double>());
                diffSP_normal_list.Add(new List<Sensorpoints.v3d>());
                diffSP_coord_list.Add(new List<Sensorpoints.p3d>());
                diffobstacleindex_list.Add(new List<int>());
                diffdomevertindex_list.Add(new List<int>());

                diffSP_beta_list_temp.Add(new List<double>());
                diffSP_psi_list_temp.Add(new List<double>());
                diffSP_normal_list_temp.Add(new List<Sensorpoints.v3d>());
                diffSP_coord_list_temp.Add(new List<Sensorpoints.p3d>());
            }




            Vector3d betaangle = new Vector3d(0, 0, 1);
            Vector3d psiangle = new Vector3d(0, 1, 0);
            Plane psiplane = new Plane(new Point3d(0, 0, 0), new Vector3d(0, 0, 1));


            //foreach SP
            Vector3d vecZ = new Vector3d(0, 0, 1);
            Parallel.For(0, SP.Length, paropts, i =>
            {
                //offset SP, otherwise there is self-intersection 
                Point3d SPoffset = CMisc.OffsetPt(SP[i], SPnormal[i], SPmesh.tolerance);

                //create a dome, rotate to SP normal
                SkyDome dome = new SkyDome(difDomeRes);
                double[,] R = CMisc.RotationMatrix(new Vector3d(0, 0, 1), SPnormal[i]);
                dome.RotateVertexVectors(R);
                Idiff_domes_temp[i] = dome;

                for (int j = 0; j < dome.VerticesHemisphere.Count; j++)
                {
                    //check for obstruction foreach vertex
                    Vector3d vertexvec = new Vector3d(dome.VertexVectorsSphere[dome.VerticesHemisphere[j]][0], dome.VertexVectorsSphere[dome.VerticesHemisphere[j]][1], dome.VertexVectorsSphere[dome.VerticesHemisphere[j]][2]);
                    Ray3d vertexray = new Ray3d(SPoffset, vertexvec);
                    Dictionary<int, double> inters_dic = new Dictionary<int, double>();
                    for (int u = 0; u < _obstacles.Count; u++)
                    {
                        double inters = Rhino.Geometry.Intersect.Intersection.MeshRay(_obstacles[u].mesh, vertexray);
                        if (inters >= 0)
                        {
                            inters_dic.Add(u, inters);
                        }
                    }
                    if (inters_dic.Count == 0) continue;

                    //sort intersections from closest to furthest.
                    var ordered = inters_dic.OrderBy(x => x.Value).ToDictionary(pair => pair.Key, pair => pair.Value);
                    //from closest intersection, check if normal of that obstacle is < 90, if yes, take that obstacle for later calculating radiation on it.
                    //                                                                          (return SP: Idiffuse_SPs[i].add(new sensorpoint) )
                    int uu = ordered.Keys.ElementAt(0);

                    if (_obstacles[uu].reflType != 0 && _obstacles[uu].reflType != 2) continue;

                    Point3d pX = vertexray.PointAt(ordered[uu]);
                    MeshPoint mshp = _obstacles[uu].mesh.ClosestMeshPoint(pX, 0.0);
                    Vector3d pNormal = _obstacles[uu].mesh.NormalAt(mshp);
                    double vAngle = Vector3d.VectorAngle(pNormal, Vector3d.Negate(vertexvec)) * (180.0 / Math.PI);
                    if (vAngle < 90)
                    {
                        double beta = Vector3d.VectorAngle(pNormal, betaangle) * (180.0 / Math.PI);
                        double psi = Vector3d.VectorAngle(pNormal, psiangle, psiplane) * (180.0 / Math.PI);
                        if (double.IsInfinity(psi)) psi = 0.0;
                        diffSP_beta_list_temp[i].Add(beta);
                        diffSP_psi_list_temp[i].Add(psi);
                        Sensorpoints.v3d _v3d;
                        _v3d.X = pNormal.X;
                        _v3d.Y = pNormal.Y;
                        _v3d.Z = pNormal.Z;
                        diffSP_normal_list_temp[i].Add(_v3d);
                        Sensorpoints.p3d _p3d;
                        _p3d.X = pNormal.X;
                        _p3d.Y = pNormal.Y;
                        _p3d.Z = pNormal.Z;
                        diffSP_coord_list_temp[i].Add(_p3d);
                        diffobstacleindex_list[i].Add(uu);
                        diffdomevertindex_list[i].Add(j);
                    }
                }
            });


            Idiff_domes = Idiff_domes_temp;
            for (int i = 0; i < SP.Length; i++)
            {
                Idiff_obstacles[i] = diffobstacleindex_list[i].ToArray(); //should be of same length as Idiffuse_SPs[i] has sensorpoints
                Idiff_domevertexindex[i] = diffdomevertindex_list[i].ToArray();

                diffSP_beta_list[i] = diffSP_beta_list_temp[i];
                diffSP_psi_list[i] = diffSP_psi_list_temp[i];
                diffSP_normal_list[i] = diffSP_normal_list_temp[i];
                diffSP_coord_list[i] = diffSP_coord_list_temp[i];
            }
        }


        /// <summary>
        /// Calc diffuse irradiation on a list of sensor points. Multi-threading version.
        /// </summary>
        /// <param name="Idiff_SP"></param>
        /// <param name="Idiff_obst"></param>
        /// <param name="Idiff_domevert"></param>
        /// <param name="Idiff_dome"></param>
        /// <param name="DOY"></param>
        /// <param name="LT"></param>
        /// <param name="weather"></param>
        /// <param name="sunvectors"></param>
        /// <param name="obstacles"></param>
        /// <param name="Idiffuse"></param>
        public static void CalcDiffuse2MT(List<List<double>> diffSP_beta_list, List<List<double>> diffSP_psi_list,
            List<List<Sensorpoints.v3d>> diffSP_normal_list, List<List<Sensorpoints.p3d>> diffSP_coord_list, int difDomeRes,
            int[][] Idiff_obst, int[][] Idiff_domevert, SkyDome[] Idiff_dome,
            int DOY, int LT, Context.cWeatherdata weather, SunVector[] sunvectors, Mesh[] obst, List<CObstacleObject> obstacles,
            double tolerance, double snow_threshold, double tilt_treshold, ParallelOptions paropts,
            out double[] Idiffuse)
        {
            ParallelOptions singleparopts = new ParallelOptions { MaxDegreeOfParallelism = 1 };

            int SPiicount = diffSP_beta_list.Count;

            int HOY = (DOY - 1) * 24 + LT;

            Idiffuse = new double[SPiicount];
            double[] Idiffuse_temp = new double[SPiicount];

            double[] Idiff_domearea = new double[SPiicount];
            for (int i = 0; i < Idiff_dome.Length; i++)
            {
                Idiff_domearea[i] = 0.0;
                for (int l = 0; l < Idiff_dome[i].Faces.Count; l++)
                {
                    Idiff_domearea[i] += Idiff_dome[i].FaceAreas[l];
                }
            }


            Parallel.For(0, SPiicount, paropts, i =>
            {
                List<bool> ShdwBeam_hour = new List<bool>();
                List<bool[]> ShdwSky = new List<bool[]>();

                Sensorpoints SPdiff = new Sensorpoints(diffSP_beta_list[i].ToArray(), diffSP_psi_list[i].ToArray(), diffSP_coord_list[i].ToArray(), diffSP_normal_list[i].ToArray(), difDomeRes);
                for (int ii = 0; ii < SPdiff.SPCount; ii++)
                {
                    Point3d orig = new Point3d(SPdiff.coord[ii].X, SPdiff.coord[ii].Y, SPdiff.coord[ii].Z);
                    Vector3d mshvrtnorm = new Vector3d(SPdiff.normal[ii].X, SPdiff.normal[ii].Y, SPdiff.normal[ii].Z);


                    /////////////////////////////////////////////////////////////////////
                    //beam for one hour only.
                    Vector3d[] vec_beam = new Vector3d[1];
                    vec_beam[0] = new Vector3d(sunvectors[HOY].udtCoordXYZ.x, sunvectors[HOY].udtCoordXYZ.y, sunvectors[HOY].udtCoordXYZ.z);
                    bool[] shdw_beam = new bool[1];
                    CShadow.CalcShadowMT(orig, mshvrtnorm, tolerance, vec_beam, obst, ref shdw_beam, singleparopts);
                    ShdwBeam_hour.Add(shdw_beam[0]);
                    /////////////////////////////////////////////////////////////////////


                    /////////////////////////////////////////////////////////////////////
                    //sky dome diffuse
                    Vector3d[] vec_sky = new Vector3d[SPdiff.sky[ii].VerticesHemisphere.Count];
                    for (int u = 0; u < vec_sky.Length; u++)
                    {
                        vec_sky[u] = new Vector3d(
                            SPdiff.sky[ii].VertexVectorsSphere[SPdiff.sky[ii].VerticesHemisphere[u]][0],
                            SPdiff.sky[ii].VertexVectorsSphere[SPdiff.sky[ii].VerticesHemisphere[u]][1],
                            SPdiff.sky[ii].VertexVectorsSphere[SPdiff.sky[ii].VerticesHemisphere[u]][2]);
                    }
                    bool[] shdw_sky = new bool[SPdiff.sky[ii].VerticesHemisphere.Count];
                    CShadow.CalcShadowMT(orig, mshvrtnorm, 0.01, vec_sky, obst, ref shdw_sky, singleparopts);
                    ShdwSky.Add(shdw_sky);
                    /////////////////////////////////////////////////////////////////////
                }



                SPdiff.SetShadows(ShdwBeam_hour, ShdwSky, HOY);
                SPdiff.SetSnowcover(snow_threshold, tilt_treshold, weather);
                SPdiff.CalcIrradiation(DOY, LT, weather, sunvectors);

                double totarea = 0.0;
                for (int f = 0; f < Idiff_dome[i].Faces.Count; f++)
                {
                    totarea += Idiff_dome[i].FaceAreas[f];
                }

                double[] domevertfilled = new double[Idiff_dome[i].VertexVectorsSphere.Count];
                for (int vi = 0; vi < Idiff_dome[i].VertexVectorsSphere.Count; vi++)
                {
                    domevertfilled[vi] = 0.0;
                }
                for (int ii = 0; ii < SPdiff.SPCount; ii++)
                {
                    int v = Idiff_domevert[i][ii];
                    domevertfilled[v] = SPdiff.I[ii][HOY] * obstacles[Idiff_obst[i][ii]].albedos[HOY];
                }

                double totI = 0.0;
                for (int vi = 0; vi < Idiff_dome[i].Faces.Count; vi++)
                {
                    int index1 = Idiff_dome[i].Faces[vi][0];
                    int index2 = Idiff_dome[i].Faces[vi][1];
                    int index3 = Idiff_dome[i].Faces[vi][2];

                    double Isum = domevertfilled[index1] + domevertfilled[index2] + domevertfilled[index3];
                    Isum /= 3;
                    Isum *= Idiff_dome[i].FaceAreas[vi];

                    totI += Isum;
                }
                totI /= totarea;
                Idiffuse_temp[i] = totI;
            });

            Idiffuse = Idiffuse_temp;
        }


        /// <summary>
        /// Calc diffuse irradiation on a list of sensor points for the entire year. Permeable objects treated as blind solids. No further interreflection. 3 days interpolation for beam and sky obstruction calculation for each secondary SP. Multi-threading version.
        /// </summary>
        /// <param name="diffSP_beta_list">List of secondary sensor point tilt angles, requested for diffuse interreflection for each main sensor point. [i][ii], i = each SP, ii = each secondary SP.</param>
        /// <param name="diffSP_psi_list">List of secondary sensor point azimuth angles, requested for diffuse interreflection for each main sensor point. [i][ii], i = each SP, ii = each secondary SP.</param>
        /// <param name="diffSP_normal_list">List of secondary sensor point normals, requested for diffuse interreflection for each main sensor point. [i][ii], i = each SP, ii = each secondary SP.</param>
        /// <param name="diffSP_coord_list">List of secondary sensor point 3d coordinates, requested for diffuse interreflection for each main sensor point. [i][ii], i = each SP, ii = each secondary SP.</param>
        /// <param name="Idiff_obst">For each main SP, an array of indices of diffuse obstacles on which secondary SP lie on. [i][ii], i = each SP, ii = each secondary SP.</param>
        /// <param name="Idiff_domevert">For each main SP, an array of indices of skydome vertices, which are used for secondary SP calculations. [i][ii], i = each SP, ii = each secondary SP.</param>
        /// <param name="Idiff_dome">For each main SP, one SkyDome. [i] = each SP.</param>
        /// <param name="SecondarydifDomeRes">Resolution of skydome spanned over secondary SP.</param>
        /// <param name="year">Year of calculation.</param>
        /// <param name="weather">Weather data.</param>
        /// <param name="sunvectors">8760 sunvectors, [t] = for each hour of the year.</param>
        /// <param name="obstacles">List of obstacle objects.</param>
        /// <param name="permeables">List of permeable obstacle objects.</param>
        /// <param name="tolerance">Tolerance for obstruction calculations, to avoid self-obstruction of obstacles.</param>
        /// <param name="snow_threshold">Snow thickness threshold, after which snow does not remain on surface.</param>
        /// <param name="tilt_treshold">Tilt angle threshold, after which snow does not remain on surface.</param>
        /// <param name="groundalbedo">Albedo of the ground. 8760 time series, [t] = for each hour of the year.</param>
        /// <param name="Idiffuse">Time series of diffuse interreflection. [i][t], i=each sensor point, t=each hour of the year 1-8760.</param>
        public static void CalcDiffuse_AnnualMT(List<List<double>> diffSP_beta_list, List<List<double>> diffSP_psi_list,
            List<List<Sensorpoints.v3d>> diffSP_normal_list, List<List<Sensorpoints.p3d>> diffSP_coord_list,
            int[][] Idiff_obst, int[][] Idiff_domevert, SkyDome[] Idiff_dome, int SecondarydifDomeRes,
            int year, Context.cWeatherdata weather, SunVector[] sunvectors, List<CObstacleObject> obstacles, List<CPermObject> permeables,
            double tolerance, double snow_threshold, double tilt_treshold, double[] groundalbedo, ParallelOptions paropts,
            out double[][] Idiffuse)
        {
            ParallelOptions singleparopts = new ParallelOptions { MaxDegreeOfParallelism = 1 };

            //using 3 day interpolation of beam radiation. and no interreflections. and no trees.
            List<Mesh> obstlist = new List<Mesh>();
            foreach (CObstacleObject obstobj in obstacles)
            {
                obstlist.Add(obstobj.mesh);
            }
            foreach (CPermObject perm in permeables)
            {
                obstlist.Add(perm.mesh);
            }
            Mesh[] obst = obstlist.ToArray();


            int SPiicount = diffSP_beta_list.Count;

            Idiffuse = new double[SPiicount][];
            double[] Idiff_domearea = new double[SPiicount];
            Parallel.For(0, Idiff_dome.Length, paropts, i =>
            {
                Idiff_domearea[i] = 0.0;
                for (int l = 0; l < Idiff_dome[i].Faces.Count; l++)
                {
                    Idiff_domearea[i] += Idiff_dome[i].FaceAreas[l];
                }
            });



            int[] equsol = SunVector.GetEquinoxSolstice(year);
            int HOYequ = (equsol[0] - 1) * 24;
            int HOYsum = (equsol[1] - 1) * 24;
            int HOYwin = (equsol[3] - 1) * 24;

            for (int i = 0; i < SPiicount; i++)
            {
                Idiffuse[i] = new double[8760];

                Sensorpoints SPdiff = new Sensorpoints(diffSP_beta_list[i].ToArray(), diffSP_psi_list[i].ToArray(), diffSP_coord_list[i].ToArray(), diffSP_normal_list[i].ToArray(), SecondarydifDomeRes);
                bool[][] ShdwSky_array = new bool[SPdiff.SPCount][];
                bool[][] ShdwBeameq_array = new bool[SPdiff.SPCount][];
                bool[][] ShdwBeamwin_array = new bool[SPdiff.SPCount][];
                bool[][] ShdwBeamsum_array = new bool[SPdiff.SPCount][];
                Parallel.For(0, SPdiff.SPCount, paropts, ii =>
                {
                    Point3d orig = new Point3d(SPdiff.coord[ii].X, SPdiff.coord[ii].Y, SPdiff.coord[ii].Z);
                    Vector3d mshvrtnorm = new Vector3d(SPdiff.normal[ii].X, SPdiff.normal[ii].Y, SPdiff.normal[ii].Z);

                    /////////////////////////////////////////////////////////////////////
                    //beam for 24 hours and 3 days.... get equinox and solstices ?   
                    // equinox:             march 20
                    // summer solstice:     june 21
                    // winter solstice:     december 21
                    Vector3d[] vec_beam_equ = new Vector3d[24];
                    Vector3d[] vec_beam_sum = new Vector3d[24];
                    Vector3d[] vec_beam_win = new Vector3d[24];
                    bool[] sunshine_equ = new bool[24];
                    bool[] sunshine_sum = new bool[24];
                    bool[] sunshine_win = new bool[24];
                    for (int t = 0; t < 24; t++)
                    {
                        if (sunvectors[HOYequ + t].Sunshine)
                            sunshine_equ[t] = true;
                        if (sunvectors[HOYsum + t].Sunshine)
                            sunshine_sum[t] = true;
                        if (sunvectors[HOYwin + t].Sunshine)
                            sunshine_win[t] = true;
                        vec_beam_equ[t] = new Vector3d(sunvectors[HOYequ + t].udtCoordXYZ.x, sunvectors[HOYequ + t].udtCoordXYZ.y, sunvectors[HOYequ + t].udtCoordXYZ.z);
                        vec_beam_sum[t] = new Vector3d(sunvectors[HOYsum + t].udtCoordXYZ.x, sunvectors[HOYsum + t].udtCoordXYZ.y, sunvectors[HOYsum + t].udtCoordXYZ.z);
                        vec_beam_win[t] = new Vector3d(sunvectors[HOYwin + t].udtCoordXYZ.x, sunvectors[HOYwin + t].udtCoordXYZ.y, sunvectors[HOYwin + t].udtCoordXYZ.z);
                    }
                    bool[] shdw_beam_equ = new bool[24];
                    bool[] shdw_beam_sum = new bool[24];
                    bool[] shdw_beam_win = new bool[24];
                    CShadow.CalcShadowMT(orig, mshvrtnorm, tolerance, vec_beam_equ, sunshine_equ, obst, ref shdw_beam_equ, singleparopts);
                    CShadow.CalcShadowMT(orig, mshvrtnorm, tolerance, vec_beam_sum, sunshine_sum, obst, ref shdw_beam_sum, singleparopts);
                    CShadow.CalcShadowMT(orig, mshvrtnorm, tolerance, vec_beam_win, sunshine_win, obst, ref shdw_beam_win, singleparopts);

                    ShdwBeameq_array[ii] = new bool[shdw_beam_equ.Length];
                    ShdwBeamwin_array[ii] = new bool[shdw_beam_win.Length];
                    ShdwBeamsum_array[ii] = new bool[shdw_beam_sum.Length];
                    shdw_beam_equ.CopyTo(ShdwBeameq_array[ii], 0);
                    shdw_beam_win.CopyTo(ShdwBeamwin_array[ii], 0);
                    shdw_beam_sum.CopyTo(ShdwBeamsum_array[ii], 0);
                    /////////////////////////////////////////////////////////////////////

                    /////////////////////////////////////////////////////////////////////
                    //sky dome diffuse
                    Vector3d[] vec_sky = new Vector3d[SPdiff.sky[ii].VerticesHemisphere.Count];
                    for (int u = 0; u < vec_sky.Length; u++)
                    {
                        vec_sky[u] = new Vector3d(
                            SPdiff.sky[ii].VertexVectorsSphere[SPdiff.sky[ii].VerticesHemisphere[u]][0],
                            SPdiff.sky[ii].VertexVectorsSphere[SPdiff.sky[ii].VerticesHemisphere[u]][1],
                            SPdiff.sky[ii].VertexVectorsSphere[SPdiff.sky[ii].VerticesHemisphere[u]][2]);
                    }
                    bool[] shdw_sky = new bool[SPdiff.sky[ii].VerticesHemisphere.Count];
                    CShadow.CalcShadowMT(orig, mshvrtnorm, tolerance, vec_sky, obst, ref shdw_sky, singleparopts);
                    // ShdwSky.Add(shdw_sky);
                    ShdwSky_array[ii] = new bool[shdw_sky.Length];
                    shdw_sky.CopyTo(ShdwSky_array[ii], 0);
                    /////////////////////////////////////////////////////////////////////
                });

                List<bool[]> ShdwBeam_equinox = new List<bool[]>(ShdwBeameq_array);
                List<bool[]> ShdwBeam_summer = new List<bool[]>(ShdwBeamsum_array);
                List<bool[]> ShdwBeam_winter = new List<bool[]>(ShdwBeamwin_array);
                List<bool[]> ShdwSky = new List<bool[]>(ShdwSky_array);
                SPdiff.SetShadowsInterpolatedMT(ShdwBeam_equinox, ShdwBeam_summer, ShdwBeam_winter, ShdwSky, paropts);
                SPdiff.SetSimpleGroundReflectionMT(diffSP_beta_list[i].ToArray(), groundalbedo, weather, sunvectors, paropts);
                SPdiff.SetSnowcoverMT(snow_threshold, tilt_treshold, weather, paropts);
                SPdiff.CalcIrradiationMT(weather, sunvectors, paropts);


                double totarea = 0.0;
                for (int f = 0; f < Idiff_dome[i].Faces.Count; f++)
                {
                    totarea += Idiff_dome[i].FaceAreas[f];
                }

                double[][] domevertfilled = new double[Idiff_dome[i].VertexVectorsSphere.Count][];
                Parallel.For(0, Idiff_dome[i].VertexVectorsSphere.Count, vi =>
                {
                    domevertfilled[vi] = new double[8760];
                    for (int t = 0; t < 8760; t++)
                    {
                        domevertfilled[vi][t] = 0.0;
                    }
                });
                Parallel.For(0, SPdiff.SPCount, ii =>
                {
                    int v = Idiff_domevert[i][ii];
                    for (int t = 0; t < 8760; t++)
                    {
                        domevertfilled[v][t] = SPdiff.I[ii][t] * obstacles[Idiff_obst[i][ii]].albedos[t];
                    }
                });
                double[] Idifftemp = new double[8760];
                Parallel.For(0, 8760, paropts, t =>
                {
                    double totI = 0.0;
                    for (int vi = 0; vi < Idiff_dome[i].Faces.Count; vi++)
                    {
                        int index1 = Idiff_dome[i].Faces[vi][0];
                        int index2 = Idiff_dome[i].Faces[vi][1];
                        int index3 = Idiff_dome[i].Faces[vi][2];
                        double Isum = domevertfilled[index1][t] + domevertfilled[index2][t] + domevertfilled[index3][t];
                        Isum /= 3;
                        Isum *= Idiff_dome[i].FaceAreas[vi];
                        totI += Isum;
                    }
                    totI /= totarea;
                    Idifftemp[t] = totI;
                });
                Idifftemp.CopyTo(Idiffuse[i], 0);
            }

        }


        /// <summary>
        /// Calc simplified diffuse irradiation on a list of sensor points for the entire year. Permeable objects treated as blind solids. No further interreflection. No beam and sky obstruction calculation for secondary sensor points, only based on tilt angle of secondary SP. Multi-threading version.
        /// </summary>
        /// <param name="diffSP_beta_list">List of secondary sensor point tilt angles, requested for diffuse interreflection for each main sensor point. [i][ii], i = each SP, ii = each secondary SP.</param>
        /// <param name="diffSP_psi_list">List of secondary sensor point azimuth angles, requested for diffuse interreflection for each main sensor point. [i][ii], i = each SP, ii = each secondary SP.</param>
        /// <param name="diffSP_normal_list">List of secondary sensor point normals, requested for diffuse interreflection for each main sensor point. [i][ii], i = each SP, ii = each secondary SP.</param>
        /// <param name="diffSP_coord_list">List of secondary sensor point 3d coordinates, requested for diffuse interreflection for each main sensor point. [i][ii], i = each SP, ii = each secondary SP.</param>
        /// <param name="Idiff_obst">For each main SP, an array of indices of diffuse obstacles on which secondary SP lie on. [i][ii], i = each SP, ii = each secondary SP.</param>
        /// <param name="Idiff_domevert">For each main SP, an array of indices of skydome vertices, which are used for secondary SP calculations. [i][ii], i = each SP, ii = each secondary SP.</param>
        /// <param name="Idiff_dome">For each main SP, one SkyDome. [i] = each SP.</param>
        /// <param name="difDomeRes">Resolution of skydome spanned over secondary SP.</param>
        /// <param name="weather">Weather data.</param>
        /// <param name="sunvectors">8760 sunvectors, [t] = for each hour of the year.</param>
        /// <param name="obstacles">List of obstacle objects.</param>
        /// <param name="snow_threshold">Snow thickness threshold, after which snow does not remain on surface.</param>
        /// <param name="tilt_treshold">Tilt angle threshold, after which snow does not remain on surface.</param>
        /// <param name="groundalbedo">Albedo of the ground. 8760 time series, [t] = for each hour of the year.</param>
        /// <param name="Idiffuse">Time series of diffuse interreflection. [i][t], i=each sensor point, t=each hour of the year 1-8760.</param>
        public static void CalcDiffuse_AnnualSimpleMT(List<List<double>> diffSP_beta_list, List<List<double>> diffSP_psi_list,
            List<List<Sensorpoints.v3d>> diffSP_normal_list, List<List<Sensorpoints.p3d>> diffSP_coord_list,
            int[][] Idiff_obst, int[][] Idiff_domevert, SkyDome[] Idiff_dome, int difDomeRes,
            Context.cWeatherdata weather, SunVector[] sunvectors, List<CObstacleObject> obstacles,
            double snow_threshold, double tilt_treshold, double[] groundalbedo, ParallelOptions paropts,
            out double[][] Idiffuse)
        {
            int SPiicount = diffSP_beta_list.Count;

            Idiffuse = new double[SPiicount][];

            double[] Idiff_domearea = new double[SPiicount];
            Parallel.For(0, Idiff_dome.Length, paropts, i =>
            {
                Idiff_domearea[i] = 0.0;
                for (int l = 0; l < Idiff_dome[i].Faces.Count; l++)
                {
                    Idiff_domearea[i] += Idiff_dome[i].FaceAreas[l];
                }
            });


            for (int i = 0; i < SPiicount; i++)
            {
                Idiffuse[i] = new double[8760];

                Sensorpoints SPdiff = new Sensorpoints(diffSP_beta_list[i].ToArray(), diffSP_psi_list[i].ToArray(), diffSP_coord_list[i].ToArray(), diffSP_normal_list[i].ToArray(), difDomeRes);
                SPdiff.SetSimpleSkyMT(diffSP_beta_list[i].ToArray(), paropts);                //Simplified version: no obstruction calculations for seconardy sensor points.
                SPdiff.SetSimpleGroundReflectionMT(diffSP_beta_list[i].ToArray(), groundalbedo, weather, sunvectors, paropts);
                SPdiff.SetSnowcoverMT(snow_threshold, tilt_treshold, weather, paropts);
                SPdiff.CalcIrradiationMT(weather, sunvectors, paropts);

                double totarea = 0.0;
                for (int f = 0; f < Idiff_dome[i].Faces.Count; f++)
                {
                    totarea += Idiff_dome[i].FaceAreas[f];
                }

                double[][] domevertfilled = new double[Idiff_dome[i].VertexVectorsSphere.Count][];
                Parallel.For(0, Idiff_dome[i].VertexVectorsSphere.Count, paropts, vi =>
                {
                    domevertfilled[vi] = new double[8760];
                    for (int t = 0; t < 8760; t++)
                    {
                        domevertfilled[vi][t] = 0.0;
                    }
                });
                Parallel.For(0, SPdiff.SPCount, paropts, ii =>
                {
                    int v = Idiff_domevert[i][ii];
                    for (int t = 0; t < 8760; t++)
                    {
                        domevertfilled[v][t] = SPdiff.I[ii][t] * obstacles[Idiff_obst[i][ii]].albedos[t];
                    }
                });
                for (int t = 0; t < 8760; t++)
                {
                    double totI = 0.0;
                    for (int vi = 0; vi < Idiff_dome[i].Faces.Count; vi++)
                    {
                        int index1 = Idiff_dome[i].Faces[vi][0];
                        int index2 = Idiff_dome[i].Faces[vi][1];
                        int index3 = Idiff_dome[i].Faces[vi][2];
                        double Isum = domevertfilled[index1][t] + domevertfilled[index2][t] + domevertfilled[index3][t];
                        Isum /= 3;
                        Isum *= Idiff_dome[i].FaceAreas[vi];
                        totI += Isum;
                    }
                    totI /= totarea;
                    Idiffuse[i][t] = totI;
                }
            }
        }

    }
}
