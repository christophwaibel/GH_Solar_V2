using Rhino.Geometry;

/*
 * cSemiPermObject.cs
 * Copyright 2017 Christoph Waibel <chwaibel@student.ethz.ch>
 * 
 * This work is licensed under the GNU GPL license version 3.
*/

namespace GHSolar
{
    /// <summary>
    /// Permeable obstacle object, e.g. a tree.
    /// </summary>
    public class CPermObject
    {
        /// <summary>
        /// Rhino mesh geometry of the obstacle
        /// </summary>
        public Mesh mesh;
        /// <summary>
        /// 8760 time series with coefficients [0.0, 1.0] to define permeability. 1.0 means it is non-permeable.
        /// </summary>
        public double[] permeability = new double[8760];

        /// <summary>
        /// Permeable object (e.g. a tree).
        /// </summary>
        /// <param name="_mesh">Geometry.</param>
        /// <param name="_permeability">8760 coefficients for each hour of the year, 1.0 : non-permeable, 0.0 : fully permeable.</param>
        public CPermObject(Mesh _mesh, double[] _permeability)
        {
            mesh = _mesh;
            _permeability.CopyTo(permeability,0);
        }
    }

}
