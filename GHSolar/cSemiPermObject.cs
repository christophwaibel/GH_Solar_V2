using Rhino.Geometry;

/*
 * cSemiPermObject.cs
 * Copyright 2017 Christoph Waibel <chwaibel@student.ethz.ch>
 * 
 * This work is licensed under the GNU GPL license version 3.
*/

namespace GHSolar
{
    public class CPermObject
    {
        public Mesh mesh;
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
