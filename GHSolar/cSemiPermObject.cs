using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using Rhino.Geometry;
using SolarModel;

/*
 * cSemiPermObject.cs
 * Copyright 2017 Christoph Waibel <chwaibel@student.ethz.ch>
 * 
 * This work is licensed under the GNU GPL license version 3.
*/

namespace GHSolar
{
    internal class cSemiPermObject
    {
        internal Mesh mesh;
        internal double[] permeability = new double[8760];

        /// <summary>
        /// Semipermeable object (e.g. a tree).
        /// </summary>
        /// <param name="_mesh">Geometry.</param>
        /// <param name="_permeability">8760 coefficients for each hour of the year, 1.0 : non-permeable, 0.0 : fully permeable.</param>
        internal cSemiPermObject(Mesh _mesh, double[] _permeability)
        {
            mesh = _mesh;
            _permeability.CopyTo(permeability,0);
        }
    }

}
