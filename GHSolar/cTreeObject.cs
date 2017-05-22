using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using Rhino.Geometry;
using SolarModel;

/*
 * cTreeObject.cs
 * Copyright 2017 Christoph Waibel <chwaibel@student.ethz.ch>
 * 
 * This work is licensed under the GNU GPL license version 3.
*/

namespace GHSolar
{
    internal class cTreeObject
    {
        internal Mesh mesh;
        internal List<double> leaves;

        internal cTreeObject(Mesh _mesh, List<double> _leaves)
        {
            mesh = _mesh;
            leaves = new List<double>(_leaves);
        }
    }

}
