using System;
using System.Drawing;
using Grasshopper.Kernel;

namespace GHSolar
{
    public class GHSolarInfo : GH_AssemblyInfo
    {
        public override string Name
        {
            get
            {
                return "GHSolar";
            }
        }
        public override Bitmap Icon
        {
            get
            {
                //Return a 24x24 pixel bitmap to represent this GHA library.
                return null;
            }
        }
        public override string Description
        {
            get
            {
                //Return a short string describing the purpose of this GHA library.
                return "";
            }
        }
        public override Guid Id
        {
            get
            {
                return new Guid("5632b14a-f5fc-4113-827a-4a97fbedd22f");
            }
        }

        public override string AuthorName
        {
            get
            {
                //Return a string identifying you or your company.
                return "Empa";
            }
        }
        public override string AuthorContact
        {
            get
            {
                //Return a string representing your preferred contact details.
                return "";
            }
        }
    }
}
