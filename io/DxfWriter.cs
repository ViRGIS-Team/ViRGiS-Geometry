using System;
using System.Collections.Generic;
using System.IO;
using System.Text;
using netDxf;
using netDxf.Entities;
using netDxf.IO;

namespace VirgisGeometry
{

    public class DxfWriter : IMeshWriter
    {
        public IOWriteResult Write(BinaryWriter writer, List<WriteMesh> vMeshes, WriteOptions options)
        { 
            return Write(writer.BaseStream, vMeshes, options);
        }

        public IOWriteResult Write(TextWriter writer, List<WriteMesh> vMeshes, WriteOptions options)
        {
            throw new NotImplementedException();
        }
        
        public IOWriteResult Write(Stream writer, List<WriteMesh> vMeshes, WriteOptions options)
        { 
            DxfDocument doc = new ();
            foreach (WriteMesh wmesh in vMeshes)
            {
                IMesh imesh = wmesh.Mesh;
                for (int i = 0; i<imesh.TriangleCount; i++)
                {
                    Index3i tri = imesh.GetTriangle(i);
                    doc.AddEntityToDocument( new Face3D(
                        imesh.GetVertex(tri.a),
                        imesh.GetVertex(tri.b),
                        imesh.GetVertex(tri.c)
                        ), true);
                }
            }

            doc.Save(writer);
            return IOWriteResult.Ok;
        }
    }
}