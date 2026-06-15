using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using netDxf;
using netDxf.Entities;
using netDxf.IO;

namespace VirgisGeometry
{
    /// <summary>
    /// Try opening with netDxf - this will only open files in autoCAD version 2000 or later.
    /// </summary>
    public class DxfReader : IMeshReader
    {
        public IOReadResult Read(TextReader reader, ReadOptions options, IMeshBuilder builder)
        {
            throw new NotImplementedException();
        }

        public IOReadResult Read(BinaryReader reader, ReadOptions options, IMeshBuilder builder)
        {
            return Read(reader.BaseStream, options, builder);
        }

        public IOReadResult Read(Stream stream, ReadOptions options, IMeshBuilder builder)
        {
            try
            {
                DxfDocument doc = DxfDocument.Load(stream);

                List<Face3D> faces = doc.Entities.Faces3D.ToList();
                IEnumerable<PolyfaceMesh> pfs = doc.Entities.PolyfaceMeshes;

                builder.AppendNewMesh(false, false, false, false);
                //
                // Add the Polyface Meshes
                //
                foreach (PolyfaceMesh pfmesh in pfs)
                {
                    List<EntityObject> pffaces = pfmesh.Explode();
                    foreach (EntityObject obj in pffaces)
                    {
                        if (obj.Type == EntityType.Face3D)
                        {
                            faces.Add(obj as Face3D);
                        }
                    }
                }

                foreach (Face3D face in faces)
                {

                    int a = builder.AppendVertex(face.FirstVertex.X, face.FirstVertex.Y, face.FirstVertex.Z);
                    int b = builder.AppendVertex(face.SecondVertex.X, face.SecondVertex.Y, face.SecondVertex.Z);
                    int c = builder.AppendVertex(face.ThirdVertex.X, face.ThirdVertex.Y, face.ThirdVertex.Z);
                    builder.AppendTriangle(a, b, c);
                    if (face.FourthVertex != face.ThirdVertex)
                    {
                        a = builder.AppendVertex(face.FirstVertex.X, face.FirstVertex.Y, face.FirstVertex.Z);
                        b = builder.AppendVertex(face.ThirdVertex.X, face.ThirdVertex.Y, face.ThirdVertex.Z);
                        c = builder.AppendVertex(face.FourthVertex.X, face.FourthVertex.Y, face.FourthVertex.Z);
                        builder.AppendTriangle(a, b, c);
                    }
                }

                return IOReadResult.Ok;
            }
            catch (DxfVersionNotSupportedException ex)
            {
                return new IOReadResult(IOCode.FormatNotSupportedError, "");
            }
            catch (Exception ex)
            {
                return new IOReadResult(IOCode.GenericReaderError, ex.ToString());
            }
        }
    }
}