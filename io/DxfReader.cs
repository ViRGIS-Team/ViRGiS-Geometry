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
            return Read(reader.BaseStream, options, builder, default);
        }

        private static int GetOrAddIndex(List<Vector3d> vertices, Vector3d v)
        {
            int index = vertices.IndexOf(v);

            if (index < 0)
            {
                vertices.Add(v);
                index = vertices.Count - 1;
            }

            return index;
        }

        public IOReadResult Read(Stream stream, ReadOptions options, IMeshBuilder builder, AxisOrder ax)
        {
            try
            {
                DxfDocument doc = DxfDocument.Load(stream);

                List<Face3D> faces = doc.Entities.Faces3D.ToList();
                IEnumerable<PolyfaceMesh> pfs = doc.Entities.PolyfaceMeshes;
                
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
                
                List<Vector3d> vertices = new ();
                List<Index3i> triangles = new ();
                
                foreach (Face3D face in faces)
                {

                    int a = GetOrAddIndex(vertices, face.FirstVertex);
                    int b = GetOrAddIndex(vertices,face.SecondVertex);
                    int c = GetOrAddIndex(vertices, face.ThirdVertex);
                    triangles.Add(new Index3i(a, b, c));
                    if (face.FourthVertex != face.ThirdVertex)
                    {
                        a = GetOrAddIndex(vertices, face.FirstVertex);
                        b = GetOrAddIndex(vertices,face.ThirdVertex);
                        c = GetOrAddIndex(vertices, face.FourthVertex);
                        triangles.Add(new Index3i(a, b, c));
                    }
                }

                builder.AppendNewMesh(DMesh3Builder.Build<Vector3d, Index3i, Vector3d>(vertices, triangles ));
                
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