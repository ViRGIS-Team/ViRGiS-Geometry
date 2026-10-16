using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using netDxf;
using netDxf.Entities;
using netDxf.IO;
using UnityEngine;


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

        public IOReadResult Read(Stream stream, ReadOptions options, IMeshBuilder builder, AxisOrder ax, bool debug=true)
        {
            try
            {
                DxfDocument doc = DxfDocument.Load(stream);

                if (debug) {
                    Debug.Log($"ENTITIES for the Active Layout = {doc.Entities.ActiveLayout}:");
                    Debug.Log($"\t{EntityType.Arc}; count: {doc.Entities.Arcs.Count()}");
                    Debug.Log($"\t{EntityType.Circle}; count: {doc.Entities.Circles.Count()}");
                    Debug.Log($"\t{EntityType.Dimension}; count: {doc.Entities.Dimensions.Count()}");
                    Debug.Log($"\t{EntityType.Ellipse}; count: {doc.Entities.Ellipses.Count()}");
                    Debug.Log($"\t{EntityType.Face3D}; count: {doc.Entities.Faces3D.Count()}");
                    Debug.Log($"\t{EntityType.Hatch}; count: {doc.Entities.Hatches.Count()}");
                    Debug.Log($"\t{EntityType.Image}; count: {doc.Entities.Images.Count()}");
                    Debug.Log($"\t{EntityType.Insert}; count: {doc.Entities.Inserts.Count()}");
                    Debug.Log($"\t{EntityType.Leader}; count: {doc.Entities.Leaders.Count()}");
                    Debug.Log($"\t{EntityType.Line}; count: {doc.Entities.Lines.Count()}");
                    Debug.Log($"\t{EntityType.Mesh}; count: {doc.Entities.Meshes.Count()}");
                    Debug.Log($"\t{EntityType.MLine}; count: {doc.Entities.MLines.Count()}");
                    Debug.Log($"\t{EntityType.MText}; count: {doc.Entities.MTexts.Count()}");
                    Debug.Log($"\t{EntityType.Point}; count: {doc.Entities.Points.Count()}");
                    Debug.Log($"\t{EntityType.PolyfaceMesh}; count: {doc.Entities.PolyfaceMeshes.Count()}");
                    Debug.Log($"\t{EntityType.PolygonMesh}; count: {doc.Entities.PolygonMeshes.Count()}");
                    Debug.Log($"\t{EntityType.Polyline2D}; count: { doc.Entities.Polylines2D.Count()}");
                    Debug.Log($"\t{EntityType.Polyline3D}; count: {doc.Entities.Polylines3D.Count()}");
                    Debug.Log($"\t{EntityType.Shape}; count: {doc.Entities.Shapes.Count()}");
                    Debug.Log($"\t{EntityType.Solid}; count: {doc.Entities.Solids.Count()}");
                    Debug.Log($"\t{EntityType.Spline}; count: { doc.Entities.Splines.Count()}");
                    Debug.Log($"\t{EntityType.Text}; count: {doc.Entities.Texts.Count()}");
                    Debug.Log($"\t{EntityType.Ray}; count: {doc.Entities.Rays.Count()}");
                    Debug.Log($"\t{EntityType.Underlay}; count: {doc.Entities.Underlays.Count()}");
                    Debug.Log($"\t{EntityType.Viewport}; count: {doc.Entities.Viewports.Count()}");
                    Debug.Log($"\t{EntityType.Wipeout}; count: {doc.Entities.Wipeouts.Count()}");
                    Debug.Log($"\t{EntityType.XLine}; count: {doc.Entities.XLines.Count()}");
                }

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