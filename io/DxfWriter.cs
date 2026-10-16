using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using netDxf;
using netDxf.Entities;
using netDxf.IO;
using UnityEngine;

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
        
        public IOWriteResult Write(Stream writer, List<WriteMesh> vMeshes, WriteOptions options, bool debug=true)
        { 
            DxfDocument doc = new ();
            foreach (WriteMesh wmesh in vMeshes)
            {
                IMesh imesh = wmesh.Mesh;
                for (int i = 0; i<imesh.TriangleCount; i++)
                {
                    Index3i tri = imesh.GetTriangle(i);
                    doc.Entities.Add( new Face3D(
                        imesh.GetVertex(tri.a),
                        imesh.GetVertex(tri.b),
                        imesh.GetVertex(tri.c)
                        ));
                }
            }

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

            doc.Save(writer);
            return IOWriteResult.Ok;
        }
    }
}