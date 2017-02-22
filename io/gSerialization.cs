﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;

namespace g3
{
    public static class gSerialization
    {

        public static void Store(Vector3f v, BinaryWriter writer)
        {
            writer.Write(v.x);
            writer.Write(v.y);
            writer.Write(v.z);
        }
        public static void Restore(ref Vector3f v, BinaryReader reader)
        {
            v.x = reader.ReadSingle();
            v.y = reader.ReadSingle();
            v.z = reader.ReadSingle();
        }

        public static void Store(Quaternionf q, BinaryWriter writer)
        {
            writer.Write(q.x);
            writer.Write(q.y);
            writer.Write(q.z);
            writer.Write(q.w);
        }
        public static void Restore(ref Quaternionf q, BinaryReader reader)
        {
            q.x = reader.ReadSingle();
            q.y = reader.ReadSingle();
            q.z = reader.ReadSingle();
            q.w = reader.ReadSingle();
        }



        public static void Store(Frame3f vFrame, BinaryWriter writer)
        {
            Store(vFrame.Origin, writer);
            Store(vFrame.Rotation, writer);
        }
        public static void Restore(ref Frame3f vFrame, BinaryReader reader)
        {
            Vector3f origin = Vector3f.Zero;
            Quaternionf orientation = Quaternionf.Identity;
            Restore(ref origin, reader);
            Restore(ref orientation, reader);
            vFrame = new Frame3f(origin, orientation);
        }





        public static void Store(List<int> values, BinaryWriter writer)
        {
            writer.Write(values.Count);
            for (int i = 0; i < values.Count; ++i)
                writer.Write(values[i]);
        }
        public static void Restore(List<int> values, BinaryReader reader)
        {
            int N = reader.ReadInt32();
            for (int i = 0; i < N; ++i)
                values.Add(reader.ReadInt32());
        }



        public static void Store(Polygon2d polygon, BinaryWriter writer)
        {
            writer.Write(polygon.VertexCount);
            for (int i = 0; i < polygon.VertexCount; ++i) {
                writer.Write(polygon[i].x);
                writer.Write(polygon[i].y);
            }
        }

        public static void Restore(Polygon2d polygon, BinaryReader reader)
        {
            int count = reader.ReadInt32();
            for ( int i = 0; i < count; ++i ) {
                double x = reader.ReadDouble();
                double y = reader.ReadDouble();
                polygon.AppendVertex(new Vector2d(x, y));
            }
        }




        public static void Store(GeneralPolygon2d polygon, BinaryWriter writer)
        {
            Store(polygon.Outer, writer);
            writer.Write(polygon.Holes.Count);
            for ( int i = 0; i < polygon.Holes.Count; ++i )
                Store(polygon.Holes[i], writer);
        }

        public static void Restore(GeneralPolygon2d polygon, BinaryReader reader)
        {
            Polygon2d outer = new Polygon2d();
            Restore(outer, reader);
            polygon.Outer = outer;

            int hole_count = reader.ReadInt32();
            for ( int i = 0; i < hole_count; ++i ) {
                Polygon2d holepoly = new Polygon2d();
                Restore(holepoly, reader);
                polygon.AddHole(holepoly, false);
            }
        }




        public static int DMesh3Version = 1;

        public static void Store(DMesh3 mesh, BinaryWriter writer)
        {
            writer.Write(DMesh3Version);

            int nComponents = (int)mesh.Components;
            writer.Write(nComponents);

            Store(mesh.VerticesBuffer, writer);
            Store(mesh.TrianglesBuffer, writer);
            Store(mesh.EdgesBuffer, writer);
            Store(mesh.EdgesRefCounts.RawRefCounts, writer);

            if ((mesh.Components & MeshComponents.VertexNormals) != 0)
                Store(mesh.NormalsBuffer, writer);
            if ((mesh.Components & MeshComponents.VertexColors) != 0)
                Store(mesh.ColorsBuffer, writer);
            if ((mesh.Components & MeshComponents.VertexUVs) != 0)
                Store(mesh.UVBuffer, writer);
            if ((mesh.Components & MeshComponents.FaceGroups) != 0)
                Store(mesh.GroupsBuffer, writer);
        }



        public static void Restore(DMesh3 mesh, BinaryReader reader)
        {
            int version = reader.ReadInt32();
            if (version != DMesh3Version)
                throw new Exception("gSerialization.Restore: Incorrect DMesh3Version!");

            MeshComponents components = (MeshComponents)reader.ReadInt32();

            Restore(mesh.VerticesBuffer, reader);
            Restore(mesh.TrianglesBuffer, reader);
            Restore(mesh.EdgesBuffer, reader);
            Restore(mesh.EdgesRefCounts.RawRefCounts, reader);

            if ((components & MeshComponents.VertexNormals) != 0) {
                mesh.EnableVertexNormals(Vector3f.AxisY);
                Restore(mesh.NormalsBuffer, reader);
            } else
                mesh.DiscardVertexNormals();

            if ((components & MeshComponents.VertexColors) != 0) {
                mesh.EnableVertexColors(Vector3f.One);
                Restore(mesh.ColorsBuffer, reader);
            } else
                mesh.DiscardVertexColors();

            if ((components & MeshComponents.VertexUVs) != 0) {
                mesh.EnableVertexUVs(Vector2f.Zero);
                Restore(mesh.UVBuffer, reader);
            } else
                mesh.DiscardVertexUVs();

            if ((components & MeshComponents.FaceGroups) != 0) {
                mesh.EnableTriangleGroups(0);
                Restore(mesh.GroupsBuffer, reader);
            } else
                mesh.DiscardTriangleGroups();

            mesh.RebuildFromEdgeRefcounts();
        }


        // [TODO] these could be a lot faster if DVector had a block-iterator...
        public static void Store(DVector<double> vec, BinaryWriter writer)
        {
            byte[] buffer = new byte[vec.BlockCount * sizeof(double)];
            int N = vec.Length;
            writer.Write(N);
            foreach ( DVector<double>.DBlock block in vec.BlockIterator() ) {
                Buffer.BlockCopy(block.data, 0, buffer, 0, block.usedCount*sizeof(double));
                writer.Write(buffer, 0, block.usedCount*sizeof(double));
            }
        }
        public static void Restore(DVector<double> vec, BinaryReader reader)
        {
            int N = reader.ReadInt32();
            byte[] bytes = reader.ReadBytes(N * sizeof(double));
            double[] buffer = new double[N];
            Buffer.BlockCopy(bytes, 0, buffer, 0, bytes.Length);
            vec.Initialize(buffer);
        }

        public static void Store(DVector<float> vec, BinaryWriter writer)
        {
            byte[] buffer = new byte[vec.BlockCount * sizeof(float)];
            int N = vec.Length;
            writer.Write(N);
            foreach ( DVector<float>.DBlock block in vec.BlockIterator() ) {
                Buffer.BlockCopy(block.data, 0, buffer, 0, block.usedCount*sizeof(float));
                writer.Write(buffer, 0, block.usedCount*sizeof(float));
            }
        }
        public static void Restore(DVector<float> vec, BinaryReader reader)
        {
            int N = reader.ReadInt32();
            byte[] bytes = reader.ReadBytes(N * sizeof(float));
            float[] buffer = new float[N];
            Buffer.BlockCopy(bytes, 0, buffer, 0, bytes.Length);
            vec.Initialize(buffer);
        }

        public static void Store(DVector<int> vec, BinaryWriter writer)
        {
            byte[] buffer = new byte[vec.BlockCount * sizeof(int)];
            int N = vec.Length;
            writer.Write(N);
            foreach ( DVector<int>.DBlock block in vec.BlockIterator() ) {
                Buffer.BlockCopy(block.data, 0, buffer, 0, block.usedCount*sizeof(int));
                writer.Write(buffer, 0, block.usedCount*sizeof(int));
            }
        }
        public static void Restore(DVector<int> vec, BinaryReader reader)
        {
            int N = reader.ReadInt32();
            byte[] bytes = reader.ReadBytes(N * sizeof(int));
            int[] buffer = new int[N];
            Buffer.BlockCopy(bytes, 0, buffer, 0, bytes.Length);
            vec.Initialize(buffer);
        }

        public static void Store(DVector<short> vec, BinaryWriter writer)
        {
            byte[] buffer = new byte[vec.BlockCount * sizeof(short)];
            int N = vec.Length;
            writer.Write(N);
            foreach ( DVector<short>.DBlock block in vec.BlockIterator() ) {
                Buffer.BlockCopy(block.data, 0, buffer, 0, block.usedCount*sizeof(short));
                writer.Write(buffer, 0, block.usedCount*sizeof(short));
            }
        }
        public static void Restore(DVector<short> vec, BinaryReader reader)
        {
            int N = reader.ReadInt32();
            byte[] bytes = reader.ReadBytes(N * sizeof(short));
            short[] buffer = new short[N];
            Buffer.BlockCopy(bytes, 0, buffer, 0, bytes.Length);
            vec.Initialize(buffer);
        }
    }
}
