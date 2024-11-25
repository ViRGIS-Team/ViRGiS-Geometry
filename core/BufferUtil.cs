using System;
using System.Collections.Generic;
using System.Linq;
using System.IO;
using System.IO.Compression;
using UnityEngine;

namespace VirgisGeometry
{
    /// <summary>
    /// Convenience functions for working with arrays. 
    ///    - Math functions on arrays of floats/doubles
    ///    - "automatic" conversion from IEnumerable<T> (via introspection)
    ///    - byte[] conversions
    ///    - zlib compress/decompress byte[] buffers
    /// </summary>
    public class BufferUtil
    {
        static public void SetVertex3(double[] v, int i, double x, double y, double z) {
            v[3 * i] = x;
            v[3 * i + 1] = y;
            v[3 * i + 2] = z;
        }
        static public void SetVertex3(float[] v, int i, float x, float y, float z)
        {
            v[3 * i] = x;
            v[3 * i + 1] = y;
            v[3 * i + 2] = z;
        }

        static public void SetVertex2(double[] v, int i, double x, double y)
        {
            v[2 * i] = x;
            v[2 * i + 1] = y;
        }
        static public void SetVertex2(float[] v, int i, float x, float y)
        {
            v[2 * i] = x;
            v[2 * i + 1] = y;
        }

        static public void SetTriangle(int[] v, int i, int a, int b, int c)
        {
            v[3 * i] = a;
            v[3 * i + 1] = b;
            v[3 * i + 2] = c;
        }


        static public double Dot(double[] a, double [] b)
        {
            double dot = 0;
            for (int i = 0; i < a.Length; ++i)
                dot += a[i] * b[i];
            return dot;
        }

        static public void MultiplyAdd(double[] dest, double multiply, double[] add)
        {
            for (int i = 0; i < dest.Length; ++i)
                dest[i] += multiply * add[i];
        }

        static public void MultiplyAdd(double[] dest, double[] multiply, double[] add)
        {
            for (int i = 0; i < dest.Length; ++i)
                dest[i] += multiply[i] * add[i];
        }

        static public double MultiplyAdd_GetSqrSum(double[] dest, double multiply, double[] add)
        {
            double sum = 0;
            for (int i = 0; i < dest.Length; ++i) {
                dest[i] += multiply * add[i];
                sum += dest[i] * dest[i];
            }
            return sum;
        }

        static public double DistanceSquared(double[] a, double[] b)
        {
            double sum = 0;
            for (int i = 0; i < a.Length; ++i)
                sum += (a[i] - b[i]) * (a[i] - b[i]);
            return sum;
        }



        static public void ParallelDot(double[] a, double[][] b, double[][] result)
        {
            int N = a.Length, count = b.Length;
            gParallel.BlockStartEnd(0, N-1, (i0, i1) => {
                for (int i = i0; i <= i1; i++) {
                    for (int j = 0; j < count; ++j)
                        result[j][i] = a[i] * b[j][i];
                }
            }, 1000);
        }


        static public double[][] AllocNxM(int N, int M)
        {
            double[][] d = new double[N][];
            for (int k = 0; k <  N; ++k)
                d[k] = new double[M];
            return d;
        }

        static public double[][] InitNxM(int N, int M, double[][] init)
        {
            double[][] d = AllocNxM(N, M);
            for (int k = 0; k < N; ++k)
                Array.Copy(init[k], d[k], M);
            return d;
        }


        /// <summary>
        /// Count number of elements in array (or up to max_i) that pass FilterF test
        /// </summary>
        static public int CountValid<T>(T[] data, Func<T, bool> FilterF, int max_i = -1) {
			int n = (max_i == -1) ? data.Length : max_i;
			int valid = 0;
			for (int i = 0; i < n; ++i) {
				if (FilterF(data[i]))
					valid++;
			}
			return valid;
		}

		/// <summary>
		/// shifts elements of array (or up to max_i) that pass FilterF to front of list,
		/// and returns number that passed
		/// </summary>
		static public int FilterInPlace<T>(T[] data, Func<T,bool> FilterF, int max_i = -1) {
			int N = (max_i == -1) ? data.Length : max_i;
			int k = 0;
			for (int i = 0; i < N; ++i) {
				if (FilterF(data[i]))
					data[k++] = data[i];
			}
			return k;
		}

		/// <summary>
		/// return a new array containing only elements (or up to max_i) that pass FilterF test
		/// </summary>
		static public T[] Filter<T>(T[] data, Func<T, bool> FilterF, int max_i = -1) {
			int n = (max_i == -1) ? data.Length : max_i;
			int valid = CountValid(data, FilterF);
			if (valid == 0)
				return null;
			T[] result = new T[valid];
			int k = 0;
			for (int i = 0; i < n; ++i) {
				if (FilterF(data[i]))
					result[k++] = data[i];
			}
			return result;
		}




        /// <summary>
        /// convert input set into Vector3d.
        /// Supports packed list of float/double tuples, list of Vector3f/Vector3d
        /// </summary>
        static public IEnumerable<Vector3d> ToVector3d<T>(IEnumerable<T> values) {
            Type t = typeof(T);

            if (t == typeof(float)) {
                IEnumerator<float> valuesf = (values as IEnumerable<float>).GetEnumerator();
                while (valuesf.MoveNext() )
                {
                    float x = valuesf.Current;
                    if (!valuesf.MoveNext()) throw new Exception("Vector Data is incomplete");
                    float y = valuesf.Current;
                    if (!valuesf.MoveNext()) throw new Exception("Vector Data is incomplete");
                    float z = valuesf.Current;
                    yield return new Vector3d(x, y, z);
                }
            } else if (t == typeof(double)) {
                IEnumerator<double> valuesd = (values as IEnumerable<double>).GetEnumerator();
                while (valuesd.MoveNext() )
                {
                    double x = valuesd.Current;
                    if (!valuesd.MoveNext()) throw new Exception("Vector Data is incomplete");
                    double y = valuesd.Current;
                    if (!valuesd.MoveNext()) throw new Exception("Vector Data is incomplete");
                    double z = valuesd.Current;
                    yield return new Vector3d(x, y, z);
                }
            } else if (t == typeof(Vector3f)) {
                IEnumerable<Vector3f> valuesvf = (values as IEnumerable<Vector3f>);
                foreach (Vector3f v in valuesvf)
                    yield return (Vector3d)v;
            }
            else if (t == typeof(Vector3))
            {
                IEnumerable<Vector3> valuesuv = (values as IEnumerable<Vector3>);
                foreach (Vector3 v in valuesuv)
                    yield return (Vector3d)v;

            }
            else if (t == typeof(Vector3d)) {
                IEnumerable<Vector3d> valuesvd = values as IEnumerable<Vector3d>;
                foreach (Vector3d v in valuesvd)
                    yield return v;
            } else
                throw new NotSupportedException("ToVector3d: unknown type " + t.ToString());
        }



        /// <summary>
        /// convert input set into Vector3f.
        /// Supports packed list of float/double tuples, list of Vector3f/Vector3d
        /// </summary>
        static public IEnumerable<Vector3f> ToVector3f<T>(IEnumerable<T> values)
        {
            Type t = typeof(T);

            if (t == typeof(float))
            {
                IEnumerator<float> valuesf = (values as IEnumerable<float>).GetEnumerator();
                while(valuesf.MoveNext())
                {
                    float x = valuesf.Current;
                    if (! valuesf.MoveNext() ) throw new Exception("Vector Data is incomplete");
                    float y = valuesf.Current;
                    if (!valuesf.MoveNext()) throw new Exception("Vector Data is incomplete");
                    float z = valuesf.Current;
                    yield return new Vector3f(x, y, z);
                }
            }
            else if (t == typeof(double))
            {
                IEnumerator<double> valuesd = (values as IEnumerable<double>).GetEnumerator();
                while (valuesd.MoveNext() )
                {
                    double x = valuesd.Current;
                    if (!valuesd.MoveNext()) throw new Exception("Vector Data is incomplete");
                    double y = valuesd.Current;
                    if (!valuesd.MoveNext()) throw new Exception("Vector Data is incomplete");
                    double z = valuesd.Current;
                    yield return new Vector3f(x, y, z);
                }
            }
            else if (t == typeof(Vector3f))
            {
                IEnumerable<Vector3f> valuesvf = (values as IEnumerable<Vector3f>);
                foreach (Vector3f v in valuesvf)
                    yield return v;

            }
            else if (t == typeof(Vector3))
            {
                IEnumerable<Vector3> valuesuv = (values as IEnumerable<Vector3>);
                foreach (Vector3 v in valuesuv)
                    yield return (Vector3f)v;

            }
            else if (t == typeof(Vector3d))
            {
                IEnumerable<Vector3d> valuesvd = values as IEnumerable<Vector3d>;
                foreach (Vector3d v in valuesvd)
                    yield return (Vector3f)v;
            }
            else
                throw new NotSupportedException("ToVector3f: unknown type " + t.ToString());
        }




        /// <summary>
        /// convert input set into Index3i.
        /// Supports packed list of int tuples, list of Vector3i/Index3i
        /// </summary>
        static public IEnumerable<Index3i> ToIndex3i<T>(IEnumerable<T> values) {

            Type t = typeof(T);
            if (t == typeof(int)) {
                IEnumerator<int> valuesi = (values as IEnumerable<int>).GetEnumerator();
                while (valuesi.MoveNext())
                {
                    int a = valuesi.Current;
                    if (!valuesi.MoveNext()) throw new Exception("Triangle is incomplete");
                    int b = valuesi.Current;
                    if (!valuesi.MoveNext()) throw new Exception("Triangle is incomplete");
                    int c = valuesi.Current;
                    yield return new Index3i(a, b, c);
                }
            } else if (t == typeof(Index3i)) {
                IEnumerable<Index3i> valuesvi = values as IEnumerable<Index3i>;
                foreach (Index3i i in valuesvi)
                    yield return i;
            } else if (t == typeof(Vector3i)) {
                IEnumerable<Vector3i> valuesvi = values as IEnumerable<Vector3i>;
                foreach (Vector3i v in valuesvi)
                    yield return (Index3i)v;

            } else
                throw new NotSupportedException("ToIndex2i: unknown type " + t.ToString());
        }

        /// <summary>
        /// convert input set into Index2i.
        /// Supports packed list of int tuples, list of Vector2i/Index2i
        /// </summary>
        static public IEnumerable<Index2i> ToIndex2i<T>(IEnumerable<T> values)
        {
            Type t = typeof(T);
            if (t == typeof(int))
            {
                IEnumerator<int> valuesi = (values as IEnumerable<int>).GetEnumerator();
                while(valuesi.MoveNext())
                {
                    int a = valuesi.Current;
                    if (! valuesi.MoveNext()) throw new Exception("Index2i data is not complete");
                    int b = valuesi.Current;
                    yield return new Index2i(a, b);
                }
            }
            else if (t == typeof(Index2i))
            {
                IEnumerable<Index2i> valuesvi = values as IEnumerable<Index2i>;
                foreach (Index2i i in valuesvi)
                    yield return i;
            }
            else if (t == typeof(Vector2i))
            {
                IEnumerable<Vector2i> valuesvi = values as IEnumerable<Vector2i>;
                foreach (Vector2i v in valuesvi)
                    yield return (Index2i)v;

            }
            else
                throw new NotSupportedException("ToIndex2i: unknown type " + t.ToString());
        }



        /// <summary>
        /// convert byte array to int array
        /// </summary>
        static public int[] ToInt(byte[] buffer)
        {
            int sz = sizeof(int);
            int Nvals = buffer.Length / sz;
            int[] v = new int[Nvals];
            for (int i = 0; i < Nvals; i++) {
                v[i] = BitConverter.ToInt32(buffer, i * sz);
            }
            return v;
        }


        /// <summary>
        /// convert byte array to short array
        /// </summary>
        static public short[] ToShort(byte[] buffer)
        {
            int sz = sizeof(short);
            int Nvals = buffer.Length / sz;
            short[] v = new short[Nvals];
            for (int i = 0; i < Nvals; i++) {
                v[i] = BitConverter.ToInt16(buffer, i * sz);
            }
            return v;
        }


        /// <summary>
        /// convert byte array to double array
        /// </summary>
        static public double[] ToDouble(byte[] buffer)
        {
            int sz = sizeof(double);
            int Nvals = buffer.Length / sz;
            double[] v = new double[Nvals];
            for (int i = 0; i < Nvals; i++) {
                v[i] = BitConverter.ToDouble(buffer, i * sz);
            }
            return v;
        }


        /// <summary>
        /// convert byte array to float array
        /// </summary>
        static public float[] ToFloat(byte[] buffer)
        {
            int sz = sizeof(float);
            int Nvals = buffer.Length / sz;
            float[] v = new float[Nvals];
            for (int i = 0; i < Nvals; i++) {
                v[i] = BitConverter.ToSingle(buffer, i * sz);
            }
            return v;
        }


        /// <summary>
        /// convert byte array to VectorArray3d
        /// </summary>
        static public VectorArray3d ToVectorArray3d(byte[] buffer)
        {
            int sz = sizeof(double);
            int Nvals = buffer.Length / sz;
            int Nvecs = Nvals / 3;
            VectorArray3d v = new VectorArray3d(Nvecs);
            for (int i = 0; i < Nvecs; i++) {
                double x = BitConverter.ToDouble(buffer, (3 * i) * sz);
                double y = BitConverter.ToDouble(buffer, (3 * i + 1) * sz);
                double z = BitConverter.ToDouble(buffer, (3 * i + 2) * sz);
                v.Set(i, x, y, z);
            }
            return v;
        }



        /// <summary>
        /// convert byte array to VectorArray2f
        /// </summary>
        static public VectorArray2f ToVectorArray2f(byte[] buffer)
        {
            int sz = sizeof(float);
            int Nvals = buffer.Length / sz;
            int Nvecs = Nvals / 2;
            VectorArray2f v = new VectorArray2f(Nvecs);
            for (int i = 0; i < Nvecs; i++) {
                float x = BitConverter.ToSingle(buffer, (2 * i) * sz);
                float y = BitConverter.ToSingle(buffer, (2 * i + 1) * sz);
                v.Set(i, x, y);
            }
            return v;
        }

        /// <summary>
        /// convert byte array to VectorArray3f
        /// </summary>
        static public VectorArray3f ToVectorArray3f(byte[] buffer)
        {
            int sz = sizeof(float);
            int Nvals = buffer.Length / sz;
            int Nvecs = Nvals / 3;
            VectorArray3f v = new VectorArray3f(Nvecs);
            for (int i = 0; i < Nvecs; i++) {
                float x = BitConverter.ToSingle(buffer, (3 * i) * sz);
                float y = BitConverter.ToSingle(buffer, (3 * i + 1) * sz);
                float z = BitConverter.ToSingle(buffer, (3 * i + 2) * sz);
                v.Set(i, x, y, z);
            }
            return v;
        }




        /// <summary>
        /// convert byte array to VectorArray3i
        /// </summary>
        static public VectorArray3i ToVectorArray3i(byte[] buffer)
        {
            int sz = sizeof(int);
            int Nvals = buffer.Length / sz;
            int Nvecs = Nvals / 3;
            VectorArray3i v = new VectorArray3i(Nvecs);
            for (int i = 0; i < Nvecs; i++) {
                int x = BitConverter.ToInt32(buffer, (3 * i) * sz);
                int y = BitConverter.ToInt32(buffer, (3 * i + 1) * sz);
                int z = BitConverter.ToInt32(buffer, (3 * i + 2) * sz);
                v.Set(i, x, y, z);
            }
            return v;
        }


        /// <summary>
        /// convert byte array to IndexArray4i
        /// </summary>
        static public IndexArray4i ToIndexArray4i(byte[] buffer)
        {
            int sz = sizeof(int);
            int Nvals = buffer.Length / sz;
            int Nvecs = Nvals / 4;
            IndexArray4i v = new IndexArray4i(Nvecs);
            for (int i = 0; i < Nvecs; i++) {
                int a = BitConverter.ToInt32(buffer, (4 * i) * sz);
                int b = BitConverter.ToInt32(buffer, (4 * i + 1) * sz);
                int c = BitConverter.ToInt32(buffer, (4 * i + 2) * sz);
                int d = BitConverter.ToInt32(buffer, (4 * i + 3) * sz);
                v.Set(i, a, b, c, d);
            }
            return v;
        }


        /// <summary>
        /// convert int array to bytes
        /// </summary>
        static public byte[] ToBytes(int[] array)
        {
            byte[] result = new byte[array.Length * sizeof(int)];
            Buffer.BlockCopy(array, 0, result, 0, result.Length);
            return result;
        }

        /// <summary>
        /// convert short array to bytes
        /// </summary>
        static public byte[] ToBytes(short[] array)
        {
            byte[] result = new byte[array.Length * sizeof(short)];
            Buffer.BlockCopy(array, 0, result, 0, result.Length);
            return result;
        }

        /// <summary>
        /// convert float array to bytes
        /// </summary>
        static public byte[] ToBytes(float[] array)
        {
            byte[] result = new byte[array.Length * sizeof(float)];
            Buffer.BlockCopy(array, 0, result, 0, result.Length);
            return result;
        }

        /// <summary>
        /// convert double array to bytes
        /// </summary>
        static public byte[] ToBytes(double[] array)
        {
            byte[] result = new byte[array.Length * sizeof(double)];
            Buffer.BlockCopy(array, 0, result, 0, result.Length);
            return result;
        }




        /// <summary>
        /// Compress a byte buffer using Deflate/ZLib compression. 
        /// </summary>
        static public byte[] CompressZLib(byte[] buffer, bool bFast)
        {
            MemoryStream ms = new MemoryStream();
#if NET_2_0 || NET_2_0_SUBSET
            DeflateStream zip = new DeflateStream(ms, CompressionMode.Compress);
#else
            DeflateStream zip = new DeflateStream(ms, (bFast) ? System.IO.Compression.CompressionLevel.Fastest : System.IO.Compression.CompressionLevel.Optimal, true);
#endif
            zip.Write(buffer, 0, buffer.Length);
            zip.Close();
            ms.Position = 0;

            byte[] compressed = new byte[ms.Length];
            ms.Read(compressed, 0, compressed.Length);

            byte[] zBuffer = new byte[compressed.Length + 4];
            Buffer.BlockCopy(compressed, 0, zBuffer, 4, compressed.Length);
            Buffer.BlockCopy(BitConverter.GetBytes(buffer.Length), 0, zBuffer, 0, 4);
            return zBuffer;
        }


        /// <summary>
        /// Decompress a byte buffer that has been compressed using Deflate/ZLib compression
        /// </summary>
        static public byte[] DecompressZLib(byte[] zBuffer)
        {
            MemoryStream ms = new MemoryStream();
            int msgLength = BitConverter.ToInt32(zBuffer, 0);
            ms.Write(zBuffer, 4, zBuffer.Length - 4);

            byte[] buffer = new byte[msgLength];

            ms.Position = 0;
            DeflateStream zip = new DeflateStream(ms, CompressionMode.Decompress);
            zip.Read(buffer, 0, buffer.Length);

            return buffer;
        }


    }




    /// <summary>
    /// utility class for porting C++ code that uses this kind of idiom:
    ///    T * ptr = &array[i];
    ///    ptr[k] = value
    /// </summary>
    public struct ArrayAlias<T>
    {
        public T[] Source;
        public int Index;

        public ArrayAlias(T[] source, int i)
        {
            Source = source;
            this.Index = i;
        }

        public T this[int i]
        {
            get { return Source[Index + i]; }
            set { Source[Index + i] = value; }
        }
    }


}
