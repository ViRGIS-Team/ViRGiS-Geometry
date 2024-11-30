using System;
using System.Collections.Generic;
using System.Linq;

namespace VirgisGeometry
{
    /// <summary>
    /// DSubmesh3 done properly - with rteferences to the base DMesh3 instead of copying and mapping
    /// 
    /// Note the main change here is that any change to the Submesh will autmatically also apply to the base mesh,
    /// this is a feature.
    /// 
    /// The DSubmesh3 is now a DMesh3 - so there is no need to fetch Submesh from the DSubmesh3.
    /// 
    /// If you want the previous behaviour - do Copy on the DSubmesh3
    /// </summary>
    public class DSubmesh3: DMesh3
    {
        protected DMesh3 BaseMesh; // This is the base mesh that the submesh references
        protected IEnumerable<int> selectedTris; // keeps the selected triangles for this 
 
        // redirect the base data structures to the base mesh
        internal override DVector<double> vertices 
        { 
            get { return BaseMesh.vertices; }
            set { BaseMesh.vertices = value; }
        }
        internal override DVector<float> normals
        {
            get { return BaseMesh.normals; }
            set { BaseMesh.normals = value; }
        }
        internal override DVector<float> colors
        {
            get { return BaseMesh.colors; }
            set { BaseMesh.colors = value; }
        }
        internal override DVector<float> uv
        {
            get { return BaseMesh.uv; }
            set { BaseMesh.uv = value; }
        }
        internal override SmallListSet vertex_edges
        {
            get { return BaseMesh.vertex_edges; }
            set { BaseMesh.vertex_edges = value; }
        }
        internal override DVector<int> triangles
        {
            get { return BaseMesh.triangles; }
            set { BaseMesh.triangles = value; }
        }

        internal override DVector<int> triangle_edges
        {
            get { return BaseMesh.triangle_edges; }
            set { BaseMesh.triangle_edges = value; }
        }
        internal override DVector<int> triangle_groups
        {
            get { return BaseMesh.triangle_groups; }
            set { BaseMesh.triangle_groups = value; }
        }
        internal override DVector<int> edges
        {
            get { return BaseMesh.edges; }
            set { BaseMesh.edges = value; }
        }
        public override bool Clockwise
        {
            get { return BaseMesh.Clockwise; }
            set { BaseMesh.Clockwise = value; }
        }
        internal override bool bFlip
        {
            get { return BaseMesh.bFlip; }
            set { BaseMesh.bFlip = value; }
        }
        internal override bool bFlipNormals
        {
            get { return BaseMesh.bFlipNormals; }
            set { BaseMesh.bFlipNormals = value; }
        }
        public override AxisOrder axisOrder
        {
            get { return BaseMesh.axisOrder; }
            set { BaseMesh.axisOrder = value; }
        }
        internal override Dictionary<string, object> Metadata
        {
            get { return BaseMesh.Metadata; }
            set { BaseMesh.Metadata = value; }
        }


        // Constructors
        public DSubmesh3(DMesh3 mesh, int[] subTriangles): base(true)
        {
            BaseMesh = mesh;
            compute(subTriangles);
        }

        public DSubmesh3(DMesh3 mesh, IEnumerable<int> subTriangles, int nTriEstimate = 0) : base(true)
        {
            BaseMesh = mesh;
            compute(subTriangles);
        }

        public DSubmesh3(DMesh3 mesh): base(true)
        {
            BaseMesh = mesh;
        }



        /// <summary>
        /// Create the submesh using the supplied triangle array
        /// </summary>
        /// <param name="subTriangles"></param>
        public void Compute(int[] subTriangles)
        {
            selectedTris = subTriangles;
            compute(subTriangles);
        }

        /// <summary>
        /// Create the submesh using the supplied triangle iterator
        /// </summary>
        /// <param name="subTriangles"></param>
        public void Compute(IEnumerable<int> subTriangles, int nTriEstimate = 0)
        {
            selectedTris = subTriangles;
            compute(subTriangles);
        }

        /// <summary>
        /// Create the submesh as an n-ring around vertex vID
        /// </summary>
        /// <param name="vID"></param>
        /// <param name="n"></param>
        public void Compute(int vID, int n)
        {
            if (!BaseMesh.IsVertex(vID))
                throw new Exception("Not a valid vertex");
            MeshFaceSelection mfs = new(BaseMesh);
            mfs.ExpandToOneRingNeighbours(n - 1);
            mfs.SelectVertexOneRing(vID);
            compute(mfs);
        }


        void compute(IEnumerable<int> triangles)
        {
            compute(triangles.GetEnumerator());
        }

        void compute(IEnumerator<int> triangles)
        {
            foreach (int tri in TriangleIndices())
                base.RemoveTriangle(tri);
            while (triangles.MoveNext()) {
                if (!BaseMesh.IsTriangle(triangles.Current)) continue; // don't get shouty over incorrect triangle
                foreach (int vert in BaseMesh.GetTriangle(triangles.Current))
                {
                    if (!IsVertex(vert))
                    {
                        vertices_refcount.allocate_at_unsafe(vert);
                    } else
                    {
                        vertices_refcount.increment(vert);
                    }
                }
                foreach (int edge in BaseMesh.GetTriEdges(triangles.Current))
                {
                    if (!IsEdge(edge))
                    {
                        edges_refcount.allocate_at_unsafe(edge);
                    } else
                    {
                        edges_refcount.increment(edge);
                    }
                }
                triangles_refcount.allocate_at_unsafe(triangles.Current);
            }
            vertices_refcount.rebuild_free_list();
            edges_refcount.rebuild_free_list();
            triangles_refcount.rebuild_free_list();
        }

        public override int AppendVertex(ref NewVertexInfo info) 
        {
            int vID = BaseMesh.AppendVertex(ref info);
            vertices_refcount.allocate_at(vID);
            updateTimeStamp(true);
            return vID;
        }

        public override int AppendVertex(DMesh3 from, int fromVID)
        {
            int vID = BaseMesh.AppendVertex(from, fromVID);
            vertices_refcount.allocate_at(vID);
            updateTimeStamp(true);
            return vID;
        }

        public override MeshResult InsertVertex(int vid, ref NewVertexInfo info, bool bUnsafe = false)
        {
            bool bOK = (bUnsafe) ? vertices_refcount.allocate_at_unsafe(vid) :
                       vertices_refcount.allocate_at(vid);
            if (bOK == false)
                return MeshResult.Failed_CannotAllocateVertex;
            MeshResult mr = BaseMesh.InsertVertex(vid, ref info, bUnsafe);
            if (mr == MeshResult.Ok)
                updateTimeStamp(true);
            return mr;
        }

        public override int AppendTriangle(Index3i tv, int gid = -1) 
        {
            int tri = BaseMesh.AppendTriangle(tv, gid);
            triangles_refcount.allocate_at(tri);
            updateTimeStamp(true);
            return tri;
        }

        public override MeshResult InsertTriangle(int tid, Index3i tv, int gid = -1, bool bUnsafe = false)
        {
            bool bOK = (bUnsafe) ? triangles_refcount.allocate_at_unsafe(tid) :
                       triangles_refcount.allocate_at(tid);
            if (bOK == false)
                return MeshResult.Failed_CannotAllocateTriangle;
            MeshResult mr = BaseMesh.InsertTriangle(tid, tv, gid, bUnsafe);
            if (mr == MeshResult.Ok)
                updateTimeStamp(true);
            return mr;
        }

        internal override int add_edge(int vA, int vB, int tA, int tB = InvalidID)
        {
            int eid = BaseMesh.add_edge(vA, vB, tA, tB);
            edges_refcount.allocate_at(eid);
            return eid;
        }

        /// <summary>
        /// Cannot Remove Triangles from SubMesh. Remove the triangle from the BaseMesh and rebuild the SubMesh
        /// </summary>
        /// <param name="tID"></param>
        /// <param name="bRemoveIsolatedVertices"></param>
        /// <param name="bPreserveManifold"></param>
        /// <returns></returns>
        public override MeshResult RemoveTriangle(int tID, bool bRemoveIsolatedVertices = true, bool bPreserveManifold = false)
        {
            MeshResult mr = BaseMesh.RemoveTriangle(tID, bRemoveIsolatedVertices, bPreserveManifold);
            if (mr == MeshResult.Ok) compute(selectedTris);
            return mr;
        }

        /// <summary>
        /// Cannot Remove Vertex from SubMesh. Remove the Vertex from the  BaseMesh and rebuild the 
        /// </summary>
        /// <param name="vID"></param>
        /// <param name="bRemoveAllTriangles"></param>
        /// <param name="bPreserveManifold"></param>
        /// <returns></returns>
        public override MeshResult RemoveVertex(int vID, bool bRemoveAllTriangles = true, bool bPreserveManifold = false)
        {
            MeshResult mr = BaseMesh.RemoveVertex(vID, bRemoveAllTriangles, bPreserveManifold);
            if (mr == MeshResult.Ok) compute(selectedTris);
            return mr;
        }

        // Static Mesh Creators for backweards compatibility

        public static DMesh3 QuickSubmesh(DMesh3 mesh, int[] triangles)
        {
            DSubmesh3 submesh = new DSubmesh3(mesh, triangles);
            return submesh;
        }
        public static DMesh3 QuickSubmesh(DMesh3 mesh, IEnumerable<int> triangles)
        {
            return QuickSubmesh(mesh, triangles.ToArray());
        }
    }


    /// <summary>
    /// Legacy DSubmesh3 class 
    /// </summary>
    [Obsolete("The Legacy DSubmesh3 class is now obsolete and will be removed in a future version")]
    public class DSubmesh3Legacy
    {
        public DMesh3 BaseMesh;
        public DMesh3 SubMesh;
        public MeshComponents WantComponents = MeshComponents.All;
        public bool ComputeTriMaps = false;
        public int OverrideGroupID = -1;
        
        public IndexFlagSet BaseSubmeshV;       // vertices in base mesh that are in submesh
                                                // (we compute this anyway, might as well hang onto it)

        public IndexMap BaseToSubV;             // vertex index map from base to submesh 
        public DVector<int> SubToBaseV;         // vertex index map from submesh to base mesh

        public IndexMap BaseToSubT;             // triangle index map from base to submesh. Only computed if ComputeTriMaps = true.
        public DVector<int> SubToBaseT;         // triangle index map from submesh to base mesh. Only computed if ComputeTriMaps = true.

        // boundary info
        public IndexHashSet BaseBorderE;        // list of internal border edge indices on base mesh. Does not include mesh boundary edges.
        public IndexHashSet BaseBoundaryE;      // list of mesh-boundary edges on base mesh that are in submesh
        public IndexHashSet BaseBorderV;        // list of border vertex indices on base mesh (ie verts of BaseBorderE - does not include mesh boundary vertices)


        public DSubmesh3Legacy(DMesh3 mesh, int[] subTriangles)
        {
            BaseMesh = mesh;
            compute(subTriangles, subTriangles.Length);
        }
        public DSubmesh3Legacy(DMesh3 mesh, IEnumerable<int> subTriangles, int nTriEstimate = 0)
        {
            BaseMesh = mesh;
            compute(subTriangles, nTriEstimate);
        }

        public DSubmesh3Legacy(DMesh3 mesh)
        {
            BaseMesh = mesh;
        }
        public void Compute(int[] subTriangles) {
            compute(subTriangles, subTriangles.Length);
        }
        public void Compute(IEnumerable<int> subTriangles, int nTriEstimate = 0) {
            compute(subTriangles, nTriEstimate);
        }

        public int MapVertexToSubmesh(int base_vID) {
            return BaseToSubV[base_vID];
        }
        public int MapVertexToBaseMesh(int sub_vID) {
            if (sub_vID < SubToBaseV.Length)
                return SubToBaseV[sub_vID];
            return DMesh3.InvalidID;
        }

        public Index2i MapVerticesToSubmesh(Index2i v) {
            return new Index2i(BaseToSubV[v.a], BaseToSubV[v.b]);
        }
        public Index2i MapVerticesToBaseMesh(Index2i v) {
            return new Index2i(MapVertexToBaseMesh(v.a), MapVertexToBaseMesh(v.b));
        }

        public void MapVerticesToSubmesh(int[] vertices)
        {
            for (int i = 0; i < vertices.Length; ++i)
                vertices[i] = BaseToSubV[vertices[i]];
        }


        public int MapEdgeToSubmesh(int base_eid)
        {
            Index2i base_ev = BaseMesh.GetEdgeV(base_eid);
            Index2i sub_ev = MapVerticesToSubmesh(base_ev);
            return SubMesh.FindEdge(sub_ev.a, sub_ev.b);
        }
        public void MapEdgesToSubmesh(int[] edges)
        {
            for (int i = 0; i < edges.Length; ++i)
                edges[i] = MapEdgeToSubmesh(edges[i]);
        }

        public int MapEdgeToBaseMesh(int sub_eid)
        {
            Index2i sub_ev = SubMesh.GetEdgeV(sub_eid);
            Index2i base_ev = MapVerticesToBaseMesh(sub_ev);
            return BaseMesh.FindEdge(base_ev.a, base_ev.b);
        }


        public int MapTriangleToSubmesh(int base_tID)
        {
            if (ComputeTriMaps == false)
                throw new InvalidOperationException("DSubmesh3.MapTriangleToSubmesh: must set ComputeTriMaps = true!");
            return BaseToSubT[base_tID];
        }
        public int MapTriangleToBaseMesh(int sub_tID)
        {
            if (ComputeTriMaps == false)
                throw new InvalidOperationException("DSubmesh3.MapTriangleToBaseMesh: must set ComputeTriMaps = true!");
            if (sub_tID < SubToBaseT.Length)
                return SubToBaseT[sub_tID];
            return DMesh3.InvalidID;
        }

        public void MapTrianglesToSubmesh(int[] triangles)
        {
            if (ComputeTriMaps == false)
                throw new InvalidOperationException("DSubmesh3.MapTrianglesToSubmesh: must set ComputeTriMaps = true!");
            for (int i = 0; i < triangles.Length; ++i)
                triangles[i] = BaseToSubT[triangles[i]];
        }

        public void ComputeBoundaryInfo(int[] subTriangles) {
            ComputeBoundaryInfo(subTriangles, subTriangles.Length);
        }
        public void ComputeBoundaryInfo(IEnumerable<int> triangles, int tri_count_est)
        {
            // set of base-mesh triangles that are in submesh
            IndexFlagSet sub_tris = new IndexFlagSet(BaseMesh.MaxTriangleID, tri_count_est);
            foreach (int ti in triangles)
                sub_tris[ti] = true;

            BaseBorderV = new IndexHashSet();
            BaseBorderE = new IndexHashSet();
            BaseBoundaryE = new IndexHashSet();

            // Iterate through edges in submesh roi on base mesh. If
            // one of the tris of the edge is not in submesh roi, then this
            // is a boundary edge.
            //
            // (edge iteration via triangle iteration processes each internal edge twice...)
            foreach (int ti in triangles) {
                Index3i tedges = BaseMesh.GetTriEdges(ti);
                for ( int j = 0; j < 3; ++j ) {
                    int eid = tedges[j];
                    Index2i tris = BaseMesh.GetEdgeT(eid);
                    if ( tris.b == DMesh3.InvalidID ) {     // this is a boundary edge
                        BaseBoundaryE[eid] = true;

                    } else if (sub_tris[tris.a] != sub_tris[tris.b]) {  // this is a border edge
                        BaseBorderE[eid] = true;
                        Index2i ve = BaseMesh.GetEdgeV(eid);
                        BaseBorderV[ve.a] = true;
                        BaseBorderV[ve.b] = true;
                    } 
                }
            }
            
        }



        // [RMS] estimate can be zero
        void compute(IEnumerable<int> triangles, int tri_count_est )
        {
            int est_verts = tri_count_est / 2;

            SubMesh = new DMesh3(BaseMesh.Components & WantComponents) { axisOrder = BaseMesh.axisOrder };

            BaseSubmeshV = new IndexFlagSet(BaseMesh.MaxVertexID, est_verts);
            BaseToSubV = new IndexMap(BaseMesh.MaxVertexID, est_verts);
            SubToBaseV = new DVector<int>();

            if ( ComputeTriMaps ) {
                BaseToSubT = new IndexMap(BaseMesh.MaxTriangleID, tri_count_est);
                SubToBaseT = new DVector<int>();
            }

            foreach ( int tid in triangles ) {
                if ( ! BaseMesh.IsTriangle(tid) )
                    throw new Exception("DSubmesh3.compute: triangle " + tid + " does not exist in BaseMesh!");
                Index3i base_t = BaseMesh.GetTriangle(tid);
                Index3i new_t = Index3i.Zero;
                int gid = BaseMesh.GetTriangleGroup(tid);

                for ( int j = 0; j < 3; ++j ) {
                    int base_v = base_t[j];
                    int sub_v = -1;
                    if (BaseSubmeshV[base_v] == false) {
                        sub_v = SubMesh.AppendVertex(BaseMesh, base_v);
                        BaseSubmeshV[base_v] = true;
                        BaseToSubV[base_v] = sub_v;
                        SubToBaseV.insert(base_v, sub_v);
                    } else
                        sub_v = BaseToSubV[base_v];
                    new_t[j] = sub_v;
                }

                if (OverrideGroupID >= 0)
                    gid = OverrideGroupID;
                int sub_tid = SubMesh.AppendTriangle(new_t, gid);

                if ( ComputeTriMaps ) {
                    BaseToSubT[tid] = sub_tid;
                    SubToBaseT.insert(tid, sub_tid);
                }
            }
        }

        /// <summary>
        /// Apply any changes to the vertices back top the base mesh
        /// </summary>
        public void ApplyToBaseMesh()
        {
            NewVertexInfo vinfo = new();
            foreach (int vID in SubMesh.VertexIndices())
            {
                if (SubMesh.GetVertex(vID, ref vinfo, true, true, true))
                    BaseMesh.SetVertex(vID, vinfo, true, true, true);
            }
        }

        // Static Mesh Creators

        public static DMesh3 QuickSubmesh(DMesh3 mesh, int[] triangles) {
            DSubmesh3Legacy submesh = new DSubmesh3Legacy(mesh, triangles);
            return submesh.SubMesh;
        }
        public static DMesh3 QuickSubmesh(DMesh3 mesh, IEnumerable<int> triangles) {
            return QuickSubmesh(mesh, triangles.ToArray());
        }

    }
}
