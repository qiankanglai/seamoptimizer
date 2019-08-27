using System.Collections;
using System.Collections.Generic;
using UnityEngine;

class SeamOptimizer
{
	public class so_texel_t
	{
		public int x, y;
		public so_texel_t(int _x, int _y)
		{
			x = _x;
			y = _y;
		}
	}

	class so_texel_cmper : IComparer<so_texel_t>
	{
		public int Compare(so_texel_t lt, so_texel_t rt)
		{
			if (lt.y < rt.y) return -1;
			if (lt.y > rt.y) return 1;
			if (lt.x < rt.x) return -1;
			if (lt.x > rt.x) return 1;
			return 0;
		}
	}
    static so_texel_cmper so_texel_cmp = new so_texel_cmper();

	public class so_bilinear_sample_t
	{
		public so_texel_t[] texels;
		public float[] weights;
		public so_bilinear_sample_t()
		{
			texels = new so_texel_t[4];
			weights = new float[4];
		}
	}

	public class so_stitching_point_t
	{
		public so_bilinear_sample_t[] sides;
		public so_stitching_point_t()
		{
			sides = new so_bilinear_sample_t[2];
			for (int i = 0; i < 2; i++)
				sides[i] = new so_bilinear_sample_t();
		}
	}

	public class so_texel_set_t
	{
		public so_texel_t[] texels;
		public int count;
		public int capacity;
	}

	static bool so_accumulate_texel(float[] sums, int x, int y, Color[] data, int w, int h, int c)
	{
		bool exists = false;
		for (int i = 0; i < c; i++)
		{
			float v = data[y * w + x][i];
			sums[i] += v;
			exists |= v > 0.0f;
		}
		return exists;
	}

	static void so_fill_with_closest(int x, int y, Color[] data, int w, int h, int c, int depth = 2)
	{
		Debug.Assert(c <= 4);
        // Kanglai: if we have data, nothing happens
		for (int i = 0; i < c; i++)
			if (data[y * w + x][i] > 0.0f)
				return;

		float[] sums = { 0, 0, 0, 0 };
		int n = 0;
        // populate neighbour pixels with data to self
		if (x > 0 && so_accumulate_texel(sums, x - 1, y, data, w, h, c)) n++;
		if (x + 1 < w && so_accumulate_texel(sums, x + 1, y, data, w, h, c)) n++;
		if (y > 0 && so_accumulate_texel(sums, x, y - 1, data, w, h, c)) n++;
		if (y + 1 < h && so_accumulate_texel(sums, x, y + 1, data, w, h, c)) n++;

		if (n == 0 && depth != 0)
		{
			--depth;
			if (x > 0)
			{
				so_fill_with_closest(x - 1, y, data, w, h, c, depth);
				if (so_accumulate_texel(sums, x - 1, y, data, w, h, c)) n++;
			}
			if (x + 1 < w)
			{
				so_fill_with_closest(x + 1, y, data, w, h, c, depth);
				if (so_accumulate_texel(sums, x + 1, y, data, w, h, c)) n++;
			}
			if (y > 0)
			{
				so_fill_with_closest(x, y - 1, data, w, h, c, depth);
				if (so_accumulate_texel(sums, x, y - 1, data, w, h, c)) n++;
			}
			if (y + 1 < h)
			{
				so_fill_with_closest(x, y + 1, data, w, h, c, depth);
				if (so_accumulate_texel(sums, x, y + 1, data, w, h, c)) n++;
			}
		}
        // average
		if (n != 0)
		{
			float ni = 1.0f / (float)n;
			for (int i = 0; i < c; i++)
				data[y * w + x][i] = sums[i] * ni;
		}
	}

	static int so_texel_hash(so_texel_t texel, int capacity)
	{
		int hash = (texel.y * 104173 + texel.x * 86813) % capacity;
		if (hash < 0)
			hash += capacity;
		return hash;
	}

	static void so_texel_set_add(so_texel_set_t set, so_texel_t[] texels, int entries, int arrayLength = 0)
	{
		if (set.count + entries > set.capacity * 3 / 4) // leave some free space to avoid having many collisions
		{
			int newCapacity = set.capacity > 64 ? set.capacity * 2 : 64;
			while (set.count + entries > newCapacity * 3 / 4)
				newCapacity *= 2;

			so_texel_t[] newTexels = new so_texel_t[newCapacity];

			for (int i = 0; i < newCapacity; i++)
				newTexels[i] = new so_texel_t(-1, -1);

			if (set.texels != null)
			{
				for (int i = 0; i < set.capacity; i++) // rehash all old texels
				{
					if (set.texels[i].x != -1)
					{
						int hash = so_texel_hash(set.texels[i], newCapacity);
						while (newTexels[hash].x != -1) // collisions
							hash = (hash + 1) % newCapacity;
						newTexels[hash] = set.texels[i];
					}
				}
			}

			set.texels = newTexels;
			set.capacity = newCapacity;
		}

		if (arrayLength == 0)
			arrayLength = entries;

		for (int i = 0; i < arrayLength; i++)
		{
			if (texels[i].x != -1)
			{
				int hash = so_texel_hash(texels[i], set.capacity);
				while (set.texels[hash].x != -1) // collisions
				{
					if (set.texels[hash].x == texels[i].x && set.texels[hash].y == texels[i].y)
						break; // texel is already in the set
					hash = (hash + 1) % set.capacity;
				}

				if (set.texels[hash].x == -1)
				{
					set.texels[hash] = texels[i];
					set.count++;
				}
			}
		}
	}

	static bool so_texel_set_contains(so_texel_set_t set, so_texel_t texel)
	{
		int hash = so_texel_hash(texel, set.capacity);
		while (set.texels[hash].x != -1) // entries with same hash
		{
			if (set.texels[hash].x == texel.x && set.texels[hash].y == texel.y)
				return true; // texel is already in the set
			hash = (hash + 1) % set.capacity;
		}
		return false;
	}

	public class so_stitching_points_t
	{
		public so_stitching_point_t[] points;
		public int count;
		public int capacity;
		public so_stitching_points_t(int _capacity)
		{
			points = new so_stitching_point_t[_capacity];
			capacity = _capacity;
			count = 0;
		}
	}

	static void so_stitching_points_add(so_stitching_points_t points, so_stitching_point_t point)
	{
		Debug.Assert(points.count < points.capacity);
		points.points[points.count++] = point;
	}

	static void so_stitching_points_append(so_stitching_points_t points, so_stitching_points_t other)
	{
		so_stitching_point_t[] newPoints = new so_stitching_point_t[points.capacity + other.capacity];
		for (int i = 0; i < points.count; i++)
			newPoints[i] = points.points[i];
		for (int i = 0; i < other.count; i++)
			newPoints[i + points.count] = other.points[i];
		points.points = newPoints;
		points.capacity = points.capacity + other.capacity;
		points.count = points.count + other.count;
	}

	public class so_seam_t
	{
		public int x_min, y_min, x_max, y_max;
		public so_texel_set_t texels;
		public so_stitching_points_t stitchingPoints;
	}

	static void so_seam_add(so_seam_t seam, so_stitching_point_t point)
	{
		for (int side = 0; side < 2; side++)
		{
			for (int texel = 0; texel < 4; texel++)
			{
				so_texel_t t = point.sides[side].texels[texel];
				seam.x_min = Mathf.Min(t.x, seam.x_min);
				seam.y_min = Mathf.Min(t.y, seam.y_min);
				seam.x_max = Mathf.Max(t.x, seam.x_min);
				seam.y_max = Mathf.Max(t.y, seam.y_min);
			}
			so_texel_set_add(seam.texels, point.sides[side].texels, 4);
		}

		so_stitching_points_add(seam.stitchingPoints, point);
	}

	static bool so_seams_intersect(so_seam_t a, so_seam_t b)
	{
#if false
		// don't join intersected seams
		return false;
#endif
		// compare bounding boxes first
		if (a.x_min > b.x_max || b.x_min >= a.x_max ||
			a.y_min > b.y_max || b.y_min >= a.y_max)
			return false;

		// bounds intersect . check each individual texel for intersection
		if (a.texels.capacity > b.texels.capacity) // swap so that we always loop over the smaller set
		{
			so_seam_t tmp = a;
			a = b;
			b = tmp;
		}

		for (int i = 0; i < a.texels.capacity; i++)
			if (a.texels.texels[i].x != -1)
				if (so_texel_set_contains(b.texels, a.texels.texels[i]))
					return true;
		return false;
	}

	static void so_seams_in_place_merge(so_seam_t dst, so_seam_t src)
	{
		// expand bounding box
		dst.x_min = Mathf.Min(src.x_min, dst.x_min);
		dst.y_min = Mathf.Min(src.y_min, dst.y_min);
		dst.x_max = Mathf.Max(src.x_max, dst.x_max);
		dst.y_max = Mathf.Max(src.y_max, dst.y_max);

		// insert src elements
		so_texel_set_add(dst.texels, src.texels.texels, src.texels.count, src.texels.capacity);
		so_stitching_points_append(dst.stitchingPoints, src.stitchingPoints);
	}

	static void so_seams_add_seam(List<so_seam_t> seams, Vector2 a0, Vector2 a1, Vector2 b0, Vector2 b1, Color[] data, int w, int h, int c)
	{
		Vector2 s = new Vector2(w-1, h-1);
        a0 = Extensions.Multiply(a0,s);
        a1 = Extensions.Multiply(a1,s);
        b0 = Extensions.Multiply(b0,s);
        b1 = Extensions.Multiply(b1,s);
		Vector2 ad = a1 - a0;
		Vector2 bd = b1 - b0;
        float l = Mathf.Max(ad.magnitude, bd.magnitude);
		int iterations = (int)(l * 5.0f); // TODO: is this the best value?
		float step = 1.0f / iterations;

		so_seam_t currentSeam = new so_seam_t();
		currentSeam.x_min = w; currentSeam.y_min = h;
		currentSeam.x_max = 0; currentSeam.y_max = 0;

		// so_seam_alloc
		currentSeam.stitchingPoints = new so_stitching_points_t(iterations + 1);
		currentSeam.texels = new so_texel_set_t();

		for (int i = 0; i <= iterations; i++)
		{
			float t = i * step;
			Vector2 a = a0 + ad * t;
			Vector2 b = b0 + bd * t;
            // Kanglai: shouldn't be Round here
            int ax = Mathf.FloorToInt(a.x), ay = Mathf.FloorToInt(a.y);
            int bx = Mathf.FloorToInt(b.x), by = Mathf.FloorToInt(b.y);
			float au = a.x - ax, av = a.y - ay, nau = 1.0f - au, nav = 1.0f - av;
			float bu = b.x - bx, bv = b.y - by, nbu = 1.0f - bu, nbv = 1.0f - bv;

			so_texel_t ta0 = new so_texel_t(ax, ay);
			so_texel_t ta1 = new so_texel_t(Mathf.Min(ax + 1, w - 1), ay);
			so_texel_t ta2 = new so_texel_t(ax, Mathf.Min(ay + 1, h - 1));
			so_texel_t ta3 = new so_texel_t(Mathf.Min(ax + 1, w - 1), Mathf.Min(ay + 1, h - 1));

			so_texel_t tb0 = new so_texel_t(bx, by);
			so_texel_t tb1 = new so_texel_t(Mathf.Min(bx + 1, w - 1), by);
			so_texel_t tb2 = new so_texel_t(bx, Mathf.Min(by + 1, h - 1));
			so_texel_t tb3 = new so_texel_t(Mathf.Min(bx + 1, w - 1), Mathf.Min(by + 1, h - 1));
			/*
			so_fill_with_closest(ta0.x, ta0.y, data, w, h, c);
			so_fill_with_closest(ta1.x, ta1.y, data, w, h, c);
			so_fill_with_closest(ta2.x, ta2.y, data, w, h, c);
			so_fill_with_closest(ta3.x, ta3.y, data, w, h, c);

			so_fill_with_closest(tb0.x, tb0.y, data, w, h, c);
			so_fill_with_closest(tb1.x, tb1.y, data, w, h, c);
			so_fill_with_closest(tb2.x, tb2.y, data, w, h, c);
			so_fill_with_closest(tb3.x, tb3.y, data, w, h, c);
			*/
			so_stitching_point_t sp = new so_stitching_point_t();
			sp.sides[0].texels[0] = ta0;
			sp.sides[0].texels[1] = ta1;
			sp.sides[0].texels[2] = ta2;
			sp.sides[0].texels[3] = ta3;

			sp.sides[0].weights[0] = nau * nav;
			sp.sides[0].weights[1] = au * nav;
			sp.sides[0].weights[2] = nau * av;
			sp.sides[0].weights[3] = au * av;

			sp.sides[1].texels[0] = tb0;
			sp.sides[1].texels[1] = tb1;
			sp.sides[1].texels[2] = tb2;
			sp.sides[1].texels[3] = tb3;

			sp.sides[1].weights[0] = nbu * nbv;
			sp.sides[1].weights[1] = bu * nbv;
			sp.sides[1].weights[2] = nbu * bv;
			sp.sides[1].weights[3] = bu * bv;

			so_seam_add(currentSeam, sp);
		}
		so_seam_t dstSeam = null;
		for (int i = 0; i < seams.Count; i++)
		{
            so_seam_t seam = seams[i];
            // Kanglai: we put seams that intersects together
			if (so_seams_intersect(currentSeam, seam))
			{
				if (dstSeam == null) // found a seam that the edge is connected to . add current edge to that seam
				{
					so_seams_in_place_merge(seam, currentSeam);
					dstSeam = seam;
				}
				else // found another seam that the edge is connected to . merge those seams
				{
					so_seams_in_place_merge(dstSeam, seam);

					// remove current seam from seams
					seams.Remove(seam);
					i--;
				}
			}
		}
		if (dstSeam == null) // did not find a seam that the edge is connected to . make a new one
        {
            seams.Add(currentSeam);
		}
	}

	static bool so_should_optimize(Vector3[] pos, int tria, int trib, float cosThreshold)
	{
		Vector3 n0 = Vector3.Cross(pos[tria + 1] - pos[tria + 0], pos[tria + 2] - pos[tria + 0]).normalized;
		Vector3 n1 = Vector3.Cross(pos[trib + 1] - pos[trib + 0], pos[trib + 2] - pos[trib + 0]).normalized;
		return Mathf.Abs(Vector3.Dot(n0, n1)) > cosThreshold;
	}

    public static List<so_seam_t> so_seams_find(Mesh mesh, Color[] data, Vector4 lightmapST,
        int w, int h, int c = 3, float cosNormalThreshold = 0.99f)
    {
        // get mesh data
        var lightmapS = new Vector2(lightmapST.x, lightmapST.y);
        var lightmapT = new Vector2(lightmapST.z, lightmapST.w);
        Vector3[] _vertices = mesh.vertices;
        Vector2[] _uv = mesh.uv2;
        for (int i = 0; i < _uv.Length; i++)
            _uv[i] = Extensions.Multiply(_uv[i], lightmapS) + lightmapT;
        int[] _triangles = mesh.triangles;

        // flat data
        int vertices = _triangles.Length;
        Vector3[] pos = new Vector3[vertices];
        Vector2[] uv = new Vector2[vertices];
        for (int i = 0; i < vertices; i++)
        {
            pos[i] = _vertices[_triangles[i]];
            uv[i] = _uv[_triangles[i]];
        }

        // mesh.bounds may be modified
        Vector3 bbmin = new Vector3(float.MaxValue, float.MaxValue, float.MaxValue);
        Vector3 bbmax = new Vector3(float.MinValue, float.MinValue, float.MinValue);
        int[] hashmap = new int[vertices * 2];
        for (int i = 0; i < vertices; i++)
        {
            bbmin = Extensions.Min(bbmin, pos[i]);
            bbmax = Extensions.Max(bbmax, pos[i]);
            hashmap[i * 2 + 0] = -1;
            hashmap[i * 2 + 1] = -1;
        }

        Vector3 bbscale = new Vector3(15.9f / bbmax.x, 15.9f / bbmax.y, 15.9f / bbmax.z);

        List<so_seam_t> seams = new List<so_seam_t>();

        for (int i0 = 0; i0 < vertices; i0++)
        {
            int tri = i0 - (i0 % 3);
            int i1 = tri + ((i0 + 1) % 3);
            //int i2 = tri + ((i0 + 2) % 3);
            Vector3 p = Extensions.Multiply(pos[i0] - bbmin, bbscale);
            int hash = (281 * (int)p.x + 569 * (int)p.y + 1447 * (int)p.z) % (vertices * 2);
            while (hashmap[hash] >= 0)
            {
                int oi0 = hashmap[hash];
                if ((pos[oi0] - pos[i0]).sqrMagnitude < 0.0000001f)
                {
                    int otri = oi0 - (oi0 % 3);
                    int oi1 = otri + ((oi0 + 1) % 3);
                    int oi2 = otri + ((oi0 + 2) % 3);
                    if ((pos[oi1] - pos[i1]).sqrMagnitude < 0.0000001f && so_should_optimize(pos, tri, otri, cosNormalThreshold))
                        so_seams_add_seam(seams, uv[i0], uv[i1], uv[oi0], uv[oi1], data, w, h, c);
                    //else if (SO_EQUAL(oi1, i2) && so_should_optimize(pos + tri, pos + otri, cosNormalThreshold)) // this will already be detected by the other side of the seam!
                    //	so_seams_add_seam(&seams, uv[i0], uv[i2], uv[oi0], uv[oi1], data, w, h, c);
                    else if ((pos[oi2] - pos[i1]).sqrMagnitude < 0.0000001f && so_should_optimize(pos, tri, otri, cosNormalThreshold))
                        so_seams_add_seam(seams, uv[i0], uv[i1], uv[oi0], uv[oi2], data, w, h, c);
                    //break;
                }
                if (++hash == vertices * 2)
                    hash = 0;
            }
            hashmap[hash] = i0;
        }

        return seams;
    }

	public class so_sparse_entry_t
	{
		public int index;
		public float value;

		public so_sparse_entry_t(int _index, int _value)
		{
			index = _index;
			value = _value;
		}
	}

	class so_sparse_entry_cmper : IComparer<so_sparse_entry_t>
	{
		public int Compare(so_sparse_entry_t ae, so_sparse_entry_t be)
		{
			return ae.index - be.index;
		}
	}
    static so_sparse_entry_cmper so_sparse_entry_cmp = new so_sparse_entry_cmper();

	public class so_sparse_entries_t
	{
		public so_sparse_entry_t[] entries;
		public int count;
		public int capacity;

		public so_sparse_entries_t(int _capacity)
		{
			entries = new so_sparse_entry_t[_capacity];
			for (int i = 0; i < _capacity; i++)
				entries[i] = new so_sparse_entry_t(-1, 0);
			capacity = _capacity;
			count = 0;
		}
	}

	static void so_sparse_matrix_add(so_sparse_entries_t matrix, int index, float value)
	{
		if (matrix.count == matrix.capacity)
		{
			int newCapacity = matrix.capacity * 2;
			if (newCapacity < 64)
				newCapacity = 64;
			so_sparse_entry_t[] newEntries = new so_sparse_entry_t[newCapacity];
			for (int i = 0; i < matrix.count; i++)
				newEntries[i] = matrix.entries[i];
			for (int i = matrix.count; i < newCapacity; i++)
				newEntries[i] = new so_sparse_entry_t(-1, 0);
			matrix.entries = newEntries;
			matrix.capacity = newCapacity;
		}

		int entryIndex = matrix.count++;
		matrix.entries[entryIndex].index = index;
		matrix.entries[entryIndex].value = value;
	}

	static void so_sparse_matrix_sort(so_sparse_entries_t matrix)
	{
        System.Array.Sort(matrix.entries, 0, matrix.count, so_sparse_entry_cmp);
	}

	static bool so_sparse_matrix_advance_to_index(so_sparse_entries_t matrix, ref int position, int index, ref float outValue)
	{
		int localPosition = position;
		while (localPosition < matrix.count && matrix.entries[localPosition].index < index)
			++localPosition;
		position = localPosition;

		if (localPosition < matrix.count && matrix.entries[localPosition].index == index)
		{
			outValue = matrix.entries[localPosition].value;
			return true;
		}

		return false;
	}

	static int so_sparse_entry_hash(int entryIndex, int capacity)
	{
		int hash = (entryIndex * 104173) % capacity;
		if (hash < 0)
			hash += capacity;
		return hash;
	}

	static so_sparse_entry_t so_sparse_entry_set_get_or_add(so_sparse_entries_t set, int index)
	{
		if (set.count + 1 > set.capacity * 3 / 4) // leave some free space to avoid having many collisions
		{
			int newCapacity = set.capacity >= 64 ? set.capacity * 2 : 64;
			so_sparse_entry_t[] newEntries = new so_sparse_entry_t[newCapacity];
			for (int i = 0; i < newCapacity; i++)
				newEntries[i] = new so_sparse_entry_t(-1, 0);

			for (int i = 0; i < set.capacity; i++) // rehash all old entries
			{
				if (set.entries[i].index != -1)
				{
					int hash_ = so_sparse_entry_hash(set.entries[i].index, newCapacity);
					while (newEntries[hash_].index != -1) // collisions
						hash_ = (hash_ + 1) % newCapacity;
					newEntries[hash_] = set.entries[i];
				}
			}
			set.entries = newEntries;
			set.capacity = newCapacity;
		}

		int hash = so_sparse_entry_hash(index, set.capacity);
		while (set.entries[hash].index != -1) // collisions
		{
			if (set.entries[hash].index == index)
				return set.entries[hash]; // entry is already in the set
			hash = (hash + 1) % set.capacity;
		}

		if (set.entries[hash].index == -1) // make new entry
		{
			set.entries[hash].index = index;
			set.entries[hash].value = 0.0f;
			set.count++;
			return set.entries[hash];
		}

		Debug.Assert(false);
		return null; // shouldn't happen
	}

	static so_sparse_entries_t so_matrix_At_times_A(float[] A, int[] sparseIndices, int maxRowIndices, int m, int n)
	{
		so_sparse_entries_t AtA = new so_sparse_entries_t((n / 16) * (n / 16));

		// compute lower left triangle only since the result is symmetric
		for (int k = 0; k < m; k++)
		{
			int srcPtr = k * maxRowIndices;
			int indexPtr = k * maxRowIndices;
			for (int i = 0; i < maxRowIndices; i++)
			{
				int index_i = sparseIndices[indexPtr + i];
				if (index_i < 0) break;
				float v = A[srcPtr + i];
				//float *dstPtr = AtA + index_i * n;
				for (int j = 0; j < maxRowIndices; j++)
				{
					int index_j = sparseIndices[indexPtr + j];
					if (index_j < 0) break;
					//dstPtr[index_j] += v * srcPtr[j];
					int index = index_i * n + index_j;

					so_sparse_entry_t entry = so_sparse_entry_set_get_or_add(AtA, index);
					entry.value += v * A[srcPtr + j];
				}
			}
		}

		// compaction step (make a compact array from the scattered hash set values)
		for (int i = 0, j = 0; i < AtA.capacity; i++)
			if (AtA.entries[i].index != -1)
				AtA.entries[j++] = AtA.entries[i];

		// sort by index . this is a sparse matrix now
		so_sparse_matrix_sort(AtA);

		return AtA;
	}

	static void so_matrix_At_times_b(float[] A, int m, int n, float[] b, float[] Atb, int[] sparseIndices, int maxRowIndices)
	{
		for (int i = 0; i < n; i++)
			Atb[i] = 0;
		for (int j = 0; j < m; j++)
		{
			int rowIndices = j * maxRowIndices;
			for (int i = 0; i < maxRowIndices; i++)
			{
				int index = sparseIndices[rowIndices + i];
				if (index < 0) break;
				Atb[index] += A[j * maxRowIndices + i] * b[j];
			}
		}
	}

	static so_sparse_entries_t so_matrix_cholesky_prepare(so_sparse_entries_t AtA, int n)
	{
		// dense
		//for (int i = 0; i < n; i++)
		//{
		//	float *a = L + i * n;
		//	for (int j = 0; j <= i; j++)
		//	{
		//		float *b = L + j * n;
		//		float sum = A[i * n + j];// + (i == j ? 0.0001 : 0.0); // some regularization
		//		for (int k = 0; k < j; k++)
		//			sum -= a[k] * b[k];
		//		if (i > j)
		//			a[j] = sum / b[j];
		//		else // i == j
		//		{
		//			if (sum <= 0.0)
		//				return false;
		//			a[i] = sqrtf(sum);
		//		}
		//	}
		//}

		// sparse
		int[] indices_i = new int[n];
		float[] row_i = new float[n];
		float[] invDiag = new float[n];

		so_sparse_entries_t L = new so_sparse_entries_t((n / 16) * (n / 16));

		int AtAindex = 0;
		for (int i = 0; i < n; i++)
		{
			int index_i_count = 0;

			int row_j_index = 0;
			for (int j = 0; j <= i; j++)
			{
				//float sum = A[i * n + j]; // + (i == j ? 0.0001 : 0.0); // regularization
				int index = i * n + j;
				float sum = 0.0f;
				so_sparse_matrix_advance_to_index(AtA, ref AtAindex, index, ref sum);

				for (int k = 0; k < index_i_count; k++)
				{
					int index_i = indices_i[k];
					float Lvalue = 0;
					if (so_sparse_matrix_advance_to_index(L, ref row_j_index, j * n + index_i, ref Lvalue))
						sum -= row_i[index_i] * Lvalue;
				}

				if (i == j)
				{
					if (sum <= 0.0f)
					{
						return L;
					}
					invDiag[i] = 1.0f / Mathf.Sqrt(sum);
				}

				if (Mathf.Abs(sum) > 0.00001f)
				{
					row_i[j] = sum * invDiag[j];
					indices_i[index_i_count++] = j;
					so_sparse_matrix_add(L, index, row_i[j]);
				}
				else
					row_i[j] = 0.0f;
			}
		}

		return L;
	}

	static void so_matrix_cholesky_solve(so_sparse_entries_t Lrows, so_sparse_entries_t Lcols, float[] x, float[] b, int n)
	{
		float[] y = new float[n];

		// L * y = b
		int Lindex = 0;
		for (int i = 0; i < n; i++)
		{
			float sum = b[i];
			while (Lindex < Lrows.count && Lrows.entries[Lindex].index < i * (n + 1))
			{
				sum -= Lrows.entries[Lindex].value * y[Lrows.entries[Lindex].index - i * n];
				++Lindex;
			}
			Debug.Assert(Lrows.entries[Lindex].index == i * (n + 1));
			y[i] = sum / Lrows.entries[Lindex].value;
			++Lindex;
		}

		// L' * x = y
		Lindex = Lcols.count - 1;
		for (int i = n - 1; i >= 0; i--)
		{
			float sum = y[i];
			while (Lindex >= 0 && Lcols.entries[Lindex].index > i * (n + 1))
			{
				sum -= Lcols.entries[Lindex].value * x[Lcols.entries[Lindex].index - i * n];
				--Lindex;
			}
			Debug.Assert(Lcols.entries[Lindex].index == i * (n + 1));
			x[i] = sum / Lcols.entries[Lindex].value;
			--Lindex;
		}
	}

	public static bool so_seam_optimize(so_seam_t seam, Color[] data, int w, int h, int c = 3, float lambda = 0.1f)
	{
		so_texel_set_t texels = seam.texels;
		so_stitching_points_t stitchingPoints = seam.stitchingPoints;

		int m = stitchingPoints.count;
		int n = texels.count;

		so_texel_t[] texelsFlat = new so_texel_t[n];
		for (int i = 0; i < n; i++)
			texelsFlat[i] = new so_texel_t(0, 0);

		float[] A = new float[(m + n) * 8];

		int[] AsparseIndices = new int[(m + n) * 8];

		float[] b = new float[m + n];

		float[] Atb = new float[n];

		float[] x = new float[n];

		for (int i = 0, j = 0; i < texels.capacity && j < n; i++)
			if (texels.texels[i].x != -1)
			{
				texelsFlat[j++] = texels.texels[i];
			}

        System.Array.Sort(texelsFlat, so_texel_cmp);

		int r = 0;
		for (int i = 0; i < m; i++)
		{
			int[] column0 = { 0, 0, 0, 0 };
			int[] column1 = { 0, 0, 0, 0 };
			bool side0valid = false, side1valid = false;
			for (int k = 0; k < 4; k++)
			{
				so_texel_t t0 = stitchingPoints.points[i].sides[0].texels[k];
				so_texel_t t1 = stitchingPoints.points[i].sides[1].texels[k];
                column0[k] = System.Array.BinarySearch(texelsFlat, t0, so_texel_cmp);
                column1[k] = System.Array.BinarySearch(texelsFlat, t1, so_texel_cmp);

				if (column0[k] == -1) { side0valid = false; break; }
				if (column1[k] == -1) { side1valid = false; break; }

				// test for validity of stitching point
				for (int ci = 0; ci < c; ci++)
				{
					side0valid |= data[t0.y * w + t0.x][ci] > 0.0f;
					side1valid |= data[t1.y * w + t1.x][ci] > 0.0f;
				}
			}

			if (side0valid && side1valid)
			{
				for (int k = 0; k < 4; k++)
				{
					A[r * 8 + k * 2 + 0] = stitchingPoints.points[i].sides[0].weights[k];
					AsparseIndices[r * 8 + k * 2 + 0] = column0[k];
					A[r * 8 + k * 2 + 1] = -stitchingPoints.points[i].sides[1].weights[k];
					AsparseIndices[r * 8 + k * 2 + 1] = column1[k];
				}
				r++;
			}
		}

		m = r;

		// add error terms for deviation from original pixel value (scaled by lambda)
		for (int i = 0; i < n; i++)
		{
			A[(m + i) * 8] = lambda;
			AsparseIndices[(m + i) * 8 + 0] = i;
			AsparseIndices[(m + i) * 8 + 1] = -1;
		}

		so_sparse_entries_t AtA = so_matrix_At_times_A(A, AsparseIndices, 8, m + n, n);
		so_sparse_entries_t L = so_matrix_cholesky_prepare(AtA, n);

		if (L.count == 0)
		{
			return false; // Cholesky decomposition failed
		}

		so_sparse_entries_t Lcols = new so_sparse_entries_t(L.count);
		for (int i = 0; i < L.count; i++)
			so_sparse_matrix_add(Lcols, (L.entries[i].index % n) * n + (L.entries[i].index / n), L.entries[i].value);
		so_sparse_matrix_sort(Lcols);

		// solve each color channel independently
		for (int ci = 0; ci < c; ci++)
		{
			for (int i = 0; i < n; i++)
				b[m + i] = lambda * data[texelsFlat[i].y * w + texelsFlat[i].x][ci];

			so_matrix_At_times_b(A, m + n, n, b, Atb, AsparseIndices, 8);
			so_matrix_cholesky_solve(L, Lcols, x, Atb, n);

			// write out results
			for (int i = 0; i < n; i++)
				data[texelsFlat[i].y * w + texelsFlat[i].x][ci] = x[i];
		}

		return true;
	}
}