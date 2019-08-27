using System.Collections;
using System.Collections.Generic;
using System.IO;
using UnityEngine;
using UnityEditor;

public partial class SeamOptimizerGUI
{
	const float cosNormalThreshold = -1f;//0.99f;
	const float lambda = 2f;
	const string optimizedSuffix = "_optimized";

    class LightmappedMesh
    {
        public Mesh mesh;
        public Vector4 lightmapST;
    }
    class LightmapMapping
    {
        public Texture2D lightmap;
        public List<LightmappedMesh> meshes;
    }

    static LightmapMapping[] GatherLightmaps()
    {
        var lightmaps = LightmapSettings.lightmaps;
        LightmapMapping[] lighmapMappings = new LightmapMapping[lightmaps.Length];
        for (int i = 0; i < lightmaps.Length; i++)
        {
            lighmapMappings[i] = new LightmapMapping();
            lighmapMappings[i].lightmap = lightmaps[i].lightmapColor;
            lighmapMappings[i].meshes = new List<LightmappedMesh>();
        }

        var meshRenderers = Object.FindObjectsOfType<MeshRenderer>();
		foreach(var meshRenderer in meshRenderers)
		{
			if (!meshRenderer.gameObject.activeInHierarchy)
				continue;
			if (meshRenderer.lightmapIndex < 0)
				continue;
			var meshFilter = meshRenderer.GetComponent<MeshFilter>();
			if (meshFilter == null)
				continue;
			var sharedMesh = meshFilter.sharedMesh;
			if (sharedMesh == null)
				continue;
			lighmapMappings[meshRenderer.lightmapIndex].meshes.Add(new LightmappedMesh()
			{
				mesh = sharedMesh,
				lightmapST = meshRenderer.lightmapScaleOffset
			});
		}
		return lighmapMappings;

	}

	static bool CheckNonOptimizedLightmaps()
	{
		var lightmaps = LightmapSettings.lightmaps;
		foreach(var lightmap in lightmaps)
		{
			var path = AssetDatabase.GetAssetPath(lightmap.lightmapColor);
			if (path.Contains(optimizedSuffix))
			{
				Debug.LogError("Optimized Lightmaps found. Please Apply Non Optimize Lightmaps first!");
				return false;
			}
		}
		return true;
	}

    static Color[] ReadPixels(Texture2D tex)
    {
		var mat = new Material(Shader.Find("Hidden/DecodeLightmap"));
        var rt = new RenderTexture(tex.width, tex.height, 0, RenderTextureFormat.ARGBFloat);
        Graphics.Blit(tex, rt, mat);

        var tex2 = new Texture2D(tex.width, tex.height, TextureFormat.RGBAFloat, false);
        RenderTexture.active = rt;
        tex2.ReadPixels(new Rect(0, 0, tex.width, tex.height), 0, 0);
        RenderTexture.active = null;

        var data = tex2.GetPixels();
        Object.DestroyImmediate(rt);
        Object.DestroyImmediate(tex2);
		Object.DestroyImmediate(mat);
        return data;
    }

	[MenuItem("SO/Optimize Lightmaps")]
	static void OptimizeLightmaps()
	{
		if (!CheckNonOptimizedLightmaps())
			return;

		var lightmapMappings = GatherLightmaps();
		foreach (var lightmapMapping in lightmapMappings)
		{
			var lightmap = lightmapMapping.lightmap;
			var pixels = ReadPixels(lightmap);

			List<SeamOptimizer.so_seam_t> seams = new List<SeamOptimizer.so_seam_t>();
			foreach (var lightmapMesh in lightmapMapping.meshes)
			{
				seams.AddRange(SeamOptimizer.so_seams_find(lightmapMesh.mesh, pixels, lightmapMesh.lightmapST, lightmap.width, lightmap.height, 3, cosNormalThreshold));
			}
#if false
			var seam = seams[4];
#else
			foreach (var seam in seams)
#endif
			{
				if (!SeamOptimizer.so_seam_optimize(seam, pixels, lightmap.width, lightmap.height, 3, lambda))
					Debug.LogWarning("Seam Optimize Failed");
			}

			var result = new Texture2D(lightmap.width, lightmap.height, TextureFormat.RGBAFloat, false, true);
			result.SetPixels(pixels);
			result.Apply();
			var filename = AssetDatabase.GetAssetPath(lightmap);
			var newfilename = Path.GetDirectoryName(filename) + "/" + Path.GetFileNameWithoutExtension(filename) + optimizedSuffix + ".exr";
			File.WriteAllBytes(newfilename, result.EncodeToEXR(Texture2D.EXRFlags.CompressZIP));

			Object.DestroyImmediate(result);
			Debug.Log("Seam Optimize Done");
			AssetDatabase.Refresh();

			bool needApply = false;
			var tiOld = AssetImporter.GetAtPath(filename) as TextureImporter;
			var tiNew = AssetImporter.GetAtPath(newfilename) as TextureImporter;
			if (tiNew.textureType != TextureImporterType.Lightmap)
			{
				tiNew.textureType = TextureImporterType.Lightmap;
				needApply = true;
			}
			if (tiNew.textureCompression != tiOld.textureCompression)
			{
				tiNew.textureCompression = tiOld.textureCompression;
				needApply = true;
			}
			if (tiNew.anisoLevel != tiOld.anisoLevel)
			{
				tiNew.anisoLevel = tiOld.anisoLevel;
				needApply = true;
			}
			if (needApply)
				tiNew.SaveAndReimport();
		}
	}

    [MenuItem("SO/Apply Optimized Lightmaps")]
    static void ApplyOptimizedLightmaps()
    {
        ApplyLightmaps(true);
    }

    [MenuItem("SO/Apply NonOptimized Lightmaps")]
    static void ApplyNonOptimizedLightmaps()
    {
        ApplyLightmaps(false);
    }

    static void ApplyLightmaps(bool optimized)
    {
        var lightmaps = LightmapSettings.lightmaps;
        foreach(var lightmap in lightmaps)
        {
            var tex = lightmap.lightmapColor;
            if (tex == null)
                continue;
            var optimized_ = tex.name.Contains(optimizedSuffix);
            if (optimized == optimized_)
                continue;
            var file = AssetDatabase.GetAssetPath(tex);
            var folder = Path.GetDirectoryName(file);
            var filename = Path.GetFileNameWithoutExtension(file).Replace(optimizedSuffix, "");
            if (optimized)
                filename += optimizedSuffix;
            file = folder + "/" + filename + ".exr";
            var tex2 = AssetDatabase.LoadAssetAtPath<Texture2D>(file);
            if (tex2 == null)
                continue;
            lightmap.lightmapColor = tex2;
        }
        LightmapSettings.lightmaps = lightmaps;
    }
}
