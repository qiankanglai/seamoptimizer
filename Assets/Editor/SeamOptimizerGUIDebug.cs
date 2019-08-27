using System.Collections;
using System.Collections.Generic;
using System.IO;
using UnityEngine;
using UnityEditor;

public partial class SeamOptimizerGUI
{
	[MenuItem("SO/Visualize Mesh 2UV")]
	static void Visualize2UV()
	{
		var lightmapMappings = GatherLightmaps();
		var mat = new Material(Shader.Find("Hidden/Draw2UV"));
		var mat2 = new Material(Shader.Find("Hidden/DecodeLightmap"));
		foreach (var lightmapMapping in lightmapMappings)
		{
			var lightmap = lightmapMapping.lightmap;
			var rt = new RenderTexture(lightmap.width, lightmap.height, 0, RenderTextureFormat.ARGBFloat, RenderTextureReadWrite.Default);

			Graphics.Blit(lightmap, rt, mat2);
			Graphics.SetRenderTarget(rt);
			GL.wireframe = true;
			mat.SetPass(0);
			foreach (var lightmapMesh in lightmapMapping.meshes)
			{
				mat.SetVector("_LightmapST", lightmapMesh.lightmapST);
				mat.SetPass(0);
				Graphics.DrawMeshNow(lightmapMesh.mesh, Vector3.zero, Quaternion.identity);
			}
			GL.wireframe = false;

			var tex = new Texture2D(lightmap.width, lightmap.height, TextureFormat.RGBAFloat, false);
			tex.ReadPixels(new Rect(0, 0, lightmap.width, lightmap.height), 0, 0);
			var filename = AssetDatabase.GetAssetPath(lightmap);
			var newfilename = Path.GetDirectoryName(filename) + "/" + Path.GetFileNameWithoutExtension(filename) + "_visualize_mesh2uv.exr";
			File.WriteAllBytes(newfilename, tex.EncodeToEXR(Texture2D.EXRFlags.CompressZIP));

			Graphics.SetRenderTarget(null);
			Object.DestroyImmediate(tex);
			Object.DestroyImmediate(rt);
		}
		Object.DestroyImmediate(mat);
		Object.DestroyImmediate(mat2);
		Debug.Log("Visualize Mesh 2UV Done");
		AssetDatabase.Refresh();
	}

	[MenuItem("SO/Visualize Seams")]
	static void VisualizeSeams()
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

			// visualize
			var debugPixels = new Color[pixels.Length];
#if false
			var seam = seams[4];
#else
			foreach (var seam in seams)
#endif
			{
				var stichingPoint = seam.stitchingPoints;
				var texels = seam.texels;
				for (int i = 0; i < stichingPoint.capacity; i++)
				{
					var point = stichingPoint.points[i];
					for (int j = 0; j < 2; j++)
					{
						var side = point.sides[j];
						for (int k = 0; k < 4; k++)
						{
							int idx = side.texels[k].y * lightmap.width + side.texels[k].x;
							//Debug.LogFormat("(" + side.texels[k].x + "," + side.texels[k].y + ")" + side.weights[k]);
							//Debug.Assert(side.weights[k] > 0);
							debugPixels[idx].r += side.weights[k];
						}
					}
				}
				// Kanglai: here texels.texels are just texels from stichingPoint.points.side
				for (int i = 0; i < texels.capacity; i++)
				{
					var texel = texels.texels[i];
					if (texel.x == -1)
						continue;
					int idx = texel.y * lightmap.width + texel.x;
					//Debug.LogFormat("(" + texel.x + "," + texel.y + ")");
					debugPixels[idx].g += 0.5f;
				}
			}

			var debug = new Texture2D(lightmap.width, lightmap.height, TextureFormat.RGBAFloat, false);
			for (int i = 0; i < debugPixels.Length; i++)
			{
				//if (debugPixels[i].r + debugPixels[i].g < 1e-6)
				//    debugPixels[i] = pixels[i];
				debugPixels[i].a = 1;
			}
			debug.SetPixels(debugPixels);
			debug.Apply();
			var filename = AssetDatabase.GetAssetPath(lightmap);
			var newfilename = Path.GetDirectoryName(filename) + "/" + Path.GetFileNameWithoutExtension(filename) + "_visualize_seams.exr";
			File.WriteAllBytes(newfilename, debug.EncodeToEXR(Texture2D.EXRFlags.CompressZIP));
			Object.DestroyImmediate(debug);
		}
		Debug.Log("Visualize Seams Done");
		AssetDatabase.Refresh();
	}
}
