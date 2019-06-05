Shader "Hidden/Draw2UV"
{
	Properties
	{
        _LightmapST ("Lightmap ST", Vector) = (1,1,0,0)
	}
	SubShader
	{
		Tags { "RenderType"="Opaque" }
        ZTest Off
        Cull Off

		Pass
		{
			CGPROGRAM
			#pragma vertex vert
			#pragma fragment frag
			
			#include "UnityCG.cginc"

			struct appdata
			{
				float4 vertex : POSITION;
				float2 uv : TEXCOORD1;
			};

			struct v2f
			{
				float4 vertex : SV_POSITION;
			};
            float4 _LightmapST;
   
			v2f vert (appdata v)
			{
				v2f o;
				float2 uv = v.uv * _LightmapST.xy + _LightmapST.zw;
#if UNITY_UV_STARTS_AT_TOP
				uv.y = 1 - uv.y;
#endif
				o.vertex = float4(uv * 2 - 1, 0.5, 1);
				return o;
			}
			
			fixed4 frag (v2f i) : SV_Target
			{
                return float4(1,1,1,1);
			}
			ENDCG
		}
	}
}
