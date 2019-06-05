Shader "Hidden/DecodeLightmap"
{
    Properties
    {
        _MainTex ("Texture", 2D) = "white" {}
    }
    SubShader
    {
        Tags { "RenderType"="Opaque" }
        LOD 100

        Pass
        {
            CGPROGRAM
            #pragma vertex vert
            #pragma fragment frag

            #include "UnityCG.cginc"

            struct appdata
            {
                float4 vertex : POSITION;
                float2 uv : TEXCOORD0;
            };

            struct v2f
            {
                float2 uv : TEXCOORD0;
                float4 vertex : SV_POSITION;
            };

            sampler2D _MainTex;
            float4 _MainTex_ST;

            v2f vert (appdata v)
            {
                v2f o;
                o.vertex = UnityObjectToClipPos(v.vertex);
                o.uv = TRANSFORM_TEX(v.uv, _MainTex);
                return o;
            }

			float4 frag (v2f i) : SV_Target
            {
				float4 ocol = tex2D(_MainTex, i.uv);
                //float3 col = (unity_Lightmap_HDR.x * pow(ocol.a, unity_Lightmap_HDR.y)) * ocol.rgb;
				float3 col = DecodeLightmap(ocol);
				col = pow(col, 2.2);
                return float4(col, 1);
            }
            ENDCG
        }
    }
}
