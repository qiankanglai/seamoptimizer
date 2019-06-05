using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public static class Extensions
{
	public static Vector3 Max(Vector3 left, Vector3 right)
	{
        return new Vector3(
            Mathf.Max(left.x, right.x),
            Mathf.Max(left.y, right.y),
            Mathf.Max(left.z, right.z)
        );
	}

	public static Vector3 Min(Vector3 left, Vector3 right)
    {
        return new Vector3(
            Mathf.Min(left.x, right.x),
            Mathf.Min(left.y, right.y),
            Mathf.Min(left.z, right.z)
        );
	}

	public static Vector2 Multiply(Vector2 left, Vector2 right)
    {
        return new Vector2(
            left.x * right.x,
            left.y * right.y
        );
	}

	public static Vector3 Multiply(Vector3 left, Vector3 right)
    {
        return new Vector3(
            left.x * right.x,
            left.y * right.y,
            left.z * right.z
        );
	}
}
