// Vertex shader
#version 150
#extension GL_ARB_explicit_attrib_location : require

out vec2 v_texcoord;

uniform mat4 viewProjection;

void main()
{
    v_texcoord = vec2(gl_VertexID % 2, gl_VertexID / 2) * 2.0;
    gl_Position = vec4(v_texcoord * 2.0 - 1.0, 0.0, 1.0);
}
