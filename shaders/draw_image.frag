// Fragment shader
#version 150
#extension GL_ARB_explicit_attrib_location : require

uniform sampler2D u_texture;
uniform bool gamma_toggle;

in vec2 v_texcoord;
out vec4 frag_color;

void main()
{
    frag_color = texture(u_texture, v_texcoord);
    frag_color.rgb /= frag_color.a;
    if (gamma_toggle) {
        frag_color.rgb = pow(frag_color.rgb, vec3(1.0 / 2.2));
    }
}
