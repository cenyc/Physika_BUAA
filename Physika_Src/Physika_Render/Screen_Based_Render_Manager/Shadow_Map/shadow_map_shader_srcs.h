/*
* @file shader_map_shader_srcs.h
* @Brief GLSL shaders for shadow map
* @author Wei Chen
*
* This file is part of Physika, a versatile physics simulation library.
* Copyright (C) 2013- Physika Group.
*
* This Source Code Form is subject to the terms of the GNU General Public License v2.0.
* If a copy of the GPL was not distributed with this file, you can obtain one at:
* http://www.gnu.org/licenses/gpl-2.0.html
*
*/

#pragma once

namespace Physika {

#define STRINGIFY(A) #A

static const char * shadow_map_vertex_shader = "#version 330 core\n" STRINGIFY(

layout(location = 0) in vec3 vert_pos;

uniform mat4 light_proj_trans;
uniform mat4 light_view_trans;
uniform mat4 model_trans;

void main()
{
    gl_Position = light_proj_trans * light_view_trans * model_trans * vec4(vert_pos, 1.0);
}

);

static const char * shadow_map_frag_shader = "#version 330 core\n" STRINGIFY(

out vec4 frag_color;

void main()
{
    frag_color = vec4(0.0, 0.0, 0.0, 1.0);
}

);

//////////////////////////////////////////////////////////////////////////////////////////////////////////
static const char * shadow_map_render_vertex_shader = STRINGIFY(

void main()
{
    gl_Position = vec4(gl_Vertex.xyz, 1.0);
    gl_TexCoord[0] = gl_MultiTexCoord0;
}

);

static const char * shadow_map_render_frag_shader = STRINGIFY(

uniform sampler2D  shadow_tex;

void main()
{
    float depth_value = texture2D(shadow_tex, gl_TexCoord[0].xy).r;
    
    float mag_factor = 10.0f;
    gl_FragColor = vec4(mag_factor*(1.0f - depth_value), 0.0f, 0.0f, 1.0f);
}

);

}//end of namespace Physika