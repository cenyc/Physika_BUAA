/*
 * @file shader_convention.h 
 * @Brief shader convention
 * @author: Wei Chen
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

/*
���ĵ�����������shader�й������ݽṹ�Ͳ��ֵ���Ӧ��Լ�����������ϸ����ش˹�Լ������ȷ��shader�������ݡ�

------------------------------------------------------------------------------------------------------
�������ԣ�

layout (location = 0) in vec3 vert_pos;          //����λ��
layout (location = 1) in vec3 vert_normal;       //���㷨��
layout (location = 2) in vec2 vert_tex_coord;    //������������
layout (location = 3) in vec3 vert_col;          //������ɫ 
layout (location = 4) in vec3 vert_vector;       //�����Զ�������

layout (location = 5) in float density;          //�����ܶ�
layout (location = 6) in int phase;              //����phase

------------------------------------------------------------------------------------------------------
�������

ͶӰ����
uniform mat4 proj_trans;

��ͼ����:
uniform mat4 view_trans;

ģ�;���:
uniform mat4 model_trans;

�۲�λ��:
uniform vec3 view_pos;

------------------------------------------------------------------------------------------------------

���պ���Ӱ��

struct DirectionalLight
{
    vec3 ambient;
    vec3 diffuse;
    vec3 specular;
    
    vec3 direction;
};

uniform int directional_light_num = 0;
uniform DirectionalLight directional_lights[5];

struct PointLight
{
    vec3 ambient;
    vec3 diffuse;
    vec3 specular;

    vec3 pos;
    float constant_atten;
    float linear_atten;
    float quadratic_atten;
};

uniform int point_light_num = 0;
uniform PointLight point_lights[5];

struct SpotLight
{
    vec3 ambient;
    vec3 diffuse;
    vec3 specular;

    vec3 pos;
    float constant_atten;
    float linear_atten;
    float quadratic_atten;

    vec3 spot_direction;
    float spot_exponent;
    float spot_cutoff;       //in radians

    bool use_spot_outer_cutoff;
    float spot_outer_cutoff; //in radians

    mat4 light_trans;
};

struct SpotLightShadowMap
{
    bool has_shadow_map;
    sampler2D shadow_map_tex;
};

uniform int spot_light_num = 0;
uniform SpotLight spot_lights[10];
uniform SpotLightShadowMap spot_light_shadow_maps[10];

struct FlexSpotLight
{
    vec3 ambient;
    vec3 diffuse;
    vec3 specular;

    vec3 pos;

    vec3 spot_direction;
    float spot_min;
    float spot_max;

    mat4 light_trans;
};

struct FlexSpotLightShadowMap
{
    bool has_shadow_map;
    sampler2D shadow_map_tex;
};

uniform int flex_spot_light_num = 0;
uniform FlexSpotLight flex_spot_lights[10];
uniform FlexSpotLightShadowMap flex_spot_light_shadow_maps[10];

uniform bool use_light;
uniform bool use_shadow_map; //ȫ�������Ƿ�ʹ��shadow_map

------------------------------------------------------------------------------------------------------

���ʣ�

struct Material
{
    vec3 Ka;
    vec3 Kd;
    vec3 Ks;
    float shininess;
    float alpha;
};

vec3 default_Ka = vec3(0.2);
vec3 default_Kd = vec3(0.8);
vec3 default_Ks = vec3(0.0);
float default_shininess = 0.0;

uniform bool use_material;
uniform Material material;

------------------------------------------------------------------------------------------------------

����

uniform bool use_tex;
uniform bool has_tex;
uniform sampler2D tex; //texture unit 0

------------------------------------------------------------------------------------------------------

��ɫ��

uniform bool use_custom_col;  //��Ӧ�������� layout (location = 3) in vec3 vert_col; 


*/