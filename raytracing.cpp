#include "raytracing.h"
#include "ray.h"
#include "hitable.h"
#include "sphere.h"
#include "triangle.h"
#include "box.h"

#include "utils2.h"  // Used for OBJ-mesh loading
#include <stdlib.h>  // Needed for drand48()

float frand(void) {
    double random;
    random = float(rand()) / float(RAND_MAX);

    while (random == 1) {
        random = float(rand()) / float(RAND_MAX);
    }
    return random;
}

float squared_length(glm::vec3 p) {
    return sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
}

glm::vec3 random_in_unit_sphere() {
    glm::vec3 p;
    do {
        p = glm::vec3(2.0*frand(), 2.0*frand(), 2.0*frand()) - glm::vec3(1, 1, 1);
    } while (squared_length(p) >= 1.0);
    return p;
}

bool refract(const glm::vec3& v, const glm::vec3& n, float ni_over_nt, glm::vec3& refracted) {
    glm::vec3 uv = glm::normalize(v);
    float dt = glm::dot(uv, n);
    float discriminant = 1.0 - ni_over_nt * ni_over_nt * (1 - dt * dt);
    if (discriminant > 0) {
        refracted = ni_over_nt * (uv - n * dt) - n * sqrt(discriminant);
        return true;
    }
    else
        return false;
}

float schlick(float cosine, float ref_idx) {
    float r0 = (1 - ref_idx) / (1 + ref_idx);
    r0 = r0 * r0;
    return r0 + (1 - r0) * pow((1 - cosine), 5);
}



namespace rt {

// Store scene (world) in a global variable for convenience
struct Scene {
    Sphere ground;
    std::vector<Sphere> spheres;
    std::vector<Box> boxes;
    std::vector<Triangle> mesh;
    Box mesh_bbox;
} g_scene;

bool hit_world(const Ray &r, float t_min, float t_max, HitRecord &rec)
{
    HitRecord temp_rec;
    bool hit_anything = false;
    float closest_so_far = t_max;

    if (g_scene.ground.hit(r, t_min, closest_so_far, temp_rec)) {
        hit_anything = true;
        closest_so_far = temp_rec.t;
        rec = temp_rec;
    }
    for (int i = 0; i < g_scene.spheres.size(); ++i) {
        if (g_scene.spheres[i].hit(r, t_min, closest_so_far, temp_rec)) {
            hit_anything = true;
            closest_so_far = temp_rec.t;
            rec = temp_rec;
        }
    }
    for (int i = 0; i < g_scene.boxes.size(); ++i) {
        if (g_scene.boxes[i].hit(r, t_min, closest_so_far, temp_rec)) {
            hit_anything = true;
            closest_so_far = temp_rec.t;
            rec = temp_rec;
        }
    }
    for (int i = 0; i < g_scene.mesh.size(); ++i) {
        if (g_scene.mesh[i].hit(r, t_min, closest_so_far, temp_rec)) {
            hit_anything = true;
            closest_so_far = temp_rec.t;
            rec = temp_rec;
        }
    }
    return hit_anything;
}

// This function should be called recursively (inside the function) for
// bouncing rays when you compute the lighting for materials, like this
//
// if (hit_world(...)) {
//     ...
//     return color(rtx, r_bounce, max_bounces - 1);
// }
//
//

// See Chapter 7 in the "Ray Tracing in a Weekend" book
glm::vec3 color(RTContext& rtx, const Ray& r, int max_bounces)
{
    if (max_bounces < 0) return glm::vec3(0.0f);

    HitRecord rec;
    if (hit_world(r, 0.001f, 9999.0f, rec)) {
        rec.normal = glm::normalize(rec.normal);
        if (rtx.show_normals) {
            return rec.normal * 0.5f + 0.5f;
        }

        //Diffuse
        if (rec.material == 0) {
            glm::vec3 target = rec.p + rec.normal + random_in_unit_sphere();
            return  rec.color * color(rtx, Ray(rec.p, target - rec.p), max_bounces - 1);// 0.5f * color(rtx, Ray(rec.p, target - rec.p), max_bounces - 1);
        }

        //Metal
        if (rec.material==1) {
            glm::vec3 reflected = glm::reflect(glm::normalize(r.direction()), rec.normal);
            Ray scattered = Ray(rec.p, reflected + rtx.fuzz* random_in_unit_sphere());
            glm::vec3 attenuation = rec.color;
            if (glm::dot(scattered.direction(), rec.normal) > 0) {
                return attenuation * color(rtx, scattered, max_bounces - 1);
            }
            else {
                return glm::vec3(0, 0, 0);
            }
        }

        //Di-electric 
        if (rec.material == 2) {
            glm::vec3 outward_normal;
            glm::vec3 reflected = glm::reflect(r.direction(), rec.normal);
            float ni_over_nt;
            glm::vec3 attenuation = glm::vec3(1.0, 1.0, 1.0);
            glm::vec3 refracted;
            float reflect_prob;
            float cosine;
            Ray scattered;


            if (glm::dot(r.direction(), rec.normal) > 0) {
                outward_normal = -rec.normal;
                ni_over_nt = rec.ref_index;
                cosine = rec.ref_index * glm::dot(r.direction(), rec.normal) / r.direction().length();
            }
            else {
                outward_normal = rec.normal;
                ni_over_nt = 1.0 / rec.ref_index;
                cosine = -glm::dot(r.direction(), rec.normal) / r.direction().length();
            }
            if (refract(r.direction(), outward_normal, ni_over_nt, refracted)) {
               scattered = Ray(rec.p, refracted);
            }
            else {
                scattered = Ray(rec.p, reflected);
                return glm::vec3(0, 0, 0);
            }
        
            return attenuation * color(rtx, scattered, max_bounces - 1);
       
        }
        return rec.normal * 0.5f + 0.5f;

    }
    // If no hit, return sky color
    else {
        glm::vec3 unit_direction = glm::normalize(r.direction());
        float t = 0.5f * (unit_direction.y + 1.0f);
        return (1.0f - t) * rtx.ground_color + t * rtx.sky_color;
    }
}


// MODIFY THIS FUNCTION!
void setupScene(RTContext &rtx, const char *filename)
{
    g_scene.ground = Sphere(glm::vec3(0.0f, -1000.5f, 0.0f), 1000.0f, 0, rtx.ground_color, 0.0);
    g_scene.spheres = {
        Sphere(glm::vec3(0.0f, 0.0f, 0.0f), 0.5f, 0, glm::vec3(0.1,0.9,0.9), 0.0),
        Sphere(glm::vec3(1.0f, 0.0f, 0.0f), 0.5f, 1, glm::vec3(0.9,0.1,0.1),0.0),
        Sphere(glm::vec3(-1.0f, 0.0f, 0.0f), 0.5f, 2, glm::vec3(0.4,0.3,0.5), 1.5),
        

       

    };

    ////Boxes
    //g_scene.boxes = {
    //    Box(glm::vec3(0.0f, -0.25f, 0.0f), glm::vec3(0.25f)),
    //    Box(glm::vec3(1.0f, -0.25f, 0.0f), glm::vec3(0.25f)),
    //    Box(glm::vec3(-1.0f, -0.25f, 0.0f), glm::vec3(0.25f)),
    //};

   // Bunny
    OBJMesh mesh;
    objMeshLoad(mesh, filename);
    g_scene.mesh.clear();
    for (int i = 0; i < mesh.indices.size(); i += 3) {
        int i0 = mesh.indices[i + 0];
        int i1 = mesh.indices[i + 1];
        int i2 = mesh.indices[i + 2];
        glm::vec3 v0 = mesh.vertices[i0] + glm::vec3(0.0f, 0.135f, 0.0f);
        glm::vec3 v1 = mesh.vertices[i1] + glm::vec3(0.0f, 0.135f, 0.0f);
        glm::vec3 v2 = mesh.vertices[i2] + glm::vec3(0.0f, 0.135f, 0.0f);
        g_scene.mesh.push_back(Triangle(v0, v1, v2));
    }
}

// MODIFY THIS FUNCTION!
void updateLine(RTContext &rtx, int y)
{
    int nx = rtx.width;
    int ny = rtx.height;
    float aspect = float(nx) / float(ny);
    glm::vec3 lower_left_corner(-1.0f * aspect, -1.0f, -1.0f);
    glm::vec3 horizontal(2.0f * aspect, 0.0f, 0.0f);
    glm::vec3 vertical(0.0f, 2.0f, 0.0f);
    glm::vec3 origin(0.0f, 0.0f, 0.0f);
    glm::mat4 world_from_view = glm::inverse(rtx.view);

    // You can try to parallelise this loop by uncommenting this line:
    #pragma omp parallel for schedule(dynamic)
    for (int x = 0; x < nx; ++x) {
        float u, v;
        if (rtx.anti_aliasing) {
             u = float(x + frand()) / float(nx);
             v = float(y + frand()) / float(ny);
        }
        else {
            u = (float(x) + 0.5f) / float(nx);
            v = (float(y) + 0.5f) / float(ny);
        }
        Ray r(origin, lower_left_corner + u * horizontal + v * vertical);
        r.A = glm::vec3(world_from_view * glm::vec4(r.A, 1.0f));
        r.B = glm::vec3(world_from_view * glm::vec4(r.B, 0.0f));

        if (rtx.current_frame <= 0) {
            // Here we make the first frame blend with the old image,
            // to smoothen the transition when resetting the accumulation
            glm::vec4 old = rtx.image[y * nx + x];
            rtx.image[y * nx + x] = glm::clamp(old / glm::max(1.0f, old.a), 0.0f, 1.0f);
        }
        glm::vec3 c = color(rtx, r, rtx.max_bounces);
        rtx.image[y * nx + x] += glm::vec4(c, 1.0f);
    }
}

void updateImage(RTContext &rtx)
{
    if (rtx.freeze) return;  // Skip update
    rtx.image.resize(rtx.width * rtx.height);  // Just in case...

    updateLine(rtx, rtx.current_line % rtx.height);

    if (rtx.current_frame < rtx.max_frames) {
        rtx.current_line += 1;
        if (rtx.current_line >= rtx.height) {
            rtx.current_frame += 1;
            rtx.current_line = rtx.current_line % rtx.height;
        }
    }
}

void resetImage(RTContext &rtx)
{
    rtx.image.clear();
    rtx.image.resize(rtx.width * rtx.height);
    rtx.current_frame = 0;
    rtx.current_line = 0;
    rtx.freeze = false;
}

void resetAccumulation(RTContext &rtx)
{
    rtx.current_frame = -1;
}

} // namespace rt
