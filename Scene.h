//
// Created by YILI WANG on 6/24/19.
//

#ifndef TRACER_FINAL_SCENE_H
#define TRACER_FINAL_SCENE_H

#include "./base/Common.h"
#include "./base/Ray.h"
#include "./base/Camera.h"
#include "./object/Material.h"
#include "./object/Object.h"
#include "kdtree/Kdtree.h"
#include "./surface/Bezier.h"

//void make_scene_for_dof(Scene &scene, Scene &lt, DofCamera &camera){





void make_scene_for_mirror(Scene &scene, Scene &lt, Camera &camera){
    Plane *floor = new Plane(Vec(0, 0, 0), Vec(0, 1, 0));
    Plane *wall = new Plane(Vec(0, 0, -100), Vec(0, 0, 1));
    Plane *back = new Plane(Vec(0,0,150),Vec(0,0,-1));


    floor->set_material(new Diffuse(new ImageTexture("../data/wood7.jpg",Vec(-200, -50, 0), Vec(1, 0, 0), Vec(0, 0, -1), 20, 20, 2.2)));
    wall->set_material(new Diffuse(new ConstantTexture(Vec(235,239,212)/255.0)));
    scene.add(floor);
    scene.add(wall);
    {
        Vec l00(-20,200,-76);
        Vec l01(-20,0,-76);
        Vec l10(30,200,-76);
        Vec l11(30,0,-76);
        Object *quad = new Quad(l00, l10, l11, l01);
        quad->set_material(new Specular(new ConstantTexture(Vec(1,1,1))));
//        quad->set_light(new Light(Ve.1c(50, 50, 50)));

        scene.add(quad);
    }


    {
        Vec l00(118.0322, 90.6049, -52.5106);
        Vec l01(118.0322, 17.6197, -52.5106);
        Vec l10(118.0322, 90.6049, 7.8514);
        Vec l11(118.0322, 17.6197, 7.8514);
        Object *quad = new Quad(l00, l01, l10, l11);
        quad->set_light(new Light(Vec(50, 50, 50)));
//        quad->set_light(new Light(Ve.1c(50, 50, 50)));
        scene.add(quad);
    }

//    {
//        Vec l00(-10, 400.6049, 1000);
//        Vec l01(-10, 350.6197, 1000);
//        Vec l10(10.0322, 400.6049, 1000);
//        Vec l11(10, 350.6197, 1000);
//        Object *quad = new Quad(l00, l01, l10, l11);
//        quad->set_light(new Light(Vec(5, 5, 5)));
////        quad->set_light(new Light(Ve.1c(50, 50, 50)));
//        scene.add(quad);
//    }

//    {
//        Vec l00(-158.0322, 90.6049, -42.5106);
//        Vec l01(-158.0322, 17.6197, -42.5106);
//        Vec l10(-158.0322, 90.6049, 17.8514);
//        Vec l11(-158.0322, 17.6197, 17.8514);
//        Object *quad = new Quad(l00, l10, l11, l01);
//        quad->set_light(new Light(Vec(90, 90, 90)));
////        quad->set_light(new Light(Ve.1c(50, 50, 50)));
//        scene.add(quad);
//    }



    {
        Sphere *sphere = new Sphere(Vec(-5,80,40),5);
        sphere->set_material(new Refractive(new ConstantTexture(Vec(1,1,1))));
        scene.add(sphere);
    }

    {
        Sphere *sphere = new Sphere(Vec(-30,60,40),3);
        sphere->set_material(new Specular(new ConstantTexture(Vec(1,1,1))));
        scene.add(sphere);
    }
    {
        Sphere *sphere = new Sphere(Vec(0,30,40),1);
        sphere->set_material(new Specular(new ConstantTexture(Vec(1,1,1))));
        scene.add(sphere);
    }



//    {
//        Vec l00(-20, 300 ,-20);
//        Vec l01(-20, 300, 20);
//        Vec l10(20, 300, -20);
//        Vec l11(20, 300, 20);
//        Object *quad = new Quad(l00, l01, l10, l11);
//        quad->set_light(new Light(Vec(200, 200, 200)));
//        scene.add(quad);
//    }

//
    {
        const double bezier_div_x = 0.05;
        const double bezier_div_y = 0.05;
        double control_x[] = {0.0/bezier_div_x,0.25/bezier_div_x,0.49/bezier_div_x,0.727/bezier_div_x};
        double control_y[] = {0.0/bezier_div_y,0.53/bezier_div_y,0.275/bezier_div_y,1.01/bezier_div_y};

        BezierCurve curve = BezierCurve(control_x, control_y, 4,4,.365);
        BezierSurface* surf = new  BezierSurface(Vec(-20,20,40), curve);
        surf->set_material(new Refractive(new ConstantTexture(Vec(1,1,1)),2.41));
        scene.add(surf);
    }
    {
        const double bezier_div_x = 0.05;
        const double bezier_div_y = 0.05;
        double control_x[] = {0.87/bezier_div_x,0.866/bezier_div_x,0./bezier_div_x,0./bezier_div_x};
        double control_y[] = {0.03/bezier_div_y,0.74/bezier_div_y,0.69/bezier_div_y,1.176/bezier_div_y};

        BezierCurve curve = BezierCurve(control_x, control_y, 4,4,.365);
        BezierSurface* surf = new  BezierSurface(Vec(-20,0,40), curve);
        surf->set_material(new Refractive(new ConstantTexture(Vec(1,1,1)),2.41));
        scene.add(surf);
    }


//    TriMesh * mirror1_= new TriMesh("../data/mirror.obj",
//                                   new Diffuse(new ConstantTexture(Vec(36./255., 18./255., 8./255.))),Vec(20,20,20),Vec(0,50,20), true, -0.10,Vec::ZAxis);
//    TriMesh * mirror2_ = new TriMesh("../data/mirror.obj",
//                                    new Refractive(new ConstantTexture(Vec(255,255,255)/255.),2.417),Vec(10,10,10),Vec(0,80,20), true, 0.10,Vec::ZAxis);
//
//
//    Kdtree *mirror1 = new Kdtree(mirror1_);
//    Kdtree *mirror2 = new Kdtree(mirror2_);
//    scene.add(mirror1);
//    scene.add(mirror2);
//



//    camera = Camera(Vec(150, 100, 100), Vec(-0.4, -0.4, -1).normalize(), Vec(0, 1, 0),25,25);
    camera = Camera(Vec(150, 100, 450), Vec(-0.37, -0.1, -1).normalize(), Vec(0, 1, 0),25,25);
}



void make_scene_for_dof(Scene &scene, Scene &lt, Camera &camera){
    Plane *floor = new Plane(Vec(0, 0, 0), Vec(0, 1, 0));
    Plane *wall = new Plane(Vec(0, 0, -76), Vec(0, 0, 1));
    floor->set_material(new Diffuse(new ImageTexture("../data/wood7.jpg",Vec(-200, -50, 0), Vec(1, 0, 0), Vec(0, 0, -1), 20, 20, 2.2)));
    wall->set_material(new Diffuse(new ConstantTexture(Vec(235,239,212)/255.0)));
    {
        Vec l00(-20,200,-76);
        Vec l01(-20,0,-76);
        Vec l10(30,200,-76);
        Vec l11(30,0,-76);
        Object *quad = new Quad(l00, l10, l11, l01);
        quad->set_material(new Specular(new ConstantTexture(Vec(1,1,1))));
//        quad->set_light(new Light(Ve.1c(50, 50, 50)));

        scene.add(quad);
    }


//    {
//        Vec l00(118.0322, 90.6049, -52.5106);
//        Vec l01(118.0322, 17.6197, -52.5106);
//        Vec l10(118.0322, 90.6049, 7.8514);
//        Vec l11(118.0322, 17.6197, 7.8514);
//        Object *quad = new Quad(l00, l01, l10, l11);
//        quad->set_light(new Light(Vec(10, 10, 10)));
////        quad->set_light(new Light(Ve.1c(50, 50, 50)));
//        scene.add(quad);
//    }

//    {
//        Vec l00(-10, 400.6049, 1000);
//        Vec l01(-10, 350.6197, 1000);
//        Vec l10(10.0322, 400.6049, 1000);
//        Vec l11(10, 350.6197, 1000);
//        Object *quad = new Quad(l00, l01, l10, l11);
//        quad->set_light(new Light(Vec(5, 5, 5)));
////        quad->set_light(new Light(Ve.1c(50, 50, 50)));
//        scene.add(quad);
//    }

    {
        Vec l00(158.0322,90.6049, -22.5106);
        Vec l01(158.0322, 17.6197, -22.5106);
        Vec l10(158.0322, 90.6049, 37.8514);
        Vec l11(158.0322, 17.6197, 37.8514);
        l00 = l00.rotate(Vec::ZAxis,0.20);
        l01 = l01.rotate(Vec::ZAxis,0.20);
        l10 = l10.rotate(Vec::ZAxis,0.20);
        l11 = l11.rotate(Vec::ZAxis,0.20);

        Object *quad = new Quad(l00, l10, l11, l01);
        quad->set_light(new Light(Vec(50, 50, 50)));
//        quad->set_light(new Light(Ve.1c(50, 50, 50)));
        scene.add(quad);
    }


    {
        Sphere *sphere = new Sphere(Vec(-5,80,40),5);
        sphere->set_material(new Refractive(new ConstantTexture(Vec(1,1,1))));
        scene.add(sphere);
    }

    {
        Sphere *sphere = new Sphere(Vec(-30,60,40),3);
        sphere->set_material(new Specular(new ConstantTexture(Vec(1,1,1))));
        scene.add(sphere);
    }
    {
        Sphere *sphere = new Sphere(Vec(0,30,40),1);
        sphere->set_material(new Specular(new ConstantTexture(Vec(1,1,1))));
        scene.add(sphere);
    }



//    {
//        Vec l00(-20, 300 ,-20);
//        Vec l01(-20, 300, 20);
//        Vec l10(20, 300, -20);
//        Vec l11(20, 300, 20);
//        Object *quad = new Quad(l00, l01, l10, l11);
//        quad->set_light(new Light(Vec(200, 200, 200)));
//        scene.add(quad);
//    }


    {
        const double bezier_div_x = 0.05;
        const double bezier_div_y = 0.05;
        double control_x[] = {0.0/bezier_div_x,0.244/bezier_div_x,0.77/bezier_div_x,0.866/bezier_div_x};
        double control_y[] = {0.0/bezier_div_y,0.508/bezier_div_y,0.364/bezier_div_y,1.01/bezier_div_y};

        BezierCurve curve = BezierCurve(control_x, control_y, 4,4,.365);
        BezierSurface* surf = new  BezierSurface(Vec(-20,20,40), curve);
        surf->set_material(new Refractive(new ConstantTexture(Vec(1,1,1))));
        scene.add(surf);
    }
    {
        const double bezier_div_x = 0.05;
        const double bezier_div_y = 0.05;
        double control_x[] = {0.87/bezier_div_x,0.866/bezier_div_x,0./bezier_div_x,0./bezier_div_x};
        double control_y[] = {0.03/bezier_div_y,0.74/bezier_div_y,0.69/bezier_div_y,1.176/bezier_div_y};

        BezierCurve curve = BezierCurve(control_x, control_y, 4,4,.365);
        BezierSurface* surf = new  BezierSurface(Vec(-20,0,40), curve);
        surf->set_material(new Refractive(new ConstantTexture(Vec(1,1,1))));
        scene.add(surf);
    }


//


//    TriMesh *leave_ = new TriMesh("./data/leave.obj",
//                                  new Diffuse(new ConstantTexture(Vec(10/255.,205./255.,100./255.))),Vec(20,20,20),Vec(60,35,-20), true,0,Vec::ZAxis);
//
//    TriMesh *bark_ = new TriMesh("./data/bark.obj",
//                                 new Diffuse(new ConstantTexture(Vec(36./255., 18./255., 8./255.))),Vec(20,20,20),Vec(60,35,-20),true,0,Vec::ZAxis);
//
//    TriMesh *chair_ = new TriMesh("./data/chair_white.obj",
//                                  new Diffuse(new ConstantTexture(Vec(85./255.,85./255.,85./255.))),Vec(1,1,1),Vec(-20,40,0),true,-0.05,Vec::ZAxis);
//    TriMesh *chair_leg_ = new TriMesh("./data/chair_brown.obj",
//                                      new Diffuse(new ConstantTexture(Vec(38./255.,38./255.,38./255.))),Vec(1,1,1),Vec(-20,40,0),true,-0.05,Vec::ZAxis);
//
//    TriMesh *cat_ = new TriMesh("./data/kitten.obj",
//                                new Specular(new ConstantTexture(Vec(245./255.,245./255.,220./255.))),Vec(0.5,0.5,0.5),Vec(-100,0,20), true, 0.15,Vec::YAxis);

//    Kdtree *chair = new Kdtree(chair_);
//    Kdtree *chair_leg = new Kdtree(chair_leg_);
//
//    Kdtree *leave = new Kdtree(leave_);
//    Kdtree *bark = new Kdtree(bark_);
//    Kdtree *cat = new Kdtree(cat_);


//    scene.add(leave);
//    scene.add(bark);
    scene.add(floor);
    scene.add(wall);
//    scene.add(cat);
//    scene.add(chair);
//    scene.add(chair_leg);
//    scene.add(sphere1);

//    camera = DofCamera(Vec(150, 100, 450), Vec(-0.37, -0.1, -1).normalize(), Vec(0, 1, 0),25,0,80);
    camera = Camera(Vec(150, 100, 450), Vec(-0.37, -0.1, -1).normalize(), Vec(0, 1, 0),25,25);
//    camera = Camera(Vec(300, 100, 450), Vec(-0.37, -0.1, -1).normalize(), Vec(0, 1, 0),25,25);
}


void make_big_scene(Scene &scene, Scene &lt,Camera &camera){
    Plane *floor = new Plane(Vec(0, 0, 0), Vec(0, 1, 0));
    Plane *wall = new Plane(Vec(0, 0, -76), Vec(0, 0, 1));
    floor->set_material(new Diffuse(new ImageTexture("../data/wood.png",Vec(-200, -50, 0), Vec(1, 0, 0), Vec(0, 0, -1), 20, 20, 2.2)));
//    wall->set_material(new Diffuse(new ImageTexture("../data/fabric.jpg",Vec(0, 0, -76), Vec(1, 0, 0), Vec(0, 1, 0), 20, 20, 2.2)));
    wall->set_material(new Diffuse(new ConstantTexture(Vec(235,239,212)/255.0)));
//    wall->set_material(new Diffuse(new ImageTexture("../data/texture3.jpg",Vec(0,0,-76),Vec(1, 0, 0), Vec(0, -1, 0), 20, 20, 2.2)));

    {
        Vec l00(158.0322, 90.6049, -22.5106);
        Vec l01(158.0322, 17.6197, -22.5106);
        Vec l10(158.0322, 90.6049, 37.8514);
        Vec l11(158.0322, 17.6197, 37.8514);

        Object *quad = new Quad(l00, l01, l10, l11);
        quad->set_light(new Light(Vec(15, 15, 15)));
//        quad->set_light(new Light(Ve.1c(50, 50, 50)));
        scene.add(quad);
    }


    {
        Vec l00(394.3732, 244.5573, -25.3815);
        Vec l01(443.0328, 161.4010, -12.3249);
        Vec l10(373.0712, 248.3344, 78.0624);
        Vec l11(421.7308, 165.1782, 91.1191);
        Object *quad = new Quad(l00, l01, l10, l11);
        quad->set_light(new Light(Vec(40, 40, 40)));
        scene.add(quad);
    }

    {
        Vec l00(-118.0322, 90.6049, -52.5106);
        Vec l01(-118.0322, 17.6197, -52.5106);
        Vec l10(-118.0322, 90.6049, 7.8514);
        Vec l11(-118.0322, 17.6197, 7.8514);
        Object *quad = new Quad(l00, l01, l10, l11);
        quad->set_light(new Light(Vec(3*253/255, 3*234/255, 3*190/255)));
//        quad->set_light(new Light(Ve.1c(50, 50, 50)));
        scene.add(quad);
    }
//    const double bezier_div_x = 0.05;
//    const double bezier_div_y = 0.05;
//    double control_x[] = {0.003/bezier_div_x,2.145/bezier_div_x,0.325/bezier_div_x,0.032/bezier_div_x};
//    double control_y[] = {0.012/bezier_div_y,0.63/bezier_div_y,0.404/bezier_div_y,1.324/bezier_div_y};
//
//    BezierCurve curve = BezierCurve(control_x, control_y, 4,4,.365);
//    BezierSurface* surf = new  BezierSurface(Vec(60,0,20), curve);
//    surf->set_material(new Specular(new ConstantTexture(Vec(0.25,0.25,0.25))));
//    scene.add(surf);


    Sphere *sphere1 = new Sphere(Vec(40,14,30),14);
    sphere1->set_material(new Refractive (new ConstantTexture(Vec(1,1,1))));

    TriMesh *leave_ = new TriMesh("../data/leave.obj",
                                  new Diffuse(new ConstantTexture(Vec(10/255.,205./255.,100./255.))),Vec(20,20,20),Vec(60,35,-20), true,0,Vec::YAxis);

    TriMesh *bark_ = new TriMesh("../data/bark.obj",
                                 new Diffuse(new ConstantTexture(Vec(36./255., 18./255., 8./255.))),Vec(20,20,20),Vec(60,35,-20),true,0,Vec::YAxis);

    TriMesh *chair_ = new TriMesh("../data/chair_white.obj",
                                  new Diffuse(new ConstantTexture(Vec(243./255.,233./255.,209./255.))),Vec(0.7,0.7,0.7),Vec(-20,0,0),true,0.15,Vec::YAxis);
    TriMesh *chair_leg_ = new TriMesh("../data/chair_brown.obj",
                                      new Diffuse(new ConstantTexture(Vec(38./255.,38./255.,38./255.))),Vec(0.7,0.7,0.7),Vec(-20,0,0),true,0.15,Vec::YAxis);

    TriMesh *cat_ = new TriMesh("../data/kitten.obj",
                                new Diffuse(new ConstantTexture(Vec(245./255.,245./255.,220./255.))),Vec(0.3,0.3,0.3),Vec(-80,0,-10), true, 0.15,Vec::YAxis);

    Kdtree *chair = new Kdtree(chair_);
    Kdtree *chair_leg = new Kdtree(chair_leg_);

    Kdtree *leave = new Kdtree(leave_);
    Kdtree *bark = new Kdtree(bark_);
    Kdtree *cat = new Kdtree(cat_);


    scene.add(leave);
    scene.add(bark);
    scene.add(floor);
    scene.add(wall);
    scene.add(cat);
    scene.add(chair);
    scene.add(chair_leg);
    scene.add(sphere1);


    camera = Camera(Vec(0, 200, 500), Vec(0, -0.5, -1).normalize(), Vec(0, 1, 0), 20, 15);

}
//
//
void make_fog_scene(Scene &scene, Scene &lt,Camera &camera){


    Object *left = new Plane(Vec(100, 0, 0), Vec(-1, 0 ,0));
    Object *right = new Plane(Vec(-100, 0, 0), Vec(1, 0 ,0));
    Object *back = new Plane(Vec(0, 0, -100), Vec(0, 0, -1));
    Object *bottom = new Plane(Vec(0, 0, 0), Vec(0, 1, 0));
    Object *top = new Plane(Vec(0, 200, 0), Vec(0, -1, 0));
    right->set_material(new Diffuse(new ConstantTexture(Vec(.25, .25, .75))));
    left->set_material(new Diffuse(new ConstantTexture(Vec(.75, .25, .25))));
    back->set_material(new Diffuse(new ConstantTexture(Vec(1, 1, 1))));
    bottom->set_material(new Diffuse(new ConstantTexture(Vec(1, 1, 1))));
    top->set_material(new Diffuse(new ConstantTexture(Vec(.75, .75, .75))));

    camera = Camera(Vec(0, 100, 100), Vec(0, 0, 1).normalize(), Vec(0, 1, 0),25,25);

    Sphere *ball_light = new Sphere(Vec(25,10,-50),10);
    ball_light->set_light(new Light(Vec(255,242,114)));
    lt.add(ball_light);
    scene.add(left);
    scene.add(right);
    scene.add(back);
    scene.add(bottom);
    scene.add(top);

//    wall->set_material(new iffuse(new ImageTexture("../data/fabric.jpg",Vec(0, 0, -76), Vec(1, 0, 0), Vec(0, 1, 0), 20, 20, 2.2)));
//    wall->set_material(new Specular(new ConstantTexture(Vec(235,239,212)/255.0)));
//    wall->set_material(new Diffuse(new ImageTexture("../data/texture3.jpg",Vec(0,0,-76),Vec(1, 0, 0), Vec(0, -1, 0), 20, 20, 2.2)));
//
//    {
//        Vec l00(158.0322, 90.6049, -52.5106);
//        Vec l01(158.0322, 17.6197, -52.5106);
//        Vec l10(158.0322, 90.6049, 7.8514);
//        Vec l11(158.0322, 17.6197, 7.8514);
//        l00 = l00.rotate(Vec::YAxis,-0.15);
//        l01 = l01.rotate(Vec::YAxis,-0.15);
//        l10 = l10.rotate(Vec::YAxis,-0.15);
//        l11 = l11.rotate(Vec::YAxis,-0.15);
//        Object *quad = new Quad(l00, l01, l10, l11);
//        quad->set_light(new Light(Vec(15, 15, 15)));
////        quad->set_light(new Light(Ve.1c(50, 50, 50)));
//        scene.add(quad);
//    }

//    {
//        Vec l00(394.3732, 244.5573, -25.3815);
//        Vec l01(443.0328, 161.4010, -12.3249);
//        Vec l10(373.0712, 248.3344, 78.0624);
//        Vec l11(421.7308, 165.1782, 91.1191);
//        Object *quad = new Quad(l00, l01, l10, l11);
//        quad->set_light(new Light(Vec(40, 40, 40)));
//        scene.add(quad);
//    }
    const double bezier_div_x = 0.05;
    const double bezier_div_y = 0.05;
    double control_x[] = {0.003/bezier_div_x,2.145/bezier_div_x,0.325/bezier_div_x,0.032/bezier_div_x};
    double control_y[] = {0.012/bezier_div_y,0.63/bezier_div_y,0.404/bezier_div_y,1.324/bezier_div_y};

    BezierCurve curve = BezierCurve(control_x, control_y, 4,4,.365);
    BezierSurface* surf = new  BezierSurface(Vec(0,0,-20), curve);
    surf->set_material(new Specular(new ConstantTexture(Vec(0.25,0.25,0.25))));
    scene.add(surf);




    Sphere* sphere = new Sphere(Vec(0,100,0),200);
    Constant_medium* medium =  new Constant_medium(sphere,0.01);
    medium->set_material(new Isotropic(new ConstantTexture(Vec(1,1,1))));
    scene.add(medium);




}


//void make_bezier_scene(Scene &scene, Scene &lt, NaiveCamera &camera){
//
//    Vec l00(-0.5, 2.49, -0.5);
//    Vec l01(-0.5, 2.49,  0.5);
//    Vec l10( 0.5, 2.49, -0.5);
//    Vec l11( 0.5, 2.49,  0.5);
//    Object *light = new Quad(l00, l10, l11, l01);
////    Object *light = new Sphere(Vec(0,3,0),1);
//    light->set_light(new Light(Vec(16, 16, 16)));
//
//    Object *left = new Plane(Vec(2.5, 0, 0), Vec(-1, 0 ,0));
//    Object *right = new Plane(Vec(-2.5, 0, 0), Vec(1, 0 ,0));
//    Object *back = new Plane(Vec(0, 0, 1.2), Vec(0, 0, -1));
//    Object *front = new Plane(Vec(0, 0, -1.2), Vec(0, 0, 1));
//    Object *bottom = new Plane(Vec(0, -1, 0), Vec(0, 1, 0));
//    Object *top = new Plane(Vec(0, 2.5, 0), Vec(0, -1, 0));
//    right->set_material(new Diffuse(new ConstantTexture(Vec(.25, .25, .75))));
//    left->set_material(new Diffuse(new ConstantTexture(Vec(.75, .25, .25))));
//    back->set_material(new Diffuse(new ConstantTexture(Vec(1, 1, 1))));
//    bottom->set_material(new Diffuse(new ConstantTexture(Vec(1, 1, 1))));
//    top->set_material(new Diffuse(new ConstantTexture(Vec(.75, .75, .75))));
//
//
//    {
//        const double bezier_div_x = 1;
//        const double bezier_div_y = 1;
//        double control_x[] = {0/bezier_div_x,1.6/bezier_div_x,0/bezier_div_x,0/bezier_div_x};
//        double control_y[] = {0/bezier_div_y,1.3/bezier_div_y,0/bezier_div_y,1.3/bezier_div_y};
//
//        BezierCurve curve = BezierCurve(control_x, control_y, 4,4,.365);
//        Object* surf = new  BezierSurface(Vec(0, -1, 0), curve);
//        surf->set_material(new Specular(new ConstantTexture(Vec(0.25,0.25,0.25))));
//        scene.add(surf);
//    }
//
////    {
////        const double bezier_div_x = 5;
////        const double bezier_div_y = 5;
////        //
////        //灯柱
////        double control_x[] = {0/bezier_div_x,0/bezier_div_x,-8/bezier_div_x,0/bezier_div_x};
////        double control_y[] = {-15/bezier_div_y,-8/bezier_div_y,-8/bezier_div_y,20/bezier_div_y};
////
////        BezierCurve curve = BezierCurve(control_x, control_y, 4,4,.365);
////        Object* surf = new  BezierSurface(Vec(0, 0, 0), curve);
////        surf->set_material(new Refractive(new ConstantTexture(Vec(0,0,0))));
////        scene.add(surf);
////    }
//
//
//    scene.add(right);
//    scene.add(left);
//    scene.add(back);
//    scene.add(bottom);
//    scene.add(top);
//
//
//    scene.add(light);
//
////    camera = NaiveCamera(Vec(0, 0.3, -5.75), Vec(0, 0, 1).normalize(), Vec(0, 1, 0), 60);
//    camera = NaiveCamera(Vec(0, 0.34, -5.75), Vec(0, 0.2, 1).normalize(), Vec(0, 1, 0), 60);
//}



//void make_scene(Scene &scene, Scene &lt,Camera &camera) {
//    Vec l00(-0.75, 2.49, -0.75);
//    Vec l01(-0.75, 2.49,  0.75);
//    Vec l10( 0.75, 2.49, -0.75);
//    Vec l11( 0.75, 2.49,  0.75);
//    Object *light = new Quad(l00, l10, l11, l01); //down right -> down left -> up left -> up right
//    Object *left = new Plane(Vec(2.5, 0, 0), Vec(-1, 0 ,0));
//    Object *right = new Plane(Vec(-2.5, 0, 0), Vec(1, 0 ,0));
//    Object *back = new Plane(Vec(0, 0, 1.2), Vec(0, 0, -1));
//    Object *bottom = new Plane(Vec(0, -1, 0), Vec(0, 1, 0));
//
//    TriMesh *mesh = new TriMesh((char*)"/Users/yiliwang/Desktop/Course-computer graphics/Raytracer/model/simplify_bunny.obj", new Diffuse(new ConstantTexture(Vec(1,1,1))),Vec(1,10,10),Vec(0,-1.5,-2), true);
//    KdTree *bunny = new KdTree(mesh);
//
//    Object *top = new Plane(Vec(0, 2.5, 0), Vec(0, -1, 0));
//    Object *sphere1 = new Sphere(Vec(1, -0.5, -0.5),0.5);
//    Object *sphere2 = new Sphere(Vec(-1.3, 1, -2.5),0.7);
//    Object *sphere3 = new Sphere(Vec(0,0.3,-1.0),0.3);
//
//    light->set_light(new Light(Vec(16, 16, 16)));
//    right->set_material(new Diffuse(new ConstantTexture(Vec(.25, .25, .75))));
//    left->set_material(new Diffuse(new ConstantTexture(Vec(.75, .25, .25))));
//    back->set_material(new Diffuse(new ConstantTexture(Vec(.25, .75, .25))));
//    bottom->set_material(new Diffuse(new ImageTexture("/Users/yiliwang/Desktop/Course-computer graphics/Raytracer/model/back.ppm",Vec(2.5,0,0),Vec(1,0,0),Vec(0,0,1),20,20,2.2)));
//    top->set_material(new Diffuse(new ConstantTexture(Vec(.75, .75, .75))));
//    sphere1->set_material(new Refractive(new ConstantTexture(Vec(1,1,1))));
//    sphere2->set_material(new Refractive(new ConstantTexture(Vec(0.2,0.5,0.9))));
//    sphere3->set_material(new Diffuse(new ConstantTexture(Vec(0.3,0.9,0.5))));
////    sphere3->set_material(new Diffuse(new ConstantTexture(Vec(0.2,0.49,0.56))));
//
//    scene.add(right);
//    scene.add(left);
//    scene.add(back);
//    scene.add(bottom);
//    scene.add(top);
//    scene.add(sphere1);
//    scene.add(sphere2);
//    scene.add(sphere3);
//    scene.add(bunny);
//
//    scene.add(light);
//    lt.add(sphere2);
//    camera = Camera(Vec(0, 1, -13.75), Vec(0, -.1, 1).normalize(), Vec(0, 1, 0), 30, 0.1, 11.25);
//}



#endif //TRACER_FINAL_SCENE_H
