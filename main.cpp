#include <iostream>



#include "./tracer/Pathtracer.h"
#include "Scene.h"
#include "./surface/Bezier.h"
#include "./kdtree/PhotonKdtree.h"
#include "./tracer/PPMtracer.h"
using namespace std;

void test_beizer(){
    const double bezier_div_x = 3;
    const double bezier_div_y = 2.5;
    double control_x[] = {20./bezier_div_x,27./bezier_div_x,30./bezier_div_x,30./bezier_div_x,30./bezier_div_x,25./bezier_div_x,20./bezier_div_x,15./bezier_div_x,30./bezier_div_x};
    double control_y[] = {0./bezier_div_y,0./bezier_div_y,10./bezier_div_y,20./bezier_div_y,30./bezier_div_y,40./bezier_div_y,60./bezier_div_y,70./bezier_div_y,80./bezier_div_y};
    BezierCurve curve(control_x, control_y, 9, 9, .365);
    BezierSurface surf = BezierSurface(Vec(0, 0, -5), curve);

////    cout<<curve.eval(curve,0.5).x<<" "<<curve.eval(curve,0.5).y<<endl;
//    cout<<curve.eval(curve,0.5).x<<" "<<curve.eval(curve,0.5).y<<endl;
////    cout<<curve.deri(curve,0.5).x<<" "<<curve.deri(curve,0.5).y<<endl;
//    cout<<curve.deri(curve,0.5).x<<" "<<curve.deri(curve,0.5).y<<endl;
    ofstream outobj("test.obj");
    for (double i = 0; i < 1; i+=0.01) {
        double xt = surf.curve.getpos(i).x;
        for (double j = 0; j <2*pi ; j+=0.01) {
            double x = xt*cos(j);
            double z  = xt*sin(j);
            double y = surf.curve.getpos(i).y;
            cout<<"v "<<x<<" "<<y<<" "<<z<<endl;
            outobj<<"v "<<x<<" "<<y<<" "<<z<<endl;
        }
    }
    outobj.close();

//    cout<<curve.getpos(0)<<endl;
//    Ray r = Ray(Vec(0,2,10),Vec(0,0.1,-1));
//    Hitpoint res = surf.intersect(r);
//    cout<<res.position<<endl;
//    cout<<res.norm<<endl;
//    cout<<res.distance<<endl;
//    cout<<curve.compute_width()<<endl;
//    cout<<curve.height<<endl;
//
//    ofstream outobj("test.obj");
//
//    for(double j = 0; j <1; j+=0.01)
//        outobj<<"v "<<curve.eval(curve,j)[0]<<" "<<curve.eval(curve,j)[1]<<" "<<0<<endl;
//
//    outobj.close();
}


//void test_photon_kd(){
//    double length = 1000;
//    vector<Photon> photonlist;
//    Vec flux = Vec(0,0,0);
//    Vec direction = Vec(0,0,0);
//    Vec center = Vec(500,500,500);
//    for (int i = 0; i < 10000; ++i) {
//        Vec pos = length*Vec(drand48(),drand48(),drand48());
//        photonlist.emplace_back(Photon(pos,flux,direction));
//    }
//    PhotonKdtree tree = PhotonKdtree(photonlist,30);
//    vector<Photon> near = tree.find_r(Photon(center,flux,direction),100);
//    for (int j = 0; j < near.size(); ++j) {
//        cout<<near.at(j)<<endl;
//    }
//};

void test_constant_medium(){
    Sphere* sphere = new Sphere(Vec(0,0,0),50);
    Constant_medium* medium =  new Constant_medium(sphere,0.1);
    medium->set_material(new Isotropic(new ConstantTexture(Vec(0,0,0))));
    for (int i = 0; i < 100000; ++i) {
        Ray r = Ray(Vec(0,0,-600),Vec(0,0,1));
        Hitpoint res =  medium->intersect(r);
        if(res.object!=NULL){
            cout<<res.distance<<endl;
            cout<<res.position<<endl;
        }

    }


}

int main(int argc, char**argv) {
    int l = -1, r = -1;
    SAMPLE = atoi(argv[1]);
    SPP_NUM = SAMPLE * 4;
    if (argc > 2) {
        l = atoi(argv[2]);
        r = atoi(argv[3]);
    }
//    test_constant_medium();
//    test_photon_kd();
//    test_beizer();


    Canvas *pt_canvas = new Canvas(h, w, 3);
    Scene scene = Scene();
    Scene light = Scene();
//    {
//        DofCamera camera = DofCamera();
//        Camera camera = Camera();
//        make_scene_for_dof(scene,light,camera);
//        PPMtracer *pm_render = new PPMtracer(&scene, &light, MAX_DEPTH);
//        pm_render->render(&camera, pt_canvas);
//    }
////
//
//
////    NaiveCamera camera = NaiveCamera();

    Camera camera = Camera();
////    make_bezier_scene(scene,light,camera);
////    make_big_scene(scene,light,camera);
    make_scene_for_mirror(scene,light,camera);
//         make_fog_scene(scene,light,camera);
    Pathtracer *pt_render = new Pathtracer(scene);
    pt_render->render(camera, pt_canvas, l, r);
    char s[20];
    sprintf(s, "result_%d_%d.bmp", l, r);
    pt_canvas->write(s);
    // pt_canvas->show("PathTracing");
//     PPMtracer *pm_render = new PPMtracer(&scene, &light, MAX_DEPTH);
//     pm_render->render(&camera, pt_canvas);
    return 0;
}