#include <time.h>
#include <iostream>
#include <fstream>
#include <string>
#include  <direct.h>
#include <vector>
#include <fstream>
#include <assert.h>
#include <random>
#include "../Cloud_Dependency/trimesh2/include/TriMesh.h"
#include "../Cloud_Dependency/SDFGen/array3.h"
#include "../Cloud_Dependency/SDFGen/sdfGen.h"
#include "../Cloud_Dependency/SDFGen/vec.h"
#include <math.h>
using namespace std;

vector<float > readSdf(string path,int &nx,int &ny,int &nz, float &dx){
    ifstream infile;
    infile.open(path.c_str());
    assert(infile.is_open());
    string s;
    vector<float > r;

    float orx,ory,orz;
    float t;
    infile >> nx>>ny>>nz;
    infile >> orx>>ory>>orz;
//    assert(orx>=-1);
//    assert(ory >= -1);
//    assert(orz >= -1);
    infile >> dx;
    float distanceMin = 10,distanceMax=-10;
    for (int i=0;i < nx*ny*nz;i++){
        infile >> t;

        r.push_back(t);
        distanceMin = distanceMin<t?distanceMin:t;
        distanceMax = distanceMax>t?distanceMax:t;
    }
    for (float &d :r){
        if(d >0){
            d = 0;
        }
        else{
            d = d/distanceMin;
        }
    }
    cout <<"nx ny nz "<<nx<<" "<<ny<<" "<<nz<<endl;
    cout << r.size()<<" "<<distanceMax<<" "<<distanceMin<<endl;
    //大于0：外部 小于0：内部
    infile.close();
    return r;
}


class SamplePoint{
public:
    float x;
    float y;
    float z;
    float size=1;
    SamplePoint(float x, float y, float z){
        this->x = x;
        this->y = y;
        this->z = z;
    }

    friend std::ostream &operator<<(std::ostream& stream, const SamplePoint& p){
        stream <<p.x<<" "<<p.y<<" "<<p.z<<" "<<p.size<<endl;
        return stream;
    }
};


void index2Cord(int index,int nx,int ny,int nz, int &ix, int &iy, int &iz){

    ix = index/(ny*nz);
    iy = (index-(ix)*(ny*nz))/(nz);
    iz = index-(ix)*(ny*nz)-iy*(nz);


    assert(ix < nx);
    assert(iy<ny);
    assert(iz < nz);

}

void sdfGen(string objPath, float dx,int padding,int &nx,int &ny,int &nz,vector<float > &distanceField){

    sdfgen::Array3f r= sdfGen(objPath,dx,padding);
    nx = r.nk;
    ny = r.nj;
    nz = r.ni;
    distanceField.clear();

    double minDis = 10;
    for(unsigned int i = 0; i < r.a.size(); ++i) {
        distanceField.push_back(r.a[i]);
        minDis = minDis>r.a[i]?r.a[i]:minDis;
    }
    for (int i=0;i <distanceField.size();i++){
        distanceField[i] = distanceField[i]>0?0:distanceField[i]/minDis;
    }
    return ;
}


vector<SamplePoint> sample(vector<float> distanceField,float dx,int nx,int ny,int nz,int sampleCount){
    float distanceSum = 0;
    int gridCount=0;
    for (float &d : distanceField){
        assert(d <= 1+1e-3);
        if(d >=0){
            distanceSum += d;
            gridCount ++;
        } else{
            assert(false);
        }
    }

    int samplePerGrid = ceil((sampleCount/(float)gridCount)/((gridCount-distanceSum)/gridCount)) ;
    cout<<"sample per grid:"<<samplePerGrid<<endl;
    //samplePerGrid = 10;
    //sample
    default_random_engine generator;
    uniform_real_distribution<float > randomer(0,1);
    vector<SamplePoint> samplevec;
    for(int i=0;i < distanceField.size();i++){
        int ix,iy,iz;
        if(distanceField[i] <=0)continue;

        index2Cord(i,nx,ny,nz,ix,iy,iz);

        for(int j=0;j < samplePerGrid;j++){
            if (distanceField[i] >randomer(generator)){
                continue;
            }
            SamplePoint samplePoint((ix+randomer(generator))*dx,(iy+randomer(generator))*dx,(iz+randomer(generator))*dx);
            samplevec.push_back(samplePoint);
        }

    }
    return  samplevec;
}

void exportToPly(string plyfilePath,vector<SamplePoint> sampleVec){
    trimesh::TriMesh *mesh = new trimesh::TriMesh;
    for(SamplePoint &samplePoint : sampleVec){
        mesh->vertices.push_back(trimesh::point(samplePoint.x,samplePoint.y,samplePoint.z));
    }
    mesh->write(plyfilePath);
}

int main(int argc, char* argv[]){
    string sdfFilepath = "D:\\study\\lab\\Physika_BUAA\\Physika_Src\\Physika_Cloud\\Cloud_Render\\input\\14.sdf";
    //string exportFilePath = "D:\\study\\lab\\Physika_BUAA\\Physika_Src\\Physika_Cloud\\Cloud_Render\\input\\14.txt";
    string objFilePath = "D:\\study\\lab\\Physika_BUAA\\Physika_Src\\Physika_Cloud\\Cloud_Render\\input\\14.obj";
    int nx,ny,nz;
    float dx=0.1;
    //[-0.5-0.5]
    vector<float > distanceField;

    sdfGen(objFilePath,dx,0,nx,ny,nz,distanceField);
    int sampleCount = 1e5;
    vector<SamplePoint> sampleVec = sample(distanceField,dx,nx,ny,nz,sampleCount);

    exportToPly( "D:\\study\\lab\\Physika_BUAA\\Physika_Src\\Physika_Cloud\\Cloud_Render\\input\\14.ply",sampleVec);



    cout<<"finish"<<endl;
    return  0;
}
