#include <iostream>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <stdlib.h>
#include <algorithm>
#include <string.h>
using std::cerr;
using std::endl;
bool globaldebugflag = false;

int globalTriangleId = 0;


////////////////////////////////

double calcNorm(double * normal)
{	
	double val = sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
	return val;
}
void normalize(double * normal)
{
	double theNorm = calcNorm(normal);
	normal[0] = normal[0]/theNorm;
	normal[1] = normal[1]/theNorm;
	normal[2] = normal[2]/theNorm;
}
double cotan(double x)
{
	return 1/tan(x);
	
}

double calcDot(double * normal1, double * normal2)
{
	double theDot = (normal1[0]*normal2[0] + normal1[1]*normal2[1] + normal1[2]*normal2[2]);
	return theDot;
}

void crossProduct(double *A, double *B, double *AB)
{
	AB[0] = A[1]*B[2]-A[2]*B[1];
	AB[1] = B[0]*A[2]-A[0]*B[2];
	AB[2] = A[0]*B[1]-A[1]*B[0];
}

class Matrix
{
  public:
    double          A[4][4];

    void            TransformPoint(const double *ptIn, double *ptOut);
    static Matrix   ComposeMatrices(const Matrix &, const Matrix &);
    void            Print(ostream &o);
};

void
Matrix::Print(ostream &o)
{
    for (int i = 0 ; i < 4 ; i++)
    {
        char str[256];
        sprintf(str, "(%.7f %.7f %.7f %.7f)\n", A[i][0], A[i][1], A[i][2], A[i][3]);
        o << str;
    }
}

Matrix
Matrix::ComposeMatrices(const Matrix &M1, const Matrix &M2)
{
    Matrix rv;
    for (int i = 0 ; i < 4 ; i++)
        for (int j = 0 ; j < 4 ; j++)
        {
            rv.A[i][j] = 0;
            for (int k = 0 ; k < 4 ; k++)
                rv.A[i][j] += M1.A[i][k]*M2.A[k][j];
        }

    return rv;
}

void
Matrix::TransformPoint(const double *ptIn, double *ptOut)
{
    ptOut[0] = ptIn[0]*A[0][0]
             + ptIn[1]*A[1][0]
             + ptIn[2]*A[2][0]
             + ptIn[3]*A[3][0];
    ptOut[1] = ptIn[0]*A[0][1]
             + ptIn[1]*A[1][1]
             + ptIn[2]*A[2][1]
             + ptIn[3]*A[3][1];
    ptOut[2] = ptIn[0]*A[0][2]
             + ptIn[1]*A[1][2]
             + ptIn[2]*A[2][2]
             + ptIn[3]*A[3][2];
    ptOut[3] = ptIn[0]*A[0][3]
             + ptIn[1]*A[1][3]
             + ptIn[2]*A[2][3]
             + ptIn[3]*A[3][3];
}

class Camera
{
  public:
    double          near, far;
    double          angle;
    double          position[3];
    double          focus[3];
    double          up[3];
    Matrix			ct;
	Matrix			vt;
	Matrix			dt;
	Matrix			tt;

    void ViewTransform(void) 
    {
      vt.A[0][0] = cotan(angle/2);
      vt.A[0][1] = 0;
      vt.A[0][2] = 0;
      vt.A[0][3] = 0;

      vt.A[1][0] = 0;
      vt.A[1][1] = cotan(angle/2);
      vt.A[1][2] = 0;
      vt.A[1][3] = 0;

      vt.A[2][0] = 0;
      vt.A[2][1] = 0;
      vt.A[2][2] = (far+near)/(far-near);
      vt.A[2][3] = -1;

      vt.A[3][0] = 0;
      vt.A[3][1] = 0;
      vt.A[3][2] = (2*far*near)/(far-near);
      vt.A[3][3] = 0;
    };


    void CameraTransform(double *v1, double *v2, double *v3, double *o) 
    {

      double t[3];
		for (int i = 0; i < 3; i++){
			ct.A[i][0] = v1[i];
			ct.A[i][1] = v2[i];
			ct.A[i][2] = v3[i];
			ct.A[i][3] = 0;
			t[i] = 0-o[i];
		}

      ct.A[3][0] = calcDot(v1,t);
      ct.A[3][1] = calcDot(v2,t);
      ct.A[3][2] = calcDot(v3,t);
      ct.A[3][3] = 1;

    };



    void DeviceTransform(double *A) 
    {

      double x = (1000*(A[0]+1))/2;
      double y = (1000*(A[1]+1))/2;
      double z = A[2];
      dt.A[0][0] = 1000/2;
      dt.A[0][1] = 0;
      dt.A[0][2] = 0;
      dt.A[0][3] = 0;

      dt.A[1][0] = 0;
      dt.A[1][1] = 1000/2;
      dt.A[1][2] = 0;
      dt.A[1][3] = 0;

      dt.A[2][0] = 0;
      dt.A[2][1] = 0;
      dt.A[2][2] = 1;
      dt.A[2][3] = 0;

      dt.A[3][0] = 1000/2;
      dt.A[3][1] = 1000/2;
      dt.A[3][2] = 0;
      dt.A[3][3] = 1;
	};
};

double SineParameterize(int curFrame, int nFrames, int ramp)
{
    int nNonRamp = nFrames-2*ramp;
    double height = 1./(nNonRamp + 4*ramp/M_PI);
    if (curFrame < ramp)
    {
        double factor = 2*height*ramp/M_PI;
        double eval = cos(M_PI/2*((double)curFrame)/ramp);
        return (1.-eval)*factor;
    }
    else if (curFrame > nFrames-ramp)
    {
        int amount_left = nFrames-curFrame;
        double factor = 2*height*ramp/M_PI;
        double eval =cos(M_PI/2*((double)amount_left/ramp));
        return 1. - (1-eval)*factor;
    }
    double amount_in_quad = ((double)curFrame-ramp);
    double quad_part = amount_in_quad*height;
    double curve_part = height*(2*ramp)/M_PI;
    return quad_part+curve_part;
}

Camera
GetCamera(int frame, int nframes)
{
    double t = SineParameterize(frame, nframes, nframes/10);
    Camera c;
    c.near = 5;
    c.far = 200;
    c.angle = M_PI/6;
    c.position[0] = 40*sin(2*M_PI*t);
    c.position[1] = 40*cos(2*M_PI*t);
    c.position[2] = 40;
    c.focus[0] = 0;
    c.focus[1] = 0;
    c.focus[2] = 0;
    c.up[0] = 0;
    c.up[1] = 1;
    c.up[2] = 0;
    
    double v1[3];
	double v2[3];
	double v3[3];
	
	for (int i = 0; i < 3; i++) 
		v3[i] = c.position[i] - c.focus[i];
	normalize(v3);
    crossProduct(c.up,v3,v1);
    normalize(v1);
   	crossProduct(v3,v1,v2);
	normalize(v2);
	normalize(v3);

    
    c.CameraTransform(v1,v2,v3,c.position);
    c.ViewTransform();
    c.DeviceTransform(c.position);
	c.tt = c.tt.ComposeMatrices(c.ct,c.vt);
	c.tt = c.tt.ComposeMatrices(c.tt,c.dt);
    
    return c;
}




///////////////////////////////

struct LightingParameters
{
    LightingParameters(void)
    {
         lightDir[0] = -0.6;
         lightDir[1] = 0;
         lightDir[2] = -0.8;
         Ka = 0.3;
         Kd = 0.7;
         Ks = 5.3;
         alpha = 7.5;
    };
  

    double lightDir[3];  // The direction of the light source
    double Ka;           // The coefficient for ambient lighting.
    double Kd;           // The coefficient for diffuse lighting.
    double Ks;           // The coefficient for specular lighting.
    double alpha;        // The exponent term for specular lighting.
};

double ceil441(double f)
{
    return ceil(f-0.00001);
}

double floor441(double f)
{
    return floor(f+0.00001);
}	


vtkImageData *
NewImage(int width, int height)
{
    vtkImageData *img = vtkImageData::New();
    img->SetDimensions(width, height, 1);
    img->AllocateScalars(VTK_UNSIGNED_CHAR, 3);
	
    return img;
}

void
WriteImage(vtkImageData *img, const char *filename)
{
   std::string full_filename = filename;
   full_filename += ".png";
   vtkPNGWriter *writer = vtkPNGWriter::New();
   writer->SetInputData(img);
   writer->SetFileName(full_filename.c_str());
   writer->Write();
   writer->Delete();
}

double lerp(double start, double end, double t)
{
	return start + t*(end-start);
}

class Triangle
{
  public:
      double         X[3];
      double         Y[3];
      double		 Z[3];
      double  		 colors[3][3];
	  double         normals[3][3];
	  double 		 shading[3];
  // would some methods for transforming the triangle in place be helpful?
};

class Screen
{
  public:
		unsigned char  *buffer;
		int width, height;
		double *zbuffer;
		Camera c;
	
  // would some methods for accessing and setting pixels be helpful?
};

//double CalcPhongShading(double * normal, Camera c)
double CalcPhongShading(Triangle * t, double*normal, Camera c, int i)
{
	LightingParameters lp;
	//cerr << "phong " << normal[0] << "," << normal[1] << "," <<normal[2]<<endl;
	double diffuse = fabs(calcDot(normal,lp.lightDir));
	double coef = 2.0*calcDot(lp.lightDir,normal);
	
	double V[3];
	/*V[0] = c.position[0] - c.focus[0];
	V[1] = c.position[1] - c.focus[1];
	V[2] = c.position[2] - c.focus[2];*/
	V[0] = t->X[i] - c.position[0]; // t.X[0] - c.position[0] , t.X[1] - c.position[0] , t.X[2] - c.position[0] 
	V[1] = t->Y[i] - c.position[1];
	V[2] = t->Z[i] - c.position[2];
	//cerr << t->X[i] <<" " << t->Y[i] << " " << t->Z[i] <<endl;
	
	normalize(V);
	
	double tempV[3];
	tempV[0] = coef*normal[0];
	tempV[1] = coef*normal[1];
	tempV[2] = coef*normal[2];
	
	double R[3];
	R[0] = tempV[0] - lp.lightDir[0];
	R[1] = tempV[1] - lp.lightDir[1];
	R[2] = tempV[2] - lp.lightDir[2];
	
	double specShade = fmax(0, pow(calcDot(R, V), lp.alpha));
	//return lp.Ka + lp.Kd*diffuse;
	return (lp.Ka + lp.Kd*diffuse + lp.Ks*specShade);
}
std::vector<Triangle>
GetTriangles(void)
{
    vtkPolyDataReader *rdr = vtkPolyDataReader::New();
    rdr->SetFileName("proj1e_geometry.vtk");
    cerr << "Reading" << endl;
    rdr->Update();
    cerr << "Done reading" << endl;
    if (rdr->GetOutput()->GetNumberOfCells() == 0)
    {
        cerr << "Unable to open file!!" << endl;
        exit(EXIT_FAILURE);
    }
    vtkPolyData *pd = rdr->GetOutput();
/*
vtkDataSetWriter *writer = vtkDataSetWriter::New();
writer->SetInput(pd);
writer->SetFileName("hrc.vtk");
writer->Write();
 */

    int numTris = pd->GetNumberOfCells();
    vtkPoints *pts = pd->GetPoints();
    vtkCellArray *cells = pd->GetPolys();
    vtkDoubleArray *var = (vtkDoubleArray *) pd->GetPointData()->GetArray("hardyglobal");
    double *color_ptr = var->GetPointer(0);
    //vtkFloatArray *var = (vtkFloatArray *) pd->GetPointData()->GetArray("hardyglobal");
    //float *color_ptr = var->GetPointer(0);
    vtkFloatArray *n = (vtkFloatArray *) pd->GetPointData()->GetNormals();
    float *normals = n->GetPointer(0);
    std::vector<Triangle> tris(numTris);
    vtkIdType npts;
    vtkIdType *ptIds;
    int idx;
    for (idx = 0, cells->InitTraversal() ; cells->GetNextCell(npts, ptIds) ; idx++)
    {
        if (npts != 3)
        {
            cerr << "Non-triangles!! ???" << endl;
            exit(EXIT_FAILURE);
        }
        double *pt = NULL;
        pt = pts->GetPoint(ptIds[0]);
        tris[idx].X[0] = pt[0];
        tris[idx].Y[0] = pt[1];
        tris[idx].Z[0] = pt[2];
        tris[idx].normals[0][0] = normals[3*ptIds[0]+0];
        tris[idx].normals[0][1] = normals[3*ptIds[0]+1];
        tris[idx].normals[0][2] = normals[3*ptIds[0]+2];
        pt = pts->GetPoint(ptIds[1]);
        tris[idx].X[1] = pt[0];
        tris[idx].Y[1] = pt[1];
        tris[idx].Z[1] = pt[2];
        tris[idx].normals[1][0] = normals[3*ptIds[1]+0];
        tris[idx].normals[1][1] = normals[3*ptIds[1]+1];
        tris[idx].normals[1][2] = normals[3*ptIds[1]+2];
        pt = pts->GetPoint(ptIds[2]);
        tris[idx].X[2] = pt[0];
        tris[idx].Y[2] = pt[1];
        tris[idx].Z[2] = pt[2];
        tris[idx].normals[2][0] = normals[3*ptIds[2]+0];
        tris[idx].normals[2][1] = normals[3*ptIds[2]+1];
        tris[idx].normals[2][2] = normals[3*ptIds[2]+2];

        // 1->2 interpolate between light blue, dark blue
        // 2->2.5 interpolate between dark blue, cyan
        // 2.5->3 interpolate between cyan, green
        // 3->3.5 interpolate between green, yellow
        // 3.5->4 interpolate between yellow, orange
        // 4->5 interpolate between orange, brick
        // 5->6 interpolate between brick, salmon
        double mins[7] = { 1, 2, 2.5, 3, 3.5, 4, 5 };
        double maxs[7] = { 2, 2.5, 3, 3.5, 4, 5, 6 };
        unsigned char RGB[8][3] = { { 71, 71, 219 },
                                    { 0, 0, 91 },
                                    { 0, 255, 255 },
                                    { 0, 128, 0 },
                                    { 255, 255, 0 },
                                    { 255, 96, 0 },
                                    { 107, 0, 0 },
                                    { 224, 76, 76 }
                                  };
        for (int j = 0 ; j < 3 ; j++)
        {
            float val = color_ptr[ptIds[j]];
            int r;
            for (r = 0 ; r < 7 ; r++)
            {
                if (mins[r] <= val && val < maxs[r])
                    break;
            }
            if (r == 7)
            {
                cerr << "Could not interpolate color for " << val << endl;
                exit(EXIT_FAILURE);
            }
            double proportion = (val-mins[r]) / (maxs[r]-mins[r]);
            tris[idx].colors[j][0] = (RGB[r][0]+proportion*(RGB[r+1][0]-RGB[r][0]))/255.0;
            tris[idx].colors[j][1] = (RGB[r][1]+proportion*(RGB[r+1][1]-RGB[r][1]))/255.0;
            tris[idx].colors[j][2] = (RGB[r][2]+proportion*(RGB[r+1][2]-RGB[r][2]))/255.0;
        }
    }

    return tris;
}

///////////
 // 0top 1left 2right
void drawLine(double coord1, double coord2, int lineY, Triangle *t, Screen *s, double ts)
{
    //fprintf(stderr,"---------draw line function-------\n");
    if(coord2 >= s->width)
        coord2= s->width;
    if(coord1<=0)
        coord1 =0;
    if(globaldebugflag){
        fprintf(stderr,"coord1 before ceiling %f, after ceiling %f\n",coord1,ceil441(coord1));
        fprintf(stderr,"coord2 before flooring %f, after flooring %f\n",coord2,floor441(coord2));
    }
    double start = coord1;
    double end = coord2;
    coord1=ceil441(coord1);
    coord2=floor441(coord2);
    double z1 = lerp(t->Z[1],t->Z[0],ts);
    double z2 = lerp(t->Z[2],t->Z[0],ts);

	double r1 = lerp(t->colors[1][0],t->colors[0][0],ts);
	double r2 = lerp(t->colors[2][0],t->colors[0][0],ts);
	
	double g1 = lerp(t->colors[1][1],t->colors[0][1],ts);
	double g2 = lerp(t->colors[2][1],t->colors[0][1],ts);
	
	double b1 = lerp(t->colors[1][2],t->colors[0][2],ts);
	double b2 = lerp(t->colors[2][2],t->colors[0][2],ts);
	
	double s1 = lerp(t->shading[1],t->shading[0],ts);
	double s2 = lerp(t->shading[2],t->shading[0],ts);
	
	//find x norms
	normalize(t->normals[0]);
	normalize(t->normals[1]);
	normalize(t->normals[2]);
	
	double normX1 = lerp(t->normals[1][0], t->normals[0][0],ts);
	double normX2 = lerp(t->normals[2][0], t->normals[0][0],ts);
	
	double normY1 = lerp(t->normals[1][1], t->normals[0][1],ts);
	double normY2 = lerp(t->normals[2][1], t->normals[0][1],ts);	
	
	double normZ1 = lerp(t->normals[1][2], t->normals[0][2],ts);
	double normZ2 = lerp(t->normals[2][2], t->normals[0][2],ts);
	
	double normal1[3];
	normal1[0] = normX1;
	normal1[1] = normY1;       
	normal1[2] = normZ1;
	
	double normal2[3];
	normal2[0] = normX2;
	normal2[1] = normY2;
	normal2[2] = normZ2;
	
	normalize(normal1);
	normalize(normal2);
	//cerr << "left: " << normal1[0] << "," << normal1[1] << "," << normal1[2] <<endl; 
    //cerr << "right: " << normal2[0] << "," << normal2[1] <<  "," << normal2[2] <<endl; 

    for(double i = coord1; i<=coord2; i++)
    {
		if(i==s->width)
			return;
		
		int j=(i+lineY*s->width)*3;
		double pVal = (i-start)/(end-start);
		if(s->zbuffer[j] <= lerp(z1,z2,pVal))
		{
			s->zbuffer[j] = lerp(z1,z2,pVal); 
			double aX = lerp(normal1[0],normal2[0],pVal);
			double aY = lerp(normal1[1],normal2[1],pVal);
			double aZ = lerp(normal1[2],normal2[2],pVal);
			double normal3[3];
			normal3[0] = aX;
			normal3[1] = aY;
			normal3[2] = aZ;
			normalize(normal3);
			double shading = lerp(s1,s2,pVal);			//~157k off
			//double shading = CalcPhongShading(normal3,s->c); // 197206 pixels off deprecated
			s->buffer[j] = std::min(255.0,ceil441(shading*lerp(r1,r2,pVal)*255.0));
			s->buffer[j+1] = std::min(255.0,ceil441(shading*lerp(g1,g2,pVal)*255.0));
			s->buffer[j+2] = std::min(255.0,ceil441(shading*lerp(b1,b2,pVal)*255.0));
		}
		else
		{
			//cerr << "went under" <<endl;
		}
       
    }
}

void FlatBottomRasterizer(Triangle *t,Screen *s)
{
    double slope1 = (t->X[0] - t->X[1]) / (t->Y[0] - t->Y[1]);
    double slope2 = (t->X[0] - t->X[2]) / (t->Y[0] - t->Y[2]);

    if(globaldebugflag)
        fprintf(stderr,"-------rasterize bottom----\n");
    double yOne = ceil441(t->Y[1]);
    double yTwo = floor441(t->Y[0]);
    if(globaldebugflag){

    }
    
    //second theory
    double coord1 = t->X[1] + slope1*(yOne - t->Y[1]);
    double coord2 = t->X[2] + slope2*(yOne - t->Y[1]);
    for (int scanlineY = yOne; scanlineY <= yTwo; scanlineY++)
    {
        if(scanlineY>=s->height)
            return;

        drawLine(coord1,coord2,scanlineY,t,s,(scanlineY-t->Y[1])/(t->Y[0]-t->Y[1]));
        coord1 += slope1;
        coord2 += slope2;
    }
}

void FlatTopRasterizer(Triangle *t,Screen *s)
{
    double slope1 = (t->X[1] - t->X[0]) / (t->Y[1] - t->Y[0]) ;
    double slope2 = (t->X[2] - t->X[0]) / (t->Y[2] - t->Y[0]) ;
    if(globaldebugflag)
        fprintf(stderr,"-------rasterize top-----\n");
    double yOne = ceil441(t->Y[0]);
    double yTwo = floor441(t->Y[1]);
    if(globaldebugflag){
        //fprintf(stderr,"y[0] before ceiling %f, after ceiling %f\n",t->Y[0],yCeil);
        //fprintf(stderr,"y[1] before flooring %f, after flooring %f\n",t->Y[1],yFloor);
        //fprintf(stderr,"coord1 %f coord2 %f\n",coord1,coord2);
    }
    double coord1 = t->X[0] + slope1*(yOne - t->Y[0]);
    double coord2 = t->X[0] + slope2*(yOne - t->Y[0]);

    for (int scanlineY = yOne; scanlineY <= yTwo; scanlineY++)
    //for (int scanlineY = ceil441(t->Y[1]); scanlineY > floor441(t->Y[0]); scanlineY--)
    {
        if(scanlineY>=s->height)
            return;
        drawLine(coord1,coord2,scanlineY,t,s,(t->Y[1]-scanlineY)/(t->Y[1]-t->Y[0]));
        coord1 += slope1;
        coord2 += slope2;
    }
}

int SplitTriangle(Triangle *t,Triangle *bottom, Triangle *top)
{
    double slope = 0;

    int tempMid = -1;
    int tempUp = -1;
    int tempDown = -1;
    int tempLeft = -1;
    int tempRight = -1;
    double leftOver = 0;


    if((t->Y[1] < t->Y[0] && t->Y[0] < t->Y[2]) || (t->Y[1] > t->Y[0] && t->Y[0] > t->Y[2])) //if 0 is middle y
    {
        tempMid = 0;
        if(t->Y[1]>t->Y[2])
        {
            tempUp = 1;
            tempDown = 2;
        }
        else
        {
            tempUp = 2;
            tempDown = 1;
        }
    }
    else if((t->Y[0] < t->Y[1] && t->Y[1] < t->Y[2]) || (t->Y[0] > t->Y[1] && t->Y[1] > t->Y[2])) // if 1 is middle y
    {
        tempMid = 1;
        if(t->Y[0]>t->Y[2])
        {
            tempUp = 0;
            tempDown = 2;
        }
        else
        {
            tempUp = 2;
            tempDown = 0;
        }
    }
    else if((t->Y[1] < t->Y[2] && t->Y[2] < t->Y[0]) || (t->Y[1] > t->Y[2] && t->Y[2] > t->Y[0])) // if 2 is middle y
    {
        tempMid = 2;
        if(t->Y[1]>t->Y[0])
        {
            tempUp = 1;
            tempDown = 0;
        }
        else
        {
            tempUp = 0;
            tempDown = 1;
        }
    }
    else //if there is a flat bottom or top
    {
        if(t->Y[0]==t->Y[1]) //2 is top/bot
        {
            if(t->Y[2]>t->Y[1])
                tempUp = 2;
            else
                tempDown = 2;
            if(t->X[0]<t->X[1]) //0 is left
            {
                tempLeft = 0;
                tempRight = 1;
            }
            else // 1 is left
            {
                tempLeft = 1;
                tempRight = 0;
            }
        }
        else if(t->Y[0]==t->Y[2]) //1 is top/bot
        {
            if(t->Y[1]>t->Y[2])
                tempUp = 1;
            else
                tempDown = 1;
            if(t->X[0]<t->X[2]) //0 is left
            {
                tempLeft = 0;
                tempRight = 2;
            }
            else
            {
                tempLeft = 2;
                tempRight = 0;
            }
        }
        else // 0 is top/bot
        {
            if(t->Y[0]>t->Y[2])
                tempUp = 0;
            else
                tempDown = 0;
            if(t->X[1]<t->X[2]) //1 is left
            {
                tempLeft = 1;
                tempRight = 2;
            }
            else
            {
                tempLeft = 2;
                tempRight = 1;
            }
        }
        if(tempUp!=-1)
        {
            top->X[0] = t->X[tempUp];
            top->Y[0] = t->Y[tempUp];
            top->Z[0] = t->Z[tempUp];
            
            top->shading[0] = t->shading[tempUp];
            
            memcpy(&top->colors[0],&t->colors[tempUp], sizeof(double)*3);
            memcpy(&top->normals[0],&t->normals[tempUp], sizeof(double)*3);
            
            top->X[1] = t->X[tempLeft];
            top->Y[1] = t->Y[tempLeft];
            top->Z[1] = t->Z[tempLeft];
            
            top->shading[1] = t->shading[tempLeft];
            
            memcpy(&top->colors[1],&t->colors[tempLeft], sizeof(double)*3);
            memcpy(&top->normals[1],&t->normals[tempLeft], sizeof(double)*3);
            top->X[2] = t->X[tempRight];
            top->Y[2] = t->Y[tempRight];
            top->Z[2] = t->Z[tempRight];
            
            top->shading[2] = t->shading[tempRight];
            
            memcpy(&top->colors[2],&t->colors[tempRight], sizeof(double)*3);
            memcpy(&top->normals[2],&t->normals[tempRight], sizeof(double)*3);
            
            
            bottom->X[0]=-10;
        }
        else
        {
            bottom->X[0] = t->X[tempDown];
            bottom->Y[0] = t->Y[tempDown];
            bottom->Z[0] = t->Z[tempDown];
            
            bottom->shading[0] = t->shading[tempDown];  
            
            memcpy(&bottom->colors[0],&t->colors[tempDown], sizeof(double)*3);
            memcpy(&bottom->normals[0],&t->normals[tempDown], sizeof(double)*3);
            
            bottom->X[1] = t->X[tempLeft];
            bottom->Y[1] = t->Y[tempLeft];
            bottom->Z[1] = t->Z[tempLeft];
            
            bottom->shading[1] = t->shading[tempLeft];           
            
            memcpy(&bottom->colors[1],&t->colors[tempLeft], sizeof(double)*3);       
            memcpy(&bottom->normals[1],&t->normals[tempLeft], sizeof(double)*3);        
                     
            bottom->X[2] = t->X[tempRight];
            bottom->Y[2] = t->Y[tempRight];
            bottom->Z[2] = t->Z[tempRight];
            
            bottom->shading[2] = t->shading[tempRight];
            
            memcpy(&bottom->colors[2],&t->colors[tempRight], sizeof(double)*3);      
            memcpy(&bottom->normals[2],&t->normals[tempRight], sizeof(double)*3);            
            top->X[0]=-10;
            		
        }
        return 1;
    }

    /////////////////////////////////////////



    if(t->X[tempUp]<t->X[tempDown])
    {
        if(t->X[tempUp]<t->X[tempMid])
            tempLeft = tempUp;
        else
            tempLeft = tempMid;
        if(t->X[tempMid] > t->X[tempDown])
            tempRight = tempMid;
        else
            tempRight = tempDown;
    }
    else
    {
        if(t->X[tempDown]<t->X[tempMid])
            tempLeft = tempDown;
        else
            tempLeft = tempMid;
        if(t->X[tempMid] > t->X[tempUp])
            tempRight = tempMid;
        else
            tempRight = tempUp;
    }

    ///////////////////////////////////////////


    //assign top and bottoms to first point
    top->X[0] = t->X[tempUp];
    top->Y[0] = t->Y[tempUp];
    top->Z[0] = t->Z[tempUp];
    
    top->shading[0] = t->shading[tempUp];
    
    memcpy(&top->colors[0],&t->colors[tempUp], sizeof(double)*3);
	memcpy(&top->normals[0],&t->normals[tempUp], sizeof(double)*3);

    bottom->X[0] = t->X[tempDown];
    bottom->Y[0] = t->Y[tempDown];
    bottom->Z[0] = t->Z[tempDown];
    
    bottom->shading[0] = t->shading[tempDown];
    
	memcpy(&bottom->colors[0],&t->colors[tempDown], sizeof(double)*3);
    memcpy(&bottom->normals[0],&t->normals[tempDown], sizeof(double)*3);	

    ///////////////////////////////////////////
    
    
    slope = (t->Y[tempUp] - t->Y[tempDown]) / (t->X[tempUp] - t->X[tempDown]);
    double b = t->Y[tempUp] - slope*t->X[tempUp];
    leftOver = (t->Y[tempMid]-b)/slope;
    if(isinf(slope)==1)
		leftOver = t->X[tempDown];
	if(isnan(slope)==1)
		leftOver = t->X[tempDown];
		
	double pVal = (t->Y[tempMid]-t->Y[tempDown])/(t->Y[tempUp]-t->Y[tempDown]);
    ////////////////////////////////////////////
    if(t->X[tempMid]<leftOver)
    {
        top->X[1] = t->X[tempMid];
        top->Y[1] = t->Y[tempMid];
        top->Z[1] = t->Z[tempMid];
        
        top->shading[1] = t->shading[tempMid];
        
        memcpy(&top->colors[1],&t->colors[tempMid], sizeof(double)*3);
        memcpy(&top->normals[1],&t->normals[tempMid], sizeof(double)*3);
        
        bottom->X[1] = t->X[tempMid];
        bottom->Y[1] = t->Y[tempMid];
        bottom->Z[1] = t->Z[tempMid];
        
        bottom->shading[1] = t->shading[tempMid];        
        
        memcpy(&bottom->colors[1],&t->colors[tempMid], sizeof(double)*3);
        memcpy(&bottom->normals[1],&t->normals[tempMid], sizeof(double)*3);
        
        top->X[2] = leftOver;
        top->Y[2] = t->Y[tempMid];
        top->Z[2] = lerp(t->Z[tempDown],t->Z[tempUp],pVal);
		
		top->shading[2] = lerp(t->shading[tempDown], t->shading[tempUp], pVal);
		
        top->colors[2][0] = lerp(t->colors[tempDown][0],t->colors[tempUp][0],pVal);
        top->colors[2][1] = lerp(t->colors[tempDown][1],t->colors[tempUp][1],pVal);
		top->colors[2][2] = lerp(t->colors[tempDown][2],t->colors[tempUp][2],pVal);   
		      
		top->normals[2][0] = lerp(t->normals[tempDown][0],t->normals[tempUp][0],pVal);
        top->normals[2][1] = lerp(t->normals[tempDown][1],t->normals[tempUp][1],pVal);
		top->normals[2][2] = lerp(t->normals[tempDown][2],t->normals[tempUp][2],pVal);         
		      
                
        bottom->X[2] = leftOver;
        bottom->Y[2] = t->Y[tempMid];
        bottom->Z[2] = lerp(t->Z[tempDown],t->Z[tempUp],pVal);

		bottom->shading[2] = lerp(t->shading[tempDown], t->shading[tempUp], pVal);

        bottom->colors[2][0] = lerp(t->colors[tempDown][0],t->colors[tempUp][0],pVal);
        bottom->colors[2][1] = lerp(t->colors[tempDown][1],t->colors[tempUp][1],pVal);
		bottom->colors[2][2] = lerp(t->colors[tempDown][2],t->colors[tempUp][2],pVal);
		
		bottom->normals[2][0] = lerp(t->normals[tempDown][0],t->normals[tempUp][0],pVal);
        bottom->normals[2][1] = lerp(t->normals[tempDown][1],t->normals[tempUp][1],pVal);
		bottom->normals[2][2] = lerp(t->normals[tempDown][2],t->normals[tempUp][2],pVal);   
        
    }
    else
    {
        top->X[2] = t->X[tempMid];
        top->Y[2] = t->Y[tempMid];
        top->Z[2] = t->Z[tempMid];
        
        top->shading[2] = t->shading[tempMid];
        
        memcpy(&top->colors[2],&t->colors[tempMid], sizeof(double)*3);
        memcpy(&top->normals[2],&t->normals[tempMid], sizeof(double)*3);
        
        bottom->X[2] = t->X[tempMid];
        bottom->Y[2] = t->Y[tempMid];
        bottom->Z[2] = t->Z[tempMid];
        
        bottom->shading[2] = t->shading[tempMid];
        
        memcpy(&bottom->colors[2],&t->colors[tempMid], sizeof(double)*3);
        memcpy(&bottom->normals[2],&t->normals[tempMid], sizeof(double)*3);
        
        top->X[1] = leftOver;
        top->Y[1] = t->Y[tempMid];
        top->Z[1] = lerp(t->Z[tempDown],t->Z[tempUp],pVal);
        
        top->shading[1] = lerp(t->shading[tempDown],t->shading[tempUp],pVal);
        
        top->colors[1][0] = lerp(t->colors[tempDown][0],t->colors[tempUp][0],pVal);
        top->colors[1][1] = lerp(t->colors[tempDown][1],t->colors[tempUp][1],pVal);
		top->colors[1][2] = lerp(t->colors[tempDown][2],t->colors[tempUp][2],pVal);
		top->normals[1][0] = lerp(t->normals[tempDown][0],t->normals[tempUp][0],pVal);
        top->normals[1][1] = lerp(t->normals[tempDown][1],t->normals[tempUp][1],pVal);
		top->normals[1][2] = lerp(t->normals[tempDown][2],t->normals[tempUp][2],pVal);
        
        bottom->X[1] = leftOver;
        bottom->Y[1] = t->Y[tempMid];
        bottom->Z[1] = lerp(t->Z[tempDown],t->Z[tempUp],pVal);
        
        bottom->shading[1] = lerp(t->shading[tempDown],t->shading[tempUp],pVal);
        
        bottom->colors[1][0] = lerp(t->colors[tempDown][0],t->colors[tempUp][0],pVal);
        bottom->colors[1][1] = lerp(t->colors[tempDown][1],t->colors[tempUp][1],pVal);
		bottom->colors[1][2] = lerp(t->colors[tempDown][2],t->colors[tempUp][2],pVal);
		bottom->normals[1][0] = lerp(t->normals[tempDown][0],t->normals[tempUp][0],pVal);
        bottom->normals[1][1] = lerp(t->normals[tempDown][1],t->normals[tempUp][1],pVal);
		bottom->normals[1][2] = lerp(t->normals[tempDown][2],t->normals[tempUp][2],pVal);
    }
    return 1;
}

void printTriangle(Triangle* t)
{
    if(globaldebugflag)
    {
        fprintf(stderr,"-----------printing triangle coords!-----------\n");
        fprintf(stderr,"X[0] %f Y[0] %f\n",t->X[0],t->Y[0]);
        fprintf(stderr,"X[1] %f Y[1] %f\n",t->X[1],t->Y[1]);
        fprintf(stderr,"X[2] %f Y[2] %f\n",t->X[2],t->Y[2]);
    }
}
void toDeviceSpace(Triangle* t, Camera c)
{
	double oldPoint1[4];
	oldPoint1[0] = t->X[0];
	oldPoint1[1] = t->Y[0];
	oldPoint1[2] = t->Z[0];
	oldPoint1[3] = 1;
	
	double oldPoint2[4];
	oldPoint2[0] = t->X[1];
	oldPoint2[1] = t->Y[1];
	oldPoint2[2] = t->Z[1];
	oldPoint2[3] = 1;
	
	double oldPoint3[4];
	oldPoint3[0] = t->X[2];
	oldPoint3[1] = t->Y[2];
	oldPoint3[2] = t->Z[2];
	oldPoint3[3] = 1;
	
	
	double point1[4];
	double point2[4];
	double point3[4];
	for(int i = 0; i < 4; i++)
	{
		point1[i] = 0;
		point2[i] = 0;
		point3[i] = 0;
	}
	
	c.tt.TransformPoint(oldPoint1,point1);
	c.tt.TransformPoint(oldPoint2,point2);
	c.tt.TransformPoint(oldPoint3,point3);
	
	t->X[0] = point1[0]/point1[3];
	t->Y[0] = point1[1]/point1[3];
	t->Z[0] = point1[2]/point1[3];
	
	t->X[1] = point2[0]/point2[3];
	t->Y[1] = point2[1]/point2[3];
	t->Z[1] = point2[2]/point2[3];
	
	t->X[2] = point3[0]/point3[3];
	t->Y[2] = point3[1]/point3[3];
	t->Z[2] = point3[2]/point3[3];
	
}
void ArbitraryRasterizer(Triangle *t,Screen *s,Camera c)
{
    Triangle bottom;
    Triangle top;
    if(globaldebugflag){
        fprintf(stderr,"\noriginal triangle\n");
        printTriangle(t);
    }
    if(globaldebugflag){
		cerr << t->normals[0][0] << "," << t->normals[0][1] << "," << t->normals[0][2] <<endl; 
		cerr << t->normals[1][0] << "," << t->normals[1][1] << "," << t->normals[1][2] <<endl; 
		cerr << t->normals[2][0] << "," << t->normals[2][1] << "," << t->normals[2][2] <<endl<<endl; 
	}
	normalize(t->normals[0]);
	normalize(t->normals[1]);
	normalize(t->normals[2]);
	t->shading[0] = CalcPhongShading(t,t->normals[0],s->c,0);
	t->shading[1] = CalcPhongShading(t,t->normals[1],s->c,1);
	t->shading[2] = CalcPhongShading(t,t->normals[2],s->c,2);
	toDeviceSpace(t,c);
    if(SplitTriangle(t, &bottom,&top)!=0){

        if(bottom.X[0]!=-10)
        {
            if(globaldebugflag)
            {
                fprintf(stderr,"bottom triangle\n");
                printTriangle(&bottom);
            }
            if(globaldebugflag){
				cerr << "flat top" << endl;
				cerr << bottom.normals[0][0] << "," << bottom.normals[0][1] << "," << bottom.normals[0][2] <<endl; 
				cerr << bottom.normals[1][0] << "," << bottom.normals[1][1] << "," << bottom.normals[1][2] <<endl; 
				cerr << bottom.normals[2][0] << "," << bottom.normals[2][1] << "," << bottom.normals[2][2] <<endl<<endl; 
			}
            FlatTopRasterizer(&bottom,s);

        }
        if(top.X[0]!=-10)
        {		
			if(globaldebugflag){
			cerr << "flat bottom" << endl;
			cerr << top.normals[0][0] << "," << top.normals[0][1] << "," << top.normals[0][2] <<endl; 
			cerr << top.normals[1][0] << "," << top.normals[1][1] << "," << top.normals[1][2] <<endl; 
			cerr << top.normals[2][0] << "," << top.normals[2][1] << "," << top.normals[2][2] <<endl<<endl; }
            FlatBottomRasterizer(&top,s);
            if(globaldebugflag){
                fprintf(stderr,"top triangle\n");
                
                printTriangle(&top);
            }
        }
    }
    else
    {
        cerr << SplitTriangle(t, &bottom,&top) <<endl;
    }
}


///////////	
void createImage(int x, int total, std::vector<Triangle> t,char const* output)
{
		
	vtkImageData *image = NewImage(1000, 1000);
	unsigned char *buffer =
	 (unsigned char *) image->GetScalarPointer(0,0,0);
	double *zbuffer = new double[1000*1000*3];
	int npixels = 1000*1000;
	for (int i = 0 ; i < npixels*3 ; i++)
	{
	   buffer[i] = 0;
	   zbuffer[i] = -1;
	}
   
	Screen screen;
	screen.buffer = buffer;
	screen.width = 1000;
	screen.height = 1000;
	screen.zbuffer = zbuffer;
	
	Camera c = GetCamera(x,total);
	
	screen.c = c;
	// YOUR CODE GOES HERE TO DEPOSIT TRIANGLES INTO PIXELS USING THE SCANLINE ALGORITHM
	//cerr << "time for debugging" <<endl;
    for (std::vector<Triangle>::iterator it = t.begin() ; it != t.end(); ++it)
	{
		ArbitraryRasterizer(&(*it),&screen,c);
		globalTriangleId++;
	}
    //for(int q = 0; q<=0;q++)
    //    ArbitraryRasterizer(&(triangles[q]),&screen);
   //ArbitraryRasterizer(&t[0],&screen,c);
   
   WriteImage(image, output);
}
int main()
{
	std::vector<Triangle> triangles = GetTriangles();
	createImage(0,1000,triangles, "allTriangles0");
	createImage(250,1000,triangles, "allTriangles250");
	createImage(500,1000,triangles, "allTriangles500");
	createImage(750,1000,triangles, "allTriangles750");
}
