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

using std::cerr;
using std::endl;
bool globaldebugflag = false;

int globalTriangleId = 0;

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

class Triangle
{
  public:
      double         X[3];
      double         Y[3];
      unsigned char color[3];

  // would some methods for transforming the triangle in place be helpful?
};

class Screen
{
  public:
      unsigned char   *buffer;
      int width, height;

  // would some methods for accessing and setting pixels be helpful?
};

std::vector<Triangle>
GetTriangles(void)
{
    vtkPolyDataReader *rdr = vtkPolyDataReader::New();
    rdr->SetFileName("proj1c_geometry.vtk");
    cerr << "Reading" << endl;
    rdr->Update();
    if (rdr->GetOutput()->GetNumberOfCells() == 0)
    {
        cerr << "Unable to open file!!" << endl;
        exit(EXIT_FAILURE);
    }
    vtkPolyData *pd = rdr->GetOutput();
    int numTris = pd->GetNumberOfCells();
    vtkPoints *pts = pd->GetPoints();
    vtkCellArray *cells = pd->GetPolys();
    vtkFloatArray *colors = (vtkFloatArray *) pd->GetPointData()->GetArray("color_nodal");
    float *color_ptr = colors->GetPointer(0);
    std::vector<Triangle> tris(numTris);
    vtkIdType npts;
    vtkIdType *ptIds;
    int idx;
    for (idx = 0, cells->InitTraversal() ; cells->GetNextCell(npts, ptIds) ; idx++)
    {
        if(idx==1000000)
            cerr<<";)"<<endl;
        if (npts != 3)
        {
            cerr << "Non-triangles!! ???" << endl;
            exit(EXIT_FAILURE);
        }
        tris[idx].X[0] = pts->GetPoint(ptIds[0])[0];
        tris[idx].X[1] = pts->GetPoint(ptIds[1])[0];
        tris[idx].X[2] = pts->GetPoint(ptIds[2])[0];
        tris[idx].Y[0] = pts->GetPoint(ptIds[0])[1];
        tris[idx].Y[1] = pts->GetPoint(ptIds[1])[1];
        tris[idx].Y[2] = pts->GetPoint(ptIds[2])[1];
        tris[idx].color[0] = (unsigned char) color_ptr[4*ptIds[0]+0];
        tris[idx].color[1] = (unsigned char) color_ptr[4*ptIds[0]+1];
        tris[idx].color[2] = (unsigned char) color_ptr[4*ptIds[0]+2];
    }
    cerr << "Done reading" << endl;

    return tris;
}

///////////
 // 0top 1left 2right
void drawLine(double coord1, double coord2, int lineY, Triangle *t, Screen *s)
{
    //fprintf(stderr,"---------draw line function-------\n");
    if(coord2 > s->width)
        coord2 = s->width-1;
    if(coord1<0)
        coord1=0;
    if(globaldebugflag){
        fprintf(stderr,"coord1 before ceiling %f, after ceiling %f\n",coord1,ceil441(coord1));
        fprintf(stderr,"coord2 before flooring %f, after flooring %f\n",coord2,floor441(coord2));
    }
    coord1=ceil441(coord1);
    coord2=floor441(coord2);
    double start = coord1;
    double end = coord2;
    for(double i = start; i<=end; i++)
    {
        if(globaldebugflag)
        {
            fprintf(stderr,"printing at x%f y%d\n",coord1,lineY);
        }
                if (i == 1776 && lineY == 1343)
                {
                    //cerr << "HRC: Getting color from " << globalTriangleId << endl;
                }
        int j=(i+lineY*s->width)*3;
        s->buffer[j] = t->color[0];
        s->buffer[j+1] = t->color[1];
        s->buffer[j+2] = t->color[2];
    }
}

void FlatBottomRasterizer(Triangle *t,Screen *s)
{
       /* cerr << "TRIANGLE = " << t->X[0] << ", " << t->Y[0]
                       <<", " << t->X[1] << ", " << t->Y[1]
                       <<", " << t->X[2] << ", " << t->Y[2]  << endl;*/
    double slope1 = (t->X[0] - t->X[1]) / (t->Y[0] - t->Y[1]);
    double slope2 = (t->X[0] - t->X[2]) / (t->Y[0] - t->Y[2]);
    if(globaldebugflag)
        fprintf(stderr,"-------rasterize bottom----\n");
    double yOne = ceil441(t->Y[1]);
    double yTwo = floor441(t->Y[0]);
    if(globaldebugflag){
        //fprintf(stderr,"y[1] before ceiling %f, after ceiling %f\n",t->Y[1],yCeil);
        //fprintf(stderr,"y[0] before flooring %f, after flooring %f\n",t->Y[0],yFloor);
        //fprintf(stderr,"coord1 %f coord2 %f\n",coord1,coord2);
    }
    if(floor441(t->Y[1])==floor441(t->Y[0]))
    {
        //fprintf(stderr,"invert fix!\n");
        //yOne=floor441(t->Y[1]);
        //yTwo=ceil441(t->Y[0]);
    }

    //first theory

    //double b = t->Y[0] - slope2*t->X[0];
    //double coord1 = (t->Y[1]-b)/slope1;

    //b = t->Y[0] - slope1*t->X[0];
    //double coord2 = (t->Y[1]-b)/slope2;

    //second theory
    double coord1 = t->X[1] + slope1*(yOne - t->Y[1]);
    double coord2 = t->X[2] + slope2*(yOne - t->Y[1]);


    //if(coord1>coord2)
        //fprintf(stderr,"left end %f right end %f\n",coord1,coord2);
        //fprintf(stderr, "SCANLINES ARE %d to %d\n", (int)yOne, (int)yTwo);
    for (int scanlineY = yOne; scanlineY <= yTwo; scanlineY++)
    {
        if(scanlineY>=1344)
            return;
//fprintf(stderr, "C1 = %f, C2 = %f\n", coord1, coord2);
        drawLine(coord1,coord2,scanlineY,t,s);
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
        if(scanlineY>=1344)
            return;
        drawLine(coord1,coord2,scanlineY,t,s);
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

    //move colors over haha
    int i;
    for(i=0;i<3;i++)
    {
        top->color[i] = t->color[i];
        bottom->color[i] = t->color[i];
    }

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
            top->X[1] = t->X[tempLeft];
            top->Y[1] = t->Y[tempLeft];
            top->X[2] = t->X[tempRight];
            top->Y[2] = t->Y[tempRight];
            bottom->X[0]=-10;
        }
        else
        {
            bottom->X[0] = t->X[tempDown];
            bottom->Y[0] = t->Y[tempDown];
            bottom->X[1] = t->X[tempLeft];
            bottom->Y[1] = t->Y[tempLeft];
            bottom->X[2] = t->X[tempRight];
            bottom->Y[2] = t->Y[tempRight];
            top->X[0]=-10;
        }
        return 1;
    }

    /////////////////////////////////////////////


    /*if(t->X[tempUp]==t->X[tempDown])
    {
        top->X[0] = t->X[tempUp];
        top->Y[0] = t->Y[tempUp];

        bottom->X[0] = t->X[tempDown];
        bottom->Y[0] = t->Y[tempDown];

        int ltemp = 1;
        int rtemp = 2;
        if(t->X[tempUp]<t->X[tempMid])
        {
            ltemp = 2;
            rtemp = 1;
        }
        top->X[rtemp] = t->X[tempMid];
        top->Y[rtemp] = t->Y[tempMid];

        bottom->X[rtemp] = t->X[tempMid];
        bottom->Y[rtemp] = t->Y[tempMid];

        ///////
        top->X[ltemp] = t->X[tempUp];
        top->Y[ltemp] = t->Y[tempMid];

        bottom->X[ltemp] = t->X[tempUp];
        bottom->Y[ltemp] = t->Y[tempMid];

        return 0;
    }*/

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

    bottom->X[0] = t->X[tempDown];
    bottom->Y[0] = t->Y[tempDown];


    ///////////////////////////////////////////
    //if left is the mid point
    slope = (t->Y[tempUp] - t->Y[tempDown]) / (t->X[tempUp] - t->X[tempDown]);
    double b = t->Y[tempUp] - slope*t->X[tempUp];
    leftOver = (t->Y[tempMid]-b)/slope;



    ////////////////////////////////////////////
    if(t->X[tempMid]<leftOver)
    {
        top->X[1] = t->X[tempMid];
        top->Y[1] = t->Y[tempMid];
        bottom->X[1] = t->X[tempMid];
        bottom->Y[1] = t->Y[tempMid];

        top->X[2] = leftOver;
        top->Y[2] = t->Y[tempMid];
        bottom->X[2] = leftOver;
        bottom->Y[2] = t->Y[tempMid];
    }
    else
    {
        top->X[2] = t->X[tempMid];
        top->Y[2] = t->Y[tempMid];
        bottom->X[2] = t->X[tempMid];
        bottom->Y[2] = t->Y[tempMid];

        top->X[1] = leftOver;
        top->Y[1] = t->Y[tempMid];
        bottom->X[1] = leftOver;
        bottom->Y[1] = t->Y[tempMid];

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

void ArbitraryRasterizer(Triangle *t,Screen *s)
{
    Triangle bottom;
    Triangle top;
    if(globaldebugflag){
        fprintf(stderr,"\noriginal triangle\n");
        printTriangle(t);
    }
    if(SplitTriangle(t, &bottom,&top)!=0){
        if(bottom.X[0]!=-10)
        {
            if(globaldebugflag)
            {
                fprintf(stderr,"bottom triangle\n");
                printTriangle(&bottom);
            }
            FlatTopRasterizer(&bottom,s);

        }
        if(top.X[0]!=-10)
        {
            FlatBottomRasterizer(&top,s);
            if(globaldebugflag){
                fprintf(stderr,"top triangle\n");
                printTriangle(&top);
            }
        }
    }
    else
    {
        cerr << "error!" <<endl;
    }
}

///////////

int main()
{
   vtkImageData *image = NewImage(1786, 1344);
   unsigned char *buffer =
     (unsigned char *) image->GetScalarPointer(0,0,0);
   int npixels = 1786*1344;
   for (int i = 0 ; i < npixels*3 ; i++)
       buffer[i] = 0;

   std::vector<Triangle> triangles = GetTriangles();
   
   Screen screen;
   screen.buffer = buffer;
   screen.width = 1786;
   screen.height = 1344;



   // YOUR CODE GOES HERE TO DEPOSIT TRIANGLES INTO PIXELS USING THE SCANLINE ALGORITHM

    for (std::vector<Triangle>::iterator it = triangles.begin() ; it != triangles.end(); ++it)
	{
	ArbitraryRasterizer(&(*it),&screen);
			globalTriangleId++;
	}
    /*for(int q = 2532853; q<=2532853;q++)
        ArbitraryRasterizer(&(triangles[q]),&screen);*/
   WriteImage(image, "allTriangles");
}
