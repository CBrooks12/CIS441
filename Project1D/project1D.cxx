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

double lerp(double start, double end, double t)
{
	//return (1-t)*start + t*end;
	return start + t*(end-start);
}

class Triangle
{
  public:
      double         X[3];
      double         Y[3];
      double		 Z[3];
      double  		 colors[3][3];

  // would some methods for transforming the triangle in place be helpful?
};

class Screen
{
  public:
      unsigned char  *buffer;
      int width, height;
      double *zbuffer;

  // would some methods for accessing and setting pixels be helpful?
};

std::vector<Triangle>
GetTriangles(void)
{
    vtkPolyDataReader *rdr = vtkPolyDataReader::New();
    rdr->SetFileName("proj1d_geometry.vtk");
    cerr << "Reading" << endl;
    rdr->Update();
    cerr << "Done reading" << endl;
    if (rdr->GetOutput()->GetNumberOfCells() == 0)
    {
        cerr << "Unable to open file!!" << endl;
        exit(EXIT_FAILURE);
    }
    vtkPolyData *pd = rdr->GetOutput();
    int numTris = pd->GetNumberOfCells();
    vtkPoints *pts = pd->GetPoints();
    vtkCellArray *cells = pd->GetPolys();
    vtkFloatArray *var = (vtkFloatArray *) pd->GetPointData()->GetArray("hardyglobal");
    float *color_ptr = var->GetPointer(0);
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
        tris[idx].X[0] = pts->GetPoint(ptIds[0])[0];
        tris[idx].X[1] = pts->GetPoint(ptIds[1])[0];
        tris[idx].X[2] = pts->GetPoint(ptIds[2])[0];
        tris[idx].Y[0] = pts->GetPoint(ptIds[0])[1];
        tris[idx].Y[1] = pts->GetPoint(ptIds[1])[1];
        tris[idx].Y[2] = pts->GetPoint(ptIds[2])[1];
        tris[idx].Z[0] = pts->GetPoint(ptIds[0])[2];
        tris[idx].Z[1] = pts->GetPoint(ptIds[1])[2];
        tris[idx].Z[2] = pts->GetPoint(ptIds[2])[2];
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

    for(double i = coord1; i<=coord2; i++)
    {
		if(i==s->width)
			return;
		int j=(i+lineY*s->width)*3;
		double pVal = (i-start)/(end-start);
		if(s->zbuffer[j] <= lerp(z1,z2,pVal))
		{
			s->zbuffer[j] = lerp(z1,z2,pVal);
			s->buffer[j] = ceil441(lerp(r1,r2,pVal)*255.0);
			s->buffer[j+1] = ceil441(lerp(g1,g2,pVal)*255.0);
			s->buffer[j+2] = ceil441(lerp(b1,b2,pVal)*255.0);
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
            memcpy(&top->colors[0],&t->colors[tempUp], sizeof(double)*3);
            
            top->X[1] = t->X[tempLeft];
            top->Y[1] = t->Y[tempLeft];
            top->Z[1] = t->Z[tempLeft];
            memcpy(&top->colors[1],&t->colors[tempLeft], sizeof(double)*3);
            top->X[2] = t->X[tempRight];
            top->Y[2] = t->Y[tempRight];
            top->Z[2] = t->Z[tempRight];
            memcpy(&top->colors[2],&t->colors[tempRight], sizeof(double)*3);
            
            
            bottom->X[0]=-10;
        }
        else
        {
            bottom->X[0] = t->X[tempDown];
            bottom->Y[0] = t->Y[tempDown];
            bottom->Z[0] = t->Z[tempDown];
            memcpy(&bottom->colors[0],&t->colors[tempDown], sizeof(double)*3);
            
            bottom->X[1] = t->X[tempLeft];
            bottom->Y[1] = t->Y[tempLeft];
            bottom->Z[1] = t->Z[tempLeft];
            memcpy(&bottom->colors[1],&t->colors[tempLeft], sizeof(double)*3);            
            bottom->X[2] = t->X[tempRight];
            bottom->Y[2] = t->Y[tempRight];
            bottom->Z[2] = t->Z[tempRight];
            memcpy(&bottom->colors[2],&t->colors[tempRight], sizeof(double)*3);      
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
    memcpy(&top->colors[0],&t->colors[tempUp], sizeof(double)*3);

    bottom->X[0] = t->X[tempDown];
    bottom->Y[0] = t->Y[tempDown];
    bottom->Z[0] = t->Z[tempDown];
	memcpy(&bottom->colors[0],&t->colors[tempDown], sizeof(double)*3);

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
        memcpy(&top->colors[1],&t->colors[tempMid], sizeof(double)*3);
        
        bottom->X[1] = t->X[tempMid];
        bottom->Y[1] = t->Y[tempMid];
        bottom->Z[1] = t->Z[tempMid];
        memcpy(&bottom->colors[1],&t->colors[tempMid], sizeof(double)*3);
        
        top->X[2] = leftOver;
        top->Y[2] = t->Y[tempMid];
        top->Z[2] = lerp(t->Z[tempDown],t->Z[tempUp],pVal);
		
        top->colors[2][0] = lerp(t->colors[tempDown][0],t->colors[tempUp][0],pVal);
        top->colors[2][1] = lerp(t->colors[tempDown][1],t->colors[tempUp][1],pVal);
		top->colors[2][2] = lerp(t->colors[tempDown][2],t->colors[tempUp][2],pVal);         
                
        bottom->X[2] = leftOver;
        bottom->Y[2] = t->Y[tempMid];
        bottom->Z[2] = lerp(t->Z[tempDown],t->Z[tempUp],pVal);

        bottom->colors[2][0] = lerp(t->colors[tempDown][0],t->colors[tempUp][0],pVal);
        bottom->colors[2][1] = lerp(t->colors[tempDown][1],t->colors[tempUp][1],pVal);
		bottom->colors[2][2] = lerp(t->colors[tempDown][2],t->colors[tempUp][2],pVal);
        
    }
    else
    {
        top->X[2] = t->X[tempMid];
        top->Y[2] = t->Y[tempMid];
        top->Z[2] = t->Z[tempMid];
        memcpy(&top->colors[2],&t->colors[tempMid], sizeof(double)*3);
        
        bottom->X[2] = t->X[tempMid];
        bottom->Y[2] = t->Y[tempMid];
        bottom->Z[2] = t->Z[tempMid];
        memcpy(&bottom->colors[2],&t->colors[tempMid], sizeof(double)*3);
        
        top->X[1] = leftOver;
        top->Y[1] = t->Y[tempMid];
        top->Z[1] = lerp(t->Z[tempDown],t->Z[tempUp],pVal);
        
        top->colors[1][0] = lerp(t->colors[tempDown][0],t->colors[tempUp][0],pVal);
        top->colors[1][1] = lerp(t->colors[tempDown][1],t->colors[tempUp][1],pVal);
		top->colors[1][2] = lerp(t->colors[tempDown][2],t->colors[tempUp][2],pVal);
        
        bottom->X[1] = leftOver;
        bottom->Y[1] = t->Y[tempMid];
        bottom->Z[1] = lerp(t->Z[tempDown],t->Z[tempUp],pVal);
        
        bottom->colors[1][0] = lerp(t->colors[tempDown][0],t->colors[tempUp][0],pVal);
        bottom->colors[1][1] = lerp(t->colors[tempDown][1],t->colors[tempUp][1],pVal);
		bottom->colors[1][2] = lerp(t->colors[tempDown][2],t->colors[tempUp][2],pVal);
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
        cerr << SplitTriangle(t, &bottom,&top) <<endl;
    }
}

///////////

int main()
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
   //cerr << "time for debugging" <<endl;
   std::vector<Triangle> triangles = GetTriangles();
   
   Screen screen;
   screen.buffer = buffer;
   screen.width = 1000;
   screen.height = 1000;
   screen.zbuffer = zbuffer;

	Triangle t;
	t.X[0] = 200;
	t.X[1] = 500;
	t.X[2] = 700;
	t.Y[0] = 600;
	t.Y[1] = 800;
	t.Y[2] = 200;
	t.Z[0] = 1;
	t.Z[1] = 1;
	t.Z[2] = 1;
	t.colors[0][0] = 1;
	t.colors[0][1] = 0;
	t.colors[0][2] = 0;
	t.colors[1][0] = 0;
	t.colors[1][1] = 0;
	t.colors[1][2] = 1;
	t.colors[2][0] = 0;
	t.colors[2][1] = 1;
	t.colors[2][2] = 0;
	//ArbitraryRasterizer(&t,&screen);
	/*t.X[0] = 200;
	t.X[1] = 500;
	t.X[2] = 700;
	t.Y[0] = 600;
	t.Y[1] = 800;
	t.Y[2] = 200;
	t.Z[0] = 1;
	t.Z[1] = 1;
	t.Z[2] = 1;
	t.colors[0][0] = 1;
	t.colors[0][1] = 0;
	t.colors[0][2] = 0;
	t.colors[1][0] = 0;
	t.colors[1][1] = 0;
	t.colors[1][2] = 1;
	t.colors[2][0] = 0;
	t.colors[2][1] = 1;
	t.colors[2][2] = 0;*/
	
   // YOUR CODE GOES HERE TO DEPOSIT TRIANGLES INTO PIXELS USING THE SCANLINE ALGORITHM
	//cerr << "time for debugging" <<endl;
    for (std::vector<Triangle>::iterator it = triangles.begin() ; it != triangles.end(); ++it)
	{
			ArbitraryRasterizer(&(*it),&screen);
			globalTriangleId++;
	}
    /*for(int q = 2532853; q<=2532853;q++)
        ArbitraryRasterizer(&(triangles[q]),&screen);*/
   //ArbitraryRasterizer(&t,&screen);
   WriteImage(image, "allTriangles");
}
