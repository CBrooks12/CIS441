#include <iostream>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>

using std::cerr;
using std::endl;

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
   std::vector<Triangle> rv(100);

   unsigned char colors[6][3] = { {255,128,0}, {255, 0, 127}, {0,204,204}, 
                                  {76,153,0}, {255, 204, 204}, {204, 204, 0}};
   for (int i = 0 ; i < 100 ; i++)
   {
       int idxI = i%10;
       int posI = idxI*100;
       int idxJ = i/10;
       int posJ = idxJ*100;
       int firstPt = (i%3);
       rv[i].X[firstPt] = posI;
       if (i == 50)
           rv[i].X[firstPt] = -10;
       rv[i].Y[firstPt] = posJ;
       rv[i].X[(firstPt+1)%3] = posI+99;
       rv[i].Y[(firstPt+1)%3] = posJ;
       rv[i].X[(firstPt+2)%3] = posI+i;
       rv[i].Y[(firstPt+2)%3] = posJ+10*(idxJ+1);
       if (i == 95)
          rv[i].Y[(firstPt+2)%3] = 1050;
       rv[i].color[0] = colors[i%6][0];
       rv[i].color[1] = colors[i%6][1];
       rv[i].color[2] = colors[i%6][2];
   }

   return rv;
}

///////////
 // 0top 1left 2right
void drawLine(double coord1, double coord2, int lineY, Triangle *t, Screen *s)
{
	if(coord2>1000)
		coord2=1000;
	if(coord1<0)
		coord1=0;
	for(int i = coord1; i<=coord2; i++)
	{
		int j=(i+lineY*1000)*3;
		s->buffer[j] = t->color[0];
		s->buffer[j+1] = t->color[1];
		s->buffer[j+2] = t->color[2];
	}
}


void swap(double * x, double * y)
{
	double temp = *x;
	*x = *y;
	*y = temp;
}

void FlatBottomRasterizer(Triangle *t,Screen *s)
{
	
	if(t->Y[0]==t->Y[1])
	{
		swap(&t->X[2],&t->X[0]);
		swap(&t->Y[2],&t->Y[0]);
		//cout<<"here";
	}
	else if(t->Y[0]==t->Y[2])
	{
		swap(&t->X[1],&t->X[0]);
		swap(&t->Y[1],&t->Y[0]);
		//cout<<"there";
	}
	else
	{
		//nothing?
	}
	if(t->X[2]<t->X[1])
	{
		swap(&t->X[1],&t->X[2]);
		swap(&t->Y[1],&t->Y[2]);
	}
	/*cout<<"y[0]"<<t->Y[0]<<endl;
	cout<<"y[1]"<<t->Y[1]<<endl;
	cout<<"y[2]"<<t->Y[2]<<endl;*/
	double slope1 = (t->X[0] - t->X[1]) / (t->Y[0] - t->Y[1]) ;
	double slope2 = (t->X[0] - t->X[2]) / (t->Y[0] - t->Y[2]) ;
	double coord1 = t->X[1];
	double coord2 = t->X[2];

	for (int scanlineY = t->Y[1]; scanlineY <= t->Y[0]; scanlineY++)
	{
		if(scanlineY>=1000)
			return;
		drawLine(ceil441(coord1),floor441(coord2),scanlineY,t,s);
		coord1 += slope1;
		coord2 += slope2;
	}
}

/*FlatTopRasterizer(Triangle *t)
{
}

ArbitraryRasterizer(Triangle *t)
{
    Triangle fb;
    Triangle ft;
    SplitTriangle(t, fb, ft);
    FlatBottomRasterizer(fb);
    FlatBottomRasterizer(ft);
}*/


///////////

int main()
{
   vtkImageData *image = NewImage(1000, 1000);
   unsigned char *buffer = 
     (unsigned char *) image->GetScalarPointer(0,0,0);
   int npixels = 1000*1000;
   for (int i = 0 ; i < npixels*3 ; i++)
       buffer[i] = 0;

   std::vector<Triangle> triangles = GetTriangles();
   
   Screen screen;
   screen.buffer = buffer;
   screen.width = 1000;
   screen.height = 1000;


	
   // YOUR CODE GOES HERE TO DEPOSIT TRIANGLES INTO PIXELS USING THE SCANLINE ALGORITHM
	for (std::vector<Triangle>::iterator it = triangles.begin() ; it != triangles.end(); ++it)
		FlatBottomRasterizer(&(*it),&screen);
		
   WriteImage(image, "allTriangles");
}
