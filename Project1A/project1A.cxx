#include <iostream>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>

using std::cerr;
using std::endl;

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

void
makeImage(unsigned char *buffer)
{
	int i;
	int j;
	int k;
	int pixelCount = 0;
	for(i = 0; i<27;i++)
	{
		for(j = 0;j< 50;j++)
		{
			for(k = 0; k<1024; k++)
			{
				if(i%3==0)
				{
					buffer[pixelCount+2] = 0;
				}
				else if(i%3==1)
				{
					buffer[pixelCount+2] = 128;
				}
				else
				{
					buffer[pixelCount+2] = 255;
				}
				if((i/3)%3==0)
				{
					buffer[pixelCount+1] = 0;
				}
				else if((i/3)%3==1)
				{
					buffer[pixelCount+1] = 128;
				}
				else
				{
					buffer[pixelCount+1] = 255;
				}
				if(i/9==0)
				{
					buffer[pixelCount] = 0;
				}
				else if(i/9==1)
				{
					buffer[pixelCount] = 128;
				}
				else
				{
					buffer[pixelCount] = 255;
				}
				pixelCount+=3;
			}
		}
	} 
}


int main()
{
   std::cerr << "In main!" << endl;
   vtkImageData *image = NewImage(1024, 1350);
   unsigned char *buffer = 
     (unsigned char *) image->GetScalarPointer(0,0,0);
   makeImage(buffer);
   WriteImage(image, "out");
}
