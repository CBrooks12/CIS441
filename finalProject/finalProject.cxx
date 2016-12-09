#include "vtkSmartPointer.h"
#include "vtkSphereSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkCallbackCommand.h"
#include "vtkCommand.h"
#include "vtkActor.h"
#include "vtkInteractorStyle.h"
#include "vtkObjectFactory.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkProperty.h"
#include "vtkCamera.h"
#include "vtkLight.h"
#include "vtkOpenGLPolyDataMapper.h"
#include "vtkJPEGReader.h"
#include "vtkImageData.h"

#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkCellArray.h>

#include <vtkTextActor.h>
#include <vtkTextProperty.h>

#include <ctime>
#include <vector>
#include <string>
class Triangle
{
  public:
      double         X[3];
      double         Y[3];
      double         Z[3];
};

class vtk441Mapper;

class vtk441Mapper : public vtkOpenGLPolyDataMapper
{
  protected:
   GLuint displayList;
   bool   initialized;
   float  size;
   int Red;
   int Green;
   int Blue;
   int direction;
   std::vector<int> aListOfStars;
  vtkImageData *img;
  vtkJPEGReader *rdr;
  GLuint texture;
  int dims[3];

  public:
   vtk441Mapper()
   {
     initialized = false;
     size = 1;
     direction = 0;
     Red = 255;
     Green = 255;
     Blue = 255;
     int i = 0;
     for(i; i<100; i++)
		aListOfStars.push_back((rand() % 500)-250);
   }

   void   IncrementSize()
   {
       size += 0.01;
       if (size > 2.0)
           size = 1.0;
   }
   
	void ChangeColor()
	{
		if(direction == 0)
		{
			if(Red!=30)
				Red--;
			else if(Green!=30)
				Green--;
			else
				Blue--;
			if(Blue==30)
				direction = 1;
		}
		else //direction ==1
		{
			if(Red!=255)
				Red++;
			else if(Green!=255)
				Green++;
			else
				Blue++;
			if(Blue==255)
				direction = 0;
		}
		
	}

   void
   RemoveVTKOpenGLStateSideEffects()
   {
     float Info[4] = { 0, 0, 0, 1 };
     glLightModelfv(GL_LIGHT_MODEL_AMBIENT, Info);
     float ambient[4] = { 1,1, 1, 1.0 };
     glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, ambient);
     float diffuse[4] = { 1, 1, 1, 1.0 };
     glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diffuse);
     float specular[4] = { 1, 1, 1, 1.0 };
     glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular);
   }


   void SetupLight(void)
   {
       glEnable(GL_LIGHTING);
       glEnable(GL_LIGHT0);
       GLfloat diffuse0[4] = { 0.6, 0.6, 0.6, 1 };
       GLfloat ambient0[4] = { 0.4, 0.4, 0.4, 1 };
       GLfloat specular0[4] = { 0.2, 0.2, 0.2, 1 };
       GLfloat pos0[4] = { 300, 200, 0, 1 };
       glLightfv(GL_LIGHT0, GL_POSITION, pos0);
       glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse0);
       glLightfv(GL_LIGHT0, GL_AMBIENT, ambient0);
       glLightfv(GL_LIGHT0, GL_SPECULAR, specular0);
       glDisable(GL_LIGHT1);
       glDisable(GL_LIGHT2);
       glDisable(GL_LIGHT3);
       glDisable(GL_LIGHT5);
       glDisable(GL_LIGHT6);
       glDisable(GL_LIGHT7);
   }
};

class vtk441MapperPart1 : public vtk441Mapper
{
 public:
   static vtk441MapperPart1 *New();
    std::vector<Triangle> SplitTriangle(std::vector<Triangle> &list)
   {
       std::vector<Triangle> output(4*list.size());
       for (unsigned int i = 0 ; i < list.size() ; i++)
       {
           double mid1[3], mid2[3], mid3[3];
           mid1[0] = (list[i].X[0]+list[i].X[1])/2;
           mid1[1] = (list[i].Y[0]+list[i].Y[1])/2;
           mid1[2] = (list[i].Z[0]+list[i].Z[1])/2;
           mid2[0] = (list[i].X[1]+list[i].X[2])/2;
           mid2[1] = (list[i].Y[1]+list[i].Y[2])/2;
           mid2[2] = (list[i].Z[1]+list[i].Z[2])/2;
           mid3[0] = (list[i].X[0]+list[i].X[2])/2;
           mid3[1] = (list[i].Y[0]+list[i].Y[2])/2;
           mid3[2] = (list[i].Z[0]+list[i].Z[2])/2;
           output[4*i+0].X[0] = list[i].X[0];
           output[4*i+0].Y[0] = list[i].Y[0];
           output[4*i+0].Z[0] = list[i].Z[0];
           output[4*i+0].X[1] = mid1[0];
           output[4*i+0].Y[1] = mid1[1];
           output[4*i+0].Z[1] = mid1[2];
           output[4*i+0].X[2] = mid3[0];
           output[4*i+0].Y[2] = mid3[1];
           output[4*i+0].Z[2] = mid3[2];
           output[4*i+1].X[0] = list[i].X[1];
           output[4*i+1].Y[0] = list[i].Y[1];
           output[4*i+1].Z[0] = list[i].Z[1];
           output[4*i+1].X[1] = mid2[0];
           output[4*i+1].Y[1] = mid2[1];
           output[4*i+1].Z[1] = mid2[2];
           output[4*i+1].X[2] = mid1[0];
           output[4*i+1].Y[2] = mid1[1];
           output[4*i+1].Z[2] = mid1[2];
           output[4*i+2].X[0] = list[i].X[2];
           output[4*i+2].Y[0] = list[i].Y[2];
           output[4*i+2].Z[0] = list[i].Z[2];
           output[4*i+2].X[1] = mid3[0];
           output[4*i+2].Y[1] = mid3[1];
           output[4*i+2].Z[1] = mid3[2];
           output[4*i+2].X[2] = mid2[0];
           output[4*i+2].Y[2] = mid2[1];
           output[4*i+2].Z[2] = mid2[2];
           output[4*i+3].X[0] = mid1[0];
           output[4*i+3].Y[0] = mid1[1];
           output[4*i+3].Z[0] = mid1[2];
           output[4*i+3].X[1] = mid2[0];
           output[4*i+3].Y[1] = mid2[1];
           output[4*i+3].Z[1] = mid2[2];
           output[4*i+3].X[2] = mid3[0];
           output[4*i+3].Y[2] = mid3[1];
           output[4*i+3].Z[2] = mid3[2];
       }
       return output;
   }
   void DrawSphere()
   {
       int recursionLevel = 3;
       Triangle t;
       t.X[0] = 1;
       t.Y[0] = 0;
       t.Z[0] = 0;
       t.X[1] = 0;
       t.Y[1] = 1;
       t.Z[1] = 0;
       t.X[2] = 0;
       t.Y[2] = 0;
       t.Z[2] = 1;
       std::vector<Triangle> list;
       list.push_back(t);
       for (int r = 0 ; r < recursionLevel ; r++)
       {
           list = SplitTriangle(list);
       }

       // really draw `
       for (int octent = 0 ; octent < 8 ; octent++)
       {
           glPushMatrix();
           glRotatef(90*(octent%4), 1, 0, 0);
           if (octent >= 4)
               glRotatef(180, 0, 0, 1);
           glBegin(GL_TRIANGLES);
           for (unsigned int i = 0 ; i < list.size() ; i++)
           {
               for (int j = 0 ; j < 3 ; j++)
               {
                   double ptMag = sqrt(list[i].X[j]*list[i].X[j]+
                                       list[i].Y[j]*list[i].Y[j]+
                                       list[i].Z[j]*list[i].Z[j]);
                   glNormal3f(list[i].X[j]/ptMag, list[i].Y[j]/ptMag, list[i].Z[j]/ptMag);
                   glTexCoord2f(list[i].X[j]/ptMag, list[i].Y[j]/ptMag);
                   glVertex3f(list[i].X[j]/ptMag, list[i].Y[j]/ptMag, list[i].Z[j]/ptMag);
               }
           }
           glEnd();
           glPopMatrix();
       }
   }
   void White(void) { glColor3ub(255, 255, 255); };
   void Custom(void) { glColor3ub(220,0, 80); };
   void CustomTwo(void) { glColor3ub(120,49, 192); };
   void Yellow(void) { glColor3ub(250, 255, 0); };
   void ColorChange(void){ glColor3ub(Red, Green ,Blue); };
   
   void GetTextures(void){
	  rdr = vtkJPEGReader::New();
	  rdr->SetFileName("texture.jpg");
      rdr->Update();
      img = rdr->GetOutput();
      img->GetDimensions(dims);
      unsigned char * buffer = (unsigned char *) img->GetScalarPointer(0,0,0);
      glGenTextures(1,&texture);
      glBindTexture(GL_TEXTURE_2D, texture);
      glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, dims[0], dims[1], 0, GL_RGB, GL_UNSIGNED_BYTE, buffer);
      glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
      glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
      glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
      initialized = true;
	};
	
   virtual void RenderPiece(vtkRenderer *ren, vtkActor *act)
   {

      RemoveVTKOpenGLStateSideEffects();
      SetupLight();
      glEnable(GL_COLOR_MATERIAL);
      glMatrixMode(GL_MODELVIEW);
      
      glPushMatrix();
      //glTranslatef(20, 0, 0);
      
      glPushMatrix();
      //White();
      ColorChange();
      if(!initialized)
		GetTextures();
		
      
      
      glScalef(60, 60, 60);
      glEnable(GL_TEXTURE_2D);
      glBindTexture(GL_TEXTURE_2D,texture);
      DrawSphere();
      glDisable(GL_TEXTURE_2D);
      glPopMatrix();
      
      glPushMatrix();
      Custom();
      glScalef(20, 20, 20);
      glTranslatef(20, 0, 20);
      DrawSphere();
      glPopMatrix();
      
      glPushMatrix();
      CustomTwo();
      glScalef(25, 25, 25);
      glTranslatef(-26, -3, -20);
      DrawSphere();
      glPopMatrix();
      
      glPopMatrix();
      
	  glPushMatrix();
	  Yellow();
	  glTranslatef(400, 175, 0);
	  glScalef(20, 20, 20);
			
	  DrawSphere();
	  glPopMatrix();
	  
	  int x = 0;
	  double rad;
	  double angle;
	  double val = 1.0/100.0;
	  for(x;x<100;x++)
	  {
		  glPushMatrix();
		  White();
		  angle = 360*val*x;
		  rad = angle/360.0*2*3.14;
		  glTranslatef(500*cos(rad),aListOfStars.at(x),500*sin(rad));
		  DrawSphere();
		  glScalef(200, 200, 200);
		  glPopMatrix();
	  }
   }
};

vtkStandardNewMacro(vtk441MapperPart1);
int globalViewAngleX = 70;
int globalViewAngleY = 30;
class vtkTimerCallback : public vtkCommand
{
  public:
    static vtkTimerCallback *New()
    {
      vtkTimerCallback *cb = new vtkTimerCallback;
      cb->TimerCount = 0;
      cb->mapper = NULL;
      cb->renWin = NULL;
      cb->cam    = NULL;
      cb->angle  = 0;
      cb->ta = NULL;
      return cb;
    }

    void   SetMapper(vtk441Mapper *m) { mapper = m; };
    void   SetRenderWindow(vtkRenderWindow *rw) { renWin = rw; };
    void   SetCamera(vtkCamera *c) { cam = c; };
	void   SetTA(vtkTextActor *txA){ ta = txA;};
    virtual void Execute(vtkObject *vtkNotUsed(caller), unsigned long eventId, void *vtkNotUsed(callData))
    {
      // THIS IS WHAT GETS CALLED EVERY TIMER
      //cout << "Got a timer!!" << this->TimerCount << endl;
      // NOW DO WHAT EVER ACTIONS YOU WANT TO DO...
		if (vtkCommand::TimerEvent == eventId)
		{
			++this->TimerCount;
		}

      // Make a call to the mapper to make it alter how it renders...
      //if (mapper != NULL)
       //     mapper->IncrementSize();
		mapper->ChangeColor();
      // Modify the camera...
      if (cam != NULL)
      {
         //cam->SetFocalPoint(-30,0,0);

         float rads = angle/360.0*2*3.14;
         float angle2 = angle + 90;
         //angle2 = angle2 % 360;
         float rads2 = angle2/360.0*2*3.14;
         cam->SetPosition(200*cos(rads),0,200*sin(rads));
         cam->SetFocalPoint(globalViewAngleX*cos(rads2),globalViewAngleY,globalViewAngleX*sin(rads2));
         
         angle+=.5;
         if (angle > 360)
            angle = 0;
         cam->SetViewUp(0,1,0);
         cam->SetClippingRange(1, 200000);
         cam->SetDistance(2000);
      }

		time_t now = time(0);
		ta->SetInput ( std::ctime(&now));  
		ta->SetPosition( 30, 30 );  
		ta->GetTextProperty()->SetFontSize ( 80 );  
		ta->GetTextProperty()->SetColor ( 1.0, 0.0, 0.0 );  
		ta->Modified();


      // Force a render...
     if (renWin != NULL)
         renWin->Render();
    }

  private:
    int TimerCount;
    vtk441Mapper *mapper;
    vtkRenderWindow *renWin;
    vtkCamera *cam;
    float angle;
    
    vtkTextActor *ta;
};


void KeypressCallbackFunction (vtkObject* caller, long unsigned int eventId, void* clientData, void* callData );

int main()
{
  // Dummy input so VTK pipeline mojo is happy.
  //
  vtkSmartPointer<vtkSphereSource> sphere = vtkSmartPointer<vtkSphereSource>::New();
  sphere->SetThetaResolution(100);
  sphere->SetPhiResolution(50);

  // The mapper is responsible for pushing the geometry into the graphics
  // library. It may also do color mapping, if scalars or other attributes
  // are defined.
  //
  vtkSmartPointer<vtk441MapperPart1> win1Mapper = vtkSmartPointer<vtk441MapperPart1>::New();
  win1Mapper->SetInputConnection(sphere->GetOutputPort());

  vtkSmartPointer<vtkActor> win1Actor = vtkSmartPointer<vtkActor>::New();
  win1Actor->SetMapper(win1Mapper);

  vtkSmartPointer<vtkRenderer> ren1 = vtkSmartPointer<vtkRenderer>::New();


  vtkSmartPointer<vtkRenderWindow> renWin = vtkSmartPointer<vtkRenderWindow>::New();
  renWin->AddRenderer(ren1);
  //ren1->SetViewport(0, 0, 0.5, 1);

  vtkSmartPointer<vtkRenderWindowInteractor> iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();
  iren->SetRenderWindow(renWin);

  // Add the actors to the renderer, set the background and size.
  //
  bool doWindow1 = true;
  if (doWindow1)
     ren1->AddActor(win1Actor);
  ren1->SetBackground(0.0, .0, .0);
  renWin->SetSize(1920, 1080);

	ren1->GetActiveCamera()->SetFocalPoint(0,0,0);
	ren1->GetActiveCamera()->SetPosition(0,0,70);
	ren1->GetActiveCamera()->SetViewUp(0,1,0);
	//ren1->GetActiveCamera()->SetClippingRange(, 2000);
	ren1->GetActiveCamera()->SetDistance(200);
	
	
	vtkSmartPointer<vtkTextActor> textActor = vtkSmartPointer<vtkTextActor>::New(); 
	time_t t;
	textActor->SetInput ( std::ctime(&t));
	//textActor->SetPosition2 ( 50, 40 );  
	textActor->GetTextProperty()->SetFontSize ( 120 );  
	textActor->GetTextProperty()->SetColor ( 1.0, 0.0, 0.0 );  
	ren1->AddActor2D(textActor);
  // This starts the event loop and invokes an initial render.
  //
  ((vtkInteractorStyle *)iren->GetInteractorStyle())->SetAutoAdjustCameraClippingRange(0);
  iren->Initialize();

  // Sign up to receive TimerEvent
  vtkSmartPointer<vtkTimerCallback> cb = vtkSmartPointer<vtkTimerCallback>::New();
  iren->AddObserver(vtkCommand::TimerEvent, cb);
  cb->SetMapper(win1Mapper);
  cb->SetRenderWindow(renWin);
  cb->SetCamera(ren1->GetActiveCamera());
  cb->SetTA(textActor);
  vtkSmartPointer<vtkCallbackCommand> keypressCallback = vtkSmartPointer<vtkCallbackCommand>::New();
  keypressCallback->SetCallback ( KeypressCallbackFunction );
  iren->AddObserver ( vtkCommand::KeyPressEvent, keypressCallback );

  int timerId = iren->CreateRepeatingTimer(5);  // repeats every 10 microseconds <--> 0.01 seconds*/
  

  
  
  iren->Start();

  return EXIT_SUCCESS;
}

void KeypressCallbackFunction ( vtkObject* caller, long unsigned int vtkNotUsed(eventId), void* vtkNotUsed(clientData), void* vtkNotUsed(callData) )
{
  std::cout << "Keypress callback" << std::endl;

  vtkRenderWindowInteractor *iren = static_cast<vtkRenderWindowInteractor*>(caller);

  std::cout << "Pressed: " << iren->GetKeySym() << std::endl;
	if(!strcmp(iren->GetKeySym(),"Up"))
		globalViewAngleY++;
	if(!strcmp(iren->GetKeySym(),"Down"))
		globalViewAngleY--;
	if(!strcmp(iren->GetKeySym(),"Left"))
		globalViewAngleX++;
	if(!strcmp(iren->GetKeySym(),"Right"))
		globalViewAngleX--;
}




