#include <iostream>
#include <algorithm>
#include <DGtal/base/Common.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/shapes/GaussDigitizer.h>
#include <DGtal/io/viewers/Viewer3D.h>
#include <DGtal/shapes/ShapeFactory.h>
#include <DGtal/io/writers/VolWriter.h>

#include <QtGui/qapplication.h>
#include <DGtal/images/ImageContainerBySTLMap.h>
#include <DGtal/io/colormaps/HueShadeColorMap.h>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>


#include "FFT.h"


using namespace std;
using namespace DGtal;
using namespace Z3i;

namespace po = boost::program_options;

/**
 * Missing parameter error message.
 *
 * @param param
 */
void missingParam ( std::string param )
{
  trace.error() <<" Parameter: "<<param<<" is required..";
  trace.info() <<std::endl;
  exit ( 1 );
}


template <typename Shape2D, typename Set2D>
void createShape(const Shape2D & aShape,
              const double resolution,
              Set2D &aSet)
{
  typedef GaussDigitizer<Z2i::Space,Shape2D> Gauss;
  Gauss dig;
  dig.attach( aShape );
  dig.init( aShape.getLowerBound()+Z2i::Vector(-1,-1),
            aShape.getUpperBound()+Z2i::Vector(1,1), resolution );
  Z2i::Domain domain = dig.getDomain();
  for(Z2i::Domain::ConstIterator it = domain.begin() ; it != domain.end();
      ++it)
  {
    Z2i::Point P = *it;
    if (dig(P))
      aSet.insert( P );
  }
}




template <typename Set2D, typename Map3D>
void addLinearMotion(const Set2D & digitalObject2D,
                     const double resolution,
                     const unsigned int startFrame,
                     const unsigned int endFrame,
                     const Z2i::RealVector initialPos,
                     const Z2i::RealVector finalPos,
                     Map3D &aMap,
                     const unsigned int id)
{
  Z2i::Point shift;
  ASSERT(startFrame<endFrame);
  Z2i::RealVector delta = (finalPos - initialPos) /  (double)(endFrame - startFrame);
 
  for(unsigned int i = startFrame; i < endFrame ; ++i)
  {
    shift = 1.0/resolution * ( initialPos +  (double)(i-startFrame)  * delta );
    for(typename Set2D::const_iterator it = digitalObject2D.begin(), itend=digitalObject2D.end();
        it != itend; ++it)
      aMap [ Point( (*it)[0] + shift[0], (*it)[1] + shift[1],  i ) ] = id;
  }
}


int main(int argc, char ** argv)
{
  // parse command line ----------------------------------------------
  po::options_description general_opt ( "Allowed options are: " );
  general_opt.add_options()
  ( "help,h", "display this message." )
  ( "resolution,r", po::value<double>(),"Spatial resolution." )
  ( "frame,f", po::value<unsigned int>(),"Temporal resolution (number of frames)." )
  ( "output,o", po::value<string>(),"Output vol filename." );
  
  QApplication application(argc,argv);
  
  bool parseOK=true;
  
  po::variables_map vm;
  try{
    po::store(po::parse_command_line(argc, argv, general_opt), vm);
  }catch(const std::exception& ex){
    parseOK=false;
    trace.info()<< "Error checking program options: "<< ex.what()<< endl;
  }
  
  po::notify ( vm );
  if ( !parseOK || vm.count ( "help" ) ||argc<=1 )
  {
    trace.info() << "Create space-time vol"<<std::endl
    << std::endl << "Basic usage: "<<std::endl
    << "\tcreateVol  --o <volOutputFileName>  ...;TODO....."<<std::endl
    << general_opt << "\n";
    return 0;
  }
  //Parse options
  if ( ! ( vm.count ( "resolution" ) ) ) missingParam ( "--resolution" );
  const double resolution = vm["resolution"].as<double>();
  if ( ! ( vm.count ( "frame" ) ) ) missingParam ( "--frame" );
  const unsigned int  nbFrame = vm["frame"].as<unsigned int>();
  if ( ! ( vm.count ( "output" ) ) ) missingParam ( "--output" );
  const std::string outputFileName = vm["output"].as<std::string>();

  
  //2D shape
  typedef Flower2D<Z2i::Space> Shape2D;
  Shape2D flower( 0.0,0.0, 20.0, 5.0, 5 , 0.0);
  
  //3D DigitalSet
  typedef std::set<Z2i::Point> Set2D;
  typedef std::map<Point, unsigned int> Map3D;
  Set2D my2Dobject;
  
  createShape<Shape2D,Set2D>(flower, resolution, my2Dobject);
  
  Map3D myFinalSet;
  addLinearMotion(my2Dobject, resolution, 0, nbFrame/2.0, Z2i::RealPoint(0.0,0.0), Z2i::RealPoint(25.0,25.0), myFinalSet, 1);
  addLinearMotion(my2Dobject, resolution, nbFrame/2.0, nbFrame, Z2i::RealPoint(25.0,25.0),Z2i::RealPoint(0.0,0.0),  myFinalSet, 30);
  addLinearMotion(my2Dobject, resolution, 0, nbFrame, Z2i::RealPoint(-25.-25,0.0), Z2i::RealPoint(25.0,25.0), myFinalSet, 50);

  trace.info() << "Scence created."<<std::endl;
  trace.info() << "Sending it to the viewer..."<<std::endl;
  
  //Colormap
  HueShadeColorMap<int, 1> cmap(0,100);
  
  //Viewer 3D
  Viewer3D<> viewer;
  viewer.show();
  int previd = -1; //to reduce the number of calls to customcolor (viewer bug)
  for(Map3D::const_iterator it= myFinalSet.begin(), itend=myFinalSet.end(); it!= itend; ++it)
  {
    int id =(*it).second;
   /* if (id != previd)
    {   viewer << CustomColors3D( cmap((*it).second) , cmap((*it).second));
    previd = id;
  }*/
  
    viewer << (*it).first;
  }
  viewer << Viewer3D<>::updateDisplay;

  
  //Map->Image->FFT
  Point  plow(0,0,0) , pup(0,0,0);
  for(Map3D::const_iterator it= myFinalSet.begin(), itend=myFinalSet.end(); it!= itend; ++it)
  {
    plow = plow.inf((*it).first);
    pup = pup.sup((*it).first);
  }
  Domain domain(plow,pup);
  trace.info() << "Overall domain = "<<domain<<std::endl;
  
  ImageContainerBySTLVector<Domain, uint32_t> image(domain);
  for(Map3D::const_iterator it= myFinalSet.begin(), itend=myFinalSet.end(); it!= itend; ++it)
    image.setValue( it->first, it->second);

  //just an export of the space/time volume
  VolWriter<ImageContainerBySTLVector<Domain, uint32_t> >::exportVol(outputFileName, image);
  
  
  //FFT
  typedef FFT< ImageContainerBySTLVector<Domain, uint32_t> > FFT3D;
  FFT3D fft(image);
 
  trace.beginBlock("Computing FFT");
  FFT3D::ComplexImage fftresult(domain);
  fft.compute(fftresult);
  trace.endBlock();
  
  return application.exec();
  
}

