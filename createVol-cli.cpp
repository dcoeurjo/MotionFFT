#include <iostream>
#include <algorithm>
#include <DGtal/base/Common.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/shapes/GaussDigitizer.h>
#include <DGtal/shapes/ShapeFactory.h>
#include <DGtal/io/writers/VolWriter.h>

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


template <typename Image>
struct RealPartToChar{
  RealPartToChar() {}
  
  unsigned char operator()(const typename FFT< Image >::Complex &complex) const
  {
    return static_cast<unsigned char >( std::abs((double)complex.real() ));
  }
};
template <typename Image>

struct ImgPartToChar{
  ImgPartToChar() {}
  unsigned char operator()(const typename FFT< Image  >::Complex &complex) const
  {
    return static_cast<unsigned char >( std::abs(complex.imag() ));
  }
};
template <typename Image>
struct MagPartToChar{
  MagPartToChar() {}
  unsigned char operator()(const typename FFT< Image >::Complex &complex) const
  {
    return static_cast<unsigned char >( std::abs(complex) );
  }
};


int main(int argc, char ** argv)
{
  // parse command line ----------------------------------------------
  po::options_description general_opt ( "Allowed options are: " );
  general_opt.add_options()
  ( "help,h", "display this message." )
  ( "resolution,r", po::value<double>(),"Spatial resolution." )
  ( "frame,f", po::value<unsigned int>(),"Temporal resolution (number of frames)." )
  ( "output,o", po::value<string>(),"Output vol filename." );
    
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
    << "\tcreateVol  --o <volOutputFileName>  ...;TODO....."<<std::endl
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
  
 
  int previd = -1; //to reduce the number of calls to customcolor (viewer bug)
 
  
  //Map->Image->FFT
  Point  plow(0,0,0) , pup(0,0,0);
  for(Map3D::const_iterator it= myFinalSet.begin(), itend=myFinalSet.end(); it!= itend; ++it)
  {
    plow = plow.inf((*it).first);
    pup = pup.sup((*it).first);
  }
  
  //Square domain
  int maxp = pup[0]-plow[0];
  if ((pup[1]-plow[1])>maxp) maxp=(pup[1]-plow[1]);
  if ((pup[2]-plow[2])>maxp) maxp=(pup[2]-plow[2]);
  pup = plow + Point(maxp,maxp,maxp);
  Domain domain(plow,pup);
  trace.info() << "Overall domain = "<<domain<<std::endl;
  
  typedef ImageContainerBySTLVector<Domain, unsigned char> Image;
  Image image(domain);
  for(Map3D::const_iterator it= myFinalSet.begin(), itend=myFinalSet.end(); it!= itend; ++it)
    image.setValue( it->first, it->second);

  //just an export of the space/time volume
  VolWriter<Image  >::exportVol(outputFileName+".vol", image);
  
  
  //FFT
  ImageContainerBySTLVector<Domain, unsigned char> image2(domain);

  typedef FFT< Image> FFT3D;
  FFT3D fft(image);
 
  trace.beginBlock("Computing FFT");
  FFT3D::ComplexImage fftresult(domain);
  fft.compute(fftresult);
  
  //Exports
  VolWriter<FFT3D::ComplexImage, RealPartToChar<Image> >::exportVol(outputFileName+"-real.vol", fftresult, RealPartToChar<Image>() );
  VolWriter<FFT3D::ComplexImage, ImgPartToChar<Image> >::exportVol(outputFileName+"-imag.vol", fftresult, ImgPartToChar<Image>() );
  VolWriter<FFT3D::ComplexImage, MagPartToChar<Image> >::exportVol(outputFileName+"-mag.vol", fftresult, MagPartToChar<Image>() );
  trace.endBlock();
 
  
  //For debuging FFFT3d
  Domain dom2(Point(0,0,0),Point(16,16,16));
  Image img2(dom2);
  img2.setValue(Point(8,8,8), 128);
  FFT3D fft2(img2);
  FFT3D::ComplexImage fftres2(dom2);
  fft2.compute(fftres2);
  VolWriter<FFT3D::ComplexImage, RealPartToChar<Image> >::exportVol(outputFileName+"-2-real.vol", fftres2, RealPartToChar<Image>() );
  VolWriter<FFT3D::ComplexImage, ImgPartToChar<Image> >::exportVol(outputFileName+"-2-imag.vol", fftres2, ImgPartToChar<Image>() );
  VolWriter<FFT3D::ComplexImage, MagPartToChar<Image> >::exportVol(outputFileName+"-2-mag.vol", fftres2, MagPartToChar<Image>() );
  
  
  //End
  
  
  return 0;
  
}

