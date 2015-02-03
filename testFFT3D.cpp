#include <iostream>
#include <algorithm>
#include <DGtal/base/Common.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/shapes/GaussDigitizer.h>
#include <DGtal/io/viewers/Viewer3D.h>
#include <DGtal/shapes/ShapeFactory.h>
#include <DGtal/io/writers/VolWriter.h>


#include "FFT.h"
#include "IFFT.h"


using namespace DGtal;
using namespace Z3i;


template <typename Image>
struct MagPartToChar{
  MagPartToChar() {}
  unsigned char operator()(const typename FFT< Image >::Complex &complex) const
  {
    return static_cast<unsigned char >( std::abs(complex) );
  }
};


/**
 * Testing DGtal/FFTW API
 *
 */
int main()
{
  
  typedef ImageContainerBySTLVector<Domain, unsigned char> Image;

  //Source
  Domain domain(Point(0,0,0), Point(512,512,512));
                                    
  Image image(domain);
  
  typedef ImplicitBall<Z3i::Space> Shape3D;
  Shape3D aShape( Point(64,64,64), 20.0);
  typedef GaussDigitizer<Z3i::Space,Shape3D> Gauss;
  Gauss dig;
  dig.attach( aShape );
  dig.init( aShape.getLowerBound()+Z3i::Vector(-1,-1,-1),
           aShape.getUpperBound()+Z3i::Vector(1,1,1), 1.0 );
  Z3i::Domain d3D = dig.getDomain();
  for(Z3i::Domain::ConstIterator it = d3D.begin() ; it != d3D.end();
      ++it)
  {
    Z3i::Point P = *it;
    if (dig(P))
      image.setValue( P , 128 );
  }
  
  
  //FFT
  typedef FFT< Image > FFT3D;
  FFT3D fft(image);
  trace.beginBlock("Computing FFT");
  FFT3D::ComplexImage fftresult(domain);
  fft.compute(fftresult);
  trace.endBlock();
  
  
  //iFFT
  typedef IFFT<FFT3D::ComplexImage> IFFT3D;
  Image imagereconstructed(domain);
  IFFT3D ifft(fftresult);
  
  trace.beginBlock("Computing IFFT");
  ifft.compute(imagereconstructed);
  trace.endBlock();
  
  
  //just an export of the reconstructed image
  VolWriter<Image  >::exportVol("original.vol", image);
  VolWriter<Image  >::exportVol("inverse.vol", imagereconstructed);
  
  VolWriter<FFT3D::ComplexImage, MagPartToChar<Image> >::exportVol("magnitude.vol", fftresult, MagPartToChar<Image>() );
    
  
  trace.info()<< "Export done."<<std::endl;

  
}