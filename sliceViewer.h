/**
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 **/
/**
 * @file sliceViewer.h
 * @author Bertrand Kerautret (\c kerautre@loria.fr )
 * LORIA (CNRS, UMR 7503), University of Nancy, France
 *
 * @date 2014/07/04
 *
 * An example file named qglViewer.
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#ifndef MAINWINDOW_H
#define MAINWINDOW_H
#ifndef Q_MOC_RUN 
   #include "DGtal/io/viewers/Viewer3D.h"
   #include <DGtal/images/ImageContainerBySTLVector.h>
   #include "DGtal/helpers/StdDefs.h"
   #include "DGtal/images/ConstImageAdapter.h"
#endif 

#include <QMainWindow>

namespace Ui {
class MainWindow;

}

class MainWindow : public QMainWindow
{
  Q_OBJECT
typedef DGtal::ImageContainerBySTLVector < DGtal::Z3i::Domain, unsigned char > Image3D;
typedef DGtal::ImageContainerBySTLVector < DGtal::Z2i::Domain, unsigned char > Image2D;
typedef DGtal::ConstImageAdapter<Image3D, Image2D::Domain, DGtal::functors::Projector< DGtal::Z3i::Space>,
                                 Image3D::Value,  DGtal::functors::Identity >  SliceImageAdapter;
  
public:
explicit MainWindow(DGtal::Viewer3D<> *viewer,  DGtal::ImageContainerBySTLVector < DGtal::Z3i::Domain, unsigned char > *myImage3D,
                      QWidget *parent = 0, Qt::WindowFlags flags=0);
  ~MainWindow();
  void setImageProjX(const QPixmap &aPixMap);
  void setImageProjY(const QPixmap &aPixMap);
  void setImageProjZ(const QPixmap &aPixMap);


  void updateSliceImageX(unsigned int sliceNumber, bool init);
  void updateSliceImageY(unsigned int sliceNumber, bool init);
  void updateSliceImageZ(unsigned int sliceNumber, bool init);

  void updateZoomImageX(unsigned int sliceNumber, double gridSize);
  void updateZoomImageY(unsigned int sliceNumber, double gridSize);
  void updateZoomImageZ(unsigned int sliceNumber, double gridSize);



public slots:
  void updateSliceImageX();
  void updateSliceImageY();
  void updateSliceImageZ();

  void updateZoomImageX();
  void updateZoomImageY();
  void updateZoomImageZ();

  void setScale1_1_ImageX();
  void setScale1_1_ImageY();
  void setScale1_1_ImageZ();

  
private:
  Ui::MainWindow *ui;
  DGtal::Viewer3D<> *myViewer;  
  Image3D *myImage3D;
};

#endif // MAINWINDOW_H
  
