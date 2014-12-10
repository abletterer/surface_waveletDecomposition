#include "surface_waveletDecomposition.h"

#include "mapHandler.h"

#define NORMALE 5

namespace CGoGN
{

namespace SCHNApps
{

bool Surface_WaveletDecomposition_Plugin::enable()
{
    m_waveletDecompositionDialog = new Dialog_Surface_WaveletDecomposition(m_schnapps);

    m_waveletDecompositionAction = new QAction("Wavelet Decomposition", this);

    m_schnapps->addMenuAction(this, "Surface;Wavelet Decomposition", m_waveletDecompositionAction);

    m_colorPerVertexShader = new CGoGN::Utils::ShaderColorPerVertex();
    registerShader(m_colorPerVertexShader);

    m_positionVBO = new Utils::VBO();
    m_colorVBO = new Utils::VBO();

    m_toDraw = false;

    connect(m_waveletDecompositionAction, SIGNAL(triggered()), this, SLOT(openWaveletDecompositionDialog()));

    connect(m_waveletDecompositionDialog->button_cancel, SIGNAL(clicked()), this, SLOT(closeWaveletDecompositionDialog()));
    connect(m_waveletDecompositionDialog->button_decompose, SIGNAL(clicked()), this, SLOT(decomposeFromDialog()));
    connect(m_waveletDecompositionDialog->button_saveImages, SIGNAL(clicked()), this, SLOT(saveImagesFromDialog()));

    return true;
}

void Surface_WaveletDecomposition_Plugin::disable()
{
    delete m_colorPerVertexShader;
    delete m_positionVBO;
    delete m_colorVBO;

    disconnect(m_waveletDecompositionAction, SIGNAL(triggered()), this, SLOT(openWaveletDecompositionDialog()));

    disconnect(m_waveletDecompositionDialog->button_cancel, SIGNAL(clicked()), this, SLOT(closeWaveletDecompositionDialog()));
    disconnect(m_waveletDecompositionDialog->button_decompose, SIGNAL(clicked()), this, SLOT(decomposeFromDialog()));
    disconnect(m_waveletDecompositionDialog->button_saveImages, SIGNAL(clicked()), this, SLOT(saveImagesFromDialog()));
}

void Surface_WaveletDecomposition_Plugin::drawMap(View* view, MapHandlerGen *map)
{
//    if(m_toDraw && m_schnapps->getSelectedView() == view)
//    {
//        //If VBO are initialized
//        glPolygonMode(GL_FRONT, GL_FILL);
//        glEnable(GL_LIGHTING);
//        glEnable(GL_POLYGON_OFFSET_FILL);
//        m_colorPerVertexShader->setAttributePosition(m_positionVBO);
//        m_colorPerVertexShader->setAttributeColor(m_colorVBO);
//        m_colorPerVertexShader->setOpacity(1.0);
//        map->draw(m_colorPerVertexShader, CGoGN::Algo::Render::GL2::TRIANGLES);
//        glDisable(GL_POLYGON_OFFSET_FILL);
//    }
}

void Surface_WaveletDecomposition_Plugin::openWaveletDecompositionDialog()
{
    m_waveletDecompositionDialog->show();
}

void Surface_WaveletDecomposition_Plugin::closeWaveletDecompositionDialog()
{
    m_waveletDecompositionDialog->close();
}

void Surface_WaveletDecomposition_Plugin::decomposeFromDialog()
{
    QList<QListWidgetItem*> currentItems = m_waveletDecompositionDialog->list_maps->selectedItems();
    if(!currentItems.empty() && m_decomposition)
    {
        const QString& mapName = currentItems[0]->text();

        MapHandler<PFP2>* mh_map = static_cast<MapHandler<PFP2>*>(m_schnapps->getMap(mapName));

        if(mh_map && m_decomposition)
        {
            PFP2::MAP* map = mh_map->getMap();

            //We search the last level of decomposition
            Decomposition* decomposition = m_decomposition;
            while(decomposition->getChild())
            {
                decomposition = decomposition->getChild();
            }

            while(decomposition->getCorrectionX()>2 && decomposition->getCorrectionY()>2)
            {
                //If a decomposition is still feasible
                decomposition = decomposition->addChild();

                Decomposition* parent = decomposition->getParent();

                QImage& image_parent = parent->getImage();

                int image_width = (image_parent.width()+image_parent.width()%2)/2;
                int image_height = (image_parent.height()+image_parent.height()%2)/2;

                //Creation of a new image composed of the even of parent image
                QImage image(image_width, image_height, image_parent.format());

                for(int i = 0; i < image_width; ++i)
                {
                    for(int j = 0; j < image_height; ++j)
                    {
                        image.setPixel(i, j, image_parent.pixel(i*2, j*2));
                    }
                }

                decomposition->setImage(image);

                //Creation of a matrix composed of the difference between odd of parent image and prediction value (linear interpolation)
                QRgb cur_pixel;
                NQRgb hori_pixel, vert_pixel;
                for(int i = 0; i < image_parent.width(); ++i)
                {
                    for(int j = 0; j < image_parent.height(); ++j)
                    {
                        cur_pixel = image_parent.pixel(i, j);

                        hori_pixel = NQRgb();
                        vert_pixel = NQRgb();

                        if(i%2==1)
                        {
                            if(i==image_parent.width()-1)
                            {
                                //Mirror effect at the border of the image
                                hori_pixel.setRed(qRed(cur_pixel) - qRed(image_parent.pixel(i-1,j)));
                                hori_pixel.setGreen(qGreen(cur_pixel) - qGreen(image_parent.pixel(i-1,j)));
                                hori_pixel.setBlue(qBlue(cur_pixel) - qBlue(image_parent.pixel(i-1,j)));
                            }
                            else
                            {
                                hori_pixel.setRed(qRed(cur_pixel) - 0.5*(qRed(image_parent.pixel(i-1,j))+qRed(image_parent.pixel(i+1,j))));
                                hori_pixel.setGreen(qGreen(cur_pixel) - 0.5*(qGreen(image_parent.pixel(i-1,j))+qGreen(image_parent.pixel(i+1,j))));
                                hori_pixel.setBlue(qBlue(cur_pixel) - 0.5*(qBlue(image_parent.pixel(i-1,j))+qBlue(image_parent.pixel(i+1,j))));
                            }
                            //Odd column number
                            if(j%2==1)
                            {
                                //Odd line number
                                if(j==image_parent.height()-1)
                                {
                                    //Mirror effect at the border of the image
                                    vert_pixel.setRed(qRed(cur_pixel) - qRed(image_parent.pixel(i,j-1)));
                                    vert_pixel.setGreen(qGreen(cur_pixel) - qGreen(image_parent.pixel(i,j-1)));
                                    vert_pixel.setBlue(qBlue(cur_pixel) - qBlue(image_parent.pixel(i,j-1)));
                                }
                                else
                                {
                                    vert_pixel.setRed(qRed(cur_pixel) - 0.5*(qRed(image_parent.pixel(i,j-1))+qRed(image_parent.pixel(i,j+1))));
                                    vert_pixel.setGreen(qGreen(cur_pixel) - 0.5*(qGreen(image_parent.pixel(i,j-1))+qGreen(image_parent.pixel(i,j+1))));
                                    vert_pixel.setBlue(qBlue(cur_pixel) - 0.5*(qBlue(image_parent.pixel(i,j-1))+qBlue(image_parent.pixel(i,j+1))));
                                }
                            }
                        }
                        else
                        {
                            if(j%2==1)
                            {
                                //Odd line number
                                if(j==image_parent.height()-1)
                                {
                                    //Mirror effect at the border of the image
                                    vert_pixel.setRed(qRed(cur_pixel) - qRed(image_parent.pixel(i,j-1)));
                                    vert_pixel.setGreen(qGreen(cur_pixel) - qGreen(image_parent.pixel(i,j-1)));
                                    vert_pixel.setBlue(qBlue(cur_pixel) - qBlue(image_parent.pixel(i,j-1)));
                                }
                                else
                                {
                                    vert_pixel.setRed(qRed(cur_pixel) - 0.5*(qRed(image_parent.pixel(i,j-1))+qRed(image_parent.pixel(i,j+1))));
                                    vert_pixel.setGreen(qGreen(cur_pixel) - 0.5*(qGreen(image_parent.pixel(i,j-1))+qGreen(image_parent.pixel(i,j+1))));
                                    vert_pixel.setBlue(qBlue(cur_pixel) - 0.5*(qBlue(image_parent.pixel(i,j-1))+qBlue(image_parent.pixel(i,j+1))));
                                }
                            }
                        }

                        parent->setHorizontalCorrection(i, j, hori_pixel);
                        parent->setVerticalCorrection(i, j, vert_pixel);
                    }
                }
                CGoGNout << "New level of decomposition created" << CGoGNendl;
            }
        }
    }
}

void Surface_WaveletDecomposition_Plugin::saveImagesFromDialog()
{
    QList<QListWidgetItem*> currentItems = m_waveletDecompositionDialog->list_maps->selectedItems();
    if(!currentItems.empty() && m_decomposition)
    {
        const QString& mapName = currentItems[0]->text();

        Decomposition* decomposition = m_decomposition;
        do
        {
            QImage image = decomposition->getImage();
            QString filename("/home/blettere/Projets/Models/Decomposition/");
            filename.append(mapName);
            filename.append("-");
            QString filename2(filename);
            filename.append(QString::number(decomposition->getLevel()));
            filename.append(".png");

            if(!image.save(filename))
            {
                CGoGNerr << "Image '" << filename.toStdString() << "' has not been saved" << CGoGNendl;
            }

            filename2.append(QString::number(decomposition->getLevel()+1));
            filename = filename2;
            filename.append("-horizontalCoherence.png");
            filename2.append("-verticalCoherence.png");

            decomposition = decomposition->getChild();

            if(decomposition)
            {
                QImage image2(decomposition->getCorrectionX(), decomposition->getCorrectionY(), image.format());
                image = QImage(decomposition->getCorrectionX(), decomposition->getCorrectionY(), image2.format());
                NQRgb correction;
                for(int i = 0; i < decomposition->getCorrectionX(); ++i)
                {
                    for(int j = 0; j < decomposition->getCorrectionY(); ++j)
                    {
                        correction = decomposition->getHorizontalCorrection(i, j);
                        image.setPixel(i, j, qRgb(qAbs(correction.getRed()), qAbs(correction.getGreen()), qAbs(correction.getBlue())));
                        correction = decomposition->getVerticalCorrection(i, j);
                        image2.setPixel(i, j, qRgb(qAbs(correction.getRed()), qAbs(correction.getGreen()), qAbs(correction.getBlue())));
                    }
                }

                if(!image.save(filename))
                {
                    CGoGNerr << "Image '" << filename2.toStdString() << "' has not been saved" << CGoGNendl;
                }

                if(!image2.save(filename2))
                {
                    CGoGNerr << "Image '" << filename2.toStdString() << "' has not been saved" << CGoGNendl;
                }
            }
        } while(decomposition);
    }
}

MapHandlerGen* Surface_WaveletDecomposition_Plugin::initializeObject(const QString& view, QString& filename, const bool multiple)
{
    if(!filename.isEmpty())
    {
        QString file = filename.left(filename.lastIndexOf('.'));
        file = file.mid(file.lastIndexOf('/')).remove(0, 1);
        MapHandlerGen* mhg_map = m_schnapps->addMap(file, 2);
        MapHandler<PFP2>* mh_map = static_cast<MapHandler<PFP2>*>(mhg_map);
        PFP2::MAP* map = mh_map->getMap();

        VertexAttribute<PFP2::VEC3, PFP2::MAP> position = mh_map->getAttribute<PFP2::VEC3, VERTEX>("position");
        if(!position.isValid())
        {
            position = mh_map->addAttribute<PFP2::VEC3, VERTEX>("position");
        }

        VertexAttribute <PFP2::VEC3, PFP2::MAP> color = mh_map->getAttribute<PFP2::VEC3, VERTEX>("color");
        if(!color.isValid())
        {
            color = mh_map->addAttribute<PFP2::VEC3, VERTEX>("color");
        }

        QString extension = filename.mid(filename.lastIndexOf('.'));
        extension.toUpper();

        QImage image;
        if(!image.load(filename, extension.toUtf8().constData()))
        {
            CGoGNout << "Image has not been loaded correctly" << CGoGNendl;
            return NULL;
        }

        image = image.convertToFormat(QImage::Format_RGB32);

        int imageX = image.width(), imageY = image.height();

        Algo::Surface::Tilings::Square::Grid<PFP2> grid(*map, imageX-1, imageY-1);
        if(multiple)
        {
            QImage image_parent;
            QString filename_big_parent = filename, filename_parent;

            filename_big_parent.replace(filename_big_parent.size()-5, 1, "0");
            filename_parent.replace(filename_parent.size()-5, 1, file[file.size()-1]);
            //CORRECTION PAIR FAIT PERDRE 1 PIXEL A DROITE

            if(!image_parent.load(filename_big_parent, extension.toUtf8().constData()))
            {
                CGoGNout << "Image has not been loaded correctly" << CGoGNendl;
                return NULL;
            }

            const int level = file[file.size()-1].digitValue();

            grid.embedIntoGrid(position, image_parent.width()-1+(image.width()%2==0?1:0)*level, image_parent.height()-1-(image.height()%2)*level);
        }
        else
        {
            grid.embedIntoGrid(position, imageX-1, imageY-1);
        }

        std::vector<Dart> vDarts = grid.getVertexDarts();

        QRgb pixel;

        for(int i = 0; i < imageX; ++i)
        {
            for(int j = 0; j < imageY; ++j)
            {
                pixel = image.pixel(i,(imageY-j)-1);
                color[vDarts[j*imageX+i]] = PFP2::VEC3(qRed(pixel)/255.f, qGreen(pixel)/255.f, qBlue(pixel)/255.f);
            }
        }

        m_decomposition = new Decomposition(imageX, imageY, 0);

        m_decomposition->setImage(image);

        mh_map->updateBB(position);
        mh_map->notifyAttributeModification(position);
        mh_map->notifyAttributeModification(color);
        mh_map->notifyConnectivityModification();
        m_positionVBO->updateData(position);
        //m_colorVBO->updateData(color);

        m_toDraw = true;

        map->enableQuickTraversal<PFP2::MAP, VERTEX>();

        if(view == m_schnapps->getSelectedView()->getName())
        {
            m_schnapps->getSelectedView()->updateGL();
        }

        return mhg_map;
    }
    return NULL;
}


#ifndef DEBUG
Q_EXPORT_PLUGIN2(Surface_WaveletDecomposition_Plugin, Surface_WaveletDecomposition_Plugin)
#else
Q_EXPORT_PLUGIN2(Surface_WaveletDecomposition_PluginD, Surface_WaveletDecomposition_Plugin)
#endif

} // namespace SCHNApps

} // namespace CGoGN
