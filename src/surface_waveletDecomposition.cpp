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
    connect(m_waveletDecompositionDialog->button_synthesize, SIGNAL(clicked()), this, SLOT(synthesizeFromDialog()));
    connect(m_waveletDecompositionDialog->button_analyze, SIGNAL(clicked()), this, SLOT(analyzeFromDialog()));

    return true;
}

void Surface_WaveletDecomposition_Plugin::disable()
{
    delete m_colorPerVertexShader;
    delete m_positionVBO;
    delete m_colorVBO;

    disconnect(m_waveletDecompositionAction, SIGNAL(triggered()), this, SLOT(openWaveletDecompositionDialog()));

    disconnect(m_waveletDecompositionDialog->button_cancel, SIGNAL(clicked()), this, SLOT(closeWaveletDecompositionDialog()));
    connect(m_waveletDecompositionDialog->button_synthesize, SIGNAL(clicked()), this, SLOT(synthesizeFromDialog()));
    connect(m_waveletDecompositionDialog->button_analyze, SIGNAL(clicked()), this, SLOT(analyzeFromDialog()));
}

void Surface_WaveletDecomposition_Plugin::drawMap(View* view, MapHandlerGen *map)
{
    if(m_toDraw && m_schnapps->getSelectedView() == view && map->getName()=="Image")
    {
        //If VBO are initialized
        glPolygonMode(GL_FRONT, GL_FILL);
        glEnable(GL_LIGHTING);
        glEnable(GL_POLYGON_OFFSET_FILL);
        m_colorPerVertexShader->setAttributePosition(m_positionVBO);
        m_colorPerVertexShader->setAttributeColor(m_colorVBO);
        m_colorPerVertexShader->setOpacity(1.0);
        map->draw(m_colorPerVertexShader, CGoGN::Algo::Render::GL2::TRIANGLES);
        glDisable(GL_POLYGON_OFFSET_FILL);
    }
}

void Surface_WaveletDecomposition_Plugin::openWaveletDecompositionDialog()
{
    m_waveletDecompositionDialog->show();
}

void Surface_WaveletDecomposition_Plugin::closeWaveletDecompositionDialog()
{
    m_waveletDecompositionDialog->close();
}

void Surface_WaveletDecomposition_Plugin::synthesizeFromDialog()
{
    QList<QListWidgetItem*> currentItems = m_waveletDecompositionDialog->list_maps->selectedItems();
    if(!currentItems.empty() && m_decomposition)
    {
        const QString& mapName = currentItems[0]->text();
        const QString& positionName = m_waveletDecompositionDialog->combo_positionAttribute->currentText();

        MapHandler<PFP2>* mh_map = static_cast<MapHandler<PFP2>*>(m_schnapps->getMap(mapName));

        if(mh_map && m_decomposition)
        {
            PFP2::MAP* map = mh_map->getMap();
            VertexAttribute<PFP2::VEC3, PFP2::MAP> positionMap = mh_map->getAttribute<PFP2::VEC3, VERTEX>(positionName);

            Decomposition* decomposition = m_decomposition;
        }
    }
}

void Surface_WaveletDecomposition_Plugin::analyzeFromDialog()
{
    QList<QListWidgetItem*> currentItems = m_waveletDecompositionDialog->list_maps->selectedItems();
    if(!currentItems.empty() && m_decomposition)
    {
        const QString& mapName = currentItems[0]->text();
        const QString& positionName = m_waveletDecompositionDialog->combo_positionAttribute->currentText();

        MapHandler<PFP2>* mh_map = static_cast<MapHandler<PFP2>*>(m_schnapps->getMap(mapName));

        if(mh_map && m_decomposition)
        {
            PFP2::MAP* map = mh_map->getMap();
            VertexAttribute<PFP2::VEC3, PFP2::MAP> positionMap = mh_map->getAttribute<PFP2::VEC3, VERTEX>(positionName);

            //We search the last level of decomposition
            Decomposition* decomposition = m_decomposition;
            while(decomposition->getChild())
            {
                decomposition = decomposition->getChild();
            }

            if(decomposition->getCorrectionX()>1 && decomposition->getCorrectionY()>1)
            {
                //If a decomposition is still feasible
                decomposition = decomposition->addChild();
                //Création d'une nouvelle image composée uniquement des pairs
                //Création d'un tableau (matrice) des coefficients d'ondelettes des impairs
            }

        }
    }
}

MapHandlerGen* Surface_WaveletDecomposition_Plugin::initializeObject(const QString& view, const QString& filename)
{
    if(!filename.isEmpty())
    {
        MapHandlerGen* mhg_map = m_schnapps->addMap("Image", 2);
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

        int imageX = image.width(), imageY = image.height();

        Algo::Surface::Tilings::Square::Grid<PFP2> grid(*map, imageX-1, imageY-1);
        grid.embedIntoGrid(position, imageX-1, imageY-1);

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

        m_decomposition = new Decomposition(imageX, imageY, 0, HORIZONTAL);

        m_decomposition->setImage(image);

        mh_map->updateBB(position);
        mh_map->notifyAttributeModification(position);
        mh_map->notifyAttributeModification(color);
        mh_map->notifyConnectivityModification();
        m_colorVBO->updateData(color);
        m_positionVBO->updateData(position);

        m_toDraw = true;

        map->enableQuickTraversal<PFP2::MAP, VERTEX>();

        if(view==m_schnapps->getSelectedView()->getName())
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
